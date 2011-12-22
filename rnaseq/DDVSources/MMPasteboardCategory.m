//  MMPasteboardCategory.m
//  VocableTrainerX

//  Created by Michael Meer on Fri Feb 20 2004.
//  Copyright (c) 2004 Michael Meer. All rights reserved.

#import "DDVDocument.h"
#import "DDVLexItem.h"
#import "DDVLexiconController.h"
#import "MMPasteboardCategory.h"

@implementation DDVWindowController ( MMPasteboardCategory )

#pragma mark PASTEBOARD

- ( IBAction ) cut: ( id ) sender
{
	NSArray *selection = [ lexiconController selectedObjects ];
	NSEnumerator *se = [ selection objectEnumerator ];
	DDVLexItem *item;
	
	[ self copy: sender ];
	
	while ( item = [ se nextObject ] )
	{
		// It is just too damned complicated to deal with undoing cut-and-paste if filtering is active
		[ [ [ self document ] undoManager ] disableUndoRegistration ];
		[ lexiconController removeObject: item ];
		[ [ [ self document ] undoManager ] enableUndoRegistration ];
	}
	
	[ lexiconTable reloadData ];
	[ lexiconTable deselectAll: self ];
	[ [ self document ] updateChangeCount: NSChangeDone ];
	[ self updateTopicsList ];
}

- ( IBAction ) copy: ( id ) sender
{
	NSPasteboard *pb = [ NSPasteboard generalPasteboard ];
	[ self writeStringToPboard: pb ];
}

- ( void ) writeStringToPboard: ( NSPasteboard * ) pb
{
	if ( [ lexiconTable numberOfSelectedRows ] < 1 )
	{
		NSBeep( );
	} 
	
	// Declare types
	[ pb declareTypes: [ NSArray arrayWithObjects: NSStringPboardType, NSTabularTextPboardType, @"MMPboardType", nil ] owner: self ];
	
	// Copy data to the pasteboard
	[ pb setString: [ self prepareNSStringForPboard:
		[ lexiconController selectedObjects ] ] forType: NSStringPboardType ];
	[ pb setString: [ self prepareNSStringForPboard:
		[ lexiconController selectedObjects ] ] forType: NSTabularTextPboardType ];
	[ pb setString: [ self prepareMMPboardType:
		[ lexiconController selectedObjects ] ] forType: @"MMPboardType" ];
}

// This type is used for pasting to other programs, like Excel, TextEdit and so on

- ( NSString * ) prepareNSStringForPboard: ( NSArray * ) lexicon
{
	NSEnumerator *e = [ lexicon objectEnumerator ];
	NSMutableString *result = [ NSMutableString string ];
	DDVLexItem *item;
	
	while ( item = [ e nextObject ] )
	{
		[ result appendString: [ NSString stringWithFormat: @"%@\t%@\t%@", [ item firstL ], [ item secondL ], [ item topic ] ] ];
		[ result appendString: [ NSString stringWithFormat: @"%c", NSCarriageReturnCharacter ] ];
	}
	
	return result;
}

// Custom pasteboard format for copying inside XQuizIt, preserves the score data for each item
// Fixes a bug from version 0.7 days which inserted a space at the beginning of some fields

- ( NSString * ) prepareMMPboardType: ( NSArray * ) lexicon
{
	NSEnumerator *e = [ lexicon objectEnumerator ];
	NSMutableString *result = [ NSMutableString string ];
	DDVLexItem *item;
	
	while ( item = [ e nextObject ] )
	{
		[ result appendString: [ NSString stringWithFormat: @"%@\t%@\t%@\t%d\t%d", [ item firstL ],
			[ item secondL ], [ item topic ], [ item numCorrect ], [ item numGuesses ] ] ];
		[ result appendString: [ NSString stringWithFormat: @"%c", NSCarriageReturnCharacter ] ];
	}
	
	return result;
}

- ( IBAction ) paste: ( id ) sender
{
	NSPasteboard *pb = [ NSPasteboard generalPasteboard ];
	BOOL readSucceeded;
	readSucceeded = [ self readStringFromPasteboard: pb ];
	if ( readSucceeded )
	{
		[ [ self document ] updateChangeCount: NSChangeDone ];
		[ self updateTopicsList ];
	}
	else
	{
		NSBeep( );
	}
}

- ( BOOL ) readStringFromPasteboard: ( NSPasteboard * ) pb
{
	NSArray *types = [ pb types ];
	NSMutableArray *result;
	NSEnumerator *e;
	DDVLexItem *item;
	int selectedRowIndex = [ lexiconTable selectedRow ];
	
	if ( [ types containsObject: @"MMPboardType" ] )
	{
		result = [ self readLexItemStringFromPboard: pb ];
	}
	else if ( [ types containsObject: NSStringPboardType ] )
	{
		result = [ self readNSStringFromPboard: pb ];
	}
	else
	{ 
		return NO; 
	}
	
	if ( selectedRowIndex == -1 )
	{
		[ [ [ self document ] undoManager ] disableUndoRegistration ];
		[ lexiconController addObjects: result ];
		[ [ [ self document ] undoManager ] enableUndoRegistration ];
	}
	else
	{
		// Insert the vocables at the selected row
		e = [ result reverseObjectEnumerator ];
		
		while ( item = [ e nextObject ] )
		{
			[ [ [ self document ] undoManager ] disableUndoRegistration ];
			[ lexiconController insertObject: item atArrangedObjectIndex: selectedRowIndex ];
			[ [ [ self document ] undoManager ] enableUndoRegistration ];
		}	
	}
	
	return YES;
}

- ( NSMutableArray * ) readNSStringFromPboard: ( NSPasteboard * ) pb
{
	NSString *value;
	NSArray *lines;
	NSEnumerator *linesEnumerator;
	NSString *line;
	NSArray *words;
	DDVLexItem *newItem;
	NSMutableArray *result = [ NSMutableArray array ];
	
	value = [ pb stringForType: NSStringPboardType ];
	lines = [ value componentsSeparatedByString: [ NSString stringWithFormat: @"%c", NSCarriageReturnCharacter ] ];
	
	linesEnumerator = [ lines objectEnumerator ];
	while ( line = [ linesEnumerator nextObject ] )
	{
		words = [ line componentsSeparatedByString: @"\t" ];
		if ( [ words count ] >= 2 )
		{
			newItem = [ [ DDVLexItem alloc ] initWithFirstLanguage: [ words objectAtIndex: 0 ]
				secondLanguage: [ words objectAtIndex: 1 ] topic: NSLocalizedString( @"None", nil ) ];
			if ( [ words count ] > 2 && [ words objectAtIndex: 2 ] != @"" )
				[ newItem setTopic: [ words objectAtIndex: 2 ] ];

			[ result addObject: [ newItem autorelease ] ];
		}
	}
	return result;
}

- ( NSMutableArray * ) readLexItemStringFromPboard: ( NSPasteboard * ) pb
{
	NSString *value;
	NSArray *lines;
	NSEnumerator *linesEnumerator;
	NSString *line;
	NSArray *words;
	DDVLexItem *newItem;
	NSMutableArray *result = [ NSMutableArray array ];
	
	value = [ pb stringForType: @"MMPboardType" ];
	lines = [ value componentsSeparatedByString:
		[ NSString stringWithFormat: @"%c", NSCarriageReturnCharacter ] ];
	
	linesEnumerator = [ lines objectEnumerator ];
	while ( line = [ linesEnumerator nextObject ] )
	{
		words = [ line componentsSeparatedByString: @"\t" ];
		if ( [ words count ] == 5 )
		{
			newItem = [ [ DDVLexItem alloc ] init ];
			[ newItem setFirstL: [ words objectAtIndex: 0 ] ];
			[ newItem setSecondL: [ words objectAtIndex: 1 ] ];
			[ newItem setTopic: [ words objectAtIndex: 2 ] ];
			[ newItem setNumCorrect: [ [ words objectAtIndex: 3 ] intValue ] ];
			[ newItem setNumGuesses: [ [ words objectAtIndex: 4 ] intValue ] ];		
			[ result addObject: [ newItem autorelease ] ];
		}
	}
	return result;
}

@end
