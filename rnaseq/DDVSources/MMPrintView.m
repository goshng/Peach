//  MMPrintView.m
//  XQuizIt

//  Created by Michael Meer on Fri Feb 27 2004.
//  Copyright (c) 2004 Michael Meer. All rights reserved.

#import "MMPrintView.h"
#import "DDVLexItem.h"

@implementation MMPrintView

- ( id ) initWithLexItems: ( NSArray * ) array andTitle: ( NSString * ) newTitle printInfo: ( NSPrintInfo * ) pi
{
	normalAttributes = [ [ NSMutableDictionary alloc ] init ];
	[ normalAttributes setObject: [ NSFont fontWithName: @"Helvetica" size: 10 ] forKey: NSFontAttributeName ];
	
	boldAttributes = [ [ NSMutableDictionary alloc ] init ];
	[ boldAttributes setObject: [ NSFont fontWithName: @"Helvetica Bold" size: 10 ] forKey: NSFontAttributeName ];
	
	boldRightJustifiedAttributes = [ boldAttributes mutableCopy ];
	NSMutableParagraphStyle *rightJustifiedStyle = [ [ NSMutableParagraphStyle alloc ] init ];
	[ rightJustifiedStyle setAlignment: NSRightTextAlignment ]; 
	[ boldRightJustifiedAttributes setObject: rightJustifiedStyle forKey: NSParagraphStyleAttributeName ];
	
	title = [ newTitle retain ];
	headerHeight = 40;
	itemHeight = 15;
	
	paperSize = [ pi paperSize ];
	
	pageHeight = paperSize.height - [ pi topMargin ] - [ pi bottomMargin ];
	pageWidth = paperSize.width - [ pi leftMargin ] - [ pi rightMargin ];
	
	contentHeight = pageHeight - headerHeight;
	
	itemsPerPage = contentHeight / itemHeight;
	
	lexicon = [ array retain ];
	
	// The operator will round down
	pages = [ lexicon count ] / itemsPerPage;
	
	// Add possible remainder
	if ( ( [ lexicon count ] % itemsPerPage ) != 0 )
		pages++;
	
	// Init with the WHOLE size of the view, all the pages!
	NSRect theFrame;
	theFrame.origin = NSMakePoint( 0, 0 );
	theFrame.size.width = pageWidth;
	theFrame.size.height = pageHeight * pages;
	
	self = [ super initWithFrame: theFrame ];
	
	return self;
}

- ( BOOL ) knowsPageRange: ( NSRange * ) rptr
{
	rptr->location = 1;
	rptr->length = pages;
	return YES;
}

// Return the drawing rectangle for a particular page number

- ( NSRect ) rectForPage: ( int ) page
{
	NSRect theResult;
	theResult.size.width = pageWidth;
	theResult.size.height = pageHeight;
	theResult.origin.x = 0;	
	theResult.origin.y =  NSMaxY( [ self bounds ] ) - page * pageHeight;	
	
	return theResult;
}

- ( void ) drawRect: ( NSRect ) rect
{
	int i;
	int itemCounter;
	NSRect headerRect;
	NSRect itemRect;
	headerRect.size.height = headerHeight;
	headerRect.size.width = rect.size.width;
	headerRect.origin.x = rect.origin.x;
	itemRect.size.height = itemHeight;
	itemRect.size.width = rect.size.width / 3;
	itemRect.origin.x = rect.origin.x;
	
	for ( i = 0; i < pages; i++ )
	{
		headerRect.origin.y = NSMaxY( [ self bounds ] ) - i * pageHeight - headerHeight;
		if ( NSIntersectsRect( headerRect, rect ) )
		{			
			// Draw the header
			[ title drawInRect: headerRect withAttributes: boldAttributes ];	    
			[ [ NSString stringWithFormat: @"%@ %d / %d", NSLocalizedString( @"Page", nil ), i + 1, pages ]
				drawInRect: headerRect withAttributes: boldRightJustifiedAttributes ];	  			
			
			[ NSBezierPath strokeLineFromPoint: NSMakePoint( headerRect.origin.x, headerRect.origin.y + itemHeight + 10 )
				toPoint: NSMakePoint( headerRect.origin.x + itemRect.size.width * 3, headerRect.origin.y + itemHeight + 10 ) ];
			
			// Draw the lexicon table (modified by DS to include the topics column)
			for ( itemCounter = 0; itemCounter < itemsPerPage; itemCounter++ )
			{
				int indexOfItem = i * itemsPerPage + itemCounter;
				
				if ( indexOfItem < [ lexicon count ] )
				{
					itemRect.origin.y = headerRect.origin.y - itemCounter * itemHeight;
					[ [ [ lexicon objectAtIndex: indexOfItem ] firstL ] drawInRect: itemRect withAttributes: normalAttributes ];	    
					
					NSRect secondWordRect = itemRect;
					secondWordRect.origin.x = secondWordRect.origin.x + secondWordRect.size.width;
					[ [ [ lexicon objectAtIndex: indexOfItem ] secondL ] drawInRect: secondWordRect withAttributes: normalAttributes ];
					
					NSRect thirdWordRect = secondWordRect;
					thirdWordRect.origin.x = thirdWordRect.origin.x + thirdWordRect.size.width;
					[ [ [ lexicon objectAtIndex: indexOfItem ] score ] drawInRect: thirdWordRect withAttributes: normalAttributes ];
				}
			}
		}
	}
}					// Atsa lotsa code blocks!

- ( void ) dealloc
{
	[ lexicon release ];
	[ super dealloc ];
}

@end
