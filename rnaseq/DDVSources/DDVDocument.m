//  DDVDocument.m
//  XQuizIt

//  Created by Daniel Stein on Wed Feb 09 2005.
//  Copyright (c) 2005 DelDotVee. All rights reserved.

#import "DDVDocument.h"
#import "DDVWindowController.h"
#import "DDVLexiconModel.h"

@implementation DDVDocument

#pragma mark еее INITIALIZATION еее

- ( id ) init
{
    if ( self = [ super init ] )
	{
		lexiconModel = [ [ DDVLexiconModel alloc ] initWithDocument: self ];
    }
    return self;
}

- ( void ) dealloc
{
	[ [ NSNotificationCenter defaultCenter ] removeObserver: self ];
	[ lexiconModel release ];
	[ super dealloc ];
}

- ( void ) makeWindowControllers
{
	DDVWindowController *controller = [ [ DDVWindowController alloc ] init ];
	[ self addWindowController: controller ];
	[ controller release ];
}

#pragma mark еее CONTROLS еее

- ( BOOL ) validateMenuItem: ( NSMenuItem * ) theItem
{
	NSSearchField *sf = [ [ [ self windowControllers ] objectAtIndex: 0 ] toolbarSearchField ];
	
	// To avoid data loss, the user should never do a simple save if the table is being filtered
	if ( [ theItem action ] == @selector( saveDocument: ) )
		return [ [ sf stringValue ] length ] == 0;
	// This might be totally unnecessary and a bit of an incovenience to the user as well
	else if ( [ theItem action ] == @selector( revertDocumentToSaved: ) )
		return [ [ sf stringValue ] length ] == 0;
	else
		return true;
}

#pragma mark еее ARCHIVING еее

- ( NSData * ) dataRepresentationOfType: ( NSString * ) type
{
	DDVWindowController *winCtl = [ [ self windowControllers ] objectAtIndex: 0 ];
	NSRect frameRect = [ [ winCtl window ] frame ];
	// Credit where credit is due -> toolbarHeightForWindow: is taken directly from Apple documentation
	float toolbarHeight = [ winCtl toolbarHeightForWindow: [ winCtl window ] ];
	NSRect saveRect = NSMakeRect( frameRect.origin.x, frameRect.origin.y + toolbarHeight,
		frameRect.size.width, frameRect.size.height - toolbarHeight );
	
	// This application subsumes all model data activity within the following single model object
	// In a more complicated app, there might be several data model objects
	[ lexiconModel setFrame: saveRect ];
	
	// Get the current data from the window controller and set it in the model object
	// This model object gathers the table data and associated document-specific labeling strings
	[ lexiconModel setLexiconArray: [ winCtl lexiconArray ] ];
	[ lexiconModel setFirstLangLabel: [ [ winCtl firstLSetting ] stringValue ] ];
	[ lexiconModel setSecondLangLabel: [ [ winCtl secondLSetting ] stringValue ] ];
	[ lexiconModel setTopicLabel: [ [ winCtl topicSetting ] stringValue ] ];
	
	if ( [ type isEqualToString: @"XQuizIt Document" ] )
	{
		return [ NSArchiver archivedDataWithRootObject: [ self lexiconModel ] ];
	}
	else
		return nil;
}

- ( BOOL ) loadDataRepresentation: ( NSData * ) data ofType: ( NSString * ) type
{	
	DDVLexiconModel *model = nil;
	
	if ( [ type isEqualToString: @"XQuizIt Document" ] )
	{
		// This autoreleased object exists just long enough to load data in the model object
		model = [ NSUnarchiver unarchiveObjectWithData: data ];

		if ( model != nil )
		{
			[ lexiconModel setFirstLangLabel: [ model firstLangLabel ] ];
			[ lexiconModel setSecondLangLabel: [ model secondLangLabel ] ];
			[ lexiconModel setTopicLabel: [ model topicLabel ] ];
			[ lexiconModel setLexiconArray: [ model lexiconArray ] ];
			[ lexiconModel setFrame: [ model frame ] ];
			return YES;
		}
		else
			return NO;
	}
	else
		return NO;
}

#pragma mark еее PRINTING еее

// Uses Michael Meer's MMPrintView class to lay out a reasonable representation of the document

- ( void ) printDocument: ( id ) sender
{
	[ [ [ self windowControllers ] objectAtIndex: 0 ] printShowingPrintPanel: TRUE ];
}

#pragma mark еее ACCESSORS еее

- ( DDVLexiconModel * ) lexiconModel
{
	return lexiconModel;
}

- ( void ) setLexiconModel: ( DDVLexiconModel * ) aModel
{
	[ lexiconModel autorelease ];
	lexiconModel = [ aModel retain ];
}

@end
