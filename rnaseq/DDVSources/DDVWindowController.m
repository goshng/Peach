//  DDVWindowController.m
//  XQuizIt

//  Created by Daniel Stein on Wed Feb 09 2005.
//  Copyright (c) 2005 DelDotVee. All rights reserved.

#import "DDVWindowController.h"
#import "DDVLexiconController.h"
#import "DDVQuizController.h"
#import "DDVDocument.h"
#import "DDVLexItem.h"
#import "DDVTopic.h"
#import "DDVLexiconModel.h"
#import "DDVSearchFieldFormatter.h"
#import "MMPrintView.h"
#import "defaults.h"

@implementation DDVWindowController

// This is the workhorse of the application, and probably contains more code than it should

#pragma mark ��� INITIALIZATION ���

- ( id ) init
{
	if ( self = [ super initWithWindowNibName: @"DDVDocument" ] )
	{
		lexiconArray = [ [ NSMutableArray alloc ] init ];
		[ self setShouldCloseDocument: YES ];   // Close entire document if main window closes
	}
	return self;
}

- ( void ) dealloc
{
	[ [ NSNotificationCenter defaultCenter ] removeObserver: self ];
	[ toolbarSearchField setDelegate: nil ];
	[ self setLexiconArray: nil ];
	[ super dealloc ];
}

- ( void ) awakeFromNib
{
	DDVLexiconModel *model = [ [ self document ] lexiconModel ];
	
	if ( [ model hasFrame ] )	// We want to restore it when user loads a documment
	{
		[ self setShouldCascadeWindows: NO ];
		[ [ self window ] setFrame: [ model frame ] display: YES ];
	}
	
	//[ topicsDrawer openOnEdge: NSMaxXEdge ];	// Could add as a setting on a per-doc basis
	
	[ self setupToolbar ];
	[ [ toolbarSearchField cell ] setWraps: NO ];
	[ [ toolbarSearchField cell ] setScrollable: YES ];
	
	// Make the search field behave nicely if the user inadvertently pastes table rows in it
	DDVSearchFieldFormatter *f = [ [ [ DDVSearchFieldFormatter alloc ] init ] autorelease ];
	[ toolbarSearchField setFormatter: f ];  // The formatter is retained by the search field

	// Manual binding prevents an object leak in Panther; is this fixed in Tiger??
	[ lexiconController bind: @"contentArray" toObject: self withKeyPath: @"lexiconArray" options: nil ];
	
	// Use slightly modified version Michael Meer's Bezier Path cell for graphical score display
	//[ [ lexiconTable tableColumnWithIdentifier: @"score" ]
	//	setDataCell: [ [ [ MMBezierPathCell alloc ] init ] autorelease ] ];
	
	// Use luscious graphical globules provided by the included image files for graphical score display
	[ [ lexiconTable tableColumnWithIdentifier: @"image" ]
		setDataCell: [ [ [ NSImageCell alloc ] init ] autorelease ] ];
	
	// Bypass undo registration when loading data from a file
	[ [ [ self document ] undoManager ] disableUndoRegistration ];
	[ lexiconController addObjects: [ model lexiconArray ] ];
	[ [ [ self document ] undoManager ] enableUndoRegistration ];
	
	// To update the interface control labels when loading data from an existing file
	[ firstLSetting setStringValue: [ [ [ self document ] lexiconModel ] firstLangLabel ] ];
	[ secondLSetting setStringValue: [ [ [ self document ] lexiconModel ] secondLangLabel ] ];
	[ topicSetting setStringValue: [ [ [ self document ] lexiconModel ] topicLabel ] ];
	
	// Pretty-up the interface after everything is loaded
	[ self updateLabelStrings ];
	[ lexiconTable deselectAll: self ];
	[ lexiconTable reloadData ];
	[ self updateTopicsList ];
	
	// Update the selection status bar
	[ self updateSelectionStatus ];

	// Shiny, happy document!
	[ [ self document ] updateChangeCount: NSChangeCleared ];
}

- ( void ) windowWillClose: ( NSNotification * ) note
{
	// Finish up by letting go of manual bindings
	[ lexiconController unbind: @"contentArray" ];
	[ lexiconController setContent: nil ];
}

- ( void ) setupToolbar
{
	NSToolbar *toolbar = [ [ NSToolbar alloc ] initWithIdentifier: @"mainToolbar" ];
	[ toolbar setDelegate: self ];
	[ toolbar setAllowsUserCustomization: YES ];	
	[ toolbar setAutosavesConfiguration: YES ];
	[ [ self window ] setToolbar: [ toolbar autorelease ] ];    
}

- ( NSString * ) windowTitleForDocumentDisplayName: ( NSString * ) displayName
{
	[ inspector setTitle: [ NSLocalizedString( @"VocableTrainerX Selection: ", nil )
		stringByAppendingString: displayName ] ];
	return [ @"VocableTrainerX: " stringByAppendingString: displayName ];
}

// This function is published in non-OO form in Apple's documentation for toolbars
// We need this for properly saving the window dimensions to restore later

- ( float ) toolbarHeightForWindow: ( NSWindow * ) window
{
	NSToolbar *toolbar;
	float toolbarHeight = 0.0;
	NSRect windowFrame;
	
	toolbar = [ window toolbar ];
	
	if ( toolbar && [ toolbar isVisible ] )
	{
		windowFrame = [ NSWindow contentRectForFrameRect: [ window frame ] styleMask: [ window styleMask ] ];
		toolbarHeight = NSHeight( windowFrame ) - NSHeight( [ [ window contentView ] frame ] );
	}
	
	return toolbarHeight;
}

// It is more than a little likely that editing a cell in the topic column will affect the topics list

- ( void ) controlTextDidEndEditing: ( NSNotification * ) note
{
	[ self updateTopicsList ];
}

// Update the status bar when the selection changes in the main table view
// A selection change results from any activity in the topics table

- ( void )  tableViewSelectionDidChange: ( NSNotification * ) note
{
	if ( [ note object ] == lexiconTable )
	{
		[ self updateSelectionStatus ];
	}
	
	if ( [ note object ] == topicTable )
	{
		// Somewhat restrictive: We cancel search field filtering if topic filtering is requested
		[ toolbarSearchField setStringValue: @"" ];
		[ lexiconController setSearchString: @"" ];
		[ lexiconController rearrangeObjects ];
		[ lexiconController setSelectedObjects: [ self lexiconSelectedByTopic ] ];
	}
}

// This creates a selection in the main table based on activity in the topics table

- ( NSArray * ) lexiconSelectedByTopic
{
	DDVLexItem *item;
	NSEnumerator *te;
	NSEnumerator *le;
	NSArray *lexicon;
	NSArray *topics;
	DDVTopic *topic;
	NSMutableArray *selection;
	
	topics = [ topicController selectedObjects ];
	lexicon = [ lexiconController arrangedObjects ];
	selection = [ NSMutableArray arrayWithCapacity: [ lexicon count ] ];
	le = [ lexicon objectEnumerator ];
	
	// Go through the main table array and find all items whose topic is selected
	while ( item = [ le nextObject ] )
	{
		// Enumerate the selected items in the topics table
		// Provide a way for explicitly selecting any items with zero-length empty topic strings
		te = [ topics objectEnumerator ];
		while ( topic = [ te nextObject ] )
		{
			if ( [ [ item topic ] isEqualToString: [ topic topicString ] ] )
			{
				[ selection addObject: item ];
			}
			if ( [ [ item topic ] isEqualToString: @"" ] && 
				[ [ topic topicString ] isEqualToString: NSLocalizedString( @"<Empty>", nil ) ] )
			{
				[ selection addObject: item ];
			}
		}
	}
	
	// This result is sent back to tableViewSelectionDidChange: as an immutable array
	return [ [ selection copy ] autorelease ];
}

// Maintain a list of all unique strings in topics column of main table to go in topics table
// Include a special provision to represent empty topic strings when inserting an item

- ( void ) updateSelectionStatus
{
	int a = [ lexiconTable numberOfRows ];
	int s = [ lexiconTable numberOfSelectedRows ];
	
	[ selectionStatus setStringValue:
		[ NSString stringWithFormat: NSLocalizedString( @"%d of %d selected", nil ), s, a ] ];
}

- ( void ) updateTopicsList
{
	DDVLexItem *item;
	NSString *aTopic;
	NSEnumerator *e;
	NSMutableArray *uniqueTopics = [ NSMutableArray array ];
	
	[ topicController removeObjects: [ topicController arrangedObjects ] ];
	
	// Interestingly, the following fails to update the topic display properly when undoing and redoing
	//e = [ [ lexiconController arrangedObjects ] objectEnumerator ];
	
	// So we use this array and it works fine
	e = [ lexiconArray objectEnumerator ];
	
	while ( item = [ e nextObject ] )
	{
		if ( ![ uniqueTopics containsObject: [ item topic ] ] && !( [ [ item topic ] isEqualToString: @"" ] ) )
			[ uniqueTopics addObject: [ item topic ] ];
		if ( ![ uniqueTopics containsObject: NSLocalizedString( @"<Empty>", nil ) ]
			&& ( [ [ item topic ] isEqualToString: @"" ] ) )
			[ uniqueTopics addObject: NSLocalizedString( @"<Empty>", nil ) ];
	}
	
	e = [ uniqueTopics objectEnumerator ];
	
	while ( aTopic = [ e nextObject ] )
		[ topicController addObject: [ [ [ DDVTopic alloc ] initWithString: aTopic ] autorelease ] ];
	
	[ topicTable reloadData ];
	[ topicTable deselectAll: self ];
}

#pragma mark ����CONTROLS ���

- ( IBAction ) toggleTopicsDrawer: ( id ) sender
{
	[ topicsDrawer toggle: self ];
}

- ( IBAction ) selectNone: ( id ) sender
{
	[ lexiconTable deselectAll: self ];
	[ topicTable deselectAll: self ];
}

- ( IBAction ) showInspector: ( id ) sender
{
	[ inspector makeKeyAndOrderFront: self ];
}

- ( IBAction ) startQuiz: ( id ) sender
{
	NSUserDefaults *defaults = [ NSUserDefaults standardUserDefaults ];
	
	NSArray *arr = [ lexiconController arrangedObjects ];
	NSArray *sel = [ lexiconController selectedObjects ];
	
	// Quiz includes only selected items or else all items in document if no prior selection exists
	if ( [ sel count ] == 0 )
		quizController = [ [ DDVQuizController alloc ] initWithWindowNibName: @"QuizWindow" quizItems: arr  ];
	else
		quizController = [ [ DDVQuizController alloc ] initWithWindowNibName: @"QuizWindow" quizItems: sel  ];
	
	[ [ self document ] addWindowController: [ quizController autorelease ] ];
	
	if ( [ defaults boolForKey: DDVQuizModeCheckBoxStateKey ] )
		[ quizController showFullScreenQuizWindow ];
	else
	{
		if ( [ inspector isVisible ] )
			[ inspector orderOut: self ];
		
		//[ quizController showWindow: self ];
		quizSession = [ NSApp beginModalSessionForWindow: [ quizController window ] ];
	}
	
	[ [ self window ] orderOut: self ];
	[ NSApp runModalSession: quizSession ];
}

- ( void ) endQuizWithFlag: ( BOOL ) answersAttempted
{
	if ( [ NSApp modalWindow ] )
		[ NSApp endModalSession: quizSession ];
	
	if ( answersAttempted )
		[ [ self document ] updateChangeCount: NSChangeDone ];
}

// These two connect a couple of menu items to control table membership

- ( IBAction ) insertLexItem: ( id ) sender
{
	[ lexiconController insert: self ];
	[ self updateSelectionStatus ];
}

- ( IBAction ) deleteLexItem: ( id ) sender
{
	[ lexiconController remove: self ];
	[ self updateSelectionStatus ];
}

- ( IBAction ) resetScores: ( id ) sender
{
	NSArray *lexicon = [ lexiconController arrangedObjects ];
	NSArray *selection = [ lexiconController selectedObjects ];
	NSEnumerator *lexEnumerator;
	DDVLexItem *item;
	
	if ( [ selection count ] == 0 )
		lexEnumerator = [ lexicon objectEnumerator ];
	else
		lexEnumerator = [ selection objectEnumerator ];
	
	// There really should be a warning before doing this
	// All scores are reset if no items have been previously selected

	while ( item = [ lexEnumerator nextObject ] )
	{
		[ item setNumGuesses: 0 ];
		[ item setNumCorrect: 0 ];
	}
	
	[ lexiconTable reloadData ];
	[ [ self document ] updateChangeCount: NSChangeDone ];
}

// Provide a handy utility to correct a commonly-occurring error in data entry

- ( IBAction ) swapLanguages: ( id ) sender
{
	DDVLexItem *item;
	NSArray *lexicon = [ lexiconController selectedObjects ];
	NSEnumerator *lexEnumerator = [ lexicon objectEnumerator ];
	while ( item = [ lexEnumerator nextObject ] )
	{
		[ item swapLanguages ];
	}
	[ [ self document ] updateChangeCount: NSChangeDone ];    
}

// Let the user reclassify multiple selections of vocabulary items all at once

- ( IBAction ) categorize: ( id ) sender
{
	[ self openCategorizeSheet ];
}

// This cleverly-hidden feature allows the user to cancel topic selections with a single mouse click

- ( void ) tableView: ( NSTableView * ) aTableView mouseDownInHeaderOfTableColumn: ( NSTableColumn * ) column
{
	if ( aTableView == topicTable )
		[ topicTable deselectAll: self ];
}

// Don't start a quiz or try to export if nothing is present in the main table yet

- ( BOOL ) validateToolbarItem: ( NSToolbarItem * ) theItem
{
	if ( [ theItem action ] == @selector( startQuiz: ) )
		return [ [ lexiconController arrangedObjects ] count ] > 0;
	if ( [ theItem action ] == @selector( exportCSV: ) )
		return [ [ lexiconController arrangedObjects ] count ] > 0;
	else
		return true;
}

// Mainly to prevent modifications of the data array if search field filtering is in progress
// Also gives reasonable menu item feedback when certain actions do not make sense

- ( BOOL ) validateMenuItem: ( NSMenuItem * ) theItem
{
	if ( [ theItem action ] == @selector( exportCSV: ) )
		return [ [ lexiconController arrangedObjects ] count ] > 0;
	if ( [ theItem action ] == @selector( importCSV: ) )
		return [ [ toolbarSearchField stringValue ] length ] == 0;
	if ( [ theItem action ] == @selector( insertLexItem: ) )
		return [ [ toolbarSearchField stringValue ] length ] == 0;
	if ( [ theItem action ] == @selector( deleteLexItem: ) )
		return ( [ [ lexiconController selectedObjects ] count ] > 0 &&
		[ [ toolbarSearchField stringValue ] length ] == 0 );
	if ( [ theItem action ] == @selector( resetScores: ) )
		return [ [ lexiconController arrangedObjects ] count ] > 0;
	if ( [ theItem action ] == @selector( swapLanguages: ) )
		return [ [ lexiconController selectedObjects ] count ] > 0;
	if ( [ theItem action ] == @selector( categorize: ) )
		return [ [ lexiconController selectedObjects ] count ] > 0;
	if ( [ theItem action ] == @selector( selectNone: ) )
		return [ [ lexiconController selectedObjects ] count ] > 0;
	else
		return TRUE;
}

#pragma mark ��� SHEETS ���

// So far the settingsPanel just allows customization of table header cells
// There are labels elsewhere in the interface that are coordinated with the header cell values

- ( IBAction ) openSettingsPanel: ( id ) sender
{
	[ NSApp beginSheet: settingsPanel modalForWindow: [ self window ] modalDelegate: nil
		didEndSelector: nil contextInfo: nil ];
	[ NSApp runModalForWindow: settingsPanel ];
	
	[ NSApp endSheet: settingsPanel ];
	[ settingsPanel orderOut: self ];
}

- ( IBAction ) dismissSettingsPanel: ( id ) sender
{
	[ self updateLabelStrings ];
	[ lexiconTable reloadData ];
	[ NSApp stopModal ];
}

- ( void ) updateLabelStrings
{
	if ( ![ [ [ firstLColumn  headerCell ] stringValue ] isEqual: [ firstLSetting stringValue ] ] )
		[ [ self document ] updateChangeCount: NSChangeDone ];
	if ( ![ [ [ secondLColumn  headerCell ] stringValue ] isEqual: [ secondLSetting stringValue ] ] )
		[ [ self document ] updateChangeCount: NSChangeDone ];
	if ( ![ [ [ topicColumn  headerCell ] stringValue ] isEqual: [ topicSetting stringValue ] ] )
		[ [ self document ] updateChangeCount: NSChangeDone ];

	[ [ firstLColumn headerCell ] setStringValue: [ firstLSetting stringValue ] ];
	[ [ secondLColumn headerCell ] setStringValue: [ secondLSetting stringValue ] ];
	[ [ topicColumn headerCell ] setStringValue: [ topicSetting stringValue ] ];
	
	[ firstLangLabel setStringValue: [ firstLSetting stringValue ] ];
	[ secondLangLabel setStringValue: [ secondLSetting stringValue ] ];
	[ topicLabel setStringValue: [ topicSetting stringValue ] ];
}

// Implements mass topic string assignment - this app is now really full of sheet

- ( void ) openCategorizeSheet
{
	[ NSApp beginSheet: categorizePanel modalForWindow: [ self window ] modalDelegate: nil
		didEndSelector: nil contextInfo: nil ];
	[ NSApp runModalForWindow: categorizePanel ];
	
	[ NSApp endSheet: categorizePanel ];
	[ categorizePanel orderOut: self ];
}

- ( IBAction ) dismissCategorizeSheet: ( id ) sender
{
	NSArray *selection = [ lexiconController selectedObjects ];
	NSEnumerator *e = [ selection objectEnumerator ];
	DDVLexItem *item;
	
	while ( item = [ e nextObject ] )
		[ item setTopic: [ categoryField stringValue ] ];
	
	[ [ self document ] updateChangeCount: NSChangeDone ];
	
	[ self updateTopicsList ];
	
	[ NSApp stopModal ];
}

- ( IBAction ) cancelCategorizeAction: ( id ) sender
{
	[ NSApp stopModal ];
}

#pragma mark ��� IMPORT / EXPORT ���

// Typical spreadsheets can handle comma-separated values (CSV)
// Credit where credit is due: This is Michael Meer's work from the original VocableTrainerX

- ( IBAction ) exportCSV: ( id ) sender
{
	NSSavePanel *panel = [ NSSavePanel savePanel ];
	
	[ panel setRequiredFileType: @"csv" ];
	
	if ( [ panel runModal ] == NSOKButton )
		[ [ self getCSVData ] writeToFile: [ panel filename ] atomically: YES ];
}

- ( NSString * ) getCSVData
{
	NSMutableString *csvString = [ NSMutableString string ];
	NSEnumerator *lexE = [ [ lexiconController arrangedObjects ] objectEnumerator ];
	
	DDVLexItem *lex;
	
	while ( lex = [ lexE nextObject ] )
	{
		[ csvString appendString: [ NSString stringWithFormat: @"\"%@\",\"%@\",\"%@\"\n",
			[ lex firstL ], [ lex secondL ], [ lex topic ] ] ];
	}
	return [ NSString stringWithString: csvString ];
}

- ( IBAction ) importCSV: ( id ) sender
{
	NSMutableString *line;
	NSString *lineBuffer;
	NSString *rawString;
	NSArray *lines;
	NSArray *components;
	NSEnumerator *lineEnum;
	DDVLexItem *lex;
	
	NSOpenPanel *panel = [ NSOpenPanel openPanel ];
	
	if ( [ panel runModalForTypes: nil ] == NSOKButton )
		rawString = [ NSString stringWithContentsOfFile: [ [ panel filenames ] objectAtIndex: 0 ] ];
	else
		return;
		
	lines = [ rawString componentsSeparatedByString: @"\n" ];
	lineEnum = [ lines objectEnumerator ];
	
	while ( lineBuffer = [ lineEnum nextObject ] )
	{
		line = [ NSMutableString stringWithString: lineBuffer ];
		[ line replaceOccurrencesOfString: @"\"" withString: @"" options: NSLiteralSearch
			range: NSMakeRange( 0, [ line length ] ) ];
		components = [ line componentsSeparatedByString: @"," ];
		
		if ( [ components count ] < 3 )
			continue;
		else
		{
			lex = [ [ DDVLexItem alloc ] initWithFirstLanguage: [ components objectAtIndex: 0 ]
				secondLanguage: [ components objectAtIndex: 1 ]
				topic: [ components objectAtIndex: 2 ] ];
			[ lexiconController addObject: [ lex autorelease ] ];
		}
	}
	[ [ self document ] updateChangeCount: NSChangeDone ];
	
	[ self updateTopicsList ];

	[ lexiconTable deselectAll: self ];
}

#pragma mark ��� ARRAY ACCESS ���

// Indexed accessors for the main data array; recommended by tutorials that introduce Bindings
// In this case, it is part of Key Value Observing (KVO) to support the use of Undo Manager

- ( unsigned int ) countOfLexiconArray
{
	return [ lexiconArray count ];
}

- ( id ) objectInLexiconArrayAtIndex: ( unsigned int ) index
{
   return [ lexiconArray objectAtIndex: index ];
}

- ( void ) insertObject: ( id ) anObject inLexiconArrayAtIndex: ( unsigned int ) index
{
	// Add the inverse of this operation to the undo stack
	NSUndoManager *undo = [ [ self document ] undoManager ];
	[ [ undo prepareWithInvocationTarget: self ] removeObjectFromLexiconArrayAtIndex: index ];
	if ( ![ undo isUndoing ] )
		[ undo setActionName: NSLocalizedString( @"Insert Item", nil ) ];
	// Add the item to the array
	[ self startObservingLexItem: anObject ];
	[ lexiconArray insertObject: anObject atIndex: index ];
	[ self updateTopicsList ];
	[ lexiconTable deselectAll: self ];
}

- ( void ) removeObjectFromLexiconArrayAtIndex: ( unsigned int ) index
{
	DDVLexItem *item = [ lexiconArray objectAtIndex: index ];
	NSUndoManager *undo = [ [ self document ] undoManager ];
	// Add the inverse of this operation to the undo stack
	[ [ undo prepareWithInvocationTarget: self ] insertObject: item inLexiconArrayAtIndex: index ];
	if ( ![ undo isUndoing ] )
		[ undo setActionName: NSLocalizedString( @"Delete Item", nil ) ];
	// Remove the item from the array
	[ self stopObservingLexItem: item ];
	[ lexiconArray removeObjectAtIndex: index ];
	[ self updateTopicsList ];
	[ lexiconTable deselectAll: self ];
}

// Not currently used, but possible utility when pasting into previous table selection

- ( void ) replaceObjectInLexiconArrayAtIndex: ( unsigned int ) index withObject: ( id ) item
{
    [ lexiconArray replaceObjectAtIndex: index withObject: item ];
	[ self updateTopicsList ];
}

#pragma mark ����KEY-VALUE OBSERVING ���

// Credit where credit is due: Aaron Hillegass presents the basics for KVO with the Undo Manager
// The original example appears in his book "Cocoa Programming for Max OS X" (Second Edition)
// The publisher is Addison-Wesley (2004)

- ( void ) startObservingLexItem: ( DDVLexItem * ) item
{
	[ item addObserver: self forKeyPath: @"firstL" options: NSKeyValueObservingOptionOld context: NULL ];
	[ item addObserver: self forKeyPath: @"secondL" options: NSKeyValueObservingOptionOld context: NULL ];
	[ item addObserver: self forKeyPath: @"topic" options: NSKeyValueObservingOptionOld context: NULL ];
}

- ( void ) stopObservingLexItem: ( DDVLexItem * ) item
{
	[ item removeObserver: self forKeyPath: @"firstL" ];
	[ item removeObserver: self forKeyPath: @"secondL" ];
	[ item removeObserver: self forKeyPath: @"topic" ];
}

- ( void ) changeKeyPath: ( NSString * ) keyPath ofObject: ( id ) obj toValue: ( id ) newValue
{
	// This method is its own inverse operation for maintaining undo stack with edits
	// The message setValue: forKeyPath: invokes KVO, and takes care of undoing
	[ obj setValue: newValue forKeyPath: keyPath ];
	[ self updateTopicsList ];
}

- ( void ) observeValueForKeyPath: ( NSString * ) kp ofObject: ( id ) obj change: ( NSDictionary * ) change context: ( void * ) con
{
	NSUndoManager *undo = [ [ self document ] undoManager ];
	id oldValue = [ change objectForKey: NSKeyValueChangeOldKey ];
	[ [ undo prepareWithInvocationTarget: self ] changeKeyPath: kp ofObject: obj toValue: oldValue ];
	[ undo setActionName: NSLocalizedString( @"Edit", nil ) ];
}

#pragma mark ��� ACCESSORS ���

- ( NSPanel * ) inspector
{
	return inspector;
}

- ( NSTableView * ) lexiconTable
{
	return lexiconTable;
}

- ( NSArray * ) lexiconArray
{
	return [ [  lexiconArray copy ] autorelease ];	// Returning an immutable array by default
}

// This accessor is crucial to ensure proper KVO for array elements

- ( void ) setLexiconArray: ( NSMutableArray * ) array
{
	NSEnumerator *e;
	DDVLexItem *item;
	
	if ( array == lexiconArray )
		return;
	
	// For editing individual cells
	e = [ lexiconArray objectEnumerator ];
	while ( item = [ e nextObject ] )
		[ self stopObservingLexItem: item ];
	
	[ lexiconArray autorelease ];
	lexiconArray = [ array retain ];
	e = [ lexiconArray objectEnumerator ];
	while ( item = [ e nextObject ] )
		[ self startObservingLexItem: item ];
}

- ( NSSearchField * ) toolbarSearchField
{
	return toolbarSearchField;
}

- ( NSTextField * ) firstLSetting
{
	return firstLSetting;
}

- ( NSTextField * ) secondLSetting
{
	return secondLSetting;
}

- ( NSTextField * ) topicSetting
{
	return topicSetting;
}

// Credit where credit is due: This uses Michael Meer's original printing view, with minor mods
// See the actual MMPrintView class for the details of laying out the print view

#pragma mark ��� PRINTING ���

- ( void ) printShowingPrintPanel: ( BOOL ) showPanels
{
	NSPrintOperation *op;
	
	// Obtain a custom view that will be printed
	NSMutableString *title = [ NSMutableString stringWithString: [ [ self document ] displayName ] ];
	
	MMPrintView *printView = [ [ MMPrintView alloc ] initWithLexItems: [ lexiconController arrangedObjects ]
		andTitle: title printInfo: [ [ self document ] printInfo ] ];
	
	// Construct the print operation and setup Print panel
	op = [ NSPrintOperation printOperationWithView: printView printInfo: [ [ self document ] printInfo ] ];
	[ op setShowPanels: showPanels ];
	
	// Run operation, which shows the Print panel if showPanels was YES
	
	//[ myDocument runModalPrintOperation: op delegate: nil didRunSelector: NULL contextInfo: NULL ];
	[ op runOperation ];
	[ printView release ];
}

@end