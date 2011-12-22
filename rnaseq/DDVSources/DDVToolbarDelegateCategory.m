//  DDVToolbarDelegateCategory.h
//  XQuizIt

//  Created by Daniel Stein on Wed Feb 9 2005.
//  Copyright (c) 2005 DelDotVee. All rights reserved.

#import "DDVToolbarDelegateCategory.h"

@implementation DDVWindowController ( DDVToolbarDelegateCategory )

#pragma mark ¥¥¥ÊTOOLBAR ¥¥¥

- ( NSArray * ) toolbarAllowedItemIdentifiers: ( NSToolbar * ) toolbar
{
	return [ NSArray arrayWithObjects:  NSToolbarSeparatorItemIdentifier,
										NSToolbarSpaceItemIdentifier,
										NSToolbarFlexibleSpaceItemIdentifier,
										NSToolbarCustomizeToolbarItemIdentifier, 
										@"InsertItem", @"DeleteItem", @"ExportItem", @"DrawerItem", 
										@"InfoItem", @"SearchItem", @"QuizItem", nil ];
}

- ( NSArray * ) toolbarDefaultItemIdentifiers: ( NSToolbar * ) toolbar
{
	return [ NSArray arrayWithObjects:  @"InsertItem",
										@"DeleteItem",
										@"ExportItem",
										@"InfoItem",
										@"DrawerItem",
										@"SearchItem",
										@"QuizItem",
										NSToolbarFlexibleSpaceItemIdentifier,
										NSToolbarCustomizeToolbarItemIdentifier, nil ];
}

- ( NSToolbarItem * ) toolbar: ( NSToolbar * ) toolbar
	itemForItemIdentifier: ( NSString * ) itemIdentifier willBeInsertedIntoToolbar: ( BOOL ) flag
{
	NSToolbarItem *item = [ [ NSToolbarItem alloc ] initWithItemIdentifier: itemIdentifier ];
	
	if ( [ itemIdentifier isEqualToString: @"InsertItem" ] )
	{
		[ item setLabel: NSLocalizedString( @"Insert Item", nil ) ];
		[ item setPaletteLabel: [ item label ] ];
		[ item setImage: [ NSImage imageNamed: @"insert.png" ] ];
		[ item setTarget: lexiconController ];
		[ item setAction: @selector( insert: ) ];
    }
	else if ( [ itemIdentifier isEqualToString: @"DeleteItem" ] )
	{
		[ item setLabel: NSLocalizedString( @"Delete Item", nil ) ];
		[ item setPaletteLabel: [ item label ] ];
		[ item setImage: [ NSImage imageNamed: @"delete.png" ] ];
		[ item setTarget: lexiconController ];
		[ item setAction: @selector( remove: ) ];
    }
	else if ( [ itemIdentifier isEqualToString: @"ExportItem" ] )
	{
		[ item setLabel: NSLocalizedString( @"Export Data", nil ) ];
		[ item setPaletteLabel: [ item label ] ];
		[ item setImage: [ NSImage imageNamed: @"echelon.png" ] ];
		[ item setTarget: self ];
		[ item setAction: @selector( exportCSV: ) ];
    }
	else if ( [ itemIdentifier isEqualToString: @"InfoItem" ] )
	{
		[ item setLabel: NSLocalizedString( @"Info", nil ) ];
		[ item setPaletteLabel: [ item label ] ];
		[ item setImage: [ NSImage imageNamed: @"info.gif" ] ];
		[ item setTarget: self ];
		[ item setAction: @selector( showInspector: ) ];
    }
	else if ( [ itemIdentifier isEqualToString: @"QuizItem" ] )
	{
		[ item setLabel: NSLocalizedString( @"Start Quiz", nil ) ];
		[ item setPaletteLabel: [ item label ] ];
		[ item setImage: [ NSImage imageNamed: @"questions.png" ] ];
		[ item setTarget: self ];
		[ item setAction: @selector( startQuiz: ) ];
    }
	else if ( [ itemIdentifier isEqualToString: @"DrawerItem" ] )
	{
		[ item setLabel: NSLocalizedString( @"Topics Drawer", nil ) ];
		[ item setPaletteLabel: [ item label ] ];
		[ item setImage: [ NSImage imageNamed: @"drawer.tiff" ] ];
		[ item setTarget: self ];
		[ item setAction: @selector( toggleTopicsDrawer: ) ];
    }
	else if ( [ itemIdentifier isEqualToString: @"SearchItem" ] )
	{
		[ item setLabel: NSLocalizedString( @"Search Data", nil ) ];
		[ item setPaletteLabel: [ item label ] ];
		[ item setView: toolbarSearchField ];
		[ item setMinSize: [ toolbarSearchField bounds ].size ];
		[ item setMaxSize: [ toolbarSearchField bounds ].size ];	
    }
	
    return [ item autorelease ];
}

@end
