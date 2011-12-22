//  DDVLexiconController.m
//  XQuizIt

//  Created by Daniel Stein on Wed Feb 16 2005.
//  Copyright (c) 2005 DelDotVee. All rights reserved.

#import "DDVLexiconController.h"
#import "DDVWindowController.h"
#import "DDVLexItem.h"
#import "DDVDocument.h"

@implementation DDVLexiconController

#pragma mark INITIALIZATION

- ( id ) init
{
	self = [ super init ];
	return self;
}

- ( void ) dealloc
{
	[ searchString release ];    
	[ [ NSNotificationCenter defaultCenter ] removeObserver: self ];
	[ super dealloc ];
}

- ( void ) awakeFromNib
{
	searchString = [ [ NSString alloc ] initWithString: @"" ];
}

#pragma mark TOOLBAR

- ( BOOL ) validateToolbarItem: ( NSToolbarItem * ) theItem
{
	// Disable add and delete buttons during active filtering; no delete if no records selected
	
	if ( [ theItem action ] == @selector( insert: ) )
		return ( [ searchString length ] == 0 );
	else if ( [ theItem action ] == @selector( remove: ) )
		return ( [ [ self selectedObjects ] count ] > 0 && [ searchString length ] == 0 );
	else
		return true;
}

// Handle search field activity

- ( void ) controlTextDidChange: ( NSNotification * ) note
{
    if ( [ [ note object ] isKindOfClass: [ NSSearchField class ] ] )
	{
		[ self setSearchString: [ [ [ note object ] stringValue ] lowercaseString ] ];
		[ self rearrangeObjects ];
		
		// Now we need to update the selection status text field in the main window
		[ windowController tableViewSelectionDidChange: nil ];
	}
}

#pragma mark CONTENTS

// Credit where credit is due: http://homepage.mac.com/mmalc/CocoaExamples/controllers.html
// See mmalcolm crawford's Bindings Tutorials (Filtering Controller)

- ( NSArray * ) arrangeObjects: ( NSArray * ) lexicon
{
	DDVLexItem *item;
	NSEnumerator *lexE;
	NSMutableArray *matchedObjects = [ NSMutableArray arrayWithCapacity: [ lexicon count ] ];
	NSString *lowerSearch = [ [ self searchString ] lowercaseString ];   // case-insensitive search
	
	if ( [ lowerSearch length ] > 0 )
	{
		lexE = [ lexicon objectEnumerator ];
		while ( item = [ lexE nextObject ] )
		{
			NSString *lowerFirstL = [ [ item firstL ] lowercaseString ];
			NSString *lowerSecondL = [ [ item secondL ] lowercaseString ];
			
			if ( [ lowerFirstL rangeOfString: lowerSearch ].location != NSNotFound 
				|| [ lowerSecondL rangeOfString: lowerSearch ].location != NSNotFound )
			{
				[ matchedObjects addObject: item ];
			}
		}
		
		return [ super arrangeObjects: matchedObjects ];
	}
	
	// No filtering, just return the existing content array
	return [ super arrangeObjects: lexicon ];
}

- ( void ) insert: ( id ) sender
{
	if ( [ [  windowController lexiconTable ] selectedRow ] == -1 )
		[ super insertObject: [ super newObject ] atArrangedObjectIndex: 0 ];
	else
		[ super insert: sender ];
	
	// Update the status bar in the main window, and the topics table in the drawer
	[ windowController updateSelectionStatus ];
	[ windowController updateTopicsList ];
}

- ( void ) remove: ( id ) sender
{
	[ super remove: sender ];
	// Update the status bar in the main window, and the topics table in the drawer
	[ windowController updateSelectionStatus ];
	[ windowController updateTopicsList ];
}

#pragma mark ACCESSORS

- ( NSString * ) searchString
{
	return searchString;
}

- ( void ) setSearchString: ( NSString * ) newSearchString
{
		[ searchString autorelease ];
		searchString = [ newSearchString copy ];
}

@end