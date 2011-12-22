//  DDVToolbarDelegateCategory.h
//  XQuizIt

//  Created by Daniel Stein on Wed Feb 9 2005.
//  Copyright (c) 2005 DelDotVee. All rights reserved.

#import <Cocoa/Cocoa.h>

#import "DDVWindowController.h"

// Oh, goodness! Not ANOTHER garish toolbar in a document-based app!

@interface DDVWindowController ( ToolbarDelegateCategory )

#pragma mark ¥¥¥ÊTOOLBAR ¥¥¥

- ( NSToolbarItem * ) toolbar: ( NSToolbar * ) toolbar
	itemForItemIdentifier: ( NSString * ) itemIdentifier willBeInsertedIntoToolbar: ( BOOL ) flag;
- ( NSArray * ) toolbarAllowedItemIdentifiers: ( NSToolbar * ) toolbar;
- ( NSArray * ) toolbarDefaultItemIdentifiers: ( NSToolbar * ) toolbar;

@end
