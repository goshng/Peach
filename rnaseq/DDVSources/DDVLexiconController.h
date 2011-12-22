//  DDVLexiconController.h
//  XQuizIt

//  Created by Daniel Stein on Wed Feb 16 2005.
//  Copyright (c) 2005 DelDotVee. All rights reserved.

#import <Cocoa/Cocoa.h>

// For the core filtering capabilities of this array controller, credit where credit is due:
// See mmalcolm crawford's Bindings Tutorials (Filtering Controller)
// http://homepage.mac.com/mmalc/CocoaExamples/controllers.html

@class DDVWindowController;

@interface DDVLexiconController: NSArrayController
{
	NSString *searchString;
	IBOutlet DDVWindowController *windowController;
}

- ( NSString * ) searchString;
- ( void ) setSearchString: ( NSString * ) newSearchString;

@end
