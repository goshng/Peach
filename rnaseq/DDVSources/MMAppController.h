// MMAppController.h
// XQuizIt

//  Created by Michael Meer on Sun Feb 15 2004.
//  Copyright (c) 2004 Michael Meer. All rights reserved.

#import <Cocoa/Cocoa.h>

@class DDVPrefsController;

@interface MMAppController: NSObject
{
	DDVPrefsController *prefsController;
}

- ( BOOL ) applicationShouldOpenUntitledFile: ( NSApplication * ) theApplication;
- ( IBAction ) openNewDocument: ( id ) sender;
- ( IBAction ) showPreferences: ( id ) sender;

@end
