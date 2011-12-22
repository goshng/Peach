//  DDVQuizController.h
//  XQuizIt

//  Created by Daniel Stein on Fri Mar 4 2005.
//  Copyright (c) 2005 DelDotVee. All rights reserved.

#import "DDVPrefsController.h"
#import "defaults.h"

@implementation DDVPrefsController

- ( id ) init
{
	if ( self = [ super initWithWindowNibName: @"Preferences" ] )
	{
	}
	
	return self;
}

// This updates the state of the checkbox in the Prefs window from the database

- ( void ) windowDidLoad
{
	NSUserDefaults *defaults;
	defaults = [ NSUserDefaults standardUserDefaults ];
	
	[ newDocCheckBox setState: [ defaults boolForKey: DDVNewDocCheckBoxStateKey ] ];
	[ quizModeCheckBox setState: [ defaults boolForKey: DDVQuizModeCheckBoxStateKey ] ];
	[ swapLangCheckBox setState: [ defaults boolForKey: DDVSwapLangCheckBoxStateKey ] ];
}

// These update the default values when the user changes checkbox states in the Preferences window

- ( IBAction ) changeNewDocCBFlag: ( id ) sender
{
	NSUserDefaults *defaults = [ NSUserDefaults standardUserDefaults ];
	[ defaults setObject: [ NSNumber numberWithBool: [ sender state ] ]
		forKey: DDVNewDocCheckBoxStateKey ];
}

// This one toggles the app between normal quiz mode and full-screen quiz mode

- ( IBAction ) changeQuizModeCBFlag: ( id ) sender
{
	NSUserDefaults *defaults = [ NSUserDefaults standardUserDefaults ];
	[ defaults setObject: [ NSNumber numberWithBool: [ sender state ] ]
		forKey: DDVQuizModeCheckBoxStateKey ];
}

- ( IBAction ) changeSwapLangCBFlag: ( id ) sender
{
	NSUserDefaults *defaults = [ NSUserDefaults standardUserDefaults ];
	[ defaults setObject: [ NSNumber numberWithBool: [ sender state ] ]
		forKey: DDVSwapLangCheckBoxStateKey ];
}

@end
