// MMAppController.m
// XQuizIt

// Created by Michael Meer on Sun Feb 15 2004.
// Copyright (c) 2004 Michael Meer. All rights reserved.

#import "MMAppController.h"
#import "DDVWindowController.h"
#import "DDVPrefsController.h"
#import "defaults.h"

@implementation MMAppController

// Michael Meer used this in the original VocableTrainerX application; it survives here
// DJS added features for establishing application preferences in defaults database

+ ( void ) initialize	// Note that this is a class method
{
	NSMutableDictionary *defaults = [ NSMutableDictionary dictionary ];
	[ [ NSUserDefaults standardUserDefaults ] registerDefaults: defaults ];
}

- ( void ) dealloc
{
	[ prefsController release ];
}

- ( BOOL ) applicationShouldOpenUntitledFile: ( NSApplication * ) theApplication
{
	NSUserDefaults *defaults = [ NSUserDefaults standardUserDefaults ];
	
	if ( [ defaults boolForKey: DDVNewDocCheckBoxStateKey ] )
		return YES;
	else
		return NO;
}

// When the user creates a new document, open the language settings and inspector panels

- ( IBAction ) openNewDocument: ( id ) sender
{
	DDVWindowController *wc;
	NSDocumentController *dc = [ NSDocumentController sharedDocumentController ];
	[ dc newDocument: self ];
	wc = [ [ [ dc currentDocument ] windowControllers ] objectAtIndex: 0 ];
	[ wc openSettingsPanel: self ];
	[ wc insertLexItem: self ];
	[ [ wc inspector ] makeKeyAndOrderFront: self ];
	
	// Another way to do this
	//[ [ NSApplication sharedApplication ]
	//	sendAction: NSSelectorFromString( @"openSettingsPanel" ) to: nil from: self ];
}

// Application-wide preferences; see also code for doc settings sheet in main window controller

- ( IBAction ) showPreferences: ( id ) sender
{
	if ( !prefsController )
		prefsController = [ [ DDVPrefsController alloc ] init ];
	
	[ prefsController showWindow: self ];
}

@end
