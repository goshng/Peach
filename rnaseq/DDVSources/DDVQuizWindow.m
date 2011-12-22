//  DDVQuizWindow.m
//  XQuizIt

//  Created by Daniel Stein on Tue Mar 8 2005.
//  Copyright (c) 2005 DelDotVee. All rights reserved.

#import "DDVQuizWindow.h"

@implementation DDVQuizWindow

// This is needed strictly for our custom full-screen borderless quiz window
// The user must be able to interact with controls in its content view via the keyboard

- ( BOOL ) canBecomeKeyWindow
{
	return YES;
}

@end
