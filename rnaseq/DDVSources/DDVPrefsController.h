//  DDVQuizController.h
//  XQuizIt

//  Created by Daniel Stein on Fri Mar 4 2005.
//  Copyright (c) 2005 DelDotVee. All rights reserved.

#import <Cocoa/Cocoa.h>

@interface DDVPrefsController: NSWindowController
{
	IBOutlet NSButton *newDocCheckBox;
	IBOutlet NSButton *quizModeCheckBox;
	IBOutlet NSButton *swapLangCheckBox;
}

- ( IBAction ) changeNewDocCBFlag: ( id ) sender;
- ( IBAction ) changeQuizModeCBFlag: ( id ) sender;
- ( IBAction ) changeSwapLangCBFlag: ( id ) sender;

@end
