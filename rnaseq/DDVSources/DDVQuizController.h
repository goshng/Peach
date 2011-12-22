//  DDVQuizController.h
//  XQuizIt

//  Created by Daniel Stein on Sat Feb 19 2005.
//  Copyright (c) 2005 DelDotVee. All rights reserved.

#import <Cocoa/Cocoa.h>

@class MMRandomEnumerator;	// Michael Meer's random enumerator class
@class DDVLexItem;
@class DDVQuizWindow;

@interface DDVQuizController: NSWindowController
{
	DDVQuizWindow *ddvQuizWindow;		// Full-screen quiz window is only created programmatically
	IBOutlet NSView *quizView;			// Layout containing the actual quiz controls for any quiz window
	
	IBOutlet NSTextField *firstLangLabel;		// Document values control these labels
	IBOutlet NSTextField *secondLangLabel;
	
	IBOutlet NSTextField *promptField;
	IBOutlet NSTextField *solutionField;
	
	IBOutlet NSTextField *feedback;				// Feedback to the user when an answer is attempted
	IBOutlet NSButton *skipButton;
	
	DDVLexItem *item;
	NSNumber *ranseq;						// These are created by the MMRandomEnumerator class
	NSArray *quizItems;
	MMRandomEnumerator *quizEnum;
	
	NSColor *rightColor;					// Color cues for the feedback string	
	NSColor *wrongColor;
	BOOL isCorrect;
	BOOL attemptedAnswer;		// Flags the document-changed status if user answers any quiz items
	
	int swapLangFlag;
	BOOL doSwapPrefsFlag;
	NSUserDefaults *defaults;
}

#pragma mark INITIALIZATION

- ( id ) initWithWindowNibName: ( NSString * ) windowNibName quizItems: ( NSArray * ) items;

#pragma mark QUIZ

- ( void ) showFullScreenQuizWindow;

- ( IBAction ) solve: ( id ) sender;		// Quiz interface solution text field action
- ( IBAction ) skip: ( id ) sender;			// Quiz interface button action

- ( void ) showNextItem;					// After an answer is given or user presses Skip button

- ( IBAction ) terminate: ( id ) sender;

// Format the correct solution for display with color cues as to whether guess is correct or not

- ( NSMutableAttributedString * ) feedbackString: ( NSString * ) string;

#pragma mark ACCESSORS

- ( void ) setQuizItems: ( NSArray * ) items;

@end
