//  DDVQuizController.m
//  XQuizIt

//  Created by Daniel Stein on Sat Feb 19 2005.
//  Copyright (c) 2005 DelDotVee. All rights reserved.

#import "DDVLexItem.h"
#import "DDVDocument.h"
#import "DDVLexiconModel.h"
#import "DDVQuizController.h"
#import "DDVQuizWindow.h"
#import "DDVWindowController.h"
#import "MMRandomEnumerator.h"
#import "defaults.h"

@implementation DDVQuizController

#pragma mark INITIALIZATION

- ( id ) init
{
	return [ self initWithWindowNibName: @"QuizWindow" quizItems: nil ];
}

- ( id ) initWithWindowNibName: ( NSString * ) windowNibName quizItems: ( NSArray * ) items
{
	self = [ super initWithWindowNibName: windowNibName ];
	
	if ( self )
	{
		quizItems = [ [ NSArray alloc ] init ];
		[ self setQuizItems: items ];
	
		rightColor = [ [ NSColor colorWithCalibratedRed: 0.0 green: 0.6 blue: 0.0 alpha: 1.0 ] retain ];
		wrongColor = [ [ NSColor colorWithCalibratedRed: 0.6 green: 0.0 blue: 0.0 alpha: 1.0 ] retain ];
		
		attemptedAnswer = NO;
	}
	
	return self;
}

- ( void ) dealloc
{
	[ [ NSNotificationCenter defaultCenter ] removeObserver: self ];
	[ rightColor release ];
	[ wrongColor release ];
	[ quizItems release ];
	[ super dealloc ];
}

- ( void ) showFullScreenQuizWindow
{
	// This is critical for allowing the full-screen quiz window to display its set of controls
	if ( quizView == nil )
		[ NSBundle loadNibNamed: @"QuizWindow.nib" owner: self ];
	
	// Capture the main display or bail out, with minimalist error-reporting
	// This looks like a very crude, brute-force way to set things up
	
	if ( CGDisplayCapture( kCGDirectMainDisplay ) != kCGErrorSuccess )
	{
		NSLog( NSLocalizedString( @"Unable to capture the main display", nil ) );
		return;
	}
	
	// Create the window and set its level to be above all other windows, even the Dock!
	
	ddvQuizWindow = [ [ DDVQuizWindow alloc ]
		initWithContentRect: [ [ NSScreen mainScreen ] frame ]
		styleMask: NSBorderlessWindowMask backing: NSBackingStoreBuffered
		defer: NO screen: [ NSScreen mainScreen ] ];
	
	[ ddvQuizWindow setFrame: [ [ NSScreen mainScreen ] frame ] display: YES ];
	[ ddvQuizWindow setContentView: quizView ];			// Window needs a view!
	[ ddvQuizWindow setLevel:  CGShieldingWindowLevel( ) ];
	[ ddvQuizWindow makeKeyAndOrderFront: self ];
	
	[ self windowDidLoad ];
}

- ( void ) windowDidLoad
{
	int index;
	NSNumber *num;
	NSMutableArray *seq = [ NSMutableArray array ];
	srandom( time( NULL ) );		// Seed the random number generator
	
	if ( ddvQuizWindow == nil )
	{
		[ [ self window ] setContentView: quizView ];
	}
	
	[ promptField setFont: [ NSFont userFontOfSize: 24 ] ];
	[ solutionField setFont: [ NSFont userFontOfSize: 24 ] ];
	
	// Create an ordered sequence of integers from zero to the number of quiz items
	for ( index = 0; index < [ quizItems count ]; index++ )
	{
		num = [ NSNumber numberWithInt: index ];
		[ seq addObject: num ];
	}
	
	// Create an instance of Michael Meer's MMRandomEnumerator to shuffle the the seq array
	// These will be used to pick items at random from the array of vocables
	
	quizEnum = [ [ MMRandomEnumerator initWithArray: seq ] retain ];
	// The enumerator retained here and is released in windowShouldClose:
	
	defaults = [ NSUserDefaults standardUserDefaults ];
	[ self showNextItem ];
}

#pragma mark QUIZ

// Run the quiz; showNextItem is triggered by a timer in solve: after an answer is entered

- ( void ) showNextItem
{
	if ( ranseq = [ quizEnum nextObject ] )		// ranseq is just an NSNumber containing an integer
	{
		if ( [ defaults boolForKey: DDVSwapLangCheckBoxStateKey ] )
			swapLangFlag = random( );
		else
			swapLangFlag = 0;
		
		if ( swapLangFlag % 2 == 0 )
		{
			// Set text field labels to the current pair of labels from the document
			[ firstLangLabel setStringValue: [ [ [ self document ] lexiconModel ] firstLangLabel ] ];
			[ secondLangLabel setStringValue: [ [ [ self document ] lexiconModel ] secondLangLabel ] ];
		}
		else
		{
			// Set text field labels to the current pair of labels from the document
			[ secondLangLabel setStringValue: [ [ [ self document ] lexiconModel ] firstLangLabel ] ];
			[ firstLangLabel setStringValue: [ [ [ self document ] lexiconModel ] secondLangLabel ] ];
		}
		
		// This DDVLexItem instance will be modified and released in the solve: method
		item = [ [ quizItems objectAtIndex: [ ranseq intValue ] ] retain ];
		
		if ( [ [ item score ] floatValue ] < 0.8 || [ item numGuesses ] < 25 )
		{
			if ( swapLangFlag % 2 == 0 )
			{
				[ promptField setStringValue: [ item firstL ] ];
				[ solutionField setStringValue: @"" ];
			}
			else
			{
				[ promptField setStringValue: [ item secondL ] ];
				[ solutionField setStringValue: @"" ];
			}
		}
		else if ( [ [ item score ] floatValue ] >= 0.8 && random( ) % 1000 > 666 )
		{
			if ( swapLangFlag % 2 == 0 )
			{
				[ promptField setStringValue: [ item firstL ] ];
				[ solutionField setStringValue: @"" ];
			}
			else
			{
				[ promptField setStringValue: [ item secondL ] ];
				[ solutionField setStringValue: @"" ];
			}
		}
		else
		{
			[ item release ];
			[ self showNextItem ];
		}
	}
	else	// We have reached the last item in the quiz
	{
		[ solutionField setEnabled: FALSE ];
		[ skipButton setEnabled: FALSE ];
	}
}

// The solve: action is triggered by a <CR> in the response text field

- ( IBAction ) solve: ( id ) sender
{
	attemptedAnswer = YES;	// This is a flag that will be used to mark doc-changed status
	
	// Even if case is grammatically significant, do a somewhat-forgiving caseInsensitveCompare:
	if ( swapLangFlag % 2 == 0 )
	{
		if ( [ [ item secondL ] caseInsensitiveCompare: [ solutionField stringValue ] ] == NSOrderedSame )
		{
			[ item setNumGuesses: [ item numGuesses ] + 1 ];
			[ item setNumCorrect: [ item numCorrect ] + 1 ];
			isCorrect = YES;
		}
		else
		{
			[ item setNumGuesses: [ item numGuesses ] + 1 ];
			isCorrect = NO;
		}
	}
	else
	{
		if ( [ [ item firstL ] caseInsensitiveCompare: [ solutionField stringValue ] ] == NSOrderedSame )
		{
			[ item setNumGuesses: [ item numGuesses ] + 1 ];
			[ item setNumCorrect: [ item numCorrect ] + 1 ];
			isCorrect = YES;
		}
		else
		{
			[ item setNumGuesses: [ item numGuesses ] + 1 ];
			isCorrect = NO;
		}
	}
	
	// The user will see the correct answer in a color that indicates whether guess is correct
	if ( swapLangFlag % 2 == 0 )
	{
		[ feedback setAttributedStringValue: [ self feedbackString: [ item secondL ] ] ];
		[ NSTimer scheduledTimerWithTimeInterval: 2 target: self
			selector: @selector( showNextItem ) userInfo: nil repeats: NO ];
	}
	else
	{
		[ feedback setAttributedStringValue: [ self feedbackString: [ item firstL ] ] ];
		[ NSTimer scheduledTimerWithTimeInterval: 2 target: self
			selector: @selector( showNextItem ) userInfo: nil repeats: NO ];
	}
	
	[ item release ];	// This DDVLexItem object was retained in the showNextItem method
}

// The skip: and terminate: actions are triggered by buttons on the quiz window interface

- ( IBAction ) skip: ( id ) sender
{
	[ self showNextItem ];
	[ feedback setStringValue: @"" ];
}

- ( IBAction ) terminate: ( id ) sender
{
	if ( ddvQuizWindow != nil )		// We have been running quiz in full-screen mode
	{
		// Release the enumerator and window instances
		[ quizEnum release ];
		[ self windowShouldClose: sender ];
		[ ddvQuizWindow release ];
		
		// Release the display
		if ( CGDisplayRelease( kCGDirectMainDisplay ) != kCGErrorSuccess )
		{
			NSLog( @"Unable to release the display(s)" );
			// If using an error dialog, make sure it appears at shield window level to be visible
		}
	}
	else
	{
		// Release the enumerator object and close the window if not full-screen
		[ quizEnum release ];
		[ self windowShouldClose: sender ];
		[ [ self window ] close ];
	}
}

// Formatting a string for feedback - green for correct, red for incorrect guess
// Feedback is given in "large, friendly letters" made popular by Douglas Adams

- ( NSMutableAttributedString * ) feedbackString: ( NSString * ) string
{
	NSMutableAttributedString *fs;
	
	fs = [ [ NSMutableAttributedString alloc ] initWithString: string ];
	[ fs addAttribute: NSFontAttributeName value: [ NSFont userFontOfSize: 24 ]
		range: NSMakeRange( 0, [ string length ] ) ];
	
	if ( isCorrect )
		[ fs addAttribute: NSForegroundColorAttributeName value: rightColor
			range: NSMakeRange( 0, [ string length ] ) ];
	else
		[ fs addAttribute: NSForegroundColorAttributeName value: wrongColor
			range: NSMakeRange( 0, [ string length ] ) ];
	
	return [ fs autorelease ];
}

- ( BOOL ) windowShouldClose: ( id ) sender
{
	// The flag will mark document as changed if any quiz item's score basis was changed
	[ [ [ [ self document ] windowControllers ] objectAtIndex: 0 ]
		endQuizWithFlag: attemptedAnswer ];
	
	// Put the main document window back in front as key window
	[ [ [ [ [ self document ] windowControllers ] objectAtIndex: 0 ] window ]
		makeKeyAndOrderFront: self ];
	
	return TRUE;
}

#pragma mark ACCESSORS

- ( void ) setQuizItems: ( NSArray * ) items
{
	[ quizItems autorelease ];
	quizItems = [ items retain ];
}

@end
