//  MMPrintView.h
//  XQuizIt

//  Created by Michael Meer on Fri Feb 27 2004.
//  Copyright (c) 2004 Michael Meer. All rights reserved.

#import <Cocoa/Cocoa.h>

// Credit where credit is due: This is Michael Meer's print view
// It is reproduced virtually verbatim from the original VocableTrainerX
// Several minor modifications have been added to include the topics column in the view

@interface MMPrintView: NSView
{
	NSArray *lexicon;
	int itemsPerPage;
	float itemHeight;
	float pageHeight;
	float pageWidth;
	float headerHeight;
	float contentHeight;
	int pages;					// Number of pages being printed
	float rectHeight;			// Amount of vertical space for each vocable
	NSString *title;
	NSSize paperSize;
	
	NSMutableDictionary *normalAttributes;
	NSMutableDictionary *boldAttributes;
	NSMutableDictionary *boldRightJustifiedAttributes;
}

- ( id ) initWithLexItems: ( NSArray * ) array andTitle: ( NSString * ) newTitle printInfo: ( NSPrintInfo * ) pi;

@end
