//  DDVLexItem.h
//  XQuizIt

//  Created by Daniel Stein on Thu Feb 10 2005.
//  Copyright (c) 2005 DelDotVee. All rights reserved.

#import <Cocoa/Cocoa.h>

// This is a pure model class for a SINGLE entry in the vocabulary table (or dictionary)
// This class also keeps track of how a student is scoring on each item in vocabulary quizzes
// A topic field is associated with each item for additional categorization within a table

@interface DDVLexItem: NSObject < NSCoding >
{
	NSString *firstL;
	NSString *secondL;
	NSString *topic;
	NSString *score;
	NSImage *image;
	
	int numGuesses;
	int numCorrect;
}

#pragma mark еее INITIALIZATION еее

- ( id ) initWithFirstLanguage: ( NSString * ) fLWord secondLanguage: ( NSString * ) sLWord
	topic: ( NSString * ) newTopic;

#pragma mark еее ACCESSORS еее

- ( void ) setFirstL: ( NSString * ) aWord;
- ( void ) setSecondL: ( NSString * ) aWord;
- ( void ) setTopic: ( NSString * ) newTopic;
- ( void ) setImage: ( NSImage * ) newImage;
- ( void ) setNumCorrect: ( int ) anInt;
- ( void ) setNumGuesses: ( int ) anInt;

- ( NSString * ) firstL;
- ( NSString * ) secondL;
- ( NSString * ) topic;
- ( NSImage * ) image;		// Image used in Score column to show item performance at a glance
- ( NSString * ) score;		// Computation for assigning a color-coded image to a student's performance
- ( int ) numCorrect;		// Raw score for student performance on quizzes
- ( int ) numGuesses;		// Number of trials in vocabulary quizzes

- ( void ) swapLanguages;

@end
