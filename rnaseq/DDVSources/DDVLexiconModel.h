//  DDVLexiconModel.h
//  XQuizIt

//  Created by Daniel Stein on Fri Feb 11 2005.
//  Copyright (c) 2005 DelDotVee. All rights reserved.

#import <Cocoa/Cocoa.h>

@class DDVDocument;

@interface DDVLexiconModel: NSObject <NSCoding>
{
	DDVDocument *document;
	
	NSRect frame;					// This is for saving the doc window position and size
	NSArray *lexiconArray;			// This is the actual table data
	
	NSString *firstLangLabel;		// These are labels for column headers and other interface objects
	NSString *secondLangLabel;
	NSString *topicLabel;
}

#pragma mark INITIALIZATION

- ( id ) initWithDocument: ( DDVDocument * ) vtxDocument;	// designated initializer

#pragma mark ACCESSORS

- ( DDVDocument * ) document;

- ( NSArray * ) lexiconArray;
- ( void ) setLexiconArray: ( NSArray * ) array;

- ( NSRect ) frame;
- ( void ) setFrame: ( NSRect ) theFrame;
- ( BOOL ) hasFrame;

- ( NSString * ) firstLangLabel;
- ( void ) setFirstLangLabel: ( NSString * ) label;
- ( NSString * ) secondLangLabel;
- ( void ) setSecondLangLabel: ( NSString * ) label;
- ( NSString * ) topicLabel;
- ( void ) setTopicLabel: ( NSString * ) label;

@end
