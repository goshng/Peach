//  DDVLexiconModel.m
//  XQuizIt

//  Created by Daniel Stein on Fri Feb 11 2005.
//  Copyright (c) 2005 DelDotVee. All rights reserved.

#import "DDVLexiconModel.h"
#import "DDVDocument.h"

@implementation DDVLexiconModel

// This object will gather all data to be archived in one place

#pragma mark INITIALIZATION

- ( id ) init
{
	return [ self initWithDocument: nil ];
}

// Designated initializer

- ( id ) initWithDocument: ( DDVDocument * ) currentDocument
{
	if ( self = [ super init ] )
	{
		document = currentDocument;
		lexiconArray = [ [ NSMutableArray alloc ] init ];
		// Strings used to label user interface elements associated with particular table columns
		firstLangLabel = [ [ NSString alloc ] initWithString:
			NSLocalizedString( @"Some Language", nil ) ];
		secondLangLabel = [ [ NSString alloc ] initWithString:
			NSLocalizedString( @"Andere Sprache", nil ) ];
		topicLabel = [ [ NSString alloc ] initWithString:
			NSLocalizedString( @"Topic", nil ) ];
	}
	return self;
}

- ( void ) dealloc
{
	[ lexiconArray release ];
	[ firstLangLabel release ];
	[ secondLangLabel release ];
	[ topicLabel release ];
	[ super dealloc ];
}

#pragma mark ACCESSORS

- ( DDVDocument * ) document
{
	return document;
}

- ( NSArray * ) lexiconArray
{
	return lexiconArray;
}

- ( void ) setLexiconArray: ( NSArray * ) array
{
	[ lexiconArray autorelease ];
	lexiconArray = [ array copy ];
}

- ( NSString * ) firstLangLabel
{
	return firstLangLabel;
}

- ( void ) setFirstLangLabel: ( NSString * ) label
{
	[ firstLangLabel autorelease ];
	firstLangLabel = [ label retain ];
}

- ( NSString * ) secondLangLabel
{
	return secondLangLabel;
}

- ( void ) setSecondLangLabel: ( NSString * ) label
{
	[ secondLangLabel autorelease ];
	secondLangLabel = [ label retain ];
}

- ( NSString * ) topicLabel
{
	return topicLabel;
}

- ( void ) setTopicLabel: ( NSString * ) label
{
	[ topicLabel autorelease ];
	topicLabel = [ label retain ];
}

- ( NSRect ) frame
{
	return frame;
}

- ( void ) setFrame: ( NSRect ) theFrame
{
	frame = theFrame;
}

- ( BOOL ) hasFrame
{
	return ( frame.size.height != 0 && frame.size.width != 0 );
}

#pragma mark STORAGE

- ( void ) encodeWithCoder: ( NSCoder * ) encoder
{
	[ encoder encodeRect: frame ];
	[ encoder encodeObject: lexiconArray ];
	[ encoder encodeObject: firstLangLabel ];
	[ encoder encodeObject: secondLangLabel ];
	[ encoder encodeObject: topicLabel ];
}

- ( id ) initWithCoder: ( NSCoder * ) decoder
{
	if ( self = [ super init ] )
	{
		[ self setFrame: [ decoder decodeRect ] ];
		[ self setLexiconArray: [ decoder decodeObject ] ];
		[ self setFirstLangLabel: [ decoder decodeObject ] ];
		[ self setSecondLangLabel: [ decoder decodeObject ] ];
		[ self setTopicLabel: [ decoder decodeObject ] ];
	}
	
	return self;
}

@end
