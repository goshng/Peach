//  DDVLexItem.m
//  XQuizIt

//  Created by Daniel Stein on Thu Feb 10 2005.
//  Copyright (c) 2005 DelDotVee. All rights reserved.

#import "DDVLexItem.h"

@implementation DDVLexItem

#pragma mark ¥¥¥ÊINITIALIZATION ¥¥¥

- ( id ) init
{
	// Strings used for default table entry upon adding a new item
	return [ self initWithFirstLanguage: NSLocalizedString( @"A Word", nil )
		secondLanguage: NSLocalizedString( @"Ein Wort", nil )
		topic: NSLocalizedString( @"None", nil ) ];
}

// Designated initializer

- ( id ) initWithFirstLanguage: ( NSString * ) fLWord secondLanguage: ( NSString * ) sLWord
	topic: ( NSString * ) newTopic
{
    if ( self = [ super init ] )
	{
		firstL = [ [ NSString stringWithString: fLWord ] retain ];
		secondL = [ [ NSString stringWithString: sLWord ] retain ];
		topic = [ [ NSString stringWithString: newTopic ] retain ];
		
		[ self setNumGuesses: 0 ];
		[ self setNumCorrect: 0 ];
		
		score = [ self score ];
		image = [ [ NSImage imageNamed: @"snow.png" ] retain ];
    }
		
    return self;
}

- ( void ) dealloc
{
	[ firstL release ];
	[ secondL release ];
	[ topic release ];
	[ image release ];
	
	[ super dealloc ];
}

#pragma mark ¥¥¥ÊARCHIVING ¥¥¥

- ( void ) encodeWithCoder: ( NSCoder * ) coder		// For saving a lexicon entry
{
	[ coder encodeObject: firstL ];
	[ coder encodeObject: secondL ];
	[ coder encodeObject: topic ];        
	[ coder encodeValueOfObjCType: @encode( int ) at: &numGuesses ];        
	[ coder encodeValueOfObjCType: @encode( int ) at: &numCorrect ];        
}

- ( id ) initWithCoder: ( NSCoder * ) coder		// To load the lexicon entry from a file
{
	if ( self = [ super init ] );
	{
		[ self setFirstL: [ coder decodeObject ] ];
		[ self setSecondL: [ coder decodeObject ] ];		
		[ self setTopic: [ coder decodeObject ] ];
		[ coder decodeValueOfObjCType: @encode( int ) at: &numGuesses ];
		[ coder decodeValueOfObjCType: @encode( int ) at: &numCorrect ];
	}
	
	return self;
}

#pragma mark ¥¥¥ÊACCESSORS ¥¥¥

- ( void ) setFirstL: ( NSString * ) aWord
{
    [ firstL autorelease ];
    firstL = [ aWord copy ];			// Use copy rather than retain in case arg is mutable
}

- ( void ) setSecondL: ( NSString * ) aWord
{
    [ secondL autorelease ];
    secondL = [ aWord copy ];			// Use copy rather than retain in case arg is mutable
}

- ( void ) setTopic: ( NSString * ) newTopic
{
    [ topic autorelease ];
    topic = [ newTopic copy ];			// Use copy rather than retain in case arg is mutable
}

- ( void ) setImage: ( NSImage * ) newImage
{
	[ image autorelease ];
	image = [ newImage retain ];
}

- ( void ) setNumGuesses: ( int ) value
{
	numGuesses = value;
}

- ( void ) setNumCorrect: ( int ) value
{
	numCorrect = value;
}

// At the very least, we will provide an actual blank object in case the user forgets to enter something

- ( NSString * ) firstL
{
	if ( firstL == nil )
		return @"";
	else
		return firstL;
}

- ( NSString * ) secondL
{
	if ( secondL == nil )
		return @"";
	else
		return secondL;
}

- ( NSString * ) topic
{
	if ( topic == nil )
		return @"";
	else
		return topic;
}

- ( NSImage * ) image
{
	return image;
}

- ( int ) numCorrect
{
	return numCorrect;
}

- ( int ) numGuesses
{
	return numGuesses;
}

- ( NSString * ) score
{
	float pct;
	
	if ( numGuesses != 0 )
		pct = ( float ) numCorrect / ( float ) numGuesses;
	else
		pct = -1.0;
	
	// Return the score as a composite string with guesses and number correct
	if ( numGuesses == numCorrect )
		score = [ NSString stringWithFormat: @"%3.1f (%d/%d)", pct, numCorrect, numGuesses ];
	else
		score = [ NSString stringWithFormat: @"%3.2f (%d/%d)", pct, numCorrect, numGuesses ];
	
	// Must return the score as a string displaying the percentage correct if the Score column is data
	//score = [ NSString stringWithFormat: @"%4.2f", pct ];

	if ( pct > 0.8 )
		[ self setImage: [ NSImage imageNamed: @"violet.png" ] ];
	else if ( pct > 0.6 )
		[ self setImage: [ NSImage imageNamed: @"green.png" ] ];
	else if ( pct > 0.4 )
		[ self setImage: [ NSImage imageNamed: @"yellow.png" ] ];
	else if ( pct > 0.2 )
		[ self setImage: [ NSImage imageNamed: @"red.png" ] ];
	else if ( pct >= 0.0 )
		[ self setImage: [ NSImage imageNamed: @"black.png" ] ];
	else if ( pct == -1.0 )
		[ self setImage: [ NSImage imageNamed: @"snow.png" ] ];

	return score;
}

- ( void ) swapLanguages
{
	NSString *temp = [ NSString stringWithString: firstL ];
	[ self setFirstL: secondL ];
	[ self setSecondL: temp ];
}

@end
