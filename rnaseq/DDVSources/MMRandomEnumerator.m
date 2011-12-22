//  MMRandomEnumerator.m
//  Shuffle the order of table entries randomly for the quiz

//  Created by Michael Meer on Thu Oct 31 2002.
//  Copyright (c) 2002 Michael Meer. All rights reserved.

#import "MMRandomEnumerator.h"

@implementation MMRandomEnumerator

#pragma mark еее INITIALIZATION еее

- ( id ) init
{
	if ( self = [ super init ] )
	{
		array =  [ [ NSMutableArray alloc ] init ];
		currentIndex = 0;
	}
	return self;
}

+ ( id ) initWithArray: ( NSArray * ) newArray
{
	MMRandomEnumerator *enumerator;
	NSMutableArray *newArrayCopy;
	NSMutableArray *resultArray;
	int arrayCount;
	int index;
	
	arrayCount = [ newArray count ];
	enumerator = [ [ MMRandomEnumerator alloc ] init ];
	newArrayCopy = [ NSMutableArray arrayWithArray: newArray ];
	resultArray = [ NSMutableArray arrayWithCapacity: arrayCount ];
	
	while ( [ resultArray count ] < arrayCount )
	{
		index = random( ) % [ newArrayCopy count ];
		[ resultArray addObject: [ newArrayCopy objectAtIndex: index ] ];
		[ newArrayCopy removeObjectAtIndex: index ];
		// Another way to do this is to use a method of NSMutableArray
		// ( void ) exchangeObjectAtIndex: ( unsigned ) idx1 withObjectAtIndex: ( unsigned ) idx2
		// Applying this on one or several passes through the array would also do nicely
		// Candidate for optimization - you could easily check this with CHUD Tools' Shark
	}
	
	[ enumerator setArray: resultArray ];
	
	return [ enumerator autorelease ];
}

- ( void ) dealloc
{
	[ array release ];
	[ super dealloc ];
}

- ( id ) nextObject
{
	if ( currentIndex >= [ array count ] )
	{
		return nil;
	}
	else
	{
		return [ array objectAtIndex: currentIndex++ ];
	}
}

- ( int ) count
{
	return [ array count ];
}

#pragma mark еее ACCESSORS еее

- ( void ) setArray: ( NSMutableArray * ) newArray
{
	[ array autorelease ];
	array = [ newArray retain ];
}

- ( NSArray * ) array
{
	return array;
}

- ( int ) currentIndex
{
	return currentIndex;
}

@end
