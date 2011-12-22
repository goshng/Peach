//  DDVTopic.m
//  XQuizIt

//  Created by Daniel Stein on Mon Feb 21 2005.
//  Copyright (c) 2005 DelDotVee. All rights reserved.

#import "DDVTopic.h"

@implementation DDVTopic

// This class supports the Topics table located in the drawer attached to the document window

- ( id ) init
{
	return [ self initWithString: [ NSString stringWithString: NSLocalizedString( @"None", nil ) ] ];
}

- ( id ) initWithString: ( NSString * ) aString
{
	if ( self = [ super init ] )
	{
		topicString = [ aString retain ];
	}
	return self;
}

- ( void ) dealloc
{
	[ topicString autorelease ];
	[ super dealloc ];
}

- ( NSString * ) topicString
{
	// Maybe we could also return something interesting if topic is blank in main table
	return topicString;
}

@end
