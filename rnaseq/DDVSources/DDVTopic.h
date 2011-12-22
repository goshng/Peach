//  DDVTopic.h
//  XQuizIt

//  Created by Daniel Stein on Mon Feb 21 2005.
//  Copyright (c) 2005 DelDotVee. All rights reserved.

#import <Cocoa/Cocoa.h>

// Bindings again: The model key to which the topic table is bound is the topic string
// Another pure model class to support the topics table

@interface DDVTopic: NSObject
{
	NSString *topicString;
}

- ( id ) initWithString: ( NSString * ) aString;
- ( NSString *) topicString;

@end
