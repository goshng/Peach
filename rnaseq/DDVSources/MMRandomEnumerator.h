//  MMRandomEnumerator.h
//  XQuizIt

//  Created by Michael Meer on Thu Oct 31 2002.
//  Copyright (c) 2002 Michael Meer. All rights reserved.

#import <Cocoa/Cocoa.h>

// Credit where credit is due: This is Michael Meer's original RandomEnumerator class
// This shuffles the order of table entries randomly when presenting a quiz

@interface MMRandomEnumerator: NSObject
{
	NSMutableArray *array;
	int currentIndex;
}

#pragma mark еее INITIALIZATION еее

+ ( id ) initWithArray: ( NSArray * ) newArray;
- ( id ) nextObject;
- ( int ) count;

#pragma mark еее ACCESSORS еее

- ( void ) setArray: ( NSMutableArray * ) newArray;
- ( NSArray * ) array;
- ( int ) currentIndex;

@end
