//  DDVTopicController.m
//  XQuizIt

//  Created by Daniel Stein on Mon Feb 21 2005.
//  Copyright (c) 2005 DelDotVee. All rights reserved.

#import "DDVTopicController.h"

@implementation DDVTopicController

// This prevents selection behavior in this table if the user adds a new topic to main table

- ( BOOL ) selectsInsertedObjects
{
	return NO;
}

@end
