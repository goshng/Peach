//  DDVTopicController.h
//  XQuizIt

//  Created by Daniel Stein on Mon Feb 21 2005.
//  Copyright (c) 2005 DelDotVee. All rights reserved.

#import <Cocoa/Cocoa.h>

// We override selectsInsertedObjects so the topic controller cannot filter on new entries

@interface DDVTopicController: NSArrayController
{
}

- ( BOOL ) selectsInsertedObjects;

@end
