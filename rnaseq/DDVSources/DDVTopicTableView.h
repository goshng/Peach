//  DDVTopicTableView.h
//  XQuizIt

//  Created by Daniel Stein on Tue Feb 22 2005.
//  Copyright (c) 2005 DelDotVee. All rights reserved.

#import <Cocoa/Cocoa.h>

// We override acceptsFirstResponder so the topic table cannot become a key view
// This lets the main table highlight its topic-based selections in the proper manner

@interface DDVTopicTableView: NSTableView
{
}

- ( BOOL ) acceptsFirstResponder;

@end
