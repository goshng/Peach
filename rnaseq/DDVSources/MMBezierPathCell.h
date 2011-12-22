//  MMBezierPathCell.h
//  XQuizIt

//  Created by Michael Meer on Mon Oct 27 2003.
//  Copyright (c) 2003 Michael Meer. All rights reserved.

#import <Cocoa/Cocoa.h>

// Michael Meer developed this utility class from a feature request that I originally made to him

@interface MMBezierPathCell: NSCell
{
	double firstQuintile;
	double secondQuintile;
	double thirdQuintile;
	double fourthQuintile;
	double fifthQuintile;
}

- ( void ) setQuintiles;

@end
