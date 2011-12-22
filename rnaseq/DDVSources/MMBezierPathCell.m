//  MMBezierPathCell.m
//  XQuizIt

//  Created by Michael Meer on Mon Oct 27 2003.
//  Copyright (c) 2003 Michael Meer. All rights reserved.

#import "MMBezierPathCell.h"

@implementation MMBezierPathCell

- ( id ) init
{
	if ( self = [ super init ] )
	{
		[ self setQuintiles ];
	}
	
	return self;
}

- ( void ) drawInteriorWithFrame: ( NSRect ) cellFrame inView: ( NSView * ) controlView 
{
	NSBezierPath *ellipse;
	int offset = 2;
	NSRect bounds = cellFrame;
	
	bounds.size.width = cellFrame.size.height;
	bounds.origin.x = cellFrame.origin.x + cellFrame.size.width / 3;
	bounds.origin.y += offset;
	bounds.size.height -= 2 * offset;
	bounds.size.width = 2 * bounds.size.height;
	
	if ( [ [ self objectValue ] intValue ] == -1 )
	{
		[ [ NSColor whiteColor ] set ];
	}
	else
	{
		if ( [ [ self objectValue ] floatValue ] >= firstQuintile )
			[ [ NSColor blueColor ] set ];
		else if ( [ [ self objectValue ] floatValue ] >= secondQuintile )
			[ [ NSColor greenColor ] set ];
		else if ( [ [ self objectValue ] floatValue ] >= thirdQuintile )
			[ [ NSColor yellowColor ] set ];
		else if ( [ [ self objectValue ] floatValue ] >= fourthQuintile )
			[ [ NSColor redColor ] set ];
		else if ( [ [ self objectValue ] floatValue ] >= fifthQuintile )
			[ [ NSColor blackColor ] set ];
	}
	
	ellipse = [ NSBezierPath bezierPathWithOvalInRect: bounds ];
	[ ellipse fill ];
	[ [ NSColor blackColor ] set ];
	[ ellipse setLineWidth: 1.5 ];
	[ ellipse stroke ];
}

- ( void ) setQuintiles
{
    firstQuintile = 0.8;
    secondQuintile = 0.6;
    thirdQuintile = 0.4;
	fourthQuintile = 0.2;
    fifthQuintile = 0.0;
 }

@end
