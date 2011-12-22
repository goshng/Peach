//  DDVSearchFieldFormatter.m
//  XQuizIt

//  Created by Daniel Stein on Fri Feb 18 2005.
//  Copyright (c) 2005 DelDotVee. All rights reserved.

#import "DDVSearchFieldFormatter.h"

@implementation DDVSearchFieldFormatter

- ( id ) init
{
	if ( self = [ super init ] )
	{
	}
	return self;
}

- ( void ) dealloc
{
	[ super dealloc ];
}

// Override methods required
// This helps the search field behave properly when user inadvertently pastes table records in it

- ( BOOL ) getObjectValue: ( id * ) obj forString: ( NSString * ) string errorDescription: ( NSString ** ) e
{
	*obj = string;
	
	if ( string != NULL )
	{
		return YES;
	}
	else
	{
		if ( e != NULL )
			*e = NSLocalizedString( @"Error in DDVSearchFieldFormatter", nil );
		return NO;
	}
}

- ( NSString * ) stringForObjectValue: ( id ) obj
{
	if ( [ obj isKindOfClass: [ NSString class ] ] )
	{
		NSMutableString *s = [ NSMutableString stringWithString: obj ];
		[ s replaceOccurrencesOfString: @"\r" withString: @" "
			options: NSCaseInsensitiveSearch range: NSMakeRange( 0, [ s length ] ) ];
		[ s replaceOccurrencesOfString: @"\n" withString: @" "
			options: NSCaseInsensitiveSearch range: NSMakeRange( 0, [ s length ] ) ];
		return [ [ s copy ] autorelease ];
	}
	else
	{
		return [ obj description ];
	}
}

@end
