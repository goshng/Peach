//  DDVDocument.h
//  XQuizIt

//  Created by Daniel Stein on Wed Feb 09 2005.
//  Copyright (c) 2005 DelDotVee. All rights reserved.

#import <Cocoa/Cocoa.h>

@class DDVLexiconModel;

@interface DDVDocument: NSDocument
{
	DDVLexiconModel *lexiconModel;		// Subsumes all data model object activity
}

- ( DDVLexiconModel * ) lexiconModel;
- ( void ) setLexiconModel: ( DDVLexiconModel * ) aModel;

@end
