//  MMPasteboardCategory.h
//  XQuizIt

//  Created by Michael Meer on Fri Feb 20 2004.
//  Copyright (c) 2004 Michael Meer. All rights reserved.

#import <Cocoa/Cocoa.h>

#import "DDVWindowController.h"

// Credit where credit is due, of course...
// Michael Meer's Pasteboard code, reproduced essentially untouched from original VocableTrainerX
// The only difference is that here it is implemented as a category on the main window controller

@interface DDVWindowController ( MMPasteboardCategory )

#pragma mark PASTEBOARD

- ( IBAction ) cut: ( id ) sender;
- ( IBAction ) copy: ( id ) sender;
- ( IBAction ) paste: ( id ) sender;
- ( void ) writeStringToPboard: ( NSPasteboard * ) pb;
- ( BOOL ) readStringFromPasteboard: ( NSPasteboard * ) pb;
- ( NSString * ) prepareNSStringForPboard: ( NSArray * ) vocables;
- ( NSString * ) prepareMMPboardType: ( NSArray * ) vocables;
- ( NSMutableArray * ) readNSStringFromPboard: ( NSPasteboard * ) pb;
- ( NSMutableArray * ) readLexItemStringFromPboard: ( NSPasteboard * ) pb;

@end
