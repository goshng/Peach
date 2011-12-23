//  DDVWindowController.h
//  XQuizIt

//  Created by Daniel Stein on Wed Feb 09 2005.
//  Copyright (c) 2005 DelDotVee. All rights reserved.

#import <Cocoa/Cocoa.h>

@class DDVLexItem;
@class DDVQuizController;
@class DDVLexiconController;

// This is the workhorse of the application, and probably contains more code than it should

@interface DDVWindowController: NSWindowController
{
	IBOutlet NSDrawer *topicsDrawer;
	IBOutlet NSSearchField *toolbarSearchField;
	IBOutlet DDVLexiconController *lexiconController;
	
	IBOutlet NSPanel *inspector;
	IBOutlet NSPanel *settingsPanel;
	
	IBOutlet NSTableColumn *firstLColumn;
	IBOutlet NSTableColumn *secondLColumn;
	IBOutlet NSTableColumn *topicColumn;
	
	IBOutlet NSTextField *firstLangLabel;
	IBOutlet NSTextField *secondLangLabel;
	IBOutlet NSTextField *topicLabel;
	
	IBOutlet NSTextField *firstLSetting;
	IBOutlet NSTextField *secondLSetting;
	IBOutlet NSTextField *topicSetting;
	
	IBOutlet NSPanel *categorizePanel;
	IBOutlet NSTextField *categoryField;
	
	IBOutlet NSTableView *lexiconTable;
	IBOutlet NSTextField *selectionStatus;
	
	IBOutlet NSTableView *topicTable;
	IBOutlet NSArrayController *topicController;
	
	NSMutableArray *lexiconArray;
	DDVQuizController *quizController;
	NSModalSession quizSession;
}

#pragma mark ��� INITIALIZATION ���

- ( void ) setupToolbar;

// This function is published in non-OO form in Apple's documentation for toolbars
- ( float ) toolbarHeightForWindow: ( NSWindow * ) window;

- ( void )  tableViewSelectionDidChange: ( NSNotification * ) note;
- ( void ) controlTextDidEndEditing: ( NSNotification * ) note;
- ( void ) updateTopicsList;
- ( void ) updateSelectionStatus;

#pragma mark ��� CONTROLS ���

- ( IBAction ) toggleTopicsDrawer: ( id ) sender;
- ( IBAction ) selectNone: ( id ) sender;
- ( IBAction ) showInspector: ( id ) sender;
- ( IBAction ) startQuiz: ( id ) sender;
- ( void ) endQuizWithFlag: ( BOOL ) answersAttempted;
- ( IBAction ) insertLexItem: ( id ) sender;
- ( IBAction ) deleteLexItem: ( id ) sender;
- ( IBAction ) resetScores: ( id ) sender;
- ( IBAction ) swapLanguages: ( id ) sender;
- ( IBAction ) categorize: ( id ) sender;

#pragma mark ��� SHEETS ���

- ( IBAction ) openSettingsPanel: ( id ) sender;
- ( IBAction ) dismissSettingsPanel: ( id ) sender;
- ( void ) updateLabelStrings;
- ( void ) openCategorizeSheet;
- ( IBAction ) dismissCategorizeSheet: ( id ) sender;
- ( IBAction ) cancelCategorizeAction: ( id ) sender;

#pragma mark ��� DE RNA-seq ���

- ( IBAction ) openSamplePathPanel: ( id ) sender;
- ( IBAction ) indexGenome: ( id ) sender;
- ( IBAction ) mapSamples: ( id ) sender;
- ( IBAction ) findDEGenes: ( id ) sender;
- ( IBAction ) initializeR: ( id ) sender;

#pragma mark ��� IMPORT / EXPORT ���

- ( IBAction ) exportCSV: ( id ) sender;
- ( NSString * ) getCSVData;
- ( IBAction ) importCSV: ( id ) sender;

#pragma mark ��� ARRAY ACCESS ���

- ( unsigned int ) countOfLexiconArray;
- ( id ) objectInLexiconArrayAtIndex: ( unsigned int ) index;
- ( void ) insertObject: ( id ) anObject inLexiconArrayAtIndex: ( unsigned int ) index;
- ( void ) removeObjectFromLexiconArrayAtIndex: ( unsigned int ) index;
- ( void ) replaceObjectInLexiconArrayAtIndex: ( unsigned int ) index withObject: ( id ) anObject;

#pragma mark ��� KEY-VALUE OBSERVING ���

- ( void ) startObservingLexItem: ( DDVLexItem * ) lex;
- ( void ) stopObservingLexItem: (DDVLexItem * ) lex;

- ( void ) changeKeyPath: ( NSString * ) keyPath ofObject: ( id ) obj toValue: ( id ) newValue;
- ( void ) observeValueForKeyPath: ( NSString * ) kp ofObject: ( id ) obj change: ( NSDictionary * ) change context: ( void * ) con;

#pragma mark ��� ACCESSORS ���

- ( NSPanel * ) inspector;

- ( NSArray * ) lexiconArray;
- ( void ) setLexiconArray: ( NSMutableArray * ) array;
- ( NSSearchField * ) toolbarSearchField;

- ( NSTextField * ) firstLSetting;
- ( NSTextField * ) secondLSetting;
- ( NSTextField * ) topicSetting;

- ( NSTableView * ) lexiconTable;

- ( NSArray * ) lexiconSelectedByTopic;

@end
