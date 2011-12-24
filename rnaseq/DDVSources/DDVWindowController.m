//  DDVWindowController.m
//  XQuizIt

//  Created by Daniel Stein on Wed Feb 09 2005.
//  Copyright (c) 2005 DelDotVee. All rights reserved.
#include <stdlib.h>
#include <Rembedded.h>
#include "../Embedding/embeddedRCall.h"
#include <R_ext/Parse.h>
#include <R.h>
#include <Rinternals.h>

#import "DDVWindowController.h"
#import "DDVLexiconController.h"
#import "DDVQuizController.h"
#import "DDVDocument.h"
#import "DDVLexItem.h"
#import "DDVTopic.h"
#import "DDVLexiconModel.h"
#import "DDVSearchFieldFormatter.h"
#import "MMPrintView.h"
#import "defaults.h"

@implementation DDVWindowController

// This is the workhorse of the application, and probably contains more code than it should

#pragma mark еее INITIALIZATION еее

- ( id ) init
{
	if ( self = [ super initWithWindowNibName: @"DDVDocument" ] )
	{
		lexiconArray = [ [ NSMutableArray alloc ] init ];
		[ self setShouldCloseDocument: YES ];   // Close entire document if main window closes
	}
  
//  // BEGIN: R Initialization
  setenv("R_TEXI2DVICMD", "/opt/local/bin/texi2dvi", 1);
  setenv("R_VERSION", "2.14.0", 1);
  setenv("R_TEXI2DVICMD", "/opt/local/bin/texi2dvi", 1);
  setenv("R_PDFVIEWER", "/usr/bin/open", 1);
  setenv("BIBINPUTS", ".:/Library/Frameworks/R.framework/Resources/share/texmf/bibtex/bib:", 1);
  setenv("BSTINPUTS", ".:/Library/Frameworks/R.framework/Resources/share/texmf/bibtex/bst:", 1);
  setenv("SED", "/usr/bin/sed", 1);
  setenv("R_INCLUDE_DIR", "/Library/Frameworks/R.framework/Resources/include", 1);
  setenv("R_PRINTCMD", "lpr", 1);
  setenv("R_RD4DVI", "ae", 1);
  setenv("R_SYSTEM_ABI", "osx,gcc,gxx,gfortran,?", 1);
  setenv("R_RD4PDF", "times,inconsolata,hyper", 1);
  setenv("R_PAPERSIZE", "letter", 1);
  setenv("R_ZIPCMD", "/usr/bin/zip", 1);
  setenv("PAGER", "/usr/bin/less", 1);
  setenv("R_GZIPCMD", "/opt/local/bin/gzip", 1);
  setenv("R_SHARE_DIR", "/Library/Frameworks/R.framework/Resources/share", 1);
  setenv("R_OSTYPE", "unix", 1);
  setenv("R_BROWSER", "/usr/bin/open", 1);
  setenv("PERL5LIB", "/opt/local/lib/perl5/site_perl/5.12.3", 1);
  setenv("R_CMD", "/Library/Frameworks/R.framework/Resources/bin/Rcmd", 1);
  setenv("TEXINPUTS", ".:/Library/Frameworks/R.framework/Resources/share/texmf/tex/latex:", 1);
  setenv("R_ARCH", "", 1);
  setenv("MAKE", "make", 1);
  setenv("R_PAPERSIZE_USER", "", 1);
  setenv("LD_LIBRARY_PATH", "/usr/lib:/usr/local/lib:/opt/local/lib:/Library/Frameworks/R.framework/Resources/lib", 1);
  setenv("DYLD_LIBRARY_PATH", "/usr/lib:/usr/local/lib:/opt/local/lib:/Library/Frameworks/R.framework/Resources/lib", 1);
  setenv("R_UNZIPCMD", "/usr/bin/unzip", 1);
  setenv("R_BZIPCMD", "/opt/local/bin/bzip2", 1);
  setenv("R_HOME", "/Library/Frameworks/R.framework/Resources", 1);
//  setenv("R_PLATFORM", "x86_64-apple-darwin10.8.0", 1);
  setenv("R_PLATFORM", "i386-apple-darwin10.8.0", 1);
  setenv("R_LIBS_USER", "~/Library/R/2.14/library", 1);
  setenv("R_DOC_DIR", "/Library/Frameworks/R.framework/Resources/doc", 1);
  
  if (!getenv("R_HOME")) {
    NSString* lastInitRError = @"R_HOME is not set. Please set all required environment variables before running this program. You have to exit the program because this program will not work with R packages.";
    NSAlert *alert = [[[NSAlert alloc] init] autorelease];
    [alert setMessageText:lastInitRError];
    [alert runModal];
    return self;
  }
  else
  {
    NSString* lastinitRError = [NSString stringWithUTF8String:getenv("R_HOME")];
    NSAlert *alert = [[[NSAlert alloc] init] autorelease];
    [alert setMessageText:lastinitRError];
    [alert runModal];
  }
  
  // This made it worked! This seems to be essential: "--no-save", 
//  char *localArgs[]={ "R", "--silent", "--no-save", "--no-restore-data", "--gui=none"};  // WORKED VERSION!
  char *localArgs[]={"R", "--no-save"};
  
//  char *localArgs[] = {"R", "--gui=none", "--silent"};
//  int stat=Rf_initialize_R(sizeof(localArgs)/sizeof(localArgs[0]), localArgs);
//  if (stat<0) {
//    NSString* commandBwaIndex = 
//    [NSString stringWithFormat:@"Error # from Rf_initialize_R: %d",stat];
//    NSAlert *alert = [[[NSAlert alloc] init] autorelease];
//    [alert setMessageText:commandBwaIndex];
//    [alert runModal];
//    return self;
//  }
////  Rf_mainloop();
  
//  char *localArgs[] = {"R", "--silent"};
  init_R(sizeof(localArgs)/sizeof(localArgs[0]), localArgs);
  
  {
    NSString* lastinitRError = [NSString stringWithUTF8String:getenv("LD_LIBRARY_PATH")];
    NSAlert *alert = [[[NSAlert alloc] init] autorelease];
    [alert setMessageText:lastinitRError];
    [alert runModal];
  }
  
  // END: R Initialization
  
	return self;
}

- ( void ) dealloc
{
	[ [ NSNotificationCenter defaultCenter ] removeObserver: self ];
	[ toolbarSearchField setDelegate: nil ];
	[ self setLexiconArray: nil ];
	[ super dealloc ];
  
  end_R();
}

- ( void ) awakeFromNib
{
	DDVLexiconModel *model = [ [ self document ] lexiconModel ];
	
	if ( [ model hasFrame ] )	// We want to restore it when user loads a documment
	{
		[ self setShouldCascadeWindows: NO ];
		[ [ self window ] setFrame: [ model frame ] display: YES ];
	}
	
	//[ topicsDrawer openOnEdge: NSMaxXEdge ];	// Could add as a setting on a per-doc basis
	
	[ self setupToolbar ];
	[ [ toolbarSearchField cell ] setWraps: NO ];
	[ [ toolbarSearchField cell ] setScrollable: YES ];
	
	// Make the search field behave nicely if the user inadvertently pastes table rows in it
	DDVSearchFieldFormatter *f = [ [ [ DDVSearchFieldFormatter alloc ] init ] autorelease ];
	[ toolbarSearchField setFormatter: f ];  // The formatter is retained by the search field

	// Manual binding prevents an object leak in Panther; is this fixed in Tiger??
	[ lexiconController bind: @"contentArray" toObject: self withKeyPath: @"lexiconArray" options: nil ];
	
	// Use slightly modified version Michael Meer's Bezier Path cell for graphical score display
	//[ [ lexiconTable tableColumnWithIdentifier: @"score" ]
	//	setDataCell: [ [ [ MMBezierPathCell alloc ] init ] autorelease ] ];
	
	// Use luscious graphical globules provided by the included image files for graphical score display
	[ [ lexiconTable tableColumnWithIdentifier: @"image" ]
		setDataCell: [ [ [ NSImageCell alloc ] init ] autorelease ] ];
	
	// Bypass undo registration when loading data from a file
	[ [ [ self document ] undoManager ] disableUndoRegistration ];
	[ lexiconController addObjects: [ model lexiconArray ] ];
	[ [ [ self document ] undoManager ] enableUndoRegistration ];
	
	// To update the interface control labels when loading data from an existing file
	[ firstLSetting setStringValue: [ [ [ self document ] lexiconModel ] firstLangLabel ] ];
	[ secondLSetting setStringValue: [ [ [ self document ] lexiconModel ] secondLangLabel ] ];
	[ topicSetting setStringValue: [ [ [ self document ] lexiconModel ] topicLabel ] ];
	
	// Pretty-up the interface after everything is loaded
	[ self updateLabelStrings ];
	[ lexiconTable deselectAll: self ];
	[ lexiconTable reloadData ];
	[ self updateTopicsList ];
	
	// Update the selection status bar
	[ self updateSelectionStatus ];

	// Shiny, happy document!
	[ [ self document ] updateChangeCount: NSChangeCleared ];
}

- ( void ) windowWillClose: ( NSNotification * ) note
{
	// Finish up by letting go of manual bindings
	[ lexiconController unbind: @"contentArray" ];
	[ lexiconController setContent: nil ];
}

- ( void ) setupToolbar
{
	NSToolbar *toolbar = [ [ NSToolbar alloc ] initWithIdentifier: @"mainToolbar" ];
	[ toolbar setDelegate: self ];
	[ toolbar setAllowsUserCustomization: YES ];	
	[ toolbar setAutosavesConfiguration: YES ];
	[ [ self window ] setToolbar: [ toolbar autorelease ] ];    
}

- ( NSString * ) windowTitleForDocumentDisplayName: ( NSString * ) displayName
{
	[ inspector setTitle: [ NSLocalizedString( @"dernaseq Selection: ", nil )
		stringByAppendingString: displayName ] ];
	return [ @"dernaseq: " stringByAppendingString: displayName ];
}

// This function is published in non-OO form in Apple's documentation for toolbars
// We need this for properly saving the window dimensions to restore later

- ( float ) toolbarHeightForWindow: ( NSWindow * ) window
{
	NSToolbar *toolbar;
	float toolbarHeight = 0.0;
	NSRect windowFrame;
	
	toolbar = [ window toolbar ];
	
	if ( toolbar && [ toolbar isVisible ] )
	{
		windowFrame = [ NSWindow contentRectForFrameRect: [ window frame ] styleMask: [ window styleMask ] ];
		toolbarHeight = NSHeight( windowFrame ) - NSHeight( [ [ window contentView ] frame ] );
	}
	
	return toolbarHeight;
}

// It is more than a little likely that editing a cell in the topic column will affect the topics list

- ( void ) controlTextDidEndEditing: ( NSNotification * ) note
{
	[ self updateTopicsList ];
}

// Update the status bar when the selection changes in the main table view
// A selection change results from any activity in the topics table

- ( void )  tableViewSelectionDidChange: ( NSNotification * ) note
{
	if ( [ note object ] == lexiconTable )
	{
		[ self updateSelectionStatus ];
	}
	
	if ( [ note object ] == topicTable )
	{
		// Somewhat restrictive: We cancel search field filtering if topic filtering is requested
		[ toolbarSearchField setStringValue: @"" ];
		[ lexiconController setSearchString: @"" ];
		[ lexiconController rearrangeObjects ];
		[ lexiconController setSelectedObjects: [ self lexiconSelectedByTopic ] ];
	}
}

// This creates a selection in the main table based on activity in the topics table

- ( NSArray * ) lexiconSelectedByTopic
{
	DDVLexItem *item;
	NSEnumerator *te;
	NSEnumerator *le;
	NSArray *lexicon;
	NSArray *topics;
	DDVTopic *topic;
	NSMutableArray *selection;
	
	topics = [ topicController selectedObjects ];
	lexicon = [ lexiconController arrangedObjects ];
	selection = [ NSMutableArray arrayWithCapacity: [ lexicon count ] ];
	le = [ lexicon objectEnumerator ];
	
	// Go through the main table array and find all items whose topic is selected
	while ( item = [ le nextObject ] )
	{
		// Enumerate the selected items in the topics table
		// Provide a way for explicitly selecting any items with zero-length empty topic strings
		te = [ topics objectEnumerator ];
		while ( topic = [ te nextObject ] )
		{
			if ( [ [ item topic ] isEqualToString: [ topic topicString ] ] )
			{
				[ selection addObject: item ];
			}
			if ( [ [ item topic ] isEqualToString: @"" ] && 
				[ [ topic topicString ] isEqualToString: NSLocalizedString( @"<Empty>", nil ) ] )
			{
				[ selection addObject: item ];
			}
		}
	}
	
	// This result is sent back to tableViewSelectionDidChange: as an immutable array
	return [ [ selection copy ] autorelease ];
}

// Maintain a list of all unique strings in topics column of main table to go in topics table
// Include a special provision to represent empty topic strings when inserting an item

- ( void ) updateSelectionStatus
{
	int a = [ lexiconTable numberOfRows ];
	int s = [ lexiconTable numberOfSelectedRows ];
	
	[ selectionStatus setStringValue:
		[ NSString stringWithFormat: NSLocalizedString( @"%d of %d selected", nil ), s, a ] ];
}

- ( void ) updateTopicsList
{
	DDVLexItem *item;
	NSString *aTopic;
	NSEnumerator *e;
	NSMutableArray *uniqueTopics = [ NSMutableArray array ];
	
	[ topicController removeObjects: [ topicController arrangedObjects ] ];
	
	// Interestingly, the following fails to update the topic display properly when undoing and redoing
	//e = [ [ lexiconController arrangedObjects ] objectEnumerator ];
	
	// So we use this array and it works fine
	e = [ lexiconArray objectEnumerator ];
	
	while ( item = [ e nextObject ] )
	{
		if ( ![ uniqueTopics containsObject: [ item topic ] ] && !( [ [ item topic ] isEqualToString: @"" ] ) )
			[ uniqueTopics addObject: [ item topic ] ];
		if ( ![ uniqueTopics containsObject: NSLocalizedString( @"<Empty>", nil ) ]
			&& ( [ [ item topic ] isEqualToString: @"" ] ) )
			[ uniqueTopics addObject: NSLocalizedString( @"<Empty>", nil ) ];
	}
	
	e = [ uniqueTopics objectEnumerator ];
	
	while ( aTopic = [ e nextObject ] )
		[ topicController addObject: [ [ [ DDVTopic alloc ] initWithString: aTopic ] autorelease ] ];
	
	[ topicTable reloadData ];
	[ topicTable deselectAll: self ];
}

#pragma mark еее╩CONTROLS еее

- ( IBAction ) toggleTopicsDrawer: ( id ) sender
{
	[ topicsDrawer toggle: self ];
}

- ( IBAction ) selectNone: ( id ) sender
{
	[ lexiconTable deselectAll: self ];
	[ topicTable deselectAll: self ];
}

- ( IBAction ) showInspector: ( id ) sender
{
	[ inspector makeKeyAndOrderFront: self ];
}

- ( IBAction ) startQuiz: ( id ) sender
{
	NSUserDefaults *defaults = [ NSUserDefaults standardUserDefaults ];
	
	NSArray *arr = [ lexiconController arrangedObjects ];
	NSArray *sel = [ lexiconController selectedObjects ];
  
  int arrayCount = [arr count];
  NSAutoreleasePool *pool =  [[NSAutoreleasePool alloc] init];
  for (int i = 0; i < arrayCount; i++) {
    NSLog([[arr objectAtIndex:i] firstL]);
    NSLog([[arr objectAtIndex:i] secondL]);
  }
  [pool release];
  return;
	
	// Quiz includes only selected items or else all items in document if no prior selection exists
	if ( [ sel count ] == 0 )
		quizController = [ [ DDVQuizController alloc ] initWithWindowNibName: @"QuizWindow" quizItems: arr  ];
	else
		quizController = [ [ DDVQuizController alloc ] initWithWindowNibName: @"QuizWindow" quizItems: sel  ];
	
	[ [ self document ] addWindowController: [ quizController autorelease ] ];
	
	if ( [ defaults boolForKey: DDVQuizModeCheckBoxStateKey ] )
		[ quizController showFullScreenQuizWindow ];
	else
	{
		if ( [ inspector isVisible ] )
			[ inspector orderOut: self ];
		
		//[ quizController showWindow: self ];
		quizSession = [ NSApp beginModalSessionForWindow: [ quizController window ] ];
	}
	
	[ [ self window ] orderOut: self ];
	[ NSApp runModalSession: quizSession ];
}

- ( void ) endQuizWithFlag: ( BOOL ) answersAttempted
{
	if ( [ NSApp modalWindow ] )
		[ NSApp endModalSession: quizSession ];
	
	if ( answersAttempted )
		[ [ self document ] updateChangeCount: NSChangeDone ];
}

// These two connect a couple of menu items to control table membership

- ( IBAction ) insertLexItem: ( id ) sender
{
	[ lexiconController insert: self ];
	[ self updateSelectionStatus ];
}

- ( IBAction ) deleteLexItem: ( id ) sender
{
	[ lexiconController remove: self ];
	[ self updateSelectionStatus ];
}

- ( IBAction ) resetScores: ( id ) sender
{
	NSArray *lexicon = [ lexiconController arrangedObjects ];
	NSArray *selection = [ lexiconController selectedObjects ];
	NSEnumerator *lexEnumerator;
	DDVLexItem *item;
	
	if ( [ selection count ] == 0 )
		lexEnumerator = [ lexicon objectEnumerator ];
	else
		lexEnumerator = [ selection objectEnumerator ];
	
	// There really should be a warning before doing this
	// All scores are reset if no items have been previously selected

	while ( item = [ lexEnumerator nextObject ] )
	{
		[ item setNumGuesses: 0 ];
		[ item setNumCorrect: 0 ];
	}
	
	[ lexiconTable reloadData ];
	[ [ self document ] updateChangeCount: NSChangeDone ];
}

// Provide a handy utility to correct a commonly-occurring error in data entry

- ( IBAction ) swapLanguages: ( id ) sender
{
	DDVLexItem *item;
	NSArray *lexicon = [ lexiconController selectedObjects ];
	NSEnumerator *lexEnumerator = [ lexicon objectEnumerator ];
	while ( item = [ lexEnumerator nextObject ] )
	{
		[ item swapLanguages ];
	}
	[ [ self document ] updateChangeCount: NSChangeDone ];    
}

// Let the user reclassify multiple selections of vocabulary items all at once

- ( IBAction ) categorize: ( id ) sender
{
	[ self openCategorizeSheet ];
}

// This cleverly-hidden feature allows the user to cancel topic selections with a single mouse click

- ( void ) tableView: ( NSTableView * ) aTableView mouseDownInHeaderOfTableColumn: ( NSTableColumn * ) column
{
	if ( aTableView == topicTable )
		[ topicTable deselectAll: self ];
}

// Don't start a quiz or try to export if nothing is present in the main table yet

- ( BOOL ) validateToolbarItem: ( NSToolbarItem * ) theItem
{
	if ( [ theItem action ] == @selector( startQuiz: ) )
		return [ [ lexiconController arrangedObjects ] count ] > 0;
	if ( [ theItem action ] == @selector( exportCSV: ) )
		return [ [ lexiconController arrangedObjects ] count ] > 0;
	else
		return true;
}

// Mainly to prevent modifications of the data array if search field filtering is in progress
// Also gives reasonable menu item feedback when certain actions do not make sense

- ( BOOL ) validateMenuItem: ( NSMenuItem * ) theItem
{
	if ( [ theItem action ] == @selector( exportCSV: ) )
		return [ [ lexiconController arrangedObjects ] count ] > 0;
	if ( [ theItem action ] == @selector( importCSV: ) )
		return [ [ toolbarSearchField stringValue ] length ] == 0;
	if ( [ theItem action ] == @selector( insertLexItem: ) )
		return [ [ toolbarSearchField stringValue ] length ] == 0;
	if ( [ theItem action ] == @selector( deleteLexItem: ) )
		return ( [ [ lexiconController selectedObjects ] count ] > 0 &&
		[ [ toolbarSearchField stringValue ] length ] == 0 );
	if ( [ theItem action ] == @selector( resetScores: ) )
		return [ [ lexiconController arrangedObjects ] count ] > 0;
	if ( [ theItem action ] == @selector( swapLanguages: ) )
		return [ [ lexiconController selectedObjects ] count ] > 0;
	if ( [ theItem action ] == @selector( categorize: ) )
		return [ [ lexiconController selectedObjects ] count ] > 0;
	if ( [ theItem action ] == @selector( selectNone: ) )
		return [ [ lexiconController selectedObjects ] count ] > 0;
	else
		return TRUE;
}

#pragma mark еее SHEETS еее

// So far the settingsPanel just allows customization of table header cells
// There are labels elsewhere in the interface that are coordinated with the header cell values

- ( IBAction ) openSettingsPanel: ( id ) sender
{
	[ NSApp beginSheet: settingsPanel modalForWindow: [ self window ] modalDelegate: nil
		didEndSelector: nil contextInfo: nil ];
	[ NSApp runModalForWindow: settingsPanel ];
	
	[ NSApp endSheet: settingsPanel ];
	[ settingsPanel orderOut: self ];
}

// Acts when users click the ok button settings.
- ( IBAction ) dismissSettingsPanel: ( id ) sender
{
	[ self updateLabelStrings ];
	[ lexiconTable reloadData ];
	[ NSApp stopModal ];
}

- ( void ) updateLabelStrings
{
//	if ( ![ [ [ firstLColumn  headerCell ] stringValue ] isEqual: [ firstLSetting stringValue ] ] )
//		[ [ self document ] updateChangeCount: NSChangeDone ];
//	if ( ![ [ [ secondLColumn  headerCell ] stringValue ] isEqual: [ secondLSetting stringValue ] ] )
//		[ [ self document ] updateChangeCount: NSChangeDone ];
//	if ( ![ [ [ topicColumn  headerCell ] stringValue ] isEqual: [ topicSetting stringValue ] ] )
//		[ [ self document ] updateChangeCount: NSChangeDone ];

  // Sets the headers of the table view.NSLocalizedString( @"Insert Item", nil )
//	[ [ firstLColumn headerCell ] setStringValue: [ firstLSetting stringValue ] ];
//	[ [ secondLColumn headerCell ] setStringValue: [ secondLSetting stringValue ] ];
//	[ [ topicColumn headerCell ] setStringValue: [ topicSetting stringValue ] ];

//	[ firstLangLabel setStringValue: [ firstLSetting stringValue ] ];
//	[ secondLangLabel setStringValue: [ secondLSetting stringValue ] ];
//	[ topicLabel setStringValue: [ topicSetting stringValue ] ];
  
  [ [ firstLColumn headerCell ] setStringValue:NSLocalizedString(@"Sample ID",nil) ];
	[ [ secondLColumn headerCell ] setStringValue:NSLocalizedString(@"Location",nil) ];
	[ [ topicColumn headerCell ] setStringValue:NSLocalizedString(@"Factor",nil) ];
	
	[ firstLangLabel setStringValue:NSLocalizedString(@"Sample ID",nil) ];
	[ secondLangLabel setStringValue:NSLocalizedString(@"Location",nil) ];
	[ topicLabel setStringValue:NSLocalizedString(@"Factor",nil) ];
}

// Implements mass topic string assignment - this app is now really full of sheet

- ( void ) openCategorizeSheet
{
	[ NSApp beginSheet: categorizePanel modalForWindow: [ self window ] modalDelegate: nil
		didEndSelector: nil contextInfo: nil ];
	[ NSApp runModalForWindow: categorizePanel ];
	
	[ NSApp endSheet: categorizePanel ];
	[ categorizePanel orderOut: self ];
}

- ( IBAction ) dismissCategorizeSheet: ( id ) sender
{
	NSArray *selection = [ lexiconController selectedObjects ];
	NSEnumerator *e = [ selection objectEnumerator ];
	DDVLexItem *item;
	
	while ( item = [ e nextObject ] )
		[ item setTopic: [ categoryField stringValue ] ];
	
	[ [ self document ] updateChangeCount: NSChangeDone ];
	
	[ self updateTopicsList ];
	
	[ NSApp stopModal ];
}

- ( IBAction ) cancelCategorizeAction: ( id ) sender
{
	[ NSApp stopModal ];
}

#pragma mark еее DE RNA-seq еее

// Takes RNA-seq samples from a dialog.
- ( IBAction ) openSamplePathPanel: ( id ) sender
{
  NSOpenPanel* openDlg = [NSOpenPanel openPanel];
  [openDlg setCanChooseFiles:YES];
  [openDlg setCanChooseDirectories:YES];
  [openDlg setPrompt:@"Select"];
//  NSString *filename = @"/Users/goshng"; 
  //  NSString *filename = [pathAsNSString lastPathComponent]; 
  //  [filename stringByDeletingPathExtension];
  if ([openDlg runModalForDirectory:nil file:nil] == NSOKButton )
  {
    NSArray* files = [openDlg filenames];
    for(NSString* filePath in [openDlg filenames])
    {
      NSLog(@"%@",filePath);
      //do something with the file at filePath
      NSAlert *alert = [[[NSAlert alloc] init] autorelease];
      [alert setMessageText:filePath];
      [alert runModal];
    }
  }
}

- ( IBAction ) initializeR: ( id ) sender
{
 
}

- ( IBAction ) indexGenome: ( id ) sender
{
  NSLog(@"Index a reference genome");
  NSString *baseDir = @"/Users/goshng/Documents/Projects/RNASeq-Analysis";
  NSString *bwaDir = 
  [NSString stringWithFormat:@"%@/output/dernaseq", baseDir];
  NSString *bwa = @"bin/bwa";
  NSString *referenceGenome = @"NC_004350.fna";

  NSString *commandBwaIndex = 
  [NSString stringWithFormat:@"cp %@/input/%@ %@/output/project", 
   baseDir, referenceGenome, bwaDir];
  system([commandBwaIndex UTF8String]);
  
  commandBwaIndex = 
  [NSString stringWithFormat:@"%@/%@ index -p %@/output/project/%@-bwa -a is %@/input/%@", 
   baseDir, bwa, bwaDir, referenceGenome, baseDir, referenceGenome];
  system([commandBwaIndex UTF8String]);
}

// Takes all of the RNA-seq samples and maps them on a reference genome.
- ( IBAction ) mapSamples: ( id ) sender
{
  NSLog(@"Run BWA");
  
  NSArray *arr = [ lexiconController arrangedObjects ];
//	NSArray *sel = [ lexiconController selectedObjects ];
  
  int arrayCount = [arr count];
  NSAutoreleasePool *pool =  [[NSAutoreleasePool alloc] init];
  for (int i = 0; i < arrayCount; i++) {
    int fastqNumber = [[[arr objectAtIndex:i] firstL] integerValue];
    NSString* gzipfile = [[arr objectAtIndex:i] secondL];
    NSLog([NSString stringWithFormat:@"%03d %@", fastqNumber, gzipfile]);
    
    // Some global preferences.
    NSString *baseDir = @"/Users/goshng/Documents/Projects/RNASeq-Analysis";
    NSString *bwaDir = 
    [NSString stringWithFormat:@"%@/output/dernaseq", baseDir];
    NSString *bwa = @"bin/bwa";
    NSString *samtools = @"bin/samtools";
    NSString *referenceGenome = @"NC_004350.fna";
    
//    NSString *gzipfile = @"/Volumes/Elements/Documents/Projects/rnaseq/data/smutans/FASTQ01.subsample-100.gz";
    
    // bwa aln
    NSString *commandBwaIndex = 
    [NSString stringWithFormat:@"%@/%@ aln -I -t 2 %@/output/project/%@-bwa %@ > %@/output/project/FASTQ%d.sai", 
     baseDir, bwa, bwaDir, referenceGenome, gzipfile, bwaDir, fastqNumber];
    system([commandBwaIndex UTF8String]);
    
    // bwa samse
    commandBwaIndex = 
    [NSString stringWithFormat:@"%@/%@ samse -n 1 -f %@/output/project/FASTQ%d.sam %@/output/project/%@-bwa %@/output/project/FASTQ%d.sai %@", 
     baseDir, bwa, bwaDir, fastqNumber, bwaDir, referenceGenome, bwaDir, fastqNumber, gzipfile];
    system([commandBwaIndex UTF8String]);
    
    // samtools view -bS
    commandBwaIndex = 
    [NSString stringWithFormat:@"%@/%@ view -bS -o %@/output/project/FASTQ%d.bam %@/output/project/FASTQ%d.sam", 
     baseDir, samtools, bwaDir, fastqNumber, bwaDir, fastqNumber];
    system([commandBwaIndex UTF8String]);
    
    // samtools sort
    commandBwaIndex = 
    [NSString stringWithFormat:@"%@/%@ sort %@/output/project/FASTQ%d.bam %@/output/project/FASTQ%d.sorted", 
     baseDir, samtools, bwaDir, fastqNumber, bwaDir, fastqNumber];
    system([commandBwaIndex UTF8String]);
    
    // samtools mpileup -q 15
    int readDepth = 300000;
    commandBwaIndex = 
    [NSString stringWithFormat:@"%@/%@ mpileup -q 15 -d %d -f %@/output/project/%@ %@/output/project/FASTQ%d.sorted.bam > %@/output/project/FASTQ%d.pileup", 
     baseDir, samtools, readDepth, bwaDir, referenceGenome, bwaDir, fastqNumber, bwaDir, fastqNumber];
    system([commandBwaIndex UTF8String]);
    
    // perl pl/samtools-pileup.pl
    commandBwaIndex = 
    [NSString stringWithFormat:@"perl -I %@ %@/pl/samtools-pileup.pl wiggle -refgenome %@/output/project/%@ -in %@/output/project/FASTQ%d.pileup -out %@/output/project/FASTQ%d.wig", 
     baseDir, baseDir, bwaDir, referenceGenome, 
     bwaDir, fastqNumber, bwaDir, fastqNumber];
    system([commandBwaIndex UTF8String]);
    
    // samtools view and pl/bwa-summary.pl
    commandBwaIndex = 
    [NSString stringWithFormat:@"%@/%@ view %@/output/project/FASTQ%d.sorted.bam | perl -I %@ %@/pl/bwa-summary.pl pos > %@/output/project/FASTQ%d-sum.pos", 
     baseDir, samtools, bwaDir, fastqNumber, baseDir, baseDir, bwaDir, fastqNumber];
    system([commandBwaIndex UTF8String]);
    
    // samtools view and pl/bwa-summary.pl
    NSString *gfffile = @"/Volumes/Elements/Documents/Projects/mauve/bacteria/Streptococcus_mutans_UA159_uid57947/NC_004350.gff";
    commandBwaIndex = 
    [NSString stringWithFormat:@"%@/%@ view %@/output/project/FASTQ%d.sorted.bam | perl -I %@ %@/pl/bwa-summary.pl rrna -gff %@ > %@/output/project/FASTQ%d-sum.rrna", 
     baseDir, samtools, bwaDir, fastqNumber, baseDir, baseDir, gfffile, bwaDir, fastqNumber];
    system([commandBwaIndex UTF8String]);
    
    // perl pl/feature-genome.pl
    NSString *pttfile = @"/Volumes/Elements/Documents/Projects/mauve/bacteria/Streptococcus_mutans_UA159_uid57947/NC_004350.ptt";
    NSString *featurefile = @"feature-genome.out-geneonly";
    commandBwaIndex = 
    [NSString stringWithFormat:@"perl -I %@ %@/pl/feature-genome.pl ptt2 -geneonly -in %@ -out %@/output/project/%@", 
     baseDir, baseDir, pttfile, bwaDir, featurefile];
    system([commandBwaIndex UTF8String]);
    
    commandBwaIndex = 
    [NSString stringWithFormat:@"perl -I %@ %@/pl/de-count.pl join -first -singlegenome -shortread %@/output/project/FASTQ%d-sum.pos -genepos %@/output/project/%@ -o %@/output/project/FASTQ%d.de", 
     baseDir, baseDir, bwaDir, fastqNumber, bwaDir, featurefile, bwaDir, fastqNumber];
    system([commandBwaIndex UTF8String]);
    
    NSLog(commandBwaIndex);
  }
  [pool release];
  return;
}



- ( IBAction ) findDEGenes: ( id ) sender
{
  NSLog(@"Run edgeR or DESeq");
  SEXP e;
  int errorOccurred;
  NSString* testRFile = 
  [NSString stringWithFormat:@"/Users/goshng/error.R"];
  
  
  PROTECT(e = lang2(install("source"), mkString([testRFile UTF8String])));
  R_tryEval(e, R_GlobalEnv, &errorOccurred);
  UNPROTECT(1);
  
  return;

  // Some global preferences.
  NSString *baseDir = @"/Users/goshng/Documents/Projects/RNASeq-Analysis";
  NSString *bwaDir = 
  [NSString stringWithFormat:@"%@/output/dernaseq", baseDir];
  NSString *bwa = @"bin/bwa";
  NSString *samtools = @"bin/samtools";
  NSString *referenceGenome = @"NC_004350.fna";
  
  //    NSString *gzipfile = @"/Volumes/Elements/Documents/Projects/rnaseq/data/smutans/FASTQ01.subsample-100.gz";
  
  // bwa aln
  NSString *commandBwaIndex = 
    [NSString stringWithFormat:@"cut -f4 %@/output/project/feature-genome.out-geneonly > %@/output/project/1.de", 
     bwaDir, bwaDir];
  system([commandBwaIndex UTF8String]);

  commandBwaIndex = 
    [NSString stringWithFormat:@"paste %@/output/project/*.de > %@/output/project/count.txt", 
     bwaDir, bwaDir];
  system([commandBwaIndex UTF8String]);

  commandBwaIndex = 
    [NSString stringWithFormat:@"%@/output/project/count.txt.index", 
     bwaDir];
  [[NSFileManager defaultManager] createFileAtPath:commandBwaIndex contents:nil attributes:nil];  
  
  NSFileHandle *fileHandle = [NSFileHandle fileHandleForWritingAtPath:commandBwaIndex];
  NSString* space = @" ";
  
  NSArray *arr = [ lexiconController arrangedObjects ];
  //	NSArray *sel = [ lexiconController selectedObjects ];
  int arrayCount = [arr count];
  NSAutoreleasePool *pool =  [[NSAutoreleasePool alloc] init];
  for (int i = 0; i < arrayCount; i++) {
//    int fastqNumber = [[[arr objectAtIndex:i] firstL] integerValue];
    NSString* factor = [[arr objectAtIndex:i] topic];
//    NSLog([NSString stringWithFormat:@"%03d %@", fastqNumber, factor]);
    [fileHandle writeData:[factor dataUsingEncoding:NSUTF8StringEncoding]];
    [fileHandle writeData:[space dataUsingEncoding:NSUTF8StringEncoding]];
  }
  space = @"\n";
  [fileHandle writeData:[space dataUsingEncoding:NSUTF8StringEncoding]];
  [fileHandle closeFile];
  [pool release];  
  
  
  
  /*
   Evaluates the two expressions:
   source("error.R")
   and then calls foo()  twice
   where foo is defined in the file error.R
   */
  commandBwaIndex = 
  [NSString stringWithFormat:@"%@/deseq-cornell.R", 
   bwaDir];
  
  
  PROTECT(e = lang2(install("source"), mkString([commandBwaIndex UTF8String])));
  R_tryEval(e, R_GlobalEnv, &errorOccurred);
  UNPROTECT(1);
  
//  PROTECT(e = lang1(install("foo")));
//  R_tryEval(e, R_GlobalEnv, &errorOccurred);
//  UNPROTECT(1);
  
  // Prepares an R script to perform DE.
  // Parses the resulting output file.  
}


#pragma mark еее IMPORT / EXPORT еее

// Typical spreadsheets can handle comma-separated values (CSV)
// Credit where credit is due: This is Michael Meer's work from the original VocableTrainerX

- ( IBAction ) exportCSV: ( id ) sender
{
	NSSavePanel *panel = [ NSSavePanel savePanel ];
	
	[ panel setRequiredFileType: @"csv" ];
	
	if ( [ panel runModal ] == NSOKButton )
		[ [ self getCSVData ] writeToFile: [ panel filename ] atomically: YES ];
}

- ( NSString * ) getCSVData
{
	NSMutableString *csvString = [ NSMutableString string ];
	NSEnumerator *lexE = [ [ lexiconController arrangedObjects ] objectEnumerator ];
	
	DDVLexItem *lex;
	
	while ( lex = [ lexE nextObject ] )
	{
		[ csvString appendString: [ NSString stringWithFormat: @"\"%@\",\"%@\",\"%@\"\n",
			[ lex firstL ], [ lex secondL ], [ lex topic ] ] ];
	}
	return [ NSString stringWithString: csvString ];
}

- ( IBAction ) importCSV: ( id ) sender
{
	NSMutableString *line;
	NSString *lineBuffer;
	NSString *rawString;
	NSArray *lines;
	NSArray *components;
	NSEnumerator *lineEnum;
	DDVLexItem *lex;
	
	NSOpenPanel *panel = [ NSOpenPanel openPanel ];
	
	if ( [ panel runModalForTypes: nil ] == NSOKButton )
		rawString = [ NSString stringWithContentsOfFile: [ [ panel filenames ] objectAtIndex: 0 ] ];
	else
		return;
		
	lines = [ rawString componentsSeparatedByString: @"\n" ];
	lineEnum = [ lines objectEnumerator ];
	
	while ( lineBuffer = [ lineEnum nextObject ] )
	{
		line = [ NSMutableString stringWithString: lineBuffer ];
		[ line replaceOccurrencesOfString: @"\"" withString: @"" options: NSLiteralSearch
			range: NSMakeRange( 0, [ line length ] ) ];
		components = [ line componentsSeparatedByString: @"," ];
		
		if ( [ components count ] < 3 )
			continue;
		else
		{
			lex = [ [ DDVLexItem alloc ] initWithFirstLanguage: [ components objectAtIndex: 0 ]
				secondLanguage: [ components objectAtIndex: 1 ]
				topic: [ components objectAtIndex: 2 ] ];
			[ lexiconController addObject: [ lex autorelease ] ];
		}
	}
	[ [ self document ] updateChangeCount: NSChangeDone ];
	
	[ self updateTopicsList ];

	[ lexiconTable deselectAll: self ];
}

#pragma mark еее ARRAY ACCESS еее

// Indexed accessors for the main data array; recommended by tutorials that introduce Bindings
// In this case, it is part of Key Value Observing (KVO) to support the use of Undo Manager

- ( unsigned int ) countOfLexiconArray
{
	return [ lexiconArray count ];
}

- ( id ) objectInLexiconArrayAtIndex: ( unsigned int ) index
{
   return [ lexiconArray objectAtIndex: index ];
}

- ( void ) insertObject: ( id ) anObject inLexiconArrayAtIndex: ( unsigned int ) index
{
	// Add the inverse of this operation to the undo stack
	NSUndoManager *undo = [ [ self document ] undoManager ];
	[ [ undo prepareWithInvocationTarget: self ] removeObjectFromLexiconArrayAtIndex: index ];
	if ( ![ undo isUndoing ] )
		[ undo setActionName: NSLocalizedString( @"Insert Item", nil ) ];
	// Add the item to the array
	[ self startObservingLexItem: anObject ];
	[ lexiconArray insertObject: anObject atIndex: index ];
	[ self updateTopicsList ];
	[ lexiconTable deselectAll: self ];
}

- ( void ) removeObjectFromLexiconArrayAtIndex: ( unsigned int ) index
{
	DDVLexItem *item = [ lexiconArray objectAtIndex: index ];
	NSUndoManager *undo = [ [ self document ] undoManager ];
	// Add the inverse of this operation to the undo stack
	[ [ undo prepareWithInvocationTarget: self ] insertObject: item inLexiconArrayAtIndex: index ];
	if ( ![ undo isUndoing ] )
		[ undo setActionName: NSLocalizedString( @"Delete Item", nil ) ];
	// Remove the item from the array
	[ self stopObservingLexItem: item ];
	[ lexiconArray removeObjectAtIndex: index ];
	[ self updateTopicsList ];
	[ lexiconTable deselectAll: self ];
}

// Not currently used, but possible utility when pasting into previous table selection

- ( void ) replaceObjectInLexiconArrayAtIndex: ( unsigned int ) index withObject: ( id ) item
{
    [ lexiconArray replaceObjectAtIndex: index withObject: item ];
	[ self updateTopicsList ];
}

#pragma mark еее╩KEY-VALUE OBSERVING еее

// Credit where credit is due: Aaron Hillegass presents the basics for KVO with the Undo Manager
// The original example appears in his book "Cocoa Programming for Max OS X" (Second Edition)
// The publisher is Addison-Wesley (2004)

- ( void ) startObservingLexItem: ( DDVLexItem * ) item
{
	[ item addObserver: self forKeyPath: @"firstL" options: NSKeyValueObservingOptionOld context: NULL ];
	[ item addObserver: self forKeyPath: @"secondL" options: NSKeyValueObservingOptionOld context: NULL ];
	[ item addObserver: self forKeyPath: @"topic" options: NSKeyValueObservingOptionOld context: NULL ];
}

- ( void ) stopObservingLexItem: ( DDVLexItem * ) item
{
	[ item removeObserver: self forKeyPath: @"firstL" ];
	[ item removeObserver: self forKeyPath: @"secondL" ];
	[ item removeObserver: self forKeyPath: @"topic" ];
}

- ( void ) changeKeyPath: ( NSString * ) keyPath ofObject: ( id ) obj toValue: ( id ) newValue
{
	// This method is its own inverse operation for maintaining undo stack with edits
	// The message setValue: forKeyPath: invokes KVO, and takes care of undoing
	[ obj setValue: newValue forKeyPath: keyPath ];
	[ self updateTopicsList ];
}

- ( void ) observeValueForKeyPath: ( NSString * ) kp ofObject: ( id ) obj change: ( NSDictionary * ) change context: ( void * ) con
{
	NSUndoManager *undo = [ [ self document ] undoManager ];
	id oldValue = [ change objectForKey: NSKeyValueChangeOldKey ];
	[ [ undo prepareWithInvocationTarget: self ] changeKeyPath: kp ofObject: obj toValue: oldValue ];
	[ undo setActionName: NSLocalizedString( @"Edit", nil ) ];
}

#pragma mark еее ACCESSORS еее

- ( NSPanel * ) inspector
{
	return inspector;
}

- ( NSTableView * ) lexiconTable
{
	return lexiconTable;
}

- ( NSArray * ) lexiconArray
{
	return [ [  lexiconArray copy ] autorelease ];	// Returning an immutable array by default
}

// This accessor is crucial to ensure proper KVO for array elements

- ( void ) setLexiconArray: ( NSMutableArray * ) array
{
	NSEnumerator *e;
	DDVLexItem *item;
	
	if ( array == lexiconArray )
		return;
	
	// For editing individual cells
	e = [ lexiconArray objectEnumerator ];
	while ( item = [ e nextObject ] )
		[ self stopObservingLexItem: item ];
	
	[ lexiconArray autorelease ];
	lexiconArray = [ array retain ];
	e = [ lexiconArray objectEnumerator ];
	while ( item = [ e nextObject ] )
		[ self startObservingLexItem: item ];
}

- ( NSSearchField * ) toolbarSearchField
{
	return toolbarSearchField;
}

- ( NSTextField * ) firstLSetting
{
	return firstLSetting;
}

- ( NSTextField * ) secondLSetting
{
	return secondLSetting;
}

- ( NSTextField * ) topicSetting
{
	return topicSetting;
}

// Credit where credit is due: This uses Michael Meer's original printing view, with minor mods
// See the actual MMPrintView class for the details of laying out the print view

#pragma mark еее PRINTING еее

- ( void ) printShowingPrintPanel: ( BOOL ) showPanels
{
	NSPrintOperation *op;
	
	// Obtain a custom view that will be printed
	NSMutableString *title = [ NSMutableString stringWithString: [ [ self document ] displayName ] ];
	
	MMPrintView *printView = [ [ MMPrintView alloc ] initWithLexItems: [ lexiconController arrangedObjects ]
		andTitle: title printInfo: [ [ self document ] printInfo ] ];
	
	// Construct the print operation and setup Print panel
	op = [ NSPrintOperation printOperationWithView: printView printInfo: [ [ self document ] printInfo ] ];
	[ op setShowPanels: showPanels ];
	
	// Run operation, which shows the Print panel if showPanels was YES
	
	//[ myDocument runModalPrintOperation: op delegate: nil didRunSelector: NULL contextInfo: NULL ];
	[ op runOperation ];
	[ printView release ];
}

@end
