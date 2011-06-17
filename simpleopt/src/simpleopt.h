/*! @file SimpleOpt.h

    @version 3.5

    @brief A cross-platform command line library which can parse almost any
    of the standard command line formats in use today. It is designed 
    explicitly to be portable to any platform and has been tested on Windows 
    and Linux. See CSimpleOptTempl for the class definition.

    @section features FEATURES

    -   MIT Licence allows free use in all software (including GPL 
        and commercial)
    -   multi-platform (Windows 95/98/ME/NT/2K/XP, Linux, Unix)
    -   supports all lengths of option names:
        <table width="60%">
            <tr><td width="30%"> - 
                <td>switch character only (e.g. use stdin for input)
            <tr><td> -o          
                <td>short (single character)
            <tr><td> -long       
                <td>long (multiple character, single switch character)
            <tr><td> --longer    
                <td>long (multiple character, multiple switch characters)
        </table>
    -   supports all types of arguments for options:
        <table width="60%">
            <tr><td width="30%"> --option        
                <td>short/long option flag (no argument)
            <tr><td> --option ARG    
                <td>short/long option with separate required argument
            <tr><td> --option=ARG    
                <td>short/long option with combined required argument
            <tr><td> --option[=ARG]  
                <td>short/long option with combined optional argument
            <tr><td> -oARG           
                <td>short option with combined required argument
            <tr><td> -o[ARG]         
                <td>short option with combined optional argument
        </table>
    -   supports options with multiple or variable numbers of arguments:
        <table width="60%">
            <tr><td width="30%"> --multi ARG1 ARG2      
                <td>Multiple arguments
            <tr><td> --multi N ARG-1 ARG-2 ... ARG-N    
                <td>Variable number of arguments
        </table>
    -   supports case-insensitive option matching on short, long and/or 
        word arguments.
    -   supports options which do not use a switch character. i.e. a special 
        word which is construed as an option. 
        e.g. "foo.exe open /directory/file.txt" 
    -   supports clumping of multiple short options (no arguments) in a string 
        e.g. "foo.exe -abcdef file1" <==> "foo.exe -a -b -c -d -e -f file1"
    -   automatic recognition of a single slash as equivalent to a single 
        hyphen on Windows, e.g. "/f FILE" is equivalent to "-f FILE".
    -   file arguments can appear anywhere in the argument list:
        "foo.exe file1.txt -a ARG file2.txt --flag file3.txt file4.txt"
        files will be returned to the application in the same order they were 
        supplied on the command line
    -   short-circuit option matching: "--man" will match "--mandate"
        invalid options can be handled while continuing to parse the command 
        line valid options list can be changed dynamically during command line
        processing, i.e. accept different options depending on an option 
        supplied earlier in the command line.
    -   implemented with only a single C++ header file
    -   optionally use no C runtime or OS functions
    -   char, wchar_t and Windows TCHAR in the same program
    -   complete working examples included
    -   compiles cleanly at warning level 4 (Windows/VC.NET 2003), warning 
        level 3 (Windows/VC6) and -Wall (Linux/gcc)

    @section usage USAGE

    The SimpleOpt class is used by following these steps:

    <ol>
    <li> Include the SimpleOpt.h header file

        <pre>
        \#include "SimpleOpt.h"
        </pre>

    <li> Define an array of valid options for your program.

<pre>
@link CSimpleOptTempl::SOption CSimpleOpt::SOption @endlink g_rgOptions[] = {
    { OPT_FLAG, _T("-a"),     SO_NONE    }, // "-a"
    { OPT_FLAG, _T("-b"),     SO_NONE    }, // "-b"
    { OPT_ARG,  _T("-f"),     SO_REQ_SEP }, // "-f ARG"
    { OPT_HELP, _T("-?"),     SO_NONE    }, // "-?"
    { OPT_HELP, _T("--help"), SO_NONE    }, // "--help"
    SO_END_OF_OPTIONS                       // END
};
</pre>

        Note that all options must start with a hyphen even if the slash will
        be accepted. This is because the slash character is automatically
        converted into a hyphen to test against the list of options. 
        For example, the following line matches both "-?" and "/?" 
        (on Windows).

        <pre>
        { OPT_HELP, _T("-?"),     SO_NONE    }, // "-?"
        </pre>

   <li> Instantiate a CSimpleOpt object supplying argc, argv and the option 
        table

<pre>
@link CSimpleOptTempl CSimpleOpt @endlink args(argc, argv, g_rgOptions);
</pre>

   <li> Process the arguments by calling Next() until it returns false. 
        On each call, first check for an error by calling LastError(), then 
        either handle the error or process the argument.

<pre>
while (args.Next()) {
    if (args.LastError() == SO_SUCCESS) {
        handle option: use OptionId(), OptionText() and OptionArg()
    }
    else {
        handle error: see ESOError enums
    }
}
</pre>

   <li> Process all non-option arguments with File(), Files() and FileCount()

<pre>
ShowFiles(args.FileCount(), args.Files());
</pre>

    </ol>

    @section notes NOTES

    -   In MBCS mode, this library is guaranteed to work correctly only when
        all option names use only ASCII characters.
    -   Note that if case-insensitive matching is being used then the first
        matching option in the argument list will be returned.

    @section licence MIT LICENCE

    The licence text below is the boilerplate "MIT Licence" used from:
    http://www.opensource.org/licenses/mit-license.php

    Copyright (c) 2006-2007, Brodie Thiesfield

    Permission is hereby granted, free of charge, to any person obtaining a
    copy of this software and associated documentation files (the "Software"),
    to deal in the Software without restriction, including without limitation
    the rights to use, copy, modify, merge, publish, distribute, sublicense,
    and/or sell copies of the Software, and to permit persons to whom the
    Software is furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included
    in all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
    OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
    IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
    CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
    TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
    SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

/*! @mainpage

    <table>
        <tr><th>Library     <td>SimpleOpt
        <tr><th>Author      <td>Brodie Thiesfield [code at jellycan dot com]
        <tr><th>Source      <td>http://code.jellycan.com/simpleopt/
    </table>

    @section SimpleOpt SimpleOpt

    A cross-platform library providing a simple method to parse almost any of
    the standard command-line formats in use today.

    See the @link SimpleOpt.h SimpleOpt @endlink documentation for full 
    details.

    @section SimpleGlob SimpleGlob

    A cross-platform file globbing library providing the ability to
    expand wildcards in command-line arguments to a list of all matching 
    files.

    See the @link SimpleGlob.h SimpleGlob @endlink documentation for full 
    details.
*/

#ifndef INCLUDED_SimpleOpt
#define INCLUDED_SimpleOpt

// Default the max arguments to a fixed value. If you want to be able to 
// handle any number of arguments, then predefine this to 0 and it will 
// use an internal dynamically allocated buffer instead.
#ifdef SO_MAX_ARGS
# define SO_STATICBUF   SO_MAX_ARGS
#else
# include <stdlib.h>    // malloc, free
# include <string.h>    // memcpy
# define SO_STATICBUF   50
#endif

//! Error values
typedef enum _ESOError
{
    //! No error
    SO_SUCCESS          =  0,   

    /*! It looks like an option (it starts with a switch character), but 
        it isn't registered in the option table. */
    SO_OPT_INVALID      = -1,   

    /*! Multiple options matched the supplied option text. 
        Only returned when NOT using SO_O_EXACT. */
    SO_OPT_MULTIPLE     = -2,   

    /*! Option doesn't take an argument, but a combined argument was 
        supplied. */
    SO_ARG_INVALID      = -3,   

    /*! SO_REQ_CMB style-argument was supplied to a SO_REQ_SEP option
        Only returned when using SO_O_PEDANTIC. */
    SO_ARG_INVALID_TYPE = -4,   

    //! Required argument was not supplied
    SO_ARG_MISSING      = -5,   

    /*! Option argument looks like another option. 
        Only returned when NOT using SO_O_NOERR. */
    SO_ARG_INVALID_DATA = -6    
} ESOError;

//! Option flags
enum _ESOFlags
{
    /*! Disallow partial matching of option names */
    SO_O_EXACT       = 0x0001, 

    /*! Disallow use of slash as an option marker on Windows. 
        Un*x only ever recognizes a hyphen. */
    SO_O_NOSLASH     = 0x0002, 

    /*! Permit arguments on single letter options with no equals sign. 
        e.g. -oARG or -o[ARG] */
    SO_O_SHORTARG    = 0x0004, 

    /*! Permit single character options to be clumped into a single 
        option string. e.g. "-a -b -c" <==> "-abc" */
    SO_O_CLUMP       = 0x0008, 

    /*! Process the entire argv array for options, including the 
        argv[0] entry. */
    SO_O_USEALL      = 0x0010, 

    /*! Do not generate an error for invalid options. errors for missing 
        arguments will still be generated. invalid options will be 
        treated as files. invalid options in clumps will be silently 
        ignored. */
    SO_O_NOERR       = 0x0020, 

    /*! Validate argument type pedantically. Return an error when a 
        separated argument "-opt arg" is supplied by the user as a 
        combined argument "-opt=arg". By default this is not considered 
        an error. */
    SO_O_PEDANTIC    = 0x0040, 

    /*! Case-insensitive comparisons for short arguments */
    SO_O_ICASE_SHORT = 0x0100, 

    /*! Case-insensitive comparisons for long arguments */
    SO_O_ICASE_LONG  = 0x0200, 

    /*! Case-insensitive comparisons for word arguments 
        i.e. arguments without any hyphens at the start. */
    SO_O_ICASE_WORD  = 0x0400, 

    /*! Case-insensitive comparisons for all arg types */
    SO_O_ICASE       = 0x0700  
};

/*! Types of arguments that options may have. Note that some of the _ESOFlags
    are not compatible with all argument types. SO_O_SHORTARG requires that
    relevant options use either SO_REQ_CMB or SO_OPT. SO_O_CLUMP requires 
    that relevant options use only SO_NONE.
 */
typedef enum _ESOArgType {
    /*! No argument. Just the option flags.
        e.g. -o         --opt */
    SO_NONE,    

    /*! Required separate argument.  
        e.g. -o ARG     --opt ARG */
    SO_REQ_SEP, 

    /*! Required combined argument.  
        e.g. -oARG      -o=ARG      --opt=ARG  */
    SO_REQ_CMB, 

    /*! Optional combined argument.  
        e.g. -o[ARG]    -o[=ARG]    --opt[=ARG] */
    SO_OPT, 

    /*! Multiple separate arguments. The actual number of arguments is
        determined programatically at the time the argument is processed.
        e.g. -o N ARG1 ARG2 ... ARGN    --opt N ARG1 ARG2 ... ARGN */
    SO_MULTI
} ESOArgType;

//! this option definition must be the last entry in the table
#define SO_END_OF_OPTIONS   { -1, NULL, SO_NONE }

#ifdef _DEBUG
# ifdef _MSC_VER
#  include <crtdbg.h>
#  define SO_ASSERT(b)  _ASSERTE(b)
# else
#  include <assert.h>
#  define SO_ASSERT(b)  assert(b)
# endif
#else
# define SO_ASSERT(b)   //!< assertion used to test input data
#endif

// ---------------------------------------------------------------------------
//                              MAIN TEMPLATE CLASS
// ---------------------------------------------------------------------------

typedef struct _SOption {
    /*! ID to return for this flag. Optional but must be >= 0 */
    int nId;        

    /*! arg string to search for, e.g.  "open", "-", "-f", "--file" 
        Note that on Windows the slash option marker will be converted
        to a hyphen so that "-f" will also match "/f". */
    const char * pszArg;

    /*! type of argument accepted by this option */
    ESOArgType nArgType;   
} SOption;

int Next();
int OptionId();
ESOError LastError(); 
const char * OptionText();
char * OptionArg();

void CSimpleOpt (
        int             argc, 
        char *        argv[], 
        const SOption * a_rgOptions, 
        int             a_nFlags
        );

void CSimpleOptFinalize (); 

#endif // INCLUDED_SimpleOpt
