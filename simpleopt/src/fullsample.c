// File:    fullSample.cpp
// Library: SimpleOpt
// Author:  Brodie Thiesfield <code@jellycan.com>
// Source:  http://code.jellycan.com/simpleopt/
//
// MIT LICENCE
// ===========
// The licence text below is the boilerplate "MIT Licence" used from:
// http://www.opensource.org/licenses/mit-license.php
//
// Copyright (c) 2006-2007, Brodie Thiesfield
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
// IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
// CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
// TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
// SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#if defined(_MSC_VER)
# include <windows.h>
# include <tchar.h>
#else
# define TCHAR		char
# define _T(x)		x
# define _tprintf	printf
# define _tmain		main
# define _ttoi      atoi
#endif

#include <stdio.h>
//#define SO_MAX_ARGS   54 // fixed size, no C lib use
#include "simpleopt.h"
#include "simpleglob.h"

static void ShowUsage()
{
    _tprintf(
        _T("Usage: fullSample [OPTIONS] [FILES]\n")
        _T("\n")
        _T("--exact         Disallow partial matching of option names\n")
        _T("--noslash       Disallow use of slash as an option marker on Windows\n")
        _T("--shortarg      Permit arguments on single letter options with no equals sign\n")
        _T("--clump         Permit single char options to be clumped as long string\n")
        _T("--noerr         Do not generate any errors for invalid options\n")
        _T("--pedantic      Generate an error for petty things\n")
        _T("--icase         Case-insensitive for all types\n")
        _T("--icase-short   Case-insensitive for short args\n")
        _T("--icase-long    Case-insensitive for long argsn")
        _T("--icase-word    Case-insensitive for word args\n")
        _T("\n")
        _T("-d  -e  -E  -f  -F  -g  -flag  --flag       Flag (no arg)\n")
        _T("-s ARG   -sep ARG  --sep ARG                Separate required arg\n")
        _T("-S ARG   -SEP ARG  --SEP ARG                Separate required arg (uppercase)\n")
        _T("-cARG    -c=ARG    -com=ARG    --com=ARG    Combined required arg\n")
        _T("-o[ARG]  -o[=ARG]  -opt[=ARG]  --opt[=ARG]  Combined optional arg\n")
        _T("-man     -mandy    -mandate                 Shortcut matching tests\n")
        _T("--man    --mandy   --mandate                Shortcut matching tests\n")
        _T("--multi0 --multi1 ARG --multi2 ARG1 ARG2    Multiple argument tests\n")
        _T("--multi N ARG-1 ARG-2 ... ARG-N             Multiple argument tests\n")
        _T("--                                          Stop argument processing\n")
        _T("open read write close zip unzip UPCASE      Special words\n")
        _T("\n")
        _T("-?  -h  -help  --help                       Output this help.\n")
        _T("\n")
        _T("If a FILE is `-', read standard input.\n")
        );
}

SOption g_rgFlags[] =
{
    { SO_O_EXACT,       _T("--exact"),          SO_NONE },
    { SO_O_NOSLASH,     _T("--noslash"),        SO_NONE },
    { SO_O_SHORTARG,    _T("--shortarg"),       SO_NONE },
    { SO_O_CLUMP,       _T("--clump"),          SO_NONE },
    { SO_O_NOERR,       _T("--noerr"),          SO_NONE },
    { SO_O_PEDANTIC,    _T("--pedantic"),       SO_NONE },
    { SO_O_ICASE,       _T("--icase"),          SO_NONE },
    { SO_O_ICASE_SHORT, _T("--icase-short"),    SO_NONE },
    { SO_O_ICASE_LONG,  _T("--icase-long"),     SO_NONE },
    { SO_O_ICASE_WORD,  _T("--icase-word"),     SO_NONE },
    SO_END_OF_OPTIONS
};

enum { OPT_HELP = 0, OPT_MULTI = 100, OPT_MULTI0, OPT_MULTI1, OPT_MULTI2, OPT_STOP  };
SOption g_rgOptions[] =
{
    { OPT_HELP,   _T("-?"),           SO_NONE    },
    { OPT_HELP,   _T("-h"),           SO_NONE    },
    { OPT_HELP,   _T("-help"),        SO_NONE    },
    { OPT_HELP,   _T("--help"),       SO_NONE    },
    {  1,         _T("-"),            SO_NONE    },
    {  2,         _T("-d"),           SO_NONE    },
    {  3,         _T("-e"),           SO_NONE    },
    {  4,         _T("-f"),           SO_NONE    },
    {  5,         _T("-g"),           SO_NONE    },
    {  6,         _T("-flag"),        SO_NONE    },
    {  7,         _T("--flag"),       SO_NONE    },
    {  8,         _T("-s"),           SO_REQ_SEP },
    {  9,         _T("-sep"),         SO_REQ_SEP },
    { 10,         _T("--sep"),        SO_REQ_SEP },
    { 11,         _T("-c"),           SO_REQ_CMB },
    { 12,         _T("-com"),         SO_REQ_CMB },
    { 13,         _T("--com"),        SO_REQ_CMB },
    { 14,         _T("-o"),           SO_OPT     },
    { 15,         _T("-opt"),         SO_OPT     },
    { 16,         _T("--opt"),        SO_OPT     },
    { 17,         _T("-man"),         SO_NONE    },
    { 18,         _T("-mandy"),       SO_NONE    },
    { 19,         _T("-mandate"),     SO_NONE    },
    { 20,         _T("--man"),        SO_NONE    },
    { 21,         _T("--mandy"),      SO_NONE    },
    { 22,         _T("--mandate"),    SO_NONE    },
    { 23,         _T("open"),         SO_NONE    },
    { 24,         _T("read"),         SO_NONE    },
    { 25,         _T("write"),        SO_NONE    },
    { 26,         _T("close"),        SO_NONE    },
    { 27,         _T("zip"),          SO_NONE    },
    { 28,         _T("unzip"),        SO_NONE    },
    { 29,         _T("-E"),           SO_NONE    },
    { 30,         _T("-F"),           SO_NONE    },
    { 31,         _T("-S"),           SO_REQ_SEP },
    { 32,         _T("-SEP"),         SO_REQ_SEP },
    { 33,         _T("--SEP"),        SO_REQ_SEP },
    { 34,         _T("UPCASE"),       SO_NONE    },
    { OPT_MULTI,  _T("--multi"),      SO_MULTI   },
    { OPT_MULTI0, _T("--multi0"),     SO_MULTI   },
    { OPT_MULTI1, _T("--multi1"),     SO_MULTI   },
    { OPT_MULTI2, _T("--multi2"),     SO_MULTI   },
    { OPT_STOP,   _T("--"),           SO_NONE    },
    SO_END_OF_OPTIONS
};

void ShowFiles(int argc, TCHAR ** argv) {
  int n;
    // glob files to catch expand wildcards
    CSimpleGlob (SG_GLOB_NODOT|SG_GLOB_NOCHECK, 0);
    if (SG_SUCCESS != Add2(argc, argv)) {
        _tprintf(_T("Error while globbing files\n"));
        return;
    }

    for (n = 0; n < FileCount(); ++n) {
        _tprintf(_T("file %2d: '%s'\n"), n, File(n));
    }
}

static const TCHAR * 
GetLastErrorText(
    int a_nError
    ) 
{
    switch (a_nError) {
    case SO_SUCCESS:            return _T("Success");
    case SO_OPT_INVALID:        return _T("Unrecognized option");
    case SO_OPT_MULTIPLE:       return _T("Option matched multiple strings");
    case SO_ARG_INVALID:        return _T("Option does not accept argument");
    case SO_ARG_INVALID_TYPE:   return _T("Invalid argument format");
    case SO_ARG_MISSING:        return _T("Required argument is missing");
    case SO_ARG_INVALID_DATA:   return _T("Invalid argument data");
    default:                    return _T("Unknown error");
    }
}

static void 
DoMultiArgs(
    int             nMultiArgs
    )
{
int n;
    TCHAR ** rgpszArg = NULL;

    // get the number of arguments if necessary
    if (nMultiArgs == -1) {
        // first arg is a count of how many we have
        rgpszArg = MultiArg(1);
        if (!rgpszArg) {
            _tprintf(
                _T("%s: '%s' (use --help to get command line help)\n"),
                GetLastErrorText(LastError()), OptionText());
            return;
        }

        nMultiArgs = _ttoi(rgpszArg[0]);
    }
    _tprintf(_T("%s: expecting %d args\n"), OptionText(), nMultiArgs);

    // get the arguments to follow
    rgpszArg = MultiArg(nMultiArgs);
    if (!rgpszArg) {
        _tprintf(
            _T("%s: '%s' (use --help to get command line help)\n"),
            GetLastErrorText(LastError()), OptionText());
        return;
    }

    for (n = 0; n < nMultiArgs; ++n) {
        _tprintf(_T("MultiArg %d: %s\n"), n, rgpszArg[n]);
    }
}


int _tmain(int argc, TCHAR * argv[]) {
    // process the command line to extract that flags for SimpleOpt 
    int nFlags = SO_O_USEALL|SO_O_SHORTARG;
/*
    CSimpleOpt (argc, argv, g_rgFlags, SO_O_NOERR|SO_O_EXACT);
    while (Next()) {
        nFlags |= OptionId();
    }
*/
    // now process the remainder of the command line with these flags
    //Init(FileCount(), Files(), g_rgOptions, nFlags);
    CSimpleOpt (argc, argv, g_rgOptions, nFlags);
    while (Next()) {
        if (LastError() != SO_SUCCESS) {
            _tprintf(
                _T("%s: '%s' (use --help to get command line help)\n"),
                GetLastErrorText(LastError()), OptionText());
            continue;
        }

        switch (OptionId()) {
        case OPT_HELP:
            ShowUsage();
            return 0;
        case OPT_MULTI:  
            DoMultiArgs(-1);
            break;
        case OPT_MULTI0: 
            DoMultiArgs(0);
            break;
        case OPT_MULTI1: 
            DoMultiArgs(1);
            break;
        case OPT_MULTI2: 
            DoMultiArgs(2);
            break;
        case OPT_STOP:
            Stop();
            break;
        default:
            if (OptionArg()) {
                _tprintf(_T("option: %2d, text: '%s', arg: '%s'\n"),
                    OptionId(), OptionText(), OptionArg());
            }
            else {
                _tprintf(_T("option: %2d, text: '%s'\n"),
                    OptionId(), OptionText());
            }
        }
    }

    /* process files from command line */
    ShowFiles(FileCount(), Files());

    CSimpleOptFinalize ();


	return 0;
}
