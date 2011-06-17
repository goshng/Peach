/*! @file SimpleGlob.h

    @version 3.5

    @brief A cross-platform file globbing library providing the ability to
    expand wildcards in command-line arguments to a list of all matching 
    files. It is designed explicitly to be portable to any platform and has 
    been tested on Windows and Linux. See CSimpleGlobTempl for the class 
    definition.

    @section features FEATURES

    -   MIT Licence allows free use in all software (including GPL and 
        commercial)
    -   multi-platform (Windows 95/98/ME/NT/2K/XP, Linux, Unix)
    -   supports most of the standard linux glob() options
    -   recognition of a forward paths as equivalent to a backward slash 
        on Windows. e.g. "c:/path/foo*" is equivalent to "c:\path\foo*".
    -   implemented with only a single C++ header file
    -   char, wchar_t and Windows TCHAR in the same program
    -   complete working examples included
    -   compiles cleanly at warning level 4 (Windows/VC.NET 2003), 
        warning level 3 (Windows/VC6) and -Wall (Linux/gcc)

    @section usage USAGE

    The SimpleGlob class is used by following these steps:

    <ol>
    <li> Include the SimpleGlob.h header file

        <pre>
        \#include "SimpleGlob.h"
        </pre>

   <li> Instantiate a CSimpleGlob object supplying the appropriate flags.

        <pre>
        @link CSimpleGlobTempl CSimpleGlob @endlink glob(FLAGS);
        </pre>

   <li> Add all file specifications to the glob class.

        <pre>
        glob.Add("file*");
        glob.Add(argc, argv);
        </pre>

   <li> Process all files with File(), Files() and FileCount()

        <pre>
        for (int n = 0; n < glob.FileCount(); ++n) {
            ProcessFile(glob.File(n));
        }
        </pre>

    </ol>

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

#ifndef INCLUDED_SimpleGlob
#define INCLUDED_SimpleGlob

/*! @brief The operation of SimpleGlob is fine-tuned via the use of a 
    combination of the following flags.

    The flags may be passed at initialization of the class and used for every
    filespec added, or alternatively they may optionally be specified in the
    call to Add() and be different for each filespec.

    @param SG_GLOB_ERR
        Return upon read error (e.g. directory does not have read permission)

    @param SG_GLOB_MARK
        Append a slash (backslash in Windows) to every path which corresponds
        to a directory

    @param SG_GLOB_NOSORT
        By default, files are returned in sorted into string order. With this
        flag, no sorting is done. This is not compatible with 
        SG_GLOB_FULLSORT.

    @param SG_GLOB_FULLSORT
        By default, files are sorted in groups belonging to each filespec that
        was added. For example if the filespec "b*" was added before the 
        filespec "a*" then the argv array will contain all b* files sorted in 
        order, followed by all a* files sorted in order. If this flag is 
        specified, the entire array will be sorted ignoring the filespec 
        groups.

    @param SG_GLOB_NOCHECK
        If the pattern doesn't match anything, return the original pattern.

    @param SG_GLOB_TILDE
        Tilde expansion is carried out (on Unix platforms)

    @param SG_GLOB_ONLYDIR
        Return only directories which match (not compatible with 
        SG_GLOB_ONLYFILE)

    @param SG_GLOB_ONLYFILE
        Return only files which match (not compatible with SG_GLOB_ONLYDIR)

    @param SG_GLOB_NODOT
        Do not return the "." or ".." special directories.
 */
enum SG_Flags {
    SG_GLOB_ERR         = 1 << 0,
    SG_GLOB_MARK        = 1 << 1,
    SG_GLOB_NOSORT      = 1 << 2,
    SG_GLOB_NOCHECK     = 1 << 3,
    SG_GLOB_TILDE       = 1 << 4,
    SG_GLOB_ONLYDIR     = 1 << 5,
    SG_GLOB_ONLYFILE    = 1 << 6,
    SG_GLOB_NODOT       = 1 << 7,
    SG_GLOB_FULLSORT    = 1 << 8
};

/*! @brief Error return codes */
enum SG_Error {
    SG_SUCCESS          =  0,
    SG_ERR_NOMATCH      =  1,
    SG_ERR_MEMORY       = -1,
    SG_ERR_FAILURE      = -2
};

// ---------------------------------------------------------------------------
// Platform dependent implementations

// if we aren't on Windows and we have ICU available, then enable ICU
// by default. Define this to 0 to intentially disable it.
#ifndef SG_HAVE_ICU
# if !defined(_WIN32) && defined(USTRING_H)
#   define SG_HAVE_ICU 1
# else
#   define SG_HAVE_ICU 0
# endif
#endif

// don't include this in documentation as it isn't relevant
#ifndef DOXYGEN

// on Windows we want to use MBCS aware string functions and mimic the
// Unix glob functionality. On Unix we just use glob.
#ifdef _WIN32
# include <mbstring.h>
# define sg_strchr          ::_mbschr
# define sg_strrchr         ::_mbsrchr
# define sg_strlen          ::_mbslen
# if __STDC_WANT_SECURE_LIB__
#  define sg_strcpy_s(a,n,b) ::_mbscpy_s(a,n,b)
# else
#  define sg_strcpy_s(a,n,b) ::_mbscpy(a,b)
# endif
# define sg_strcmp          ::_mbscmp
# define sg_strcasecmp      ::_mbsicmp
# define char_T           unsigned char
#else
# include <sys/types.h>
# include <sys/stat.h>
# include <glob.h>
# include <limits.h>
# define MAX_PATH           PATH_MAX
# define sg_strchr          strchr
# define sg_strrchr         strrchr
# define sg_strlen          strlen
# define sg_strcpy_s(a,n,b) strcpy(a,b)
# define sg_strcmp          strcmp
# define sg_strcasecmp      strcasecmp
# define char_T           char
#endif

#include <stdlib.h>
#include <string.h>
#include <wchar.h>

// use assertions to test the input data
#ifdef _DEBUG
# ifdef _MSC_VER
#  include <crtdbg.h>
#  define SG_ASSERT(b)    _ASSERTE(b)
# else
#  include <assert.h>
#  define SG_ASSERT(b)    assert(b)
# endif
#else
# define SG_ASSERT(b)
#endif

typedef enum _SG_FileType {
    SG_FILETYPE_INVALID,
    SG_FILETYPE_FILE,
    SG_FILETYPE_DIR
} SG_FileType;

#endif // DOXYGEN

// ---------------------------------------------------------------------------
//                              MAIN TEMPLATE CLASS
// ---------------------------------------------------------------------------

void CSimpleGlob(unsigned int a_uiFlags, int a_nReservedSlots);
void CSimpleGlobFinalize();
typedef enum _ARG_ARRAY_TYPE { OFFSETS, POINTERS } ARG_ARRAY_TYPE;
int FileCount(); 
char * File(int n); 
char ** Files(); 
int Add(const char *a_pszFileSpec);
int Add2(int a_nCount, const char * const * a_rgpszFileSpec);





#endif // INCLUDED_SimpleGlob
