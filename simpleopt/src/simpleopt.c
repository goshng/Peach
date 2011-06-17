#include "simpleopt.h"

static    const SOption * m_rgOptions;     //!< pointer to options table 
static    int             m_nFlags;        //!< flags 
static    int             m_nOptionIdx;    //!< current argv option index
static    int             m_nOptionId;     //!< id of current option (-1 = invalid)
static    int             m_nNextOption;   //!< index of next option 
static    int             m_nLastArg;      //!< last argument, after this are files
static    int             m_argc;          //!< argc to process
static    char **       m_argv;          //!< argv
static    const char *  m_pszOptionText; //!< curr option text, e.g. "-f"
static    char *        m_pszOptionArg;  //!< curr option arg, e.g. "c:\file.txt"
static    char *        m_pszClump;      //!< clumped single character options
static    char          m_szShort[3];    //!< temp for clump and combined args
static    ESOError        m_nLastError;    //!< error status from the last call
static    char **       m_rgShuffleBuf;  //!< shuffle buffer for large argc

static void Copy(char ** ppDst, char ** ppSrc, int nCount); 
static char PrepareArg(char * a_pszString);
static int NextClumped();
static void ShuffleArg(int a_nStartIdx, int a_nCount);
static int LookupOption(const char * a_pszOption);
static int CalcMatch(const char *a_pszSource, const char *a_pszTest);
static int IsEqual(char a_cLeft, char a_cRight, int a_nArgType);

static char * FindEquals(char *s) {
        while (*s && *s != (char)'=') ++s;
        return *s ? s : NULL;
    }





static void SetOptions(const SOption * a_rgOptions) { 
        m_rgOptions = a_rgOptions; 
    }

static void SetFlags(int a_nFlags) { m_nFlags = a_nFlags; }

static int HasFlag(int a_nFlag) { 
        return (m_nFlags & a_nFlag) == a_nFlag; 
    }
ESOError LastError() { return m_nLastError; }
int OptionId() { return m_nOptionId; }
const char * OptionText() { return m_pszOptionText; }
char * OptionArg() { return m_pszOptionArg; }
static int FileCount() { return m_argc - m_nLastArg; }
static char * File(int n) {
        SO_ASSERT(n >= 0 && n < FileCount());
        return m_argv[m_nLastArg + n];
    }
static char ** Files() { return &m_argv[m_nLastArg]; }

int Init(
    int             a_argc,
    char *        a_argv[],
    const SOption * a_rgOptions,
    int             a_nFlags
    )
{
    m_argc           = a_argc;
    m_nLastArg       = a_argc;
    m_argv           = a_argv;
    m_rgOptions      = a_rgOptions;
    m_nLastError     = SO_SUCCESS;
    m_nOptionIdx     = 0;
    m_nOptionId      = -1;
    m_pszOptionText  = NULL;
    m_pszOptionArg   = NULL;
    m_nNextOption    = (a_nFlags & SO_O_USEALL) ? 0 : 1;
    m_szShort[0]     = (char)'-';
    m_szShort[2]     = (char)'\0';
    m_nFlags         = a_nFlags;
    m_pszClump       = NULL;

#ifdef SO_MAX_ARGS
	if (m_argc > SO_MAX_ARGS) {
        m_nLastError = SO_ARG_INVALID_DATA;
        m_nLastArg = 0;
		return 0;
	}
#else
    if (m_rgShuffleBuf) {
        free(m_rgShuffleBuf);
    }
    if (m_argc > SO_STATICBUF) {
        m_rgShuffleBuf = (char**) malloc(sizeof(char*) * m_argc);
        if (!m_rgShuffleBuf) {
            return 0;
        }
    }
#endif

    return 1;
}

int Next()
{
#ifdef SO_MAX_ARGS
    if (m_argc > SO_MAX_ARGS) {
        SO_ASSERT(!"Too many args! Check the return value of Init()!");
        return 0;
    }
#endif

    // process a clumped option string if appropriate
    if (m_pszClump && *m_pszClump) {
        // silently discard invalid clumped option
        int bIsValid = NextClumped();
        while (*m_pszClump && !bIsValid && HasFlag(SO_O_NOERR)) {
            bIsValid = NextClumped();
        }

        // return this option if valid or we are returning errors
        if (bIsValid || !HasFlag(SO_O_NOERR)) {
            return 1;
        }
    }
    SO_ASSERT(!m_pszClump || !*m_pszClump);
    m_pszClump = NULL;

    // init for the next option
    m_nOptionIdx    = m_nNextOption;
    m_nOptionId     = -1;
    m_pszOptionText = NULL;
    m_pszOptionArg  = NULL;
    m_nLastError    = SO_SUCCESS;

    // find the next option
    char cFirst;
    int nTableIdx = -1;
    int nOptIdx = m_nOptionIdx;
    while (nTableIdx < 0 && nOptIdx < m_nLastArg) {
        char * pszArg = m_argv[nOptIdx];
        m_pszOptionArg  = NULL;

        // find this option in the options table
        cFirst = PrepareArg(pszArg);
        if (pszArg[0] == (char)'-') {
            // find any combined argument string and remove equals sign
            m_pszOptionArg = FindEquals(pszArg);
            if (m_pszOptionArg) {
                *m_pszOptionArg++ = (char)'\0';
            }
        }
        nTableIdx = LookupOption(pszArg);

        // if we didn't find this option but if it is a short form
        // option then we try the alternative forms
        if (nTableIdx < 0
            && !m_pszOptionArg
            && pszArg[0] == (char)'-'
            && pszArg[1]
            && pszArg[1] != (char)'-'
            && pszArg[2])
        {
            // test for a short-form with argument if appropriate
            if (HasFlag(SO_O_SHORTARG)) {
                m_szShort[1] = pszArg[1];
                int nIdx = LookupOption(m_szShort);
                if (nIdx >= 0
                    && (m_rgOptions[nIdx].nArgType == SO_REQ_CMB
                        || m_rgOptions[nIdx].nArgType == SO_OPT))
                {
                    m_pszOptionArg = &pszArg[2];
                    pszArg         = m_szShort;
                    nTableIdx      = nIdx;
                }
            }

            // test for a clumped short-form option string and we didn't
            // match on the short-form argument above
            if (nTableIdx < 0 && HasFlag(SO_O_CLUMP))  {
                m_pszClump = &pszArg[1];
                ++m_nNextOption;
                if (nOptIdx > m_nOptionIdx) {
                    ShuffleArg(m_nOptionIdx, nOptIdx - m_nOptionIdx);
                }
                return Next();
            }
        }

        // The option wasn't found. If it starts with a switch character
        // and we are not suppressing errors for invalid options then it
        // is reported as an error, otherwise it is data.
        if (nTableIdx < 0) {
            if (!HasFlag(SO_O_NOERR) && pszArg[0] == (char)'-') {
                m_pszOptionText = pszArg;
                break;
            }
            
            pszArg[0] = cFirst;
            ++nOptIdx;
            if (m_pszOptionArg) {
                *(--m_pszOptionArg) = (char)'=';
            }
        }
    }

    // end of options
    if (nOptIdx >= m_nLastArg) {
        if (nOptIdx > m_nOptionIdx) {
            ShuffleArg(m_nOptionIdx, nOptIdx - m_nOptionIdx);
        }
        return 0;
    }
    ++m_nNextOption;

    // get the option id
    ESOArgType nArgType = SO_NONE;
    if (nTableIdx < 0) {
        m_nLastError    = (ESOError) nTableIdx; // error code
    }
    else {
        m_nOptionId     = m_rgOptions[nTableIdx].nId;
        m_pszOptionText = m_rgOptions[nTableIdx].pszArg;

        // ensure that the arg type is valid
        nArgType = m_rgOptions[nTableIdx].nArgType;
        switch (nArgType) {
        case SO_NONE:
            if (m_pszOptionArg) {
                m_nLastError = SO_ARG_INVALID;
            }
            break;

        case SO_REQ_SEP:
            if (m_pszOptionArg) {
                // they wanted separate args, but we got a combined one, 
                // unless we are pedantic, just accept it.
                if (HasFlag(SO_O_PEDANTIC)) {
                    m_nLastError = SO_ARG_INVALID_TYPE;
                }
            }
            // more processing after we shuffle
            break;

        case SO_REQ_CMB:
            if (!m_pszOptionArg) {
                m_nLastError = SO_ARG_MISSING;
            }
            break;

        case SO_OPT:
            // nothing to do
            break;

        case SO_MULTI:
            // nothing to do. Caller must now check for valid arguments
            // using GetMultiArg()
            break;
        }
    }

    // shuffle the files out of the way
    if (nOptIdx > m_nOptionIdx) {
        ShuffleArg(m_nOptionIdx, nOptIdx - m_nOptionIdx);
    }

    // we need to return the separate arg if required, just re-use the
    // multi-arg code because it all does the same thing
    if (   nArgType == SO_REQ_SEP 
        && !m_pszOptionArg 
        && m_nLastError == SO_SUCCESS) 
    {
        char ** ppArgs = MultiArg(1);
        if (ppArgs) {
            m_pszOptionArg = *ppArgs;
        }
    }

    return 1;
}

void Stop()
{
    if (m_nNextOption < m_nLastArg) {
        ShuffleArg(m_nNextOption, m_nLastArg - m_nNextOption);
    }
}

static char PrepareArg(
    char * a_pszString
    ) 
{
#ifdef _WIN32
    // On Windows we can accept the forward slash as a single character
    // option delimiter, but it cannot replace the '-' option used to
    // denote stdin. On Un*x paths may start with slash so it may not
    // be used to start an option.
    if (!HasFlag(SO_O_NOSLASH)
        && a_pszString[0] == (char)'/'
        && a_pszString[1]
        && a_pszString[1] != (char)'-')
    {
        a_pszString[0] = (char)'-';
        return (char)'/';
    }
#endif
    return a_pszString[0];
}

static int NextClumped()
{
    // prepare for the next clumped option
    m_szShort[1]    = *m_pszClump++;
    m_nOptionId     = -1;
    m_pszOptionText = NULL;
    m_pszOptionArg  = NULL;
    m_nLastError    = SO_SUCCESS;

    // lookup this option, ensure that we are using exact matching
    int nSavedFlags = m_nFlags;
    m_nFlags = SO_O_EXACT;
    int nTableIdx = LookupOption(m_szShort);
    m_nFlags = nSavedFlags;

    // unknown option
    if (nTableIdx < 0) {
        m_nLastError = (ESOError) nTableIdx; // error code
        return 0;
    }

    // valid option
    m_pszOptionText = m_rgOptions[nTableIdx].pszArg;
    ESOArgType nArgType = m_rgOptions[nTableIdx].nArgType;
    if (nArgType == SO_NONE) {
        m_nOptionId = m_rgOptions[nTableIdx].nId;
        return 1;
    }

    if (nArgType == SO_REQ_CMB && *m_pszClump) {
        m_nOptionId = m_rgOptions[nTableIdx].nId;
        m_pszOptionArg = m_pszClump;
        while (*m_pszClump) ++m_pszClump; // must point to an empty string
        return 1;
    }

    // invalid option as it requires an argument
    m_nLastError = SO_ARG_MISSING;
    return 1;
}

// Shuffle arguments to the end of the argv array.
//
// For example:
//      argv[] = { "0", "1", "2", "3", "4", "5", "6", "7", "8" };
//
//  ShuffleArg(1, 1) = { "0", "2", "3", "4", "5", "6", "7", "8", "1" };
//  ShuffleArg(5, 2) = { "0", "1", "2", "3", "4", "7", "8", "5", "6" };
//  ShuffleArg(2, 4) = { "0", "1", "6", "7", "8", "2", "3", "4", "5" };
static void ShuffleArg(
    int a_nStartIdx,
    int a_nCount
    )
{
    char * staticBuf[SO_STATICBUF];
    char ** buf = m_rgShuffleBuf ? m_rgShuffleBuf : staticBuf;
    int nTail = m_argc - a_nStartIdx - a_nCount;

    // make a copy of the elements to be moved
    Copy(buf, m_argv + a_nStartIdx, a_nCount);

    // move the tail down
    Copy(m_argv + a_nStartIdx, m_argv + a_nStartIdx + a_nCount, nTail);

    // append the moved elements to the tail
    Copy(m_argv + a_nStartIdx + nTail, buf, a_nCount);

    // update the index of the last unshuffled arg
    m_nLastArg -= a_nCount;
}

static int LookupOption(
    const char * a_pszOption
    ) 
{
    int nBestMatch = -1;    // index of best match so far
    int nBestMatchLen = 0;  // matching characters of best match
    int nLastMatchLen = 0;  // matching characters of last best match
    int n;

    for (n = 0; m_rgOptions[n].nId >= 0; ++n) {
        // the option table must use hyphens as the option character,
        // the slash character is converted to a hyphen for testing.
        SO_ASSERT(m_rgOptions[n].pszArg[0] != (char)'/');

        int nMatchLen = CalcMatch(m_rgOptions[n].pszArg, a_pszOption);
        if (nMatchLen == -1) {
            return n;
        }
        if (nMatchLen > 0 && nMatchLen >= nBestMatchLen) {
            nLastMatchLen = nBestMatchLen;
            nBestMatchLen = nMatchLen;
            nBestMatch = n;
        }
    }

    // only partial matches or no match gets to here, ensure that we
    // don't return a partial match unless it is a clear winner
    if (HasFlag(SO_O_EXACT) || nBestMatch == -1) {
        return SO_OPT_INVALID;
    }
    return (nBestMatchLen > nLastMatchLen) ? nBestMatch : SO_OPT_MULTIPLE;
}

// calculate the number of characters that match (case-sensitive)
// 0 = no match, > 0 == number of characters, -1 == perfect match
static int CalcMatch(
    const char *  a_pszSource,
    const char *  a_pszTest
    ) 
{
    if (!a_pszSource || !a_pszTest) {
        return 0;
    }

    // determine the argument type
    int nArgType = SO_O_ICASE_LONG;
    if (a_pszSource[0] != '-') {
        nArgType = SO_O_ICASE_WORD;
    }
    else if (a_pszSource[1] != '-' && !a_pszSource[2]) {
        nArgType = SO_O_ICASE_SHORT;
    }

    // match and skip leading hyphens
    while (*a_pszSource == (char)'-' && *a_pszSource == *a_pszTest) {
        ++a_pszSource; 
        ++a_pszTest;
    }
    if (*a_pszSource == (char)'-' || *a_pszTest == (char)'-') {
        return 0;
    }

    // find matching number of characters in the strings
    int nLen = 0;
    while (*a_pszSource && IsEqual(*a_pszSource, *a_pszTest, nArgType)) {
        ++a_pszSource; 
        ++a_pszTest; 
        ++nLen;
    }

    // if we have exhausted the source...
    if (!*a_pszSource) {
        // and the test strings, then it's a perfect match
        if (!*a_pszTest) {
            return -1;
        }

        // otherwise the match failed as the test is longer than
        // the source. i.e. "--mant" will not match the option "--man".
        return 0;
    }

    // if we haven't exhausted the test string then it is not a match
    // i.e. "--mantle" will not best-fit match to "--mandate" at all.
    if (*a_pszTest) {
        return 0;
    }

    // partial match to the current length of the test string
    return nLen;
}

static int IsEqual(
    char  a_cLeft,
    char  a_cRight,
    int     a_nArgType
    ) 
{
    // if this matches then we are doing case-insensitive matching
    if (m_nFlags & a_nArgType) {
        if (a_cLeft  >= 'A' && a_cLeft  <= 'Z') a_cLeft  += 'a' - 'A';
        if (a_cRight >= 'A' && a_cRight <= 'Z') a_cRight += 'a' - 'A';
    }
    return a_cLeft == a_cRight;
}

// calculate the number of characters that match (case-sensitive)
// 0 = no match, > 0 == number of characters, -1 == perfect match
char ** MultiArg(
    int a_nCount
    )
{
int n;
    // ensure we have enough arguments
    if (m_nNextOption + a_nCount > m_nLastArg) {
        m_nLastError = SO_ARG_MISSING;
        return NULL;
    }

    // our argument array
    char ** rgpszArg = &m_argv[m_nNextOption];

    // Ensure that each of the following don't start with an switch character.
    // Only make this check if we are returning errors for unknown arguments.
    if (!HasFlag(SO_O_NOERR)) {
        for (n = 0; n < a_nCount; ++n) {
            char ch = PrepareArg(rgpszArg[n]);
            if (rgpszArg[n][0] == (char)'-') {
                rgpszArg[n][0] = ch;
                m_nLastError = SO_ARG_INVALID_DATA;
                return NULL;
            }
            rgpszArg[n][0] = ch;
        }
    }

    // all good
    m_nNextOption += a_nCount;
    return rgpszArg;
}

static void Copy(char ** ppDst, char ** ppSrc, int nCount) {
#ifdef SO_MAX_ARGS
        // keep our promise of no CLIB usage
        while (nCount-- > 0) *ppDst++ = *ppSrc++;
#else
        memcpy(ppDst, ppSrc, nCount * sizeof(char*));
#endif
    }


void CSimpleOpt (
        int             argc, 
        char *        argv[], 
        const SOption * a_rgOptions, 
        int             a_nFlags
        ) 
{ 
  m_rgShuffleBuf = NULL;
  Init(argc, argv, a_rgOptions, a_nFlags); 
}


void CSimpleOptFinalize () { if (m_rgShuffleBuf) free(m_rgShuffleBuf); }

