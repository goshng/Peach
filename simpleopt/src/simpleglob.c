#include "simpleglob.h"

    static unsigned int        m_uiFlags;
    static ARG_ARRAY_TYPE      m_nArgArrayType;    //!< argv is indexes or pointers
    static char **           m_rgpArgs;          //!< argv 
    static int                 m_nReservedSlots;   //!< # client slots in argv array
    static int                 m_nArgsSize;        //!< allocated size of array
    static int                 m_nArgsLen;         //!< used length
    static char *            m_pBuffer;          //!< argv string buffer
    static size_t              m_uiBufferSize;     //!< allocated size of buffer
    static size_t              m_uiBufferLen;      //!< used length of buffer
    static char              m_szPathPrefix[MAX_PATH]; //!< wildcard path prefix

/* class SimpleGlobUtil */
    static const char * SimpleGlobUtil_strchr(const char *s, char c) {
        return (char *) sg_strchr((const char_T *)s, c);
    }
    static const wchar_t * SimpleGlobUtilW_strchr(const wchar_t *s, wchar_t c) {
        return wcschr(s, c);
    }

    static const char * SimpleGlobUtil_strrchr(const char *s, char c) {
        return (char *) sg_strrchr((const char_T *)s, c);
    }
    static const wchar_t * SimpleGlobUtilW_strrchr(const wchar_t *s, wchar_t c) {
        return wcsrchr(s, c);
    }

    // Note: char strlen returns number of bytes, not characters
    static size_t SimpleGlobUtil_strlen(const char *s) { return strlen(s); }
    static size_t SimpleGlobUtilW_strlen(const wchar_t *s) { return wcslen(s); }

    static void SimpleGlobUtil_strcpy_s(char *dst, size_t n, const char *src)  {
        (void) n;
        sg_strcpy_s((char_T *)dst, n, (const char_T *)src);
    }
    static void SimpleGlobUtilW_strcpy_s(wchar_t *dst, size_t n, const wchar_t *src) {
# if __STDC_WANT_SECURE_LIB__
        wcscpy_s(dst, n, src);
#else
        (void) n;
        wcscpy(dst, src);
#endif
    }

    static int SimpleGlobUtil_strcmp(const char *s1, const char *s2) {
        return sg_strcmp((const char_T *)s1, (const char_T *)s2);
    }
    static int SimpleGlobUtilW_strcmp(const wchar_t *s1, const wchar_t *s2) {
        return wcscmp(s1, s2);
    }

    static int SimpleGlobUtil_strcasecmp(const char *s1, const char *s2) {
        return sg_strcasecmp((const char_T *)s1, (const char_T *)s2);
    }
#if _WIN32
    static int SimpleGlobUtilW_strcasecmp(const wchar_t *s1, const wchar_t *s2) {
        return _wcsicmp(s1, s2);
    }
#endif // _WIN32


static glob_t  m_glob;
static size_t  m_uiCurr;
static int    m_bIsDir;


#ifdef _WIN32

#ifndef INVALID_FILE_ATTRIBUTES
# define INVALID_FILE_ATTRIBUTES    ((DWORD)-1)
#endif

#define SG_PATH_CHAR    '\\'

/*! @brief Windows glob implementation. */
template<class char>
struct SimpleGlobBase
{
    SimpleGlobBase() : m_hFind(INVALID_HANDLE_VALUE) { }

    int FindFirstFileS(const char * a_pszFileSpec, unsigned int) {
        m_hFind = FindFirstFileA(a_pszFileSpec, &m_oFindDataA);
        if (m_hFind != INVALID_HANDLE_VALUE) {
            return SG_SUCCESS;
        }
        DWORD dwErr = GetLastError();
        if (dwErr == ERROR_FILE_NOT_FOUND) {
            return SG_ERR_NOMATCH;
        }
        return SG_ERR_FAILURE;
    }
    int FindFirstFileS(const wchar_t * a_pszFileSpec, unsigned int) {
        m_hFind = FindFirstFileW(a_pszFileSpec, &m_oFindDataW);
        if (m_hFind != INVALID_HANDLE_VALUE) {
            return SG_SUCCESS;
        }
        DWORD dwErr = GetLastError();
        if (dwErr == ERROR_FILE_NOT_FOUND) {
            return SG_ERR_NOMATCH;
        }
        return SG_ERR_FAILURE;
    }

    int FindNextFileS(char) {
        return FindNextFileA(m_hFind, &m_oFindDataA) != FALSE;
    }
    int FindNextFileS(wchar_t) {
        return FindNextFileW(m_hFind, &m_oFindDataW) != FALSE;
    }

    void FindDone() {
        FindClose(m_hFind);
    }

    const char * GetFileNameS(char) const {
        return m_oFindDataA.cFileName;
    }
    const wchar_t * GetFileNameS(wchar_t) const {
        return m_oFindDataW.cFileName;
    }

    int IsDirS(char) const {
        return GetFileTypeS(m_oFindDataA.dwFileAttributes) == SG_FILETYPE_DIR;
    }
    int IsDirS(wchar_t) const {
        return GetFileTypeS(m_oFindDataW.dwFileAttributes) == SG_FILETYPE_DIR;
    }

    SG_FileType GetFileTypeS(const char * a_pszPath) {
        return GetFileTypeS(GetFileAttributesA(a_pszPath));
    }
    SG_FileType GetFileTypeS(const wchar_t * a_pszPath)  {
        return GetFileTypeS(GetFileAttributesW(a_pszPath));
    }
    SG_FileType GetFileTypeS(DWORD a_dwAttribs) const {
        if (a_dwAttribs == INVALID_FILE_ATTRIBUTES) {
            return SG_FILETYPE_INVALID;
        }
        if (a_dwAttribs & FILE_ATTRIBUTE_DIRECTORY) {
            return SG_FILETYPE_DIR;
        }
        return SG_FILETYPE_FILE;
    }

private:
    HANDLE              m_hFind;
    WIN32_FIND_DATAA    m_oFindDataA;
    WIN32_FIND_DATAW    m_oFindDataW;
};

#else // !_WIN32

#define SG_PATH_CHAR    '/'

static    void FilePrep() {
        m_bIsDir = 0;
        size_t len = strlen(m_glob.gl_pathv[m_uiCurr]);
        if (m_glob.gl_pathv[m_uiCurr][len-1] == '/') {
            m_bIsDir = 1;
            m_glob.gl_pathv[m_uiCurr][len-1] = 0;
        }
    }

static    int FindFirstFileS(const char * a_pszFileSpec, unsigned int a_uiFlags) {
        int nFlags = GLOB_MARK | GLOB_NOSORT;
        if (a_uiFlags & SG_GLOB_ERR)    nFlags |= GLOB_ERR;
        if (a_uiFlags & SG_GLOB_TILDE)  nFlags |= GLOB_TILDE;
        int rc = glob(a_pszFileSpec, nFlags, NULL, &m_glob);
        if (rc == GLOB_NOSPACE) return SG_ERR_MEMORY;
        if (rc == GLOB_ABORTED) return SG_ERR_FAILURE;
        if (rc == GLOB_NOMATCH) return SG_ERR_NOMATCH;
        m_uiCurr = 0;
        FilePrep();
        return SG_SUCCESS;
    }

static    int FindNextFileS() {
        SG_ASSERT(m_uiCurr != (size_t)-1);
        if (++m_uiCurr >= m_glob.gl_pathc) {
            return 0;
        }
        FilePrep();
        return 1;
    }

static    void FindDone() {
        globfree(&m_glob);
        memset(&m_glob, 0, sizeof(m_glob));
        m_uiCurr = (size_t)-1;
    }

static    const char * GetFileNameS() {
        SG_ASSERT(m_uiCurr != (size_t)-1);
        return m_glob.gl_pathv[m_uiCurr];
    }

static    int IsDirS() {
        SG_ASSERT(m_uiCurr != (size_t)-1);
        return m_bIsDir;
    }

static    SG_FileType GetFileTypeS(const char * a_pszPath) {
        struct stat sb;
        if (0 != stat(a_pszPath, &sb)) {
            return SG_FILETYPE_INVALID;
        }
        if (S_ISDIR(sb.st_mode)) {
            return SG_FILETYPE_DIR;
        }
        if (S_ISREG(sb.st_mode)) {
            return SG_FILETYPE_FILE;
        }
        return SG_FILETYPE_INVALID;
    }

#endif // _WIN32










static    void SetArgvArrayType(ARG_ARRAY_TYPE a_nNewType);

static int Init(unsigned int a_uiFlags, int a_nReservedSlots);
int FileCount() { return m_nArgsLen; }
char ** Files() {
        SetArgvArrayType(POINTERS);
        return m_rgpArgs;
    }

    char * File(int n) {
        SG_ASSERT(n >= 0 && n < m_nArgsLen);
        return Files()[n];
    }
    static int AppendName(const char *a_pszFileName, int a_bIsDir);
    static int GrowArgvArray(int a_nNewLen);
    static int GrowStringBuffer(size_t a_uiMinSize);
    static int fileSortCompare(const void *a1, const void *a2);



void CSimpleGlob(
    unsigned int    a_uiFlags,
    int             a_nReservedSlots
    )
{
   /* From base class */
        memset(&m_glob, 0, sizeof(m_glob));
        m_uiCurr = (size_t)-1;

    m_rgpArgs           = NULL;
    m_nArgsSize         = 0;
    m_pBuffer           = NULL;
    m_uiBufferSize      = 0;

    Init(a_uiFlags, a_nReservedSlots);
}

void CSimpleGlobFinalize()
{
    if (m_rgpArgs) free(m_rgpArgs);
    if (m_pBuffer) free(m_pBuffer);

/* From base class */
        globfree(&m_glob);
}


static int Init(
    unsigned int    a_uiFlags,
    int             a_nReservedSlots
    )
{
int n;
    m_nArgArrayType     = POINTERS;
    m_uiFlags           = a_uiFlags;
    m_nArgsLen          = a_nReservedSlots;
    m_nReservedSlots    = a_nReservedSlots;
    m_uiBufferLen       = 0;

    if (m_nReservedSlots > 0) {
        if (!GrowArgvArray(m_nReservedSlots)) {
            return SG_ERR_MEMORY;
        }
        for (n = 0; n < m_nReservedSlots; ++n) {
            m_rgpArgs[n] = NULL;
        }
    }

    return SG_SUCCESS;
}

int Add(
    const char *a_pszFileSpec
    )
{
#ifdef _WIN32
    // Windows FindFirst/FindNext recognizes forward slash as the same as 
    // backward slash and follows the directories. We need to do the same 
    // when calculating the prefix and when we have no wildcards.
    char szFileSpec[MAX_PATH];
    SimpleGlobUtil_strcpy_s(szFileSpec, MAX_PATH, a_pszFileSpec);
    const char * pszPath = SimpleGlobUtil_strchr(szFileSpec, '/');
    while (pszPath) {
        szFileSpec[pszPath - szFileSpec] = SG_PATH_CHAR;
        pszPath = SimpleGlobUtil_strchr(pszPath + 1, '/');
    }
    a_pszFileSpec = szFileSpec;
#endif

    // if this doesn't contain wildcards then we can just add it directly
    m_szPathPrefix[0] = 0;
    if (!SimpleGlobUtil_strchr(a_pszFileSpec, '*') &&
        !SimpleGlobUtil_strchr(a_pszFileSpec, '?'))
    {
        SG_FileType nType = GetFileTypeS(a_pszFileSpec);
        if (nType == SG_FILETYPE_INVALID) {
            if (m_uiFlags & SG_GLOB_NOCHECK) {
                return AppendName(a_pszFileSpec, 0);
            }
            return SG_ERR_NOMATCH;
        }
        return AppendName(a_pszFileSpec, nType == SG_FILETYPE_DIR);
    }

#ifdef _WIN32
    // Windows doesn't return the directory with the filename, so we need to 
    // extract the path from the search string ourselves and prefix it to the 
    // filename we get back.
    const char * pszFilename = 
        SimpleGlobUtil_strrchr(a_pszFileSpec, SG_PATH_CHAR);
    if (pszFilename) {
        SimpleGlobUtil_strcpy_s(m_szPathPrefix, MAX_PATH, a_pszFileSpec);
        m_szPathPrefix[pszFilename - a_pszFileSpec + 1] = 0;
    }
#endif

    // search for the first match on the file
    int rc = FindFirstFileS(a_pszFileSpec, m_uiFlags);
    if (rc != SG_SUCCESS) {
        if (rc == SG_ERR_NOMATCH && (m_uiFlags & SG_GLOB_NOCHECK)) {
            int ok = AppendName(a_pszFileSpec, 0);
            if (ok != SG_SUCCESS) rc = ok;
        }
        return rc;
    }

    // add it and find all subsequent matches
    int nError, nStartLen = m_nArgsLen;
    int bSuccess;
    do {
        nError = AppendName(GetFileNameS(), IsDirS());
        bSuccess = FindNextFileS();
    }
    while (nError == SG_SUCCESS && bSuccess);
    FindDone();

    // sort these files if required
    if (m_nArgsLen > nStartLen && !(m_uiFlags & SG_GLOB_NOSORT)) {
        if (m_uiFlags & SG_GLOB_FULLSORT) {
            nStartLen = m_nReservedSlots;
        }
        SetArgvArrayType(POINTERS);
        qsort(
            m_rgpArgs + nStartLen,
            m_nArgsLen - nStartLen,
            sizeof(m_rgpArgs[0]), fileSortCompare);
    }

    return nError;
}

int
Add2(
    int                     a_nCount,
    const char * const *  a_rgpszFileSpec
    )
{
 int n;
    int nResult;
    for (n = 0; n < a_nCount; ++n) {
        nResult = Add(a_rgpszFileSpec[n]);
        if (nResult != SG_SUCCESS) {
            return nResult;
        }
    }
    return SG_SUCCESS;
}

int
AppendName(
    const char *  a_pszFileName,
    int            a_bIsDir
    )
{
    // we need the argv array as offsets in case we resize it
    SetArgvArrayType(OFFSETS);

    // check for special cases which cause us to ignore this entry
    if ((m_uiFlags & SG_GLOB_ONLYDIR) && !a_bIsDir) {
        return SG_SUCCESS;
    }
    if ((m_uiFlags & SG_GLOB_ONLYFILE) && a_bIsDir) {
        return SG_SUCCESS;
    }
    if ((m_uiFlags & SG_GLOB_NODOT) && a_bIsDir) {
        if (a_pszFileName[0] == '.') {
            if (a_pszFileName[1] == '\0') {
                return SG_SUCCESS;
            }
            if (a_pszFileName[1] == '.' && a_pszFileName[2] == '\0') {
                return SG_SUCCESS;
            }
        }
    }

    // ensure that we have enough room in the argv array
    if (!GrowArgvArray(m_nArgsLen + 1)) {
        return SG_ERR_MEMORY;
    }

    // ensure that we have enough room in the string buffer (+1 for null)
    size_t uiPrefixLen = SimpleGlobUtil_strlen(m_szPathPrefix);
    size_t uiLen = uiPrefixLen + SimpleGlobUtil_strlen(a_pszFileName) + 1; 
    if (a_bIsDir && (m_uiFlags & SG_GLOB_MARK) == SG_GLOB_MARK) {
        ++uiLen;    // need space for the backslash
    }
    if (!GrowStringBuffer(m_uiBufferLen + uiLen)) {
        return SG_ERR_MEMORY;
    }

    // add this entry. m_uiBufferLen is offset from beginning of buffer.
    m_rgpArgs[m_nArgsLen++] = (char*)m_uiBufferLen;
    SimpleGlobUtil_strcpy_s(m_pBuffer + m_uiBufferLen,
        m_uiBufferSize - m_uiBufferLen, m_szPathPrefix);
    SimpleGlobUtil_strcpy_s(m_pBuffer + m_uiBufferLen + uiPrefixLen,
        m_uiBufferSize - m_uiBufferLen - uiPrefixLen, a_pszFileName);
    m_uiBufferLen += uiLen;

    // add the directory slash if desired
    if (a_bIsDir && (m_uiFlags & SG_GLOB_MARK) == SG_GLOB_MARK) {
        const static char szDirSlash[] = { SG_PATH_CHAR, 0 };
        SimpleGlobUtil_strcpy_s(m_pBuffer + m_uiBufferLen - 2,
            m_uiBufferSize - (m_uiBufferLen - 2), szDirSlash);
    }

    return SG_SUCCESS;
}

void
SetArgvArrayType(
    ARG_ARRAY_TYPE  a_nNewType
    )
{
int n;
    if (m_nArgArrayType == a_nNewType) return;
    if (a_nNewType == POINTERS) {
        SG_ASSERT(m_nArgArrayType == OFFSETS);
        for (n = 0; n < m_nArgsLen; ++n) {
            m_rgpArgs[n] = (m_rgpArgs[n] == (char*)-1) ?
                NULL : m_pBuffer + (size_t) m_rgpArgs[n];
        }
    }
    else {
        SG_ASSERT(a_nNewType == OFFSETS);
        SG_ASSERT(m_nArgArrayType == POINTERS);
        for (n = 0; n < m_nArgsLen; ++n) {
            m_rgpArgs[n] = (m_rgpArgs[n] == NULL) ?
                (char*) -1 : (char*) (m_rgpArgs[n] - m_pBuffer);
        }
    }
    m_nArgArrayType = a_nNewType;
}

int
GrowArgvArray(
    int a_nNewLen
    )
{
    if (a_nNewLen >= m_nArgsSize) {
        static const int SG_ARGV_INITIAL_SIZE = 32;
        int nNewSize = (m_nArgsSize > 0) ? 
            m_nArgsSize * 2 : SG_ARGV_INITIAL_SIZE;
        while (a_nNewLen >= nNewSize) {
            nNewSize *= 2;
        }
        void * pNewBuffer = realloc(m_rgpArgs, nNewSize * sizeof(char*));
        if (!pNewBuffer) return 0;
        m_nArgsSize = nNewSize;
        m_rgpArgs = (char**) pNewBuffer;
    }
    return 1;
}

int
GrowStringBuffer(
    size_t a_uiMinSize
    )
{
    if (a_uiMinSize >= m_uiBufferSize) {
        static const int SG_BUFFER_INITIAL_SIZE = 1024;
        size_t uiNewSize = (m_uiBufferSize > 0) ? 
            m_uiBufferSize * 2 : SG_BUFFER_INITIAL_SIZE;
        while (a_uiMinSize >= uiNewSize) {
            uiNewSize *= 2;
        }
        void * pNewBuffer = realloc(m_pBuffer, uiNewSize * sizeof(char));
        if (!pNewBuffer) return 0;
        m_uiBufferSize = uiNewSize;
        m_pBuffer = (char*) pNewBuffer;
    }
    return 1;
}

static int
fileSortCompare(
    const void *a1,
    const void *a2
    )
{
    const char * s1 = *(const char **)a1;
    const char * s2 = *(const char **)a2;
    if (s1 && s2) {
        return SimpleGlobUtil_strcasecmp(s1, s2);
    }
    // NULL sorts first
    return s1 == s2 ? 0 : (s1 ? 1 : -1);
}

