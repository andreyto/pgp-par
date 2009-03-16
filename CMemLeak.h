// CMemLeak.c and CMemLeak.h are taken from the public domain.  If the
// build flag DEBUG_MEMORY_LEAKS is set, then malloc is redefined,
// to assist in tracking down memory leaks.  Using Purify or Valgrind
// is better, though.
#ifndef CMEMLEAK_H
#define CMEMLEAK_H

#include <stdlib.h>
#include <string.h>

// Used for tracking allocations 
extern void* XWBMalloc(unsigned int iSize, const char* iFile, const unsigned int iLine);
extern void* XWBCalloc(unsigned int iNum, unsigned int iSize, const char* iFile,
    const unsigned int iLine);
extern char* XWBStrDup(const char* iOrig, const char* iFile, const unsigned int iLine);

// Used for tracking reallocations 
extern void* XWBRealloc(void* iPrev, unsigned int iSize, const char* iFile, const unsigned int iLine);

// Used for tracking deallocations 
extern void XWBFree(void* iPtr, const char* iDesc, const char* iFile, const unsigned int iLine);

// Used for reporting 
extern void XWBReport(const char* iTag);
extern void XWBReportFinal(void);

// Used for detecting FMW 
extern void XWBNoFree(void);
extern void XWBPreallocate(const int iInitialAllocations);

//#define DEBUG_MEMORY_LEAKS

// Change this ifdef, in order to redefine malloc(etc.) and track memory leaks:
#ifdef DEBUG_MEMORY_LEAKS
#define malloc(x) XWBMalloc((x), __FILE__, __LINE__)
#define realloc(x,size) XWBRealloc(x,(size),__FILE__,__LINE__)
#define free(x)   XWBFree(x, #x, __FILE__, __LINE__)
#define strdup(x) XWBStrDup(x, __FILE__, __LINE__)
#define calloc(num,size) XWBCalloc((num), (size), __FILE__, __LINE__)
#endif

#endif // CMEMLEAK_H

