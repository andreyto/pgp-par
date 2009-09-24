//Copyright 2007, The Regents of the University of California
//All Rights Reserved
//
//Permission to use, copy, modify and distribute any part of this 
//program for educational, research and non-profit purposes, without fee, 
//and without a written agreement is hereby granted, provided that the 
//above copyright notice, this paragraph and the following three paragraphs 
//appear in all copies.
//
//Those desiring to incorporate this work into commercial 
//products or use for commercial purposes should contact the Technology 
//Transfer & Intellectual Property Services, University of California, 
//San Diego, 9500 Gilman Drive, Mail Code 0910, La Jolla, CA 92093-0910, 
//Ph: (858) 534-5815, FAX: (858) 534-7345, E-MAIL:invent@ucsd.edu.
//
//IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY 
//FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, 
//INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN 
//IF THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY 
//OF SUCH DAMAGE.
//
//THE SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE UNIVERSITY 
//OF CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, 
//ENHANCEMENTS, OR MODIFICATIONS.  THE UNIVERSITY OF CALIFORNIA MAKES NO 
//REPRESENTATIONS AND EXTENDS NO WARRANTIES OF ANY KIND, EITHER IMPLIED OR 
//EXPRESS, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
//MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE USE OF 
//THE SOFTWARE WILL NOT INFRINGE ANY PATENT, TRADEMARK OR OTHER RIGHTS.

#ifndef UTILS_H
#define UTILS_H

#include <memory.h>
#include <string.h>
#include <stdio.h>

#ifndef _WIN32
// linux lacks these:
#define max(X,Y) (((X)>(Y)) ? (X) : (Y))
#define min(X,Y) (((X)>(Y)) ? (Y) : (X))
#endif

#ifdef _WIN32
#define SEPARATOR '\\'
#define SEPARATOR_STRING "\\"
#else
#define SEPARATOR '/'
#define SEPARATOR_STRING "/"
#endif

// We don't like compiler warnings.  Therefore, we cast all our 
// qsort comparison callbacks, in order to avoid this:
// Warning: "passing arg 4 of `qsort' from incompatible pointer type"
typedef int (*QSortCompare)(const void*, const void*);

// For tokenizing strings:
#define WHITESPACE " \r\n\t"

// It seems that tolower is not defined on OSX, so we define our own:
#define ConvertToLower(c) (((c)<'A' || (c)>'Z') ? (c) : ((c)-'A' + 'a'))

#define DIGEST_TYPE_UNKNOWN 0
#define DIGEST_TYPE_TRYPSIN 1

#define MIN_VALID_PEPTIDE_LENGTH 6

// For debugging tag generation.  (Comment it out to disable)
//#define DEBUG_TAG_GENERATION

// PaH!  sizeof() is returning 4 for sizeof(long long), which is INCORRECT.
#define LONG_LONG_SIZE 8

// There are many places where we index into an array by an amino acid's index (A == 0, C == 2, and so on).
// These arrays waste a little space because there is no entry for B, J, O, U, X, Z.  But they are
// very time-efficient.
#define AMINO_ACIDS 26

#define SCORE_TOP_MATCH_VERBOSELY 1
#define BUFFER_SIZE 1024
#define MAX_FILENAME_LEN 1024
#define MAX_MATCHES 10

// We divide the full m/z range into PRM BINS.  The width of each PRM bin
// is 0.1Da, or BIN_SIZE / DALTON.  Masses are stored as ints, where one dalton equals a mass 
// of 1000.  Mass bins are 0.1Da (10 units) in width.  There
// are still some places where bin size is hard-coded.
#define PRM_BIN_SIZE 100

// Number of padding entries added to the end of the PRMScores array.  Useful if we want the score for a PRM that's
// outside our mass range, but just barely.
#define PRM_ARRAY_SLACK 10

// Maximum allowed length for a line in line-based data files.  (If the line is longer
// than this, we report an error)
#define MAX_LINE_LENGTH 2048

#define TAG_EDGE_SCORE_MULTIPLIER 20

#define FAST_ROUND(Float, Int)\
{\
    Int = (int)((Float) + 0.5);\
}

#define ROUND_MASS(FloatMass, IntMass)\
{\
    (IntMass) = (int)(FloatMass * 1000 + 0.5);\
}

#define ROUND_MASS_TO_DELTA_BIN(x, bin) \
{\
(bin) = (int)(((x) + 200000) / 1000.0 + 0.5);\
}

#define MAX_FILENAME_LEN 1024


// Instrument type LTQ is the default.  
#define INSTRUMENT_TYPE_LTQ 0
// QTOF: Fragmentation properties are different (b series is weaker, y series is stronger).
// Parent masses are quite accurate, so parent mass correction is NOT performed
#define INSTRUMENT_TYPE_QTOF 1
// FT hybrid: Parent masses are extremely accurate, so parent mass correction is NOT performed.
// The fragment masses can still be a bit inaccurate however.
#define INSTRUMENT_TYPE_FT_HYBRID 2

#define RUN_MODE_DEFAULT 0
#define RUN_MODE_TAGS_ONLY 1
#define RUN_MODE_MUTATION 2
#define RUN_MODE_BLIND 4
#define RUN_MODE_VERBOSE 8
#define RUN_MODE_PMC_ONLY 16
#define RUN_MODE_PREP_MS2DB 32
#define RUN_MODE_BLINDTAG 64

#define PMC_FEATURE_RAW 1
#define PMC_FEATURE_AVG_RATIO 2
#define PMC_FEATURE_AVG_DIFF 4

typedef enum DatabaseType
{
    evDBTypeTrie = 0,
    evDBTypeMS2DB,
    evDBTypeSpliceDB
} DatabaseType;

typedef struct DatabaseFile
{
    char FileName[MAX_FILENAME_LEN + 1];
    char IndexFileName[MAX_FILENAME_LEN + 1];
    int Type;
    struct DatabaseFile* Next;
    FILE* DBFile;
    FILE* IndexFile;
} DatabaseFile;

typedef struct StringNode
{
    struct StringNode* Next;
    char* String;
} StringNode;



// Global options.  (Set on command-line or in config-file)
typedef struct Options
{
    // RunMode is a set of flags describing which overall code path to take.
    int RunMode; 

    // maximum number of post-translational mods to allow in a match
    int MaxPTMods; 

    // maximum allowed mass error for prefix/suffix peptides
    int Epsilon; 
    int PeakPPM;

    // maximum allowed mass error for prefix/suffix masses
    int FlankingMassEpsilon;

    // return at most this many matches in a search
    int MaxMatches; 

    // -v provides extended debugging info
    int VerboseFlag; 

    // amino acid input-file
    char AminoFileName[MAX_FILENAME_LEN];

    // -o output file (if not set, print matches to stdout)
    char FinalOutputFileName[MAX_FILENAME_LEN]; 
    char OutputFileName[MAX_FILENAME_LEN];  // Intermediate output, before p-value computation
    char ErrorFileName[MAX_FILENAME_LEN];
    int ErrorCount;
    int WarningCount;
    DatabaseFile* FirstDatabase;
    DatabaseFile* LastDatabase;

    // -m file listing legal post-translational modifications
    char PTModFileName[MAX_FILENAME_LEN]; 

    // -i input file name
    char InputFileName[MAX_FILENAME_LEN];
    char ResourceDir[MAX_FILENAME_LEN];

    // either stdout, or opened OutputFileName)
    FILE* OutputFile; 
    // either stderr, or opened ErrorFileName:
    FILE* ErrorFile;

    // -t requests unit tests
    int TestingFlag; 

    // if true, we remember *all* the occurrences of matched peptides.
    int ReportAllMatches; 

    // How far we're allowed to tweak the parent mass of the spectrum.  (Parent masses are often off
    // by one or two amu)
    int ParentMassEpsilon;
    int ParentMassPPM;

    struct Peptide* TruePeptide;

    char MandatoryModName[256];

    int MandatoryModIndex;

    // How many matches to report.  Defaults to 5.
    int ReportMatchCount;

    // How many matches to store for detailed scoring.  Defaults to 100.
    int StoreMatchCount;

    // How many tags shall we generate, and how long shall they be?
    int GenerateTagCount;
    int GenerateTagLength;

    // Nonzero if this is this a trypsin digest, or some other type of specific digest.
    // If DigestType != 0, then we can give a penalty for missed cleavages, and a bonus for matching termini
    int DigestType; 

    // Linked list of SpectrumNodes:
    struct SpectrumNode* FirstSpectrum;
    struct SpectrumNode* LastSpectrum;

    // Linked list of InputFiles:
    struct InputFileNode* FirstFile;
    struct InputFileNode* LastFile;

    int SpectrumCount;
    int DynamicRangeMin;
    int DynamicRangeMax;
    int TaglessSearchFlag;

    // If PhosphorylationFlag, then attempt to interpret phosphorylated peptides.  This has implications
    // for tag-generation, as well as candidate scoring.
    int PhosphorylationFlag;

    int TagPTMMode; // 0 is free, 1 is forbidden, and 2 is penalized

    int MultiChargeMode; // if 1, try multiple parent charge states.

    int TrieBlockSize;
    int InstrumentType;
    // Options for unrestrictive PTM search:
    // DeltaBinCount is the number of mass bins in the range [MinPTMDelta, MaxPTMDelta], 
    // by default it equals 400 * 10 = 4000. 
    int MinPTMDelta;
    int MaxPTMDelta;
    int DeltaBinCount;
    int DeltasPerAA; // == max(DeltaBinCount*2, 512)
    // If TRUE, then use PepNovo for tag generation (assumed to live in working directory!)
    int ExternalTagger; 

    // Options for producing an .ms2db file from .gff files:
    StringNode* FirstGFFFileName;
    StringNode* LastGFFFileName;
    char GenomeFileName[MAX_FILENAME_LEN + 1];
    char ChromosomeName[256 + 1];

    // If XMLStrictFlag is set, then we'll complain about any unexpected
    // tags or attributes.  This is useful when debugging .ms2db file
    // generation.  In production, this flag won't generally be set, 
    // because it is officially Allowable to add new tags and 
    // attributes to an .ms2db file.
    int XMLStrictFlag;

    // if RequireTermini is 1 or 2, then we accept only semi-tryptic or fully-tryptic matches.
    int RequireTermini;

    int NewScoring; //temporary flag while we work on a new code path for scoring
} Options;

extern Options* GlobalOptions;

int CopyBufferLine(char* Source, int BufferPos, int BufferEnd, char* LineBuffer, int StripWhitespace);
int CompareFloats(const float* a, const float* b);
int CompareInts(const int* a, const int* b);
int CompareStrings(const char* StringA, const char* StringB);
char TranslateCodon(char* DNA);
void WriteReverseComplement(char* Source, char* Destination);
void ReverseString(char* String);

#ifdef __ppc__
size_t ReadBinary(void* Buffer, size_t ItemSize, size_t ItemCount, FILE* stream);
size_t WriteBinary(void* Buffer, size_t ItemSize, size_t ItemCount, FILE* stream);
#define BYTEORDER_BIG_ENDIAN
#else
#define ReadBinary fread
#define WriteBinary fwrite
#define BYTEORDER_LITTLE_ENDIAN
#endif

void AssertionFailed(char* Assertion, char* FileName, int LineNumber);

#define INSPECT_ASSERT(expr) \
    if (!(expr)) \
    AssertionFailed(#expr, __FILE__, __LINE__)

#define SafeFree(Pointer)\
    if (Pointer) \
    {\
    free(Pointer);\
    }

// a FileLineParser is called once per line as a callback from ParseFileByLines()
typedef int (*FileLineParser)(int FilePos, int LineNumber, char* LineBuffer, void* ParseData);

void ParseFileByLines(FILE* File, FileLineParser Parser, void* ParseData, int ProcessCommentLines);
float GetMedian(float* Values, int ValueCount);
//#define PMC_USE_SVM
#define MQSCORE_USE_SVM

#define MQ_FEATURE_COUNT 7

#endif //UTILS_H
