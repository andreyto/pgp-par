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

// ParseInput.c is responsible for parsing the Inspect input file.

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <locale.h>
#include <ctype.h>

#include "Trie.h"
#include "Utils.h"
#include "Spectrum.h"
#include "Mods.h"
#include "Score.h"
#include "Tagger.h"
#include "FreeMod.h"
#include "CMemLeak.h"
#include "SVM.h"
#include "BN.h"
#include "Run.h"
#include "SNP.h"
#include "SpliceDB.h"
#include "ChargeState.h"
#include "Scorpion.h"
#include "ParseXML.h"
#include "SpliceScan.h"
#include "Errors.h"
#include "IonScoring.h"

// If the input file specifies a directory full of spectra, we must iterate over the files in that directory.
// That works a bit differently on Windows and on Unix.
#ifdef _WIN32
#include <windows.h>
#include <sys/types.h>
#include <sys/stat.h>
#else
#include <dirent.h>
#include <sys/stat.h>
#endif

// Global variables:
extern Options* GlobalOptions;
extern MSSpectrum* Spectrum;

// Array of spectra to be searched.  We put them into an array so that we can qsort
// them.  (Not crucial, but it's nice to get output in order)
SpectrumNode* g_BigNodeArray = NULL;

extern StringNode* FirstTagCheckNode;
extern StringNode* LastTagCheckNode;

void AddDatabase(DatabaseFile* Database)
{
    if (!Database)
    {
        return;
    }
    if (!GlobalOptions->FirstDatabase)
    {
        GlobalOptions->FirstDatabase = Database;
        GlobalOptions->LastDatabase = Database;
    }
    else
    {
        GlobalOptions->LastDatabase->Next = Database;
        GlobalOptions->LastDatabase = Database;
    }
}

// Parse a FASTA file, and convert it into a .trie file.  This is the same
// thing as PrepDB.py.
void PrepareSecondarySequenceFile(char* FileName)
{
    FILE* FastaFile;
    FILE* DBFile;
    FILE* IndexFile;
    int Dummy = 0;
    char Char;
    int BytesRead;
    int NameLength = 0;
    int ReadingName = 0;
    int TargetDBPos = 0;
    int SourceFilePos = 0;
    char* StarChar = "*";
    char* NullChar = "\0";
    DatabaseFile* Database;
    //
    FastaFile = fopen(FileName, "rb");
    if (!FastaFile)
    {
        REPORT_ERROR_S(8, FileName);
        return;
    }
    DBFile = fopen("AdditionalDB.trie", "wb");
    IndexFile = fopen("AdditionalDB.index", "wb");
    if (!DBFile || !IndexFile)
    {
        printf("Unable to write out processed secondary database!  Skipping.\n");
        return;
    }
    while (1)
    {
        BytesRead = ReadBinary(&Char, sizeof(char), 1, FastaFile);
        if (!BytesRead)
        {
            break;
        }
        if (Char == '>')
        {
            ReadingName = 1;
            if (TargetDBPos)
            {
                WriteBinary(StarChar, sizeof(char), 1, DBFile);
                TargetDBPos++;
            }
            WriteBinary(&SourceFilePos, sizeof(int), 1, IndexFile);
            // Source file pos is a long long; assume we must write another 4 bytes:
            WriteBinary(&Dummy, sizeof(int), 1, IndexFile);
            WriteBinary(&TargetDBPos, sizeof(int), 1, IndexFile);
            NameLength = 0;
            continue;
        }
        if (Char == '\r' || Char == '\n')
        {
            if (ReadingName)
            {
                // Pad the protein name out:
                while (NameLength < 80)
                {
                    WriteBinary(NullChar, sizeof(char), 1, IndexFile);
                    NameLength++;
                }

            }
            ReadingName = 0;
            continue;
        }
        if (ReadingName)
        {
            if (NameLength < 79)
            {
                WriteBinary(&Char, sizeof(char), 1, IndexFile);
                NameLength++;
            }

            continue;
        }
        if (Char == ' ' || Char == '\t' || Char == '*')
        {
            continue;
        }
        WriteBinary(&Char, sizeof(char), 1, DBFile);
        TargetDBPos++;
    }
    fclose(DBFile);
    fclose(IndexFile);
    Database = (DatabaseFile*)calloc(1, sizeof(DatabaseFile));
    strcpy(Database->FileName, "AdditionalDB.trie");
    Database->Type = evDBTypeTrie;
    AddDatabase(Database);
}

// Free the linked list of TagCheckNodes
void FreeTagCheckNodes()
{
    StringNode* Node;
    StringNode* Prev = NULL;
    for (Node = FirstTagCheckNode; Node; Node = Node->Next)
    {
        if (Prev)
        {
            SafeFree(Prev->String);
            SafeFree(Prev);
        }
        Prev = Node;
    }
    if (Prev)
    {
        SafeFree(Prev->String);
        SafeFree(Prev);
    }
}

void FreeInputFileNodes()
{
    InputFileNode* Node;
    InputFileNode* Prev = NULL;

    // At this point, freeing doesn't matter much, since this function's called
    // just before program exit.
    for (Node = GlobalOptions->FirstFile; Node; Node = Node->Next)
    {
        SafeFree(Prev);
        Prev = Node;
    }
    SafeFree(Prev);
}

char* ProteaseNames[] = {"none", "trypsin", "chymotrypsin", "lysc", "aspn", "gluc"};

typedef struct ParseSpectraFromFileInfo
{
    int FirstScan;
    int LastScan;
    int ScanNumber;
    InputFileNode* InputFile;
} ParseSpectraFromFileInfo;

// Callback for ParseSpectraFromMS2File: Handle one line of an .ms2 spectrum file.
int ParseSpectraFromMS2FileCallback(int LineNumber, int FilePos, char* LineBuffer, void* UserData)
{
    ParseSpectraFromFileInfo* Info;
    //
    Info = (ParseSpectraFromFileInfo*)UserData;
    if (LineBuffer[0] == ':')
    {
        Info->InputFile->Format = SPECTRUM_FORMAT_MS2_COLONS;
        Info->ScanNumber = atoi(LineBuffer + 1);
        if (Info->ScanNumber >= Info->FirstScan && (Info->LastScan < 0 || Info->ScanNumber < Info->LastScan))
        {
            AddSpectrumToList(Info->InputFile, FilePos, Info->ScanNumber);
        }
    }
    if (LineBuffer[0] == 'S' && (LineBuffer[1] == ' ' || LineBuffer[1] == '\t'))
    {
		// Start of a spectrum:
        Info->ScanNumber = atoi(LineBuffer + 1);
		Info->InputFile->Format = SPECTRUM_FORMAT_MS2;
        if (Info->ScanNumber >= Info->FirstScan && (Info->LastScan < 0 || Info->ScanNumber < Info->LastScan))
        {
            AddSpectrumToList(Info->InputFile, FilePos, Info->ScanNumber);
        }
    }
    if (LineBuffer[0] == 'Z' && (LineBuffer[1] == ' ' || LineBuffer[1] == '\t'))
    {
        //// This is the start of a spectrum:
        //Info->InputFile->Format = SPECTRUM_FORMAT_MS2;
        //if (Info->ScanNumber >= Info->FirstScan && (Info->LastScan < 0 || Info->ScanNumber < Info->LastScan))
        //{
        //    AddSpectrumToList(Info->InputFile, FilePos, Info->ScanNumber);
        //}
    }
    return 1;
}

// Iterate over lines in the MS2 file.
// If you reach a line of the form "Z [charge] [mass]", that's the beginning of a spectrum record.
// The spectrum code knows that it should process the first Z it sees, skip any others, then process peaks
// until it sees something that is not a peak.
void ParseSpectraFromMS2File(char* FileName, InputFileNode* InputFile, int FirstScan, int LastScan)
{
    FILE* MS2File;
    ParseSpectraFromFileInfo Info;
    //
    MS2File = fopen(FileName, "rb");
    if (!MS2File)
    {
        REPORT_ERROR_S(8, FileName);
        return;
    }
    printf("Count spectra from '%s'...\n", FileName);
    Info.FirstScan = FirstScan;
    Info.LastScan = LastScan;
    Info.InputFile = InputFile;
    ParseFileByLines(MS2File, ParseSpectraFromMS2FileCallback, &Info, 0);
    fclose(MS2File);
}

// Callback for ParseSpectraFromMGFFile: Handle one line of an .mgf spectrum file.
int ParseSpectraFromMGFFileCallback(int LineNumber, int FilePos, char* LineBuffer, void* UserData)
{
    ParseSpectraFromFileInfo* Info;
    //
    Info = (ParseSpectraFromFileInfo*)UserData;
    if (!strncmp(LineBuffer, "BEGIN IONS", 10))
    {
        if (Info->ScanNumber >= Info->FirstScan && (Info->LastScan < 0 || Info->ScanNumber <= Info->LastScan))
        {
            AddSpectrumToList(Info->InputFile, FilePos, Info->ScanNumber);
        }
        Info->ScanNumber++;
    }
    return 1;
}

// Callback for ParseSpectraFromPKLFile: Handle one line of a .pkl spectrum file.
// Assume: If there are three numbers, then this line is the header of a spectrum.
int ParseSpectraFromPKLFileCallback(int LineNumber, int FilePos, char* LineBuffer, void* UserData)
{
    ParseSpectraFromFileInfo* Info;
    char* StrA;
    char* StrB;
    char* StrC;
    float FloatValue;
    int IntValue;
    //
    Info = (ParseSpectraFromFileInfo*)UserData;

    // First, check to see that there are three fields on this line of the file:
    StrA = strtok(LineBuffer, WHITESPACE);
    StrB = strtok(NULL, WHITESPACE);
    if (!StrB)
    {
        return 1;
    }
    StrC = strtok(NULL, WHITESPACE);
    if (!StrC)
    {
        return 1;
    }
    // Now, check to see that the three fields are valid numbers:
    FloatValue = (float)atof(StrA);
    if (!FloatValue)
    {
        return 1;
    }
    FloatValue = (float)atof(StrB);
    if (!FloatValue)
    {
        return 1;
    }
    IntValue = atoi(StrC);
    if (Info->ScanNumber >= Info->FirstScan && (Info->LastScan < 0 || Info->ScanNumber <= Info->LastScan))
    {
        AddSpectrumToList(Info->InputFile, FilePos, Info->ScanNumber);
    }
    Info->ScanNumber++;
    return 1;
}

void ParseSpectraFromPKLFile(char* FileName, InputFileNode* InputFile,
    int FirstScan, int LastScan)
{
    FILE* SpectrumFile;
    ParseSpectraFromFileInfo Info;
    //
    SpectrumFile = fopen(FileName, "rb");
    if (!SpectrumFile)
    {
        REPORT_ERROR_S(8, FileName);
        return;
    }
    printf("Count spectra from '%s'...\n", FileName);
    Info.FirstScan = FirstScan;
    Info.LastScan = LastScan;
    Info.InputFile = InputFile;
    Info.ScanNumber = 0;
    ParseFileByLines(SpectrumFile, ParseSpectraFromPKLFileCallback, &Info, 0);
    fclose(SpectrumFile);
}

// Iterate over lines in the MGF file.
// If you reach a line of the form "BEGIN IONS", that starts a spectrum record.
// The spectrum parse code knows that it should process a CHARGE line, a PEPMASS line, then peak lines.
// Note: Scan numbers from the MGF file are *ignored*!  The first scan we see is number 0,
// the next is number 1, etc.
void ParseSpectraFromMGFFile(char* FileName, InputFileNode* InputFile,
    int FirstScan, int LastScan)
{
    FILE* MGFFile;
    ParseSpectraFromFileInfo Info;
    //
    MGFFile = fopen(FileName, "rb");
    if (!MGFFile)
    {
        REPORT_ERROR_S(8, FileName);
        return;
    }
    printf("Count spectra from '%s'...\n", FileName);
    Info.FirstScan = FirstScan;
    Info.LastScan = LastScan;
    Info.InputFile = InputFile;
    Info.ScanNumber = 0;
    ParseFileByLines(MGFFile, ParseSpectraFromMGFFileCallback, &Info, 0);
    fclose(MGFFile);
}

void AddSpectrumNodesForFile(char* FileName, InputFileNode* InputFile,
    int FirstScan, int LastScan)
{
    int Format;
    //
    // Based upon the file extension, decide whether and how to parse the input file
    Format = GuessSpectrumFormatFromExtension(FileName);
    InputFile->Format = Format;
    switch (Format)
    {
    case SPECTRUM_FORMAT_MS2:
    case SPECTRUM_FORMAT_MS2_COLONS:
        ParseSpectraFromMS2File(FileName, InputFile, FirstScan, LastScan);
        break;
    case SPECTRUM_FORMAT_MZXML:
        ParseSpectraFromMZXML(FileName, InputFile, FirstScan, LastScan);
        break;
    case SPECTRUM_FORMAT_MZDATA:
        ParseSpectraFromMZData(FileName, InputFile, FirstScan, LastScan);
        break;
    case SPECTRUM_FORMAT_MGF:
        ParseSpectraFromMGFFile(FileName, InputFile, FirstScan, LastScan);
        break;
    case SPECTRUM_FORMAT_PKL:
        ParseSpectraFromPKLFile(FileName, InputFile, FirstScan, LastScan);
        break;
    case SPECTRUM_FORMAT_DTA:
        // Let's assume that we can treat it as a .dta file.
        AddSpectrumToList(InputFile, 0, 0);
        break;
    default:
        printf("Not parsing unknown spectrum file format:%s\n", FileName);
        break;
    }
}

// Add a spectrum file to our input list.  If the file contains multiple spectra,
// then we'll create several SpectrumNode instances.  If FirstScan is set and >0,
// skip all scans with scan number < FirstScan.  If LastScan is set and >-1, then skip
// all scans with scan number > LastScan.  (INCLUSIVE ends)
void AddSpectraToList(char* FileName, int FirstScan, int LastScan)
{
    InputFileNode* NewFile;
    //
    NewFile = (InputFileNode*)calloc(1, sizeof(InputFileNode));
    strncpy(NewFile->FileName, FileName, MAX_FILENAME_LEN);
    if (GlobalOptions->LastFile)
    {
        GlobalOptions->LastFile->Next = NewFile;
    }
    else
    {
        GlobalOptions->FirstFile = NewFile;
    }
    GlobalOptions->LastFile = NewFile;
    AddSpectrumNodesForFile(FileName, NewFile, FirstScan, LastScan);
    //strncpy(NewNode->FileName, FileName, MAX_FILENAME_LEN);

}

#ifdef _WIN32
// The WINDOWS way to iterate over a directory:
void ProcessInputCommandSpectra(char* FileName, int FirstScan, int LastScan)
{
    char DirBuffer[1024];
    char StarBuffer[1024];
    char FileNameBuffer[1024];
    struct stat StatBuffer;
    int StatResult;
    int Len;
    int Result;
    int SkipFile;
    if (!FileName || !FileName[0])
    {
        printf("* Error: null filename specified in 'spectra' command\n");
        return;
    }
    StatResult = stat(FileName, &StatBuffer);
    if (StatResult < 0)
    {
        REPORT_ERROR_S(8, FileName);
        //printf("Unable to stat '%s' - skipping.\n", FileName);
        return;
    }
    if (StatBuffer.st_mode & _S_IFDIR)
    {
        HANDLE hFindFile;
        WIN32_FIND_DATA wFileFindData;
        sprintf(DirBuffer, FileName);
        Len = strlen(FileName);
        if (DirBuffer[Len-1] != '\\')
        {
            strcat(DirBuffer, "\\");
        }
        sprintf(StarBuffer, "%s*.*", DirBuffer);
        hFindFile = FindFirstFile(StarBuffer, &wFileFindData);
        while (hFindFile != INVALID_HANDLE_VALUE)
        {
            SkipFile = 0;
            if (wFileFindData.cFileName[0]=='\0')
            {
                SkipFile = 1;
            }
            if (wFileFindData.cFileName[0]=='.' && wFileFindData.cFileName[1]=='\0')
            {
                SkipFile = 1;
            }
            if (wFileFindData.cFileName[0]=='.' && wFileFindData.cFileName[1]=='.' && wFileFindData.cFileName[2]=='\0')
            {
                SkipFile = 1;
            }
            if (!SkipFile)
            {
                sprintf(FileNameBuffer, "%s%s", DirBuffer, wFileFindData.cFileName);
                StatResult = stat(FileNameBuffer, &StatBuffer);
                if (StatBuffer.st_mode & _S_IFREG)
                {
                    //printf("Adding file to list: '%s'\n", FileNameBuffer);
                    AddSpectraToList(FileNameBuffer, FirstScan, LastScan);
                }
            }
            Result = FindNextFile(hFindFile, &wFileFindData);
            if (!Result)
            {
                break;
            }
        }
    }
    else
    {
        AddSpectraToList(FileName, FirstScan, LastScan);
    }
}
#else
// The UNIX way to iterate over a directory:
void ProcessInputCommandSpectra(char* FileName, int FirstScan, int LastScan)
{
    char DirBuffer[1024];
    char FileNameBuffer[1024];
    struct stat StatBuffer;
    DIR *dirp;
    struct dirent *ep;
    int StatResult;
    int Len;
    if (!FileName || !FileName[0])
    {
        printf("* Error: null filename specified in 'spectra' command\n");
        return;
    }

    StatResult = stat(FileName, &StatBuffer);
    if (S_ISDIR(StatBuffer.st_mode))
    {
        sprintf(DirBuffer, FileName);
        Len = strlen(FileName);
        if (DirBuffer[Len-1] != '/')
        {
            strcat(DirBuffer, "/");
        }
        dirp = opendir(DirBuffer);
        while ((ep = readdir(dirp)) != NULL)
        {
            Len = strlen(ep->d_name);
            if (ep->d_name[0]=='\0')
            {
                continue;
            }
            if (ep->d_name[0]=='.' && ep->d_name[1]=='\0')
            {
                continue;
            }
            if (ep->d_name[0]=='.' && ep->d_name[1]=='.' && ep->d_name[2]=='\0')
            {
                continue;
            }
            sprintf(FileNameBuffer, "%s%s", DirBuffer, ep->d_name);
            StatResult = stat(FileNameBuffer, &StatBuffer);
            if (StatBuffer.st_mode & S_IFREG)
            {
                AddSpectraToList(FileNameBuffer, FirstScan, LastScan);
            }
        }
        closedir(dirp);
    }
    else
    {
        AddSpectraToList(FileName, FirstScan, LastScan);
    }
}
#endif

#define REJECT_NULL_VALUE(name)\
{\
if (!CommandValue || !CommandValue[0]) \
{\
    printf("* Error: Null value for '(name)'\n");\
    return;\
}\
}

typedef int (*InputParameterParser)(char* CommandValue);

#define INPUT_VALUE_TYPE_NONE 0
#define INPUT_VALUE_TYPE_INT 1
#define INPUT_VALUE_TYPE_STRING 2
typedef struct InputParameter
{
    char* Name;
    InputParameterParser ParseFunction;
    int ValueType;
} InputParameter;

int ParseInputTagCheck(char* Value)
{
    StringNode* StrNode;
    //
    StrNode = (StringNode*)calloc(1, sizeof(StringNode));
    StrNode->String = strdup(Value);
    if (FirstTagCheckNode)
    {
        LastTagCheckNode->Next = StrNode;
    }
    else
    {
        FirstTagCheckNode = StrNode;
    }
    LastTagCheckNode = StrNode;
    return 1;
}

int ParseInputTagsOnly(char* Value)
{
    GlobalOptions->RunMode |= RUN_MODE_TAGS_ONLY;
    return 1;
}

int ParseInputExternalTagger(char* Value)
{
    GlobalOptions->ExternalTagger = 1;
    return 1;
}

int ParseInputSpectra(char* Value)
{
    char* ScanStr;
    int FirstScan;
    int LastScan;
    // Spectrum file:
    // Note: If the file is a directory, we will iterate over all the files in
    // the directory.  We don't recurse into subdirectories.
    //ScanStr = strtok(NULL, ",");
    FirstScan = 0; //default
    LastScan = -1; //default
    ScanStr = strtok(NULL, ",");
    if (ScanStr)
    {
        FirstScan = atoi(ScanStr);
        ScanStr = strtok(NULL, ",");
        if (ScanStr)
        {
            LastScan = atoi(ScanStr);
            // LastScan of -1 means no upper limit...but otherwise, LastScan should
            // not be below FirstScan!
            if (LastScan < FirstScan && LastScan >= 0)
            {
                REPORT_WARNING_II(9, FirstScan, LastScan);
                return 0;
            }
        }
    }
    ProcessInputCommandSpectra(Value, FirstScan, LastScan);
    return 1;
}

int ParseInputInstrument(char* CommandValue)
{
    // Instrument name is LCQ or QTOF or FT-Hybrid.  If QTOF, use a different
    // scoring model, and don't perform parent-mass correction.
    if (!CompareStrings(CommandValue, "ESI-ION-TRAP"))
    {
        GlobalOptions->InstrumentType = INSTRUMENT_TYPE_LTQ;
    }
    else if (!CompareStrings(CommandValue, "QTOF"))
    {
        GlobalOptions->InstrumentType = INSTRUMENT_TYPE_QTOF;
        GlobalOptions->ParentMassPPM = 100;
    }
    else if (!CompareStrings(CommandValue, "FT-HYBRID"))
    {
        GlobalOptions->InstrumentType = INSTRUMENT_TYPE_FT_HYBRID;
        GlobalOptions->ParentMassPPM = 100;
    }
    else
    {
        printf("** Warning: unknown instrument type '%s'\n", CommandValue);
        return 0;
    }
    return 1;
}

int ParseInputProtease(char* CommandValue)
{
    int ProteaseIndex;
    for (ProteaseIndex = 0; ProteaseIndex < sizeof(ProteaseNames)/sizeof(char*); ProteaseIndex++)
    {
        if (!CompareStrings(ProteaseNames[ProteaseIndex], CommandValue))
        {
            GlobalOptions->DigestType = ProteaseIndex;
            return 1;
        }
    }
    printf("* Error: Protease '%s' not understood\n", CommandValue);
    return 0;
}

int GuessDBTypeFromExtension(char* FileName)
{
    char* Extension;
    Extension = FileName + strlen(FileName);
    while (Extension > FileName)
    {
        Extension--;
        if (*Extension == '.')
        {
            if (!CompareStrings(Extension, ".ms2db"))
            {
                return evDBTypeMS2DB;
            }
            if (!CompareStrings(Extension, ".dat"))
            {
                return evDBTypeSpliceDB;
            }
        }
    }
    return evDBTypeTrie; // default guess
}

int ParseInputDB(char* CommandValue)
{
    DatabaseFile* Database;
    FILE* TempFile;
    char DBFileName[MAX_FILENAME_LEN + 1];

    /// If what we got looks like a complete path, then keep the path:
    if (CommandValue[0]=='/' || CommandValue[0]=='.' || CommandValue[1]==':')
    {
        strncpy(DBFileName, CommandValue, MAX_FILENAME_LEN);
    }
    else
    {
        // Otherwise, go to $resourcedir\database
        sprintf(DBFileName, "%sDatabase%c%s", GlobalOptions->ResourceDir, SEPARATOR, CommandValue);
    }
    //To-ju: Putting protein databases in a subfolder of the inspect executable will cause unnatural coupling.

    //strncpy(DBFileName, CommandValue, MAX_FILENAME_LEN);


    TempFile = fopen(DBFileName, "rb");
    if (!TempFile)
    {
        REPORT_ERROR_S(8, DBFileName);
        return 0;
    }
    else
    {
        fclose(TempFile);
    }
    Database = (DatabaseFile*)calloc(1, sizeof(DatabaseFile));
    strcpy(Database->FileName, DBFileName);
    Database->Type = GuessDBTypeFromExtension(Database->FileName);
    AddDatabase(Database);
    return 1;
}

int ParseInputPMTolerance(char* CommandValue)
{
    GlobalOptions->ParentMassEpsilon = rint(atof(CommandValue) * MASS_SCALE);
    return 1;
}

int ParseInputReportMatches(char* CommandValue)
{
    GlobalOptions->ReportMatchCount = atoi(CommandValue);
    GlobalOptions->ReportMatchCount = min(100, max(1, GlobalOptions->ReportMatchCount));
    return 1;
}

int ParseInputRequireTermini(char* CommandValue)
{
    int RequireTerminiCount;
    //
    RequireTerminiCount = atoi(CommandValue);
    if (RequireTerminiCount < 0 || RequireTerminiCount > 2)
    {
        REPORT_ERROR_I(47, RequireTerminiCount);
    }

    GlobalOptions->RequireTermini = RequireTerminiCount;
    return 1;
}

int ParseInputRequiredMod(char* CommandValue)
{
    strncpy(GlobalOptions->MandatoryModName, CommandValue, 256);
    return 1;
}
int ParseInputTagCount(char* CommandValue)
{
    GlobalOptions->GenerateTagCount = atoi(CommandValue);
    return 1;
}

int ParseInputTagLength(char* CommandValue)
{
    GlobalOptions->GenerateTagLength = atoi(CommandValue);
    if (GlobalOptions->GenerateTagLength <= 0 || GlobalOptions->GenerateTagLength > 6)
    {
        REPORT_ERROR_I(38, GlobalOptions->GenerateTagLength);
        GlobalOptions->GenerateTagLength = DEFAULT_TAG_LENGTH;
        return 0;
    }
    return 1;
}

int ParseInputIonTolerance(char* CommandValue)
{
    //if (!CompareStrings(CommandName, "IonTolerance") || !CompareStrings(CommandName, "Ion_Tolerance"))
    GlobalOptions->Epsilon = rint(atof(CommandValue) * MASS_SCALE);
    return 1;
}

int ParseInputMods(char* CommandValue)
{
    GlobalOptions->MaxPTMods = atoi(CommandValue);
    return 1;
}

int ParseInputFreeMods(char* CommandValue)
{
    char Path[MAX_FILENAME_LEN];

    GlobalOptions->MaxPTMods = atoi(CommandValue);
    // "freemods,1" or "freemods,2" allows mutations plus a rich PTM set.
    if (GlobalOptions->MaxPTMods && !(GlobalOptions->RunMode & RUN_MODE_BLIND))
    {
        GlobalOptions->RunMode |= RUN_MODE_MUTATION;
        GlobalOptions->PhosphorylationFlag = 1;
        sprintf(Path, "%s%s", GlobalOptions->ResourceDir, FILENAME_MASS_DELTAS);
        LoadMassDeltas(Path,  GlobalOptions->RunMode & RUN_MODE_MUTATION);
    }
    return 1;
}

int ParseInputPMCOnly(char* CommandValue)
{
    if (atoi(CommandValue))
    {
        GlobalOptions->RunMode |= RUN_MODE_PMC_ONLY;
    }
    return 1;
}
int ParseInputTagless(char* CommandValue)
{
    GlobalOptions->TaglessSearchFlag = atoi(CommandValue);
    return 1;
}
int ParseInputBlind(char* CommandValue)
{
    if (atoi(CommandValue))
    {
        GlobalOptions->RunMode |= RUN_MODE_BLIND;
    }
    return 1;
}
int ParseInputBlindTagging(char* CommandValue)
{
    if (atoi(CommandValue))
    {
        GlobalOptions->RunMode |= RUN_MODE_BLINDTAG;
    }
    return 1;
}

// Maximum size, in daltons, of PTMs to consider.  (Blind search only)
int ParseInputMaxPTMSize(char* CommandValue)
{
    GlobalOptions->MaxPTMDelta = atoi(CommandValue);
    if (GlobalOptions->MaxPTMDelta < 1 || GlobalOptions->MaxPTMDelta >= 2000)
    {
        printf("** Error: Invalid maxptmsize '%s' - please select a value between 10 and 2000Da\n", CommandValue);
        GlobalOptions->MaxPTMDelta = 200;
        return 0;
    }
    GlobalOptions->DeltaBinCount = (GlobalOptions->MaxPTMDelta - GlobalOptions->MinPTMDelta) * 10 + 1;
    GlobalOptions->DeltasPerAA = max(512, GlobalOptions->DeltaBinCount * 2);
    return 1;
}

// If multicharge flag is set, then ALWAYS try charge correction on spectra.  (Otherwise, do it only
// if the source file provides no charge, or says the charge is zero)
int ParseInputMultiCharge(char* CommandValue)
{
    GlobalOptions->MultiChargeMode = atoi(CommandValue);
    return 1;
}
int ParseInputXMLStrict(char* CommandValue)
{
    GlobalOptions->XMLStrictFlag = atoi(CommandValue);
    return 1;
}

int ParseInputPRMModel(char* CommandValue)
{
    char* StrCharge = NULL;
    char* FileName;
    int Charge;
    //
    StrCharge = CommandValue;
    FileName = strtok(CommandValue, ",");
    Charge = atoi(StrCharge);
    if (Charge < 2 || Charge > 3)
    {
        REPORT_ERROR(46);
        return 0;
    }
    return ReplacePRMScoringModel(Charge, FileName);
}

int ParseInputTAGModel(char* CommandValue)
{
    char* StrCharge = NULL;
    char* FileName;
    int Charge;
    //
    StrCharge = CommandValue;
    FileName = strtok(CommandValue, ",");
    Charge = atoi(StrCharge);
    if (Charge < 2 || Charge > 3)
    {
        REPORT_ERROR(46);
        return 0;
    }
    return ReplaceTAGScoringModel(Charge, FileName);
}

int ParseInputMod(char* CommandValue)
{
    int ModFlags;
    char* StrMass = NULL;
    char* StrAminos = NULL;
    char* StrType = NULL;
    char* StrName = NULL;
    float MassDelta;
    char* Amino;
    int AminoIndex;
    int AminoFoundFlag;
    int Bin;
    int ModIndex;
    char ModNameBuffer[64];
    //
    if (!MassDeltas)
    {
        LoadMassDeltas(NULL, 0);
    }
    if (AllPTModCount == MAX_PT_MODTYPE)
    {
        // Too many!
        REPORT_ERROR_S(35, CommandValue);
        return 0;
    }
    ModFlags = DELTA_FLAG_VALID;
    StrMass = CommandValue;
    StrAminos = strtok(NULL, ","); // required, can be "*" for no specificity
    if (!StrAminos || !*StrAminos)
    {
        printf("* Error: Modification must have amino acids specified!\n");
        return 0;
    }
    StrType = strtok(NULL, ","); // optional: fix/opt/cterminal/nterminal
    if (StrType)
    {
        StrName = strtok(NULL, ","); // optional: name
    }
    if (!StrMass || !StrAminos || !StrAminos[0])
    {
        printf("** Error: invalid modification in input file.  Skipping!\n");
        return 0;
    }
    if (strstr(StrAminos, "*"))
    {
        StrAminos = "ACDEFGHIKLMNPQRSTVWY";
    }
    MassDelta = (float)atof(StrMass);
    if (MassDelta == 0 || MassDelta > 1000 || MassDelta < -200)
    {
        printf("** Error: invalid modification in input file; mass is %.2f.  Skipping!\n", MassDelta);
        return 0;
    }
    // Default modification type is OPTIONAL.
    if (!StrType)
    {
        StrType = "opt";
    }
    // Default name is the mass (rounded to integer, with sign indicated)
    if (!StrName)
    {
        sprintf(ModNameBuffer, "%+d", (int)rint(MassDelta));
        StrName = ModNameBuffer;
    }
    // If it's a fixed modification, then adjust the amino acid mass:
    if (!CompareStrings(StrType, "fix") || !CompareStrings(StrType, "fixed"))
    {
        for (Amino = StrAminos; *Amino; Amino++)
        {
            AminoIndex = *Amino - 'A';
            if (AminoIndex >= 0 && AminoIndex < TRIE_CHILD_COUNT)
            {
                PeptideMass[Amino[0]] += rint(MassDelta * MASS_SCALE);
                // We haven't yet called PopulateJumpingHash(), so that's all we need to do
            }
        }
        return 1;
    }
    else if (!CompareStrings(StrType, "cterminal") || !CompareStrings(StrType, "c-terminal"))
    {
        ModFlags |= DELTA_FLAG_C_TERMINAL;
    }
    else if (!CompareStrings(StrType, "nterminal") || !CompareStrings(StrType, "n-terminal"))
    {
        ModFlags |= DELTA_FLAG_N_TERMINAL;
    }
    else if (!CompareStrings(StrType, "opt") || !CompareStrings(StrType, "optional"))
    {
        ; // pass
    }
    else
    {
        REPORT_ERROR_S(36, StrType);
    }

    if (!CompareStrings(StrName, "phosphorylation"))
    {
        g_PhosphorylationMod = AllPTModCount;
        GlobalOptions->PhosphorylationFlag = 1;
        ModFlags |= DELTA_FLAG_PHOSPHORYLATION;
    }
    AllKnownPTMods[AllPTModCount].Flags = ModFlags;
    strncpy(AllKnownPTMods[AllPTModCount].Name, StrName, 40);
    // Add another modification to each amino acid's mod-array:
    AminoFoundFlag = 0;
    for (Amino = StrAminos; *Amino; Amino++)
    {
        AminoIndex = *Amino - 'A';
        if (AminoIndex >= 0 && AminoIndex < TRIE_CHILD_COUNT)
        {
            AminoFoundFlag = 1;
            AllKnownPTMods[AllPTModCount].Allowed[AminoIndex] = 1;
            // Add to the first still-available slot:
            for (ModIndex = 0; ModIndex < GlobalOptions->DeltasPerAA; ModIndex++)
            {
                if (!MassDeltas[AminoIndex][ModIndex].Flags)
                {
                    strncpy(MassDeltas[AminoIndex][ModIndex].Name, StrName, 40);
                    MassDeltas[AminoIndex][ModIndex].RealDelta = rint(MassDelta * MASS_SCALE);
                    ROUND_MASS_TO_DELTA_BIN(MassDelta, Bin);
                    MassDeltas[AminoIndex][ModIndex].Delta = Bin;
                    MassDeltas[AminoIndex][ModIndex].Index = AllPTModCount;
                    MassDeltaByIndex[AminoIndex * MAX_PT_MODTYPE + AllPTModCount] = &MassDeltas[AminoIndex][ModIndex];
                    MassDeltaByIndex[MDBI_ALL_MODS * MAX_PT_MODTYPE + AllPTModCount] = &MassDeltas[AminoIndex][ModIndex];
                    MassDeltas[AminoIndex][ModIndex].Flags = ModFlags;
                    break;
                }
            }
        }
    }
    if (!AminoFoundFlag)
    {
        REPORT_ERROR_S(37, StrAminos);
        return 0;
    }
    AllKnownPTMods[AllPTModCount].Mass = rint(MassDelta * MASS_SCALE);
    g_PTMLimit[AllPTModCount] = 2; // allow 2 per peptide by default
    // But, only allow ONE c-terminal one:
    if ((ModFlags & DELTA_FLAG_C_TERMINAL) || (ModFlags & DELTA_FLAG_N_TERMINAL))
    {
        g_PTMLimit[AllPTModCount] = 1;
    }
    AllPTModCount++;
    return 1;
}
int ParseInputPTM(char* CommandValue)
{
    printf("*** The 'ptm' input command is no longer supported - please use 'mod' instead.\n");
    printf(" (Refer to the documentation for details)\n");
    return 0;
}
int ParseInputSequenceFile(char* CommandValue)
{
    PrepareSecondarySequenceFile(CommandValue);
    return 1;
}

int ParseInputReadGFF(char* CommandValue)
{
    FILE* GFFFile;
    StringNode* Node;
    // Check to be sure we can read the file:
    GFFFile = fopen(CommandValue, "rb");
    if (!GFFFile)
    {
        REPORT_ERROR_S(8, CommandValue);
    }
    else
    {
        // File is ok - add it to the GFF file list.
        fclose(GFFFile);
        Node = (StringNode*)calloc(1, sizeof(StringNode));
        Node->String = strdup(CommandValue);
        if (GlobalOptions->LastGFFFileName)
        {
            GlobalOptions->LastGFFFileName->Next = Node;
        }
        else
        {
            GlobalOptions->FirstGFFFileName = Node;
        }
        GlobalOptions->LastGFFFileName = Node;
    }
    GlobalOptions->RunMode = RUN_MODE_PREP_MS2DB;
    return 1;
}

int ParseInputGenomeFile(char* CommandValue)
{
    strncpy(GlobalOptions->GenomeFileName, CommandValue, MAX_FILENAME_LEN);
    GlobalOptions->RunMode = RUN_MODE_PREP_MS2DB;
    return 1;
}

int ParseInputChromosomeName(char* CommandValue)
{
    strncpy(GlobalOptions->ChromosomeName, CommandValue, 256);
    return 1;
}

int ParseInputParentPPM(char* ValueString)
{
    int CommandValue = atoi(ValueString);
    if (CommandValue < 1 || CommandValue > 4000)
    {
        REPORT_ERROR_I(44, CommandValue);
        return 0;
    }
    GlobalOptions->ParentMassPPM = CommandValue;
    return 1;
}

int ParseInputPeakPPM(char* ValueString)
{
    int CommandValue = atoi(ValueString);
    if (CommandValue < 1 || CommandValue > 1000)
    {
        REPORT_ERROR_I(44, CommandValue);
        return 0;
    }
    GlobalOptions->PeakPPM = CommandValue;
    return 1;
}

int ParseInputNewScoring(char* Value)
{
    GlobalOptions->NewScoring = 1;
    return 1;
}
static const InputParameter InputParameters[] =
{
    {"Blind", ParseInputBlind, INPUT_VALUE_TYPE_INT},
    {"Unrestrictive", ParseInputBlind, INPUT_VALUE_TYPE_INT},
    {"BlindTagging", ParseInputBlindTagging, INPUT_VALUE_TYPE_INT},
    {"Database", ParseInputDB, INPUT_VALUE_TYPE_STRING},
    {"DB", ParseInputDB, INPUT_VALUE_TYPE_STRING},
    {"ExternalTagger", ParseInputExternalTagger, INPUT_VALUE_TYPE_NONE},
    {"FreeMods", ParseInputFreeMods, INPUT_VALUE_TYPE_INT},
    {"Instrument", ParseInputInstrument, INPUT_VALUE_TYPE_STRING},
    {"IonTolerance", ParseInputIonTolerance, INPUT_VALUE_TYPE_STRING},
    {"MaxPTMSize", ParseInputMaxPTMSize, INPUT_VALUE_TYPE_INT},
    {"Mod", ParseInputMod, INPUT_VALUE_TYPE_STRING},
    {"Mods", ParseInputMods, INPUT_VALUE_TYPE_INT},
    {"MultiCharge", ParseInputMultiCharge, INPUT_VALUE_TYPE_INT},
    {"PMCOnly", ParseInputPMCOnly, INPUT_VALUE_TYPE_INT},
    {"PMTolerance", ParseInputPMTolerance, INPUT_VALUE_TYPE_STRING},
    {"PM_Tolerance", ParseInputPMTolerance, INPUT_VALUE_TYPE_STRING},  // deprecated
    {"PRMModel", ParseInputPRMModel, INPUT_VALUE_TYPE_STRING},
    {"Protease", ParseInputProtease, INPUT_VALUE_TYPE_STRING},
    {"ReportMatches", ParseInputReportMatches, INPUT_VALUE_TYPE_INT},
    {"RequireTermini", ParseInputRequireTermini, INPUT_VALUE_TYPE_INT},
    {"RequiredMod", ParseInputRequiredMod, INPUT_VALUE_TYPE_STRING},
    {"SequenceFile", ParseInputSequenceFile, INPUT_VALUE_TYPE_STRING},
    {"Spectra", ParseInputSpectra, INPUT_VALUE_TYPE_STRING},
    {"TagCheck", ParseInputTagCheck, INPUT_VALUE_TYPE_STRING},
    {"TagCount", ParseInputTagCount, INPUT_VALUE_TYPE_INT},
    {"TagCountB", ParseInputTagCount, INPUT_VALUE_TYPE_INT}, // deprecated
    {"TagLength", ParseInputTagLength, INPUT_VALUE_TYPE_INT},
    {"TAGModel", ParseInputTAGModel, INPUT_VALUE_TYPE_STRING},
    {"Tagless", ParseInputTagless, INPUT_VALUE_TYPE_INT},
    {"TagsOnly", ParseInputTagsOnly, INPUT_VALUE_TYPE_NONE},
    {"XMLStrict", ParseInputXMLStrict, INPUT_VALUE_TYPE_INT},


    // Commands for preparing MS2DB files:
    {"ReadGFF", ParseInputReadGFF, INPUT_VALUE_TYPE_STRING},
    {"GenomeFile", ParseInputGenomeFile, INPUT_VALUE_TYPE_STRING},
    {"ChromosomeName", ParseInputChromosomeName, INPUT_VALUE_TYPE_STRING},
    {"ParentPPM", ParseInputParentPPM, INPUT_VALUE_TYPE_INT},
    {"PeakPPM", ParseInputPeakPPM, INPUT_VALUE_TYPE_INT},
    {"NewScoring",ParseInputNewScoring,INPUT_VALUE_TYPE_NONE},

    // Sentinel:
    {NULL}
};

// Process one line from the inspect input file.  Lines have the form "command,value".
int ProcessInputCommand(int LineNumber, int FilePos, char* LineBuffer, void* UserData)
{
    const InputParameter* Parameter;
    int CommandMatched = 0;
    int ValueOK = 1;
    char* CheckChar;
    char* CommandName;
    char* Value;
    //
    CommandName = strtok(LineBuffer, ",");
    Value = strtok(NULL, ",");
    for (Parameter = InputParameters; Parameter->Name; Parameter++)
    {
        if (CompareStrings(CommandName, Parameter->Name))
        {
            continue;
        }
        CommandMatched = 1;
        // Validate the value:
        switch (Parameter->ValueType)
        {
        case INPUT_VALUE_TYPE_NONE:
            if (Value && *Value)
            {
                REPORT_ERROR_S(39, CommandName);
                ValueOK = 0;
            }
            break;
        case INPUT_VALUE_TYPE_STRING:
            if (!Value || !*Value)
            {
                REPORT_ERROR_S(40, CommandName);
                ValueOK = 0;
            }
            break;
        case INPUT_VALUE_TYPE_INT:
            if (!Value || !*Value)
            {
                REPORT_ERROR_S(41, CommandName);
                ValueOK = 0;
                break;
            }
            for (CheckChar = Value; *CheckChar; CheckChar++)
            {
                if (!isdigit(*CheckChar))
                {
                    REPORT_ERROR_S(41, CommandName);
                    ValueOK = 0;
                    break;
                }
            }
            break;
        }
        if (ValueOK)
        {
            Parameter->ParseFunction(Value);
        }
    }
    if (!CommandMatched)
    {
        REPORT_ERROR_S(13, CommandName);
    }
    return 1;
}

// Parse the input file; return TRUE if successful.
int ParseInputFile()
{
    FILE* InputFile;
    int ModIndex;

    ///////////////////
    InputFile = fopen(GlobalOptions->InputFileName, "rb");
    if (!InputFile)
    {
        REPORT_ERROR_S(8, GlobalOptions->InputFileName);
        return 0;
    }
    ParseFileByLines(InputFile, ProcessInputCommand, NULL, 0);
    fclose(InputFile);

    // PTM processing:
    if (AllPTModCount && !GlobalOptions->MaxPTMods)
    {
        // This is worrisome - the user has defined modifications, but matches are not
        // permitted to USE modifications.  That is reasonable only under weird circumstances.
        if (GlobalOptions->RunMode & RUN_MODE_TAGS_ONLY && GlobalOptions->ExternalTagger)
        {
            //
        }
        else
        {
            REPORT_ERROR(34);
        }
    }
    for (ModIndex = 0; ModIndex < AllPTModCount; ModIndex++)
    {
        g_PTMLimit[ModIndex] = min(g_PTMLimit[ModIndex], GlobalOptions->MaxPTMods);
    }
    if (GlobalOptions->MaxPTMods > 2)
    {
        if (GlobalOptions->RunMode & (RUN_MODE_MUTATION | RUN_MODE_BLIND))
        {
            printf("** Warning: Unrestrictive search with more than two mods is NOT recommended.\n");
        }
    }
    // Set the flanking mass tolerance: Equal to parent mass tolerance plus ion tolerance
    // plus 0.1
    GlobalOptions->FlankingMassEpsilon = GlobalOptions->ParentMassEpsilon + GlobalOptions->Epsilon + 10;

    if (GlobalOptions->ErrorCount)
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

int CompareSpectrumNodes(const SpectrumNode* NodeA, const SpectrumNode* NodeB)
{
    int NameResult;
    NameResult = strcmp(NodeA->InputFile->FileName, NodeB->InputFile->FileName);
    if (NameResult)
    {
        return NameResult;
    }
    //return (strcmp(NodeA->FileName, NodeB->FileName));
    return (NodeA->FilePosition - NodeB->FilePosition);
}

// Sort spectra by filename.
void SortSpectra()
{
    SpectrumNode* Node;
    SpectrumNode* Prev;
    int NodeIndex;
    int NodeCount;
    //
    if (!GlobalOptions->FirstSpectrum)
    {
        return;
    }
    g_BigNodeArray = (SpectrumNode*)calloc(GlobalOptions->SpectrumCount, sizeof(SpectrumNode));
    NodeIndex = 0;
    for (Node = GlobalOptions->FirstSpectrum; Node; Node = Node->Next)
    {
        memcpy(g_BigNodeArray + NodeIndex, Node, sizeof(SpectrumNode));
        NodeIndex++;
    }
    NodeCount = NodeIndex;

    // Free old list:
    Prev = NULL;
    for (Node = GlobalOptions->FirstSpectrum; Node; Node = Node->Next)
    {
        SafeFree(Prev);
        Prev = Node;
    }
    SafeFree(Prev);

    // Sort array:
    qsort(g_BigNodeArray, NodeCount, sizeof(SpectrumNode), (QSortCompare)CompareSpectrumNodes);
    for (NodeIndex = 0; NodeIndex < NodeCount; NodeIndex++)
    {
        if (NodeIndex < NodeCount-1)
        {
            g_BigNodeArray[NodeIndex].Next = g_BigNodeArray + NodeIndex + 1;
        }
        else
        {
            g_BigNodeArray[NodeIndex].Next = NULL;
        }
    }
    GlobalOptions->FirstSpectrum = g_BigNodeArray;
    GlobalOptions->LastSpectrum = g_BigNodeArray + (NodeCount - 1);
}


