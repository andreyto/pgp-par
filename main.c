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


// Inspect - a tool for efficient peptide MS/MS interpretation in the
// presence of post-translational modifications.

// Inspect can use partial de novo to generate tags, then use the tags 
// to search a protein database for matching peptides.  A tag
// has a prefix mass, a sequence of peptides, and a suffix mass.  Typically, 
// tags are tripeptides.  We use a trie structure (Aho-Corasic algorithm) to find 
// occurrences of our tag strings in the database, then examine the flanking masses 
// to be sure they match.  The flanking mass comparison is a d.p. 'hit extension' 
// algorithm.
//
// Inspect requires a database-file in the correct format.  The file
// should contain peptides concatenated together, separated by asterisks.
// No whitespace or newlines.  Like this:  PANTS*GWWYTT*GAAH
// The PrepDB.py script compresses a Swiss-prot or FASTA database into
// concatenated format.  An accompanying .index file is produced, so that the
// name of a matched protein can be reported.

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
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
#include "LDA.h"
#include "Run.h"
#include "SNP.h"
#include "SpliceDB.h"
#include "ChargeState.h"
#include "Scorpion.h"
#include "ParseXML.h"
#include "SpliceScan.h"
#include "ParseInput.h"
#include "PValue.h"
#include "Errors.h"
#include "BuildMS2DB.h"
#include "IonScoring.h"

// Global variables, shared between main.c and Trie.c:
extern Options* GlobalOptions;
extern MSSpectrum* Spectrum;

// Array of spectra to be searched.  We put them into an array so that we can qsort
// them.  (Not crucial, but it's nice to get output in order)
extern SpectrumNode* g_BigNodeArray;

extern StringNode* FirstTagCheckNode;
extern StringNode* LastTagCheckNode;


void PrintUsageInfo()
{
    printf("\nSample command-line:\n");
    printf("Inspect.exe -i Foo.in -o Foo.txt -e ErrorsFoo.txt\n");
    printf("Command-line arguments:\n");
    printf(" -i InputFileName: Path to a config file specifying search parameters.\n");
    printf(" -o OutputFileName: Output file for match results.  If not\n");
    printf("          specified, output goes to stdout.\n");
    printf(" -e ErrorFileName: Output file for errors and warnings, if any.  If not\n");
    printf("          specified, any errors go to Inspect.err; if there are no errors.\n");
    printf("          or warnings reported, this file will be erased at end of run.\n");
    printf(" -r ResourceDir: Directory for resource files (such \n");
    printf("     as AminoAcidMasses.txt).  Defaults to current directory. \n");
    printf("  Consult the documentation (Inspect.html) for further details.\n");
}

void FreeSpectra()
{
    SpectrumNode* Node;
    SpectrumNode* Prev = NULL;
    //
    for (Node = GlobalOptions->FirstSpectrum; Node; Node = Node->Next)
    {
        if (Prev)
        {
            FreeSpectrum(Prev->Spectrum);
            Prev->Spectrum = NULL;
            // Important: don't free spectrum nodes, because they all come from one big array!
            //FreeSpectrumNode(Prev);
        }
        Prev = Node;
        //FreeSpectrum(Node->Spectrum);
    }
    if (Prev)
    {
        FreeSpectrum(Prev->Spectrum);
        Prev->Spectrum = NULL;
    }
    GlobalOptions->FirstSpectrum = NULL;
    GlobalOptions->LastSpectrum = NULL;
}


void FreeGlobalOptions()
{

    StringNode* Prev;
    StringNode* GFFNode;
    DatabaseFile* PrevDB;
    DatabaseFile* DatabaseNode;

    if (!GlobalOptions)
    {
        return;
    }

    // Free the list FirstGFFFileName...LastGFFFileName
    Prev = NULL;
    for (GFFNode = GlobalOptions->FirstGFFFileName; GFFNode; GFFNode = GFFNode->Next)
    {
        if (Prev)
        {
            SafeFree(Prev->String);
            SafeFree(Prev);
        }
        Prev = GFFNode;
    }
    if (Prev)
    {
        SafeFree(Prev->String);
        SafeFree(Prev);
    }

    // Free the DatabaseFile list:
    PrevDB = NULL;
    for (DatabaseNode = GlobalOptions->FirstDatabase; DatabaseNode; DatabaseNode = DatabaseNode->Next)
    {
        SafeFree(PrevDB);
        PrevDB = DatabaseNode;
    }

    // Save the overall struct:
    SafeFree(GlobalOptions);
    GlobalOptions = NULL;
}

// Free various structs that we built up.  (This isn't strictly necessary, since we're about
// to exit process anyway, but it's good practice)
// NOTE: After calling Cleanup(), you can't call Log() any more, because GlobalOptions no longer
// points at a log file.
void Cleanup()
{
    //printf("Cleaning up...\n");
    FreeMassDeltaByMass();
    FreeMassDeltas();
    FreeIsSubDecoration();
    //FreeTaggingModel();
    FreeJumpingHash();
    FreeSVMModels();
    FreeBayesianModels();
    FreeTagCheckNodes();
    FreeInputFileNodes();
    FreeLDAModels();
    FreeCCModelSVM();
    FreeTagSkewScores();
    if (GlobalOptions)
    {
        FreeSpectra();
        // Close our error file.  And if we never wrote errors or warnings, erase it!
        if (GlobalOptions->ErrorFile)
        {
            fclose(GlobalOptions->ErrorFile);
            GlobalOptions->ErrorFile = NULL;
            if (!GlobalOptions->ErrorCount && !GlobalOptions->WarningCount)
            {
                unlink(GlobalOptions->ErrorFileName);
            }
        }
        FreeGlobalOptions();
    }
    SafeFree(g_BigNodeArray);
    g_BigNodeArray = NULL;
    SafeFree(GlobalStats);
    GlobalStats = NULL;
}


// Parse the command-line arguments, and populate GlobalOptions.  
// Returns true on success, 0 if the args are invalid.
int ReadCommandLineArgs(int argc, char** argv)
{
    int Index = 1;
    int MoreArgs;
    int Result;
    char PeptideFilePath[2048];
    if (argc<2)
    {
        return 0;
    }
    while (Index < argc)
    {
        if (argv[Index][0] != '-')
        {
            REPORT_ERROR_S(18, argv[Index]);
            return 0;
        }
        // Are there args after this one?
        if (Index < argc-1)
        {
            MoreArgs = 1;
        }
        else
        {
            MoreArgs = 0;
        }
        switch (ConvertToLower(argv[Index][1]))
        {
        case 'i': // Input options file name
            if (!MoreArgs)
            {
                REPORT_ERROR_S(19, "-i");
                return 0;
            }
            strncpy(GlobalOptions->InputFileName, argv[Index + 1], MAX_FILENAME_LEN);
            Index += 2;
            break;
        case 'o':
            if (!MoreArgs)
            {
                REPORT_ERROR_S(19, "-o");
                return 0;
            }
            strncpy(GlobalOptions->FinalOutputFileName, argv[Index + 1], MAX_FILENAME_LEN);
            Index += 2;
            break;        
        case 'e':
            if (!MoreArgs)
            {
                REPORT_ERROR_S(19, "-e");
                return 0;
            }
            strncpy(GlobalOptions->ErrorFileName, argv[Index + 1], MAX_FILENAME_LEN);
            Index += 2;
            break;        

        case 'r':
            if (!MoreArgs)
            {
                REPORT_ERROR_S(19, "-r");
                return 0;
            }
            strcpy(GlobalOptions->ResourceDir, argv[Index + 1]);
            printf("Setting resource directory: '%s'\n", argv[Index + 1]);
            if (*(GlobalOptions->ResourceDir + strlen(GlobalOptions->ResourceDir) - 1) != SEPARATOR)
            {
                strcat(GlobalOptions->ResourceDir, SEPARATOR_STRING);
            }
            printf("Resource directory is: '%s'\n", GlobalOptions->ResourceDir);
            Index += 2;
            break;
        case 'v':
            GlobalOptions->VerboseFlag = 1;
            Index++;
            break;
        default:
            printf("Error: I don't understand this argument '%s'.\n", argv[Index]);
            return 0;
        }
    }

    // Read the table of amino acid masses:
    sprintf(PeptideFilePath, "%s%s", GlobalOptions->ResourceDir, GlobalOptions->AminoFileName);
    Result = LoadPeptideMasses(PeptideFilePath);
    if (!Result)
    {
        printf("Error - couldn't load amino acid masses!\n");
        return 1;
    }
    // If -r argument wasn't passed, then use the current working directory:
    if (!GlobalOptions->ResourceDir[0])
    {
        sprintf(GlobalOptions->ResourceDir, ".%c", SEPARATOR);
    }
    if (GlobalOptions->InputFileName)
    {
        //printf("Parse input file:\n");
        Result = ParseInputFile();
        //printf("Input file parse result %d\n", Result);
        if (!Result)
        {
            return 0;
        }
        SortSpectra();
    }

    // If no spectra were specified, then error out - unless we're running a 
    // mode that requires no spectra.
    if (!GlobalOptions->FirstSpectrum)
    {
        if (!GlobalOptions->RunMode & (RUN_MODE_PREP_MS2DB))
        {
            REPORT_ERROR(11);
            return 0;
        }
    }

    if (!(*GlobalOptions->FinalOutputFileName))
    {
        sprintf(GlobalOptions->FinalOutputFileName, "Inspect.txt");
    }

    return 1;
}

// Perform miscellaneous chores *after* reading the input script and *before* starting to search.
int Initialize()
{
    char Path[2048];

    if (!MassDeltas)
    {
        sprintf(Path, "%s%s", GlobalOptions->ResourceDir, FILENAME_MASS_DELTAS);
        if (GlobalOptions->RunMode & RUN_MODE_BLIND)
        {
            //LoadMassDeltas(Path, 0);
        }
        else
        {
            LoadMassDeltas(Path, GlobalOptions->RunMode & RUN_MODE_MUTATION);
        }
    }

    InitBayesianModels();
    SetTagSkewScores();
    if (GlobalOptions->RunMode & (RUN_MODE_BLINDTAG | RUN_MODE_BLIND))
    {
        FreeMassDeltas();
        LoadMassDeltas(NULL, 0);
        AddBlindMods();
    }
    else
    {
        InitMassDeltaByMass();
    }

    PopulateJumpingHash();
    //LoadFlankingAminoEffects();
    //sprintf(Path, "%s%s", GlobalOptions->ResourceDir, FILENAME_SCORING_MODEL);
    //Result = InitScoringModel(Path);
    //if (!Result)
    //{
    //    printf("Error loading scoring model from file '%s'\n", Path);
    //    return 0;
    //}

#ifdef MQSCORE_USE_SVM
    InitPValueSVM();
#else
    InitPValueLDA();
#endif

    return 1;
}

// Offshoot of main() for handling spliced-database creation and maintenance:
void MainSpliceDB(int argc, char** argv)
{
    int ChromosomeNumber;
    int ReverseFlag;
    char* GeneName;
    char* CustomFileName;
    int IntervalStart = -1;
    int IntervalEnd = -1;
    int MinORFLength;
    char SNPFileName[256];
    //
    // inspect <chromosome> <reverseflag> [ GeneName, OutputFileName, IntervalStart, IntervalEnd ]
    
    ChromosomeNumber = atoi(argv[1]);
    ReverseFlag = atoi(argv[2]);
    if (argc > 3)
    {
        MinORFLength = atoi(argv[3]);
    }
    else
    {
        MinORFLength = 50;//DEFAULT_MINIMUM_ORF_LENGTH;
    }
    
    if (MinORFLength == 0)
    {
        MinORFLength = -1;
    }
    printf("MainSpliceDB() chrom %d reverse %d minorf %d\n", ChromosomeNumber, ReverseFlag, MinORFLength);
    // Read a linked-list of all the polymorphisms we'd like to account for:
    sprintf(SNPFileName, "SNP\\%d.snp", ChromosomeNumber);
    ParsePolyNodes(SNPFileName); // %%% ARABIDOPSIS: No polynodes available
    printf("PolyNodes parsed\n");
    if (argc > 4)
    {
        GeneName = argv[4];
        CustomFileName = argv[5];
        IntervalStart = atoi(argv[6]);
        IntervalEnd = atoi(argv[7]);
        PrepareOneGeneSpliceDB(ChromosomeNumber, ReverseFlag, IntervalStart, IntervalEnd, CustomFileName, GeneName, MinORFLength);
    }
    else
    {
        printf("PrepareSpliceDB...\n");
        PrepareSpliceDB(ChromosomeNumber, ReverseFlag, MinORFLength);
    }
    FreePolyNodes();
    
}

// MainTraining() is called if the first command-line argument is "train".
// Syntax is:
// inspect.exe train [model] [OracleFile] [SpectrumDir] [extra]
// Example:
// inspect.exe train pmc c:\ms\TrainingSet.txt c:\ms\TrainingSet
// 
// Output format depends on the particular model, but generally we spew out a delimited text file
// which can be processed by a wrapper-script.
int MainTraining(int argc, char** argv)
{
    char* ModelName;
    char* OracleFile;
    char OracleDir[1024];
    int Len;
    //
    if (argc < 5)
    {
        printf("Error: Not enough arguments to train!\n");
        printf("Please provide model name, oracle file, and spectrum directory.\n");
        printf("Sample: inspect.exe train pmc c:\\ms\\TrainingSet.txt c:\\ms\\TrainingSet\n");
        return -1;
    }
    InitOptions();
    ModelName = argv[2];
    OracleFile = argv[3];
    // Guarantee that OracleDir ends with a delimiter:
    strcpy(OracleDir, argv[4]);
    Len = strlen(OracleDir);

    if (Len && OracleDir[Len] != SEPARATOR)
    {
        OracleDir[Len] = SEPARATOR;
        OracleDir[Len + 1] = '\0';
    }
    // Various trainings are available:
    if (!CompareStrings(ModelName, "pmc"))
    {
        //TrainPMC(OracleFile, OracleDir);
    }
    else if (!CompareStrings(ModelName, "cc"))
    {
        //TrainCC(OracleFile, OracleDir);
    }
    else if (!CompareStrings(ModelName, "pepprm"))
    {
        LoadPeptideMasses("AminoAcidMasses.txt");
        PeptideMass['C'] += CAM_MASS; // ASSUMED: All cysteines in the training set carry the +57 modification.
        //GlobalOptions->InstrumentType = INSTRUMENT_TYPE_QTOF; 
        TrainPepPRM(OracleFile, OracleDir);
    }
    else if (!CompareStrings(ModelName, "tag"))
    {
        LoadPeptideMasses("AminoAcidMasses.txt");
        PeptideMass['C'] += CAM_MASS; // ASSUMED: All cysteines in the training set carry the +57 modification.
        LoadMassDeltas(NULL, 0);
        InitMassDeltaByMass();
        PopulateJumpingHash();
        TrainTagging(OracleFile, OracleDir);
    }

    else
    {
        printf("Unknown model name '%s' - no training performed.\n", ModelName);
    }
    return 0;
}

int MainTesting(int argc, char** argv)
{
    char* ModelName;
    char* OracleFile;
    char OracleDir[1024];
    int Len;
    //

    InitOptions();
    ModelName = argv[2];
    OracleFile = argv[3];
    // Guarantee that OracleDir ends with a delimiter:
    if (argc > 4)
    {
        strcpy(OracleDir, argv[4]);
        Len = strlen(OracleDir);
        if (Len && OracleDir[Len] != SEPARATOR)
        {
            OracleDir[Len] = SEPARATOR;
            OracleDir[Len + 1] = '\0';
        }
    }
    // Various tests are available:
    if (!CompareStrings(ModelName, "pmc"))
    {
        //TestPMC(OracleFile, OracleDir);
    }
    else if (!CompareStrings(ModelName, "splicedbug"))
    {
        TestSpliceDB(argc, argv);
    }
    else if (!CompareStrings(ModelName, "cc"))
    {
        LoadPeptideMasses("AminoAcidMasses.txt");
        PeptideMass['C'] += CAM_MASS; // ASSUMED: All cysteines in the training set carry the +57 modification.
        //TestCC(OracleFile, OracleDir);
    }
    else if (!CompareStrings(ModelName, "prmq"))
    {
        LoadPeptideMasses("AminoAcidMasses.txt");
        PeptideMass['C'] += CAM_MASS; // ASSUMED: All cysteines in the training set carry the +57 modification.
        // The oracle file contains the true match for a spectrum, followed by many false matches.
        // Compute the total (average) PRM score for each, sort them, and report the position of the 
        // true peptide within the list.  (Report a histogram of these positions)
        TestPRMQuickScoring(OracleFile, OracleDir);
    }
    else if (!CompareStrings(ModelName, "pepprm"))
    {
        LoadPeptideMasses("AminoAcidMasses.txt");
        PeptideMass['C'] += CAM_MASS; // ASSUMED: All cysteines in the training set carry the +57 modification.
        TestPepPRM(OracleFile, OracleDir);
    }
    else if (!CompareStrings(ModelName, "lda"))
    {
        LoadPeptideMasses("AminoAcidMasses.txt");
        PeptideMass['C'] += CAM_MASS; // ASSUMED: All cysteines in the training set carry the +57 modification.
        TestLDA(OracleFile, OracleDir);
    }
    else if (!CompareStrings(ModelName, "tag"))
    {
        LoadPeptideMasses("AminoAcidMasses.txt");
        PeptideMass['C'] += CAM_MASS; // ASSUMED: All cysteines in the training set carry the +57 modification.
        LoadMassDeltas(NULL, 0);
        InitMassDeltaByMass();
        PopulateJumpingHash();
        // The oracle file contains the true match for a spectrum, followed by many false matches.
        // Compute the total (average) PRM score for each, sort them, and report the position of the 
        // true peptide within the list.  (Report a histogram of these positions)
        TestTagging(OracleFile, OracleDir);
    }
    else if (!CompareStrings(ModelName, "pvalue"))
    {
        // Read in positive and negative feature-vectors, and produce a histogram:
        TestPValue(OracleFile);
    }
    else
    {
        printf("Unknown model name '%s' - no testing performed.\n", ModelName);
    }

    return 0;

}

// SpliceFind: Arguments are a genomic database, and "string-table" protein database.
// We also get a range of protein numbers.  We then look through the genome to find the 
// best (approximate!) match for each protein.
int MainSpliceFind(int argc, char* argv[])
{
    int FirstRecord;
    int LastRecord; // inclusive!
    char IndexFileName[512];
    char* Temp;
    //
    if (argc < 6)
    {
        printf("** Not enough args for splice find.  Sample run:\n");
        printf("inspect splicefind database\\ipiv313.trie ESTSpliceDB\\Genome.dat 0 1000\n");
        return -1;
    }
    FirstRecord = atoi(argv[4]);
    LastRecord = atoi(argv[5]);
    if (LastRecord <= FirstRecord && LastRecord > -1)
    {
        printf("** Bad record numbers: %s to %s\n", argv[4], argv[5]);
        return -1;
    }
    strcpy(IndexFileName, argv[2]);
    for (Temp = IndexFileName + strlen(IndexFileName); Temp >= IndexFileName; Temp--)
    {
        if (*Temp == '.')
        {
            *Temp = '\0';
            break;
        }
    }
    strcat(IndexFileName, ".index");
    SSDatabaseScan(argv[2], IndexFileName, argv[3], FirstRecord, LastRecord);
    return 1;
}

int LoadAndScoreSpectrum()
{
    //char* FilePath = "PTMScore\\HEKMerged\\Spectra\\H\\R.HIADLAGNSEVILPVPAFNVINGGS+244HAG.N.2.dta";
    //char* Annotation = "R.HIADLAGNSEVILPVPAFNVINGGS+244HAG.N";
    char* FilePath = "SystemTest\\TestSpectrum.dta";
    char* Annotation = "VKEAMAPK";
    MSSpectrum* Spectrum;
    int FilePosition = 0; // Default: byte offset 0
    SpectrumNode* Node;
    FILE* SpectrumFile;
    //
    Node = (SpectrumNode*)calloc(1, sizeof(SpectrumNode));
    Node->FilePosition = FilePosition;
    Node->ScanNumber = 0;
    Node->InputFile = (InputFileNode*)calloc(1, sizeof(InputFileNode));
    strncpy(Node->InputFile->FileName, FilePath, MAX_FILENAME_LEN);
    // Guess the file format:
    Node->InputFile->Format = GuessSpectrumFormatFromExtension(FilePath);
    SpectrumFile = fopen(FilePath, "rb");
    fseek(SpectrumFile, Node->FilePosition, 0);
    Node->Spectrum = (MSSpectrum*)calloc(1, sizeof(MSSpectrum));
    Spectrum = Node->Spectrum;
    Node->Spectrum->Node = Node;
    SpectrumLoadFromFile(Node->Spectrum, SpectrumFile);
    fclose(SpectrumFile);
    WindowFilterPeaks(Node->Spectrum, 0, 0);
    IntensityRankPeaks(Node->Spectrum);
    //SpectrumComputeNoiseDistributions(Node);
    //SpectrumComputeBinnedIntensities(Node);
    printf("Tweak and score...\n");
    TweakSpectrum(Node);
    ////////////////////////////////////
    // Score:
    ////////////////////////////////////
    // Free:
    // The PySpectrum object wraps a Spectrum object, but also a SpectrumNode and an InputFileNode.  
    // So, free those as well:
    if (Spectrum->Node->InputFile)
    {
        free(Spectrum->Node->InputFile);
        Spectrum->Node->InputFile = NULL;
    }
    if (Spectrum->Node)
    {
        FreeSpectrumNode(Spectrum->Node);
    }
    else
    {
        FreeSpectrum(Spectrum);
    }
    return 0;
}

int TestMain(int argc, char* argv[])
{
    char Buffer[2048];
    // For temp test scaffolding
    InitOptions();
    InitErrors();
    InitStats();
    Initialize();
    printf(">>> Start <<<\n");
    Cleanup();
    printf(">>> End <<<\n");
    ReadBinary(Buffer, sizeof(char), 1, stdin);
    return 1;
}

// Program entry point.  Parses arguments, does initialization of global data, 
// then either runs unit tests or calls RunTrieSearch.
int main(int argc, char** argv)
{
    int Result;
    clock_t StartTime;
    clock_t EndTime;
    float ElapsedTime;
    int ChromosomeNumber;
    //

    //return TestMain(argc, argv);
    // Jump into the training/testing code, maybe:
    if (argc > 1 && !CompareStrings(argv[1], "train"))
    {
        return MainTraining(argc, argv);
    }
    if (argc > 1 && !CompareStrings(argv[1], "test"))
    {
        return MainTesting(argc, argv);
    }
    if (argc > 1 && !CompareStrings(argv[1], "splicefind"))
    {
        return MainSpliceFind(argc, argv);
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////

    printf("\nInsPecT vesrion %s\n  Interpretation of Peptides with Post-translational Modifications.\n", INSPECT_VERSION_NUMBER);
    printf("  Copyright 2007, The Regents of the University of California\n");
    printf("  [See Docs directory for usage manual and copyright information]\n\n");
    fflush(stdout);
    // Allocate stuff:
    AllocMassDeltaByIndex();

    // Slightly hacky behavior: If the first argument is an integer, then
    // jump to the splice-db code:
    if (argc > 1)
    {
        ChromosomeNumber = atoi(argv[1]);
        if (ChromosomeNumber)
        {
            MainSpliceDB(argc, argv);
            goto cleanup;
        }
    }

    // Set the (global) default options:
    InitOptions();
    InitErrors();
    InitStats();

    // Parse arguments.  If ReadCommandLineArgs returns false, we didn't get
    // valid arguments, so we print usage info and quit.
    Result = ReadCommandLineArgs(argc, argv);
    if (!Result)
    {
        PrintUsageInfo();
        goto cleanup;
    }

    // Open the error file *after* parsing the command-line:
    GlobalOptions->ErrorFile = fopen(GlobalOptions->ErrorFileName, "wb");
    if (!GlobalOptions->ErrorFile)
    {
        GlobalOptions->ErrorFile = stderr;
    }

    printf("Initialize:\n");
    Result = Initialize();
    if (!Result)
    {
        printf("Initialization FAILED - aborting search.\n");
        goto cleanup;
    }
    printf("Initialization successful\n");

    ///////////////////////////////////////////////////
    // Main function: Run the search!
    StartTime = clock();

    // Set an intermediate output file name, if we're performing a search.
    // (We write to the intermediate file, then perform p-value computation)
    if (!(GlobalOptions->RunMode & (RUN_MODE_TAGS_ONLY | RUN_MODE_PMC_ONLY | RUN_MODE_PREP_MS2DB)))
    {
        sprintf(GlobalOptions->OutputFileName, "%s.tmp", GlobalOptions->FinalOutputFileName);
    }
    else
    {
        sprintf(GlobalOptions->OutputFileName, "%s", GlobalOptions->FinalOutputFileName);
    }
    GlobalOptions->OutputFile = fopen(GlobalOptions->OutputFileName, "wb");
    if (!GlobalOptions->OutputFile)
    {
        REPORT_ERROR_S(8, GlobalOptions->OutputFileName);
        goto cleanup;
    }

    if (GlobalOptions->RunMode & RUN_MODE_PREP_MS2DB)
    {
        BuildMS2DB();
    }
    else if (GlobalOptions->RunMode & RUN_MODE_PMC_ONLY)
    {
        // Just correct charges and parent masses, don't search anything:
        PerformSpectrumTweakage();
    }
    else if ((GlobalOptions->RunMode & RUN_MODE_TAGS_ONLY) && !GlobalOptions->ExternalTagger)
    {
        PerformTagGeneration();
    }
    else
    {
        RunSearch();
    }
    
    EndTime = clock();
    ElapsedTime = (float)((EndTime - StartTime) / (float)CLOCKS_PER_SEC);
    printf("Elapsed time: %.4f seconds.\n", ElapsedTime);
    printf("Inspect run complete.\n");
 
cleanup:
    Cleanup(); 
    return 0;
}
