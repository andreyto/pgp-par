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
#include "CMemLeak.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include "Trie.h"
#include "Utils.h"
#include "Run.h"
#include "Tagger.h"
#include "Score.h"
#include "FreeMod.h"
#include "Spliced.h"
#include "BN.h"
#include "SVM.h"
#include "Scorpion.h"
#include "ChargeState.h"
#include "PValue.h"
#include "MS2DB.h"
#include "IonScoring.h"

extern float g_CutScores[];
extern PRMBayesianModel* PRMModelCharge2;

// Forward Declaration
void DebugPrintBlindTagExtensions(SearchInfo* Info);
void AttemptParentMassPeakRemoval(MSSpectrum* Spectrum);
void RestoreParentMassPeakRemoval(MSSpectrum* Spectrum);


TrieTag* TagGraphGenerateTags(TagGraph* Graph, MSSpectrum* Spectrum, int* TagCount, 
    int MaximumTagCount, SpectrumTweak* Tweak, float TagEdgeScoreMultiplier,
    PRMBayesianModel* Model);

void OutputMatchesForSpectrum(SpectrumNode* Node, FILE* OutputFile)
{
    char MatchedPeptideVerbose[256];
    char PeptideName[256];
    Peptide* Match;
    //Peptide* NextMatch;
    int MatchNumber = 0;
    int FeatureIndex;
    static int FirstCall = 1;
    PeptideSpliceNode* SpliceNode;
    //
    if (!OutputFile)
    {
        return;
    }
    PeptideName[0] = '\0'; 
    SetMatchDeltaCN(Node);

    //SetMatchDeltaCN(Node);
    Match = Node->FirstMatch;
    while (Match)
    {
        GetProteinID(Match->RecordNumber, Match->DB, PeptideName);

        // Write a header line:
        if (FirstCall)
        {
            FirstCall = 0;
            fprintf(OutputFile, "#SpectrumFile\tScan#\tAnnotation\tProtein\tCharge\t");
            fprintf(OutputFile, "MQScore\tLength\tTotalPRMScore\tMedianPRMScore\tFractionY\tFractionB\tIntensity\tNTT\t");
            fprintf(OutputFile, "p-value\tF-Score\t");
            fprintf(OutputFile, "DeltaScore\tDeltaScoreOther\tRecordNumber\tDBFilePos\tSpecFilePos\t");
            if (GlobalOptions->FirstDatabase->Type != evDBTypeTrie)
            {
                fprintf(OutputFile, "Chromosome\tStrand\tGenomicPos\tSplicedSequence\tSplices\t");
            }
            fprintf(OutputFile, "\n");
        }
        
        //GetProteinID(Match->RecordNumber, IndexFile, PeptideName, 1);
        WriteMatchToString(Match, MatchedPeptideVerbose, 1);

        // Which spectrum?
        fprintf(OutputFile, "%s\t%d\t", Node->InputFile->FileName, Node->ScanNumber);
        // What's the match?
        fprintf(OutputFile, "%s\t%s\t%d\t%.3f\t", MatchedPeptideVerbose, PeptideName, Match->Tweak->Charge, Match->MatchQualityScore);
        // How good is the match?
        for (FeatureIndex = 0; FeatureIndex < MQ_FEATURE_COUNT; FeatureIndex++)
        {
            fprintf(OutputFile, "%.3f\t", Match->ScoreFeatures[FeatureIndex]);
        }
        //fprintf(OutputFile, "%.3f\t", Match->InitialScore / 1000.0);
        fprintf(OutputFile, "%.5f\t", Match->FScore);
        fprintf(OutputFile, "%.5f\t", Match->PValue);
        fprintf(OutputFile, "%.3f\t", Match->DeltaCN);
        fprintf(OutputFile, "%.3f\t", Match->DeltaCNOther);
        //fprintf(OutputFile, "%.3f\t", Match->ParentMassError / 100.0); // Temp: Parent mass error (for FT)
        // Extra fields, for debugging:
        fprintf(OutputFile, "%d\t%d\t%d\t", Match->RecordNumber, Match->FilePos, Node->FilePosition);
        ////////////////////////////////////////////////////////////
        // If it's a splice-tolerant search, then output some information about the match:
        if (Match->DB->Type != evDBTypeTrie)
        {
            if (Match->ChromosomeNumber >= 0)
            {
                fprintf(OutputFile, "%d\t", Match->ChromosomeNumber);
                fprintf(OutputFile, "%d\t", Match->ChromosomeForwardFlag);
            }
            else
            {
                fprintf(OutputFile, "\t");
                fprintf(OutputFile, "\t");
            }
            if (Match->GenomicLocationStart >= 0)
            {
                fprintf(OutputFile, "%d-%d\t", 
                    min(Match->GenomicLocationStart, Match->GenomicLocationEnd),
                    max(Match->GenomicLocationStart, Match->GenomicLocationEnd));
            }
            else
            {
                fprintf(OutputFile, "\t");
            }
            fprintf(OutputFile, "%s\t", Match->SplicedBases);
            for (SpliceNode = Match->SpliceHead; SpliceNode; SpliceNode = SpliceNode->Next)
            {
                fprintf(OutputFile, "%d-%d ",  SpliceNode->DonorPos, SpliceNode->AcceptorPos);
            }
            fprintf(OutputFile, "\t");
        }
        fprintf(OutputFile, "\n");
        Match = Match->Next;
        MatchNumber++;
        if (MatchNumber >= GlobalOptions->ReportMatchCount)
        {
            break;
        }
    }
    //printf("Wrote out %d matches for '%s'.\n", MatchNumber, Node->FileName);
    fflush(OutputFile);
}

#define vprintf(x) if (VerboseFlag) printf(x)

void MutationModeSearch(SearchInfo* Info)
{
    Peptide* FirstMatch = NULL;
    Peptide* LastMatch = NULL;
    Peptide* NextOldMatchNode;
    Peptide* OldMatchNode;
    Peptide* MatchNode;
    Peptide* NextMatchNode;
    Peptide* FreeNode;
    Peptide* FreePrev;
    int MatchCount = 0;
    int VerboseFlag = 0;
    int TweakIndex;
    MSSpectrum* Spectrum = Info->Spectrum;
    SpectrumNode* Node = Info->Spectrum->Node;

    if (Spectrum->PeakCount < 10) // Demand AT LEAST ten peaks (even that many is a bit silly; 50 is more like it)
    {
        return;
    }
    for (TweakIndex = 0; TweakIndex < TWEAK_COUNT; TweakIndex++)
    {
        if (!Node->Tweaks[TweakIndex].Charge)
        {
            continue;
        }
        fseek(Info->DB->DBFile, 0, 0);
        // *** PRM scores now *** 
        Spectrum->Charge = Node->Tweaks[TweakIndex].Charge;
        Spectrum->ParentMass = Node->Tweaks[TweakIndex].ParentMass;
        //vprintf("[V] Assign isotope neighbors\n");
        //SpectrumAssignIsotopeNeighbors(Node->Spectrum);
        //vprintf("[V] Find isotopic peaks\n");
        //SpectrumFindIsotopicPeaks(Node->Spectrum);
        FreeTagGraph(Node->Spectrum->Graph);
        vprintf("[V] Construct tag graph\n");
        Node->Spectrum->Graph = ConstructTagGraph(Node->Spectrum);
        vprintf("[V] Add nodes\n");
        TagGraphAddNodes(Node->Spectrum->Graph, Node->Spectrum);
        vprintf("[V] Score nodes\n");
        TagGraphScorePRMNodes(NULL, Node->Spectrum->Graph, Node->Spectrum, Node->Tweaks + TweakIndex);
        vprintf("[V] Populate back edges\n");
        if (GlobalOptions->MaxPTMods > 1)
        {
            TagGraphPopulateBackEdges(Node->Spectrum->Graph);
        }
        vprintf("[V] Set PRM scores\n");
        SetSpectrumPRMScores(Node->Spectrum, Node->Tweaks + TweakIndex); 
        vprintf("[V] Tagless search 1:\n");
        SearchDatabaseTagless(Info, GlobalOptions->MaxPTMods, VerboseFlag, Node->Tweaks + TweakIndex);
        ////////////
        vprintf("[V] Score matches:\n");
        vprintf("[V] merge multi-charge list:\n");
        OldMatchNode = FirstMatch;
        if (FirstMatch)
        {
            NextOldMatchNode = FirstMatch->Next;
        }
        else
        {
            NextOldMatchNode = NULL;
        }
        MatchNode = Node->FirstMatch;
        if (MatchNode)
        {
            NextMatchNode = MatchNode->Next;
        }
        else
        {
            NextMatchNode = NULL;
        }
        MatchCount = 0;
        FirstMatch = NULL;
        LastMatch = NULL;
        while (MatchNode || OldMatchNode)
        {
            if (!MatchNode || (OldMatchNode && MatchNode->InitialScore < OldMatchNode->InitialScore))
            {
                // Add one of the old matches to the master-list:
                if (FirstMatch)
                {
                    LastMatch->Next = OldMatchNode;
                    OldMatchNode->Prev = LastMatch;
                    LastMatch = OldMatchNode;
                    LastMatch->Next = NULL;
                }
                else
                {
                    FirstMatch = OldMatchNode;
                    LastMatch = OldMatchNode;
                    OldMatchNode->Prev = NULL;
                    OldMatchNode->Next = NULL;
                }
                OldMatchNode = NextOldMatchNode;
                if (OldMatchNode)
                {
                    NextOldMatchNode = OldMatchNode->Next;
                }
                else
                {
                    NextOldMatchNode = NULL;
                }

            }
            else
            {
                // Add one of the new matches to the master-list:
                if (FirstMatch)
                {
                    LastMatch->Next = MatchNode;
                    MatchNode->Prev = LastMatch;
                    LastMatch = MatchNode;
                    LastMatch->Next = NULL;
                }
                else
                {
                    FirstMatch = MatchNode;
                    LastMatch = MatchNode;
                    MatchNode->Prev = NULL;
                    MatchNode->Next = NULL;
                }
                MatchNode = NextMatchNode;
                if (MatchNode)
                {
                    NextMatchNode = MatchNode->Next;
                }
                else
                {
                    NextMatchNode = NULL;
                }
            }
            MatchCount++;
            if (MatchCount >= GlobalOptions->StoreMatchCount)
            {
                break;
            }
        }
        // Now we can free any remaining matches from these lists:
        FreeNode = MatchNode;
        FreePrev = NULL;
        while (FreeNode)
        {
            if (FreePrev)
            {
                FreePeptideNode(FreePrev);
            }
            FreePrev = FreeNode;
            FreeNode = FreeNode->Next;
        }
        if (FreePrev)
        {
            FreePeptideNode(FreePrev);
        }
        FreeNode = OldMatchNode;
        FreePrev = NULL;
        while (FreeNode)
        {
            if (FreePrev)
            {
                FreePeptideNode(FreePrev);
            }
            FreePrev = FreeNode;
            FreeNode = FreeNode->Next;
        }
        if (FreePrev)
        {
            FreePeptideNode(FreePrev);
        }
        Node->FirstMatch = NULL;
        Node->LastMatch = NULL;
        Node->MatchCount = 0;
        // Check the master-list for duplicates:
        for (MatchNode = FirstMatch; MatchNode; MatchNode = MatchNode->Next)
        {
            for (OldMatchNode = MatchNode->Next; OldMatchNode; OldMatchNode = OldMatchNode->Next)
            {
                if (OldMatchNode->RecordNumber == MatchNode->RecordNumber && !strcmp(OldMatchNode->Bases, MatchNode->Bases) && 
                    !memcmp(MatchNode->AminoIndex, OldMatchNode->AminoIndex, sizeof(int) * MAX_PT_MODS) && 
                    !memcmp(MatchNode->ModType, OldMatchNode->ModType, sizeof(int) * MAX_PT_MODS))
                {
                    // Free OldMatchNode, it's a duplicate!
                    if (OldMatchNode->Prev)
                    {
                        OldMatchNode->Prev->Next = OldMatchNode->Next;
                    }
                    if (OldMatchNode->Next)
                    {
                        OldMatchNode->Next->Prev = OldMatchNode->Prev;
                    }
                    if (LastMatch == OldMatchNode)
                    {
                        LastMatch = OldMatchNode->Prev;
                    }
                    FreePeptideNode(OldMatchNode);
                    OldMatchNode = MatchNode->Next;
                    if (!OldMatchNode)
                    {
                        break;
                    }
                }
            }
        }
    } // tweak loop
    Node->FirstMatch = FirstMatch;
    Node->LastMatch = LastMatch;
    Node->MatchCount = MatchCount;
    vprintf("[V] Complete.\n");
}

TrieNode* ConstructTagsForSpectrum(TrieNode* Root, SpectrumNode* Node, int TagCount)
{
    int TweakIndex;
    MSSpectrum* Spectrum;
    //
    Spectrum = Node->Spectrum;
    for (TweakIndex = 0; TweakIndex < TWEAK_COUNT; TweakIndex++)
    {
        if (!Node->Tweaks[TweakIndex].Charge)
        {
            continue;
        }
        Spectrum->Charge = Node->Tweaks[TweakIndex].Charge;
        Spectrum->ParentMass = Node->Tweaks[TweakIndex].ParentMass;
        //SpectrumAssignIsotopeNeighbors(Spectrum);
        //SpectrumFindIsotopicPeaks(Spectrum);
		//sam Temp Insert
		AttemptParentMassPeakRemoval( Spectrum);
        Root = GenerateTagsFromSpectrum(Spectrum, Root, TagCount, Node->Tweaks + TweakIndex);
		//Sam Temp Insert
		RestoreParentMassPeakRemoval(Spectrum);

    }
    return Root;
}

void OutputTagsToFile(FILE* OutputFile, char* SpectrumFileName, int ScanNumber, int SpectrumFilePos, TrieTag* TagArray, int TagCount)
{
    int TagIndex;
    //TrieTagHanger* Tag;
    TrieTag* Tag;
    TagCount = min(TagCount, GlobalOptions->GenerateTagCount);
    for (TagIndex = 0; TagIndex < TagCount; TagIndex++)
    {
        Tag = TagArray + TagIndex;
        fprintf(OutputFile, "%s\t", SpectrumFileName);
        fprintf(OutputFile, "%d\t", ScanNumber);
        fprintf(OutputFile, "%d\t", SpectrumFilePos);
        fprintf(OutputFile, "%d\t", TagIndex);
        fprintf(OutputFile, "%.2f\t", Tag->PrefixMass / (float)MASS_SCALE);
        fprintf(OutputFile, "%s\t", Tag->Tag);
        fprintf(OutputFile, "%.2f\t", Tag->SuffixMass / (float)MASS_SCALE);
        fprintf(OutputFile, "%.2f\t", Tag->Score);
        fprintf(OutputFile, "\n");
    }
}

int MergeIdenticalTags(TrieTag* TagArray, int TagCount)
{
    int TagIndexA;
    int TagIndexB;
    TrieTag* TagA;
    TrieTag* TagB;
    int Diff;
    //
    for (TagIndexA = 0; TagIndexA < TagCount; TagIndexA++)
    {
        TagA = TagArray + TagIndexA;
        TagIndexB = TagIndexA + 1;
        while (TagIndexB < TagCount)
        {
            TagB = TagArray + TagIndexB;
            if (strcmp(TagA->Tag, TagB->Tag))
            {
                TagIndexB++;
                continue;
            }
            Diff = abs(TagA->PrefixMass - TagB->PrefixMass);
            if (Diff > GlobalOptions->Epsilon)
            {
                TagIndexB++;
                continue;
            }
            Diff = abs(TagA->SuffixMass - TagB->SuffixMass);
            if (Diff > GlobalOptions->Epsilon)
            {
                TagIndexB++;
                continue;
            }
            // These tags are essentially identical!  Remove B.
            memmove(TagArray + TagIndexB, TagArray + TagIndexB + 1, sizeof(TrieTag) * (TagCount - TagIndexB));
            TagCount--;
            // TagIndexB is unchanged.
        }
    }
    return TagCount;
}

static TrieTag* _TagsOnlyTagList = NULL;

// Perform ONLY tag generation...and output the resulting tags.
void PerformTagGeneration()
{
    SpectrumNode* SNode;
    FILE* SpectrumFile;
    int Result;
    int TagCount = GlobalOptions->GenerateTagCount;
    int TweakIndex;
    int TotalTagCount;
    int TagIndex;
    int TagsGenerated;
    SpectrumTweak* Tweak;
    TrieTag* Tags;
    //
    // Write a HEADER to the output file:
    fprintf(GlobalOptions->OutputFile, "#File\tScan\tFilePos\tTagIndex\tPrefix\tTag\tSuffix\tTagscore\t\n");
    if (!_TagsOnlyTagList)
    {
        _TagsOnlyTagList = (TrieTag*)calloc(TWEAK_COUNT * TagCount + 1, sizeof(TrieTag));
    }
    BuildDecorations();
    for (SNode = GlobalOptions->FirstSpectrum; SNode; SNode = SNode->Next)
    {
        SpectrumFile = fopen(SNode->InputFile->FileName, "rb");
        if (SpectrumFile)
        {
            SNode->Spectrum = (MSSpectrum*)calloc(1, sizeof(MSSpectrum));
            SNode->Spectrum->Node = SNode;
            fseek(SpectrumFile, SNode->FilePosition, 0);
            Result = SpectrumLoadFromFile(SNode->Spectrum, SpectrumFile);
            fclose(SpectrumFile);
            if (!Result)
            {
                SafeFree(SNode->Spectrum);
                SNode->Spectrum = NULL;
                continue;
            }
            else
            {
                WindowFilterPeaks(SNode->Spectrum, 0, 0);
                IntensityRankPeaks(SNode->Spectrum);
            }
            if (!SNode->PMCFlag)
            {
                TweakSpectrum(SNode);
            }
            TotalTagCount = 0;
            for (TweakIndex = 0; TweakIndex < TWEAK_COUNT; TweakIndex++)
            {
                if (!SNode->Tweaks[TweakIndex].Charge)
                {
                    continue;
                }
                Tweak = SNode->Tweaks + TweakIndex;
                SNode->Spectrum->Charge = Tweak->Charge;
                SNode->Spectrum->ParentMass = Tweak->ParentMass;
                //SpectrumAssignIsotopeNeighbors(SNode->Spectrum);
                //SpectrumFindIsotopicPeaks(SNode->Spectrum);
                SNode->Spectrum->Graph = ConstructTagGraph(SNode->Spectrum);
                TagGraphAddNodes(SNode->Spectrum->Graph, SNode->Spectrum);
                TagGraphScorePRMNodes(NULL, SNode->Spectrum->Graph, SNode->Spectrum, Tweak);
                TagGraphPopulateEdges(SNode->Spectrum->Graph);
                Tags = TagGraphGenerateTags(SNode->Spectrum->Graph, SNode->Spectrum, &TagsGenerated, TagCount, Tweak, TAG_EDGE_SCORE_MULTIPLIER, NULL);

                for (TagIndex = 0; TagIndex < min(TagCount, TagsGenerated); TagIndex++)
                {
                    memcpy(_TagsOnlyTagList + TotalTagCount, Tags + TagIndex, sizeof(TrieTag));
                    TotalTagCount++;
                }
                FreeTagGraph(SNode->Spectrum->Graph);
                SNode->Spectrum->Graph = NULL;
            } // Tweak list
            qsort(_TagsOnlyTagList, TotalTagCount, sizeof(TrieTag), (QSortCompare)CompareTagScores);
            TotalTagCount = MergeIdenticalTags(_TagsOnlyTagList, TotalTagCount);
            OutputTagsToFile(GlobalOptions->OutputFile, SNode->InputFile->FileName, SNode->ScanNumber, SNode->FilePosition, _TagsOnlyTagList, TotalTagCount);
            // Clean up the spectrum:
            FreeSpectrum(SNode->Spectrum);
            SNode->Spectrum = NULL;
        }        
    }
    SafeFree(_TagsOnlyTagList);
    _TagsOnlyTagList = NULL;
}

#define TEMP_TAGGING_INPUT "TempTagging.dta"
#define TEMP_TAGGING_OUTPUT "TempTags.txt"

// Call upon PepNovo to generate some tags for us.
// To do so, we need to write out a temporary .dta file!
TrieNode* ConstructTagsExternalTagger(TrieNode* Root, SpectrumNode* Node, int TagCount)
{
    FILE* TempDTAFile;
    FILE* TempTagOutputFile;
    int TweakIndex;
    int PeakIndex;
    SpectralPeak* Peak;
    char CommandLine[2048];
    char LineBuffer[MAX_LINE_LENGTH];
    int BytesToRead;
    int BufferPos = 0;
    int BytesRead;
    int BufferEnd = 0;
    int LineNumber = 0;
    int PrevLineFilePos = 0;
    int LineFilePos = 0;
    char TextBuffer[BUFFER_SIZE * 2];
    char* BitA;
    char* BitB;
    char* BitC;
    int WithinTagsFlag = 0;
    float PrefixMass;
    float Probability;
    int DuplicateFlag;
    TrieTag* NewTag;
    char* TempAA;
    char AnnotationBuffer[256];
    char ModBuffer[256];
    int TagIndex = 0;
    int AminoIndex;
    MassDelta* Delta;
    int ModBufferPos;
    int ModIndex;
    int TotalTagCount = 0;
    //
    if (!Root)
    {
        Root = NewTrieNode();
        Root->FailureNode = Root;
    }
    // Initialization for tags-only:
    if (GlobalOptions->RunMode & RUN_MODE_TAGS_ONLY)
    {
        if (!_TagsOnlyTagList)
        {
            _TagsOnlyTagList = (TrieTag*)calloc(TWEAK_COUNT * TagCount + 1, sizeof(TrieTag));
        }
    }

    for (TweakIndex = 0; TweakIndex < TWEAK_COUNT; TweakIndex++)
    {
        // Skip this mass-tweak, if it's not a valid charge/mass combo
        if (!Node->Tweaks[TweakIndex].Charge)
        {
            continue;
        }
        unlink(TEMP_TAGGING_INPUT);
        TempDTAFile = fopen(TEMP_TAGGING_INPUT, "wb");
        if (!TempDTAFile)
        {
            printf("** Error opening tag input file %s for writing!", TEMP_TAGGING_INPUT);
            return Root;
        }
        fprintf(TempDTAFile, "%.3f %d\n", Node->Tweaks[TweakIndex].ParentMass / (float)MASS_SCALE, Node->Tweaks[TweakIndex].Charge);
        for (PeakIndex = 0; PeakIndex < Node->Spectrum->PeakCount; PeakIndex++)
        {
            Peak = Node->Spectrum->Peaks + PeakIndex;
            fprintf(TempDTAFile, "%.3f %.3f\n", Peak->Mass / (float)MASS_SCALE, Peak->Intensity);
        }
        fclose(TempDTAFile);
        // Call out to pepnovo:
        unlink(TEMP_TAGGING_OUTPUT);
        sprintf(CommandLine, "pepnovo.exe -dta %s -model tryp_model.txt -num_tags %d > %s", TEMP_TAGGING_INPUT, TagCount, TEMP_TAGGING_OUTPUT);
        system(CommandLine);
        TempTagOutputFile = fopen(TEMP_TAGGING_OUTPUT, "rb");
        if (!TempTagOutputFile)
        {
            printf("** Error: Unable to open tag output file '%s'\n", TEMP_TAGGING_OUTPUT);
            return Root;
        }
        WithinTagsFlag = 0;
        while (1)
        {
            BytesToRead = BUFFER_SIZE - BufferEnd;
            BytesRead = ReadBinary(TextBuffer + BufferEnd, sizeof(char), BytesToRead, TempTagOutputFile);
            BufferEnd += BytesRead;
            TextBuffer[BufferEnd] = '\0';
            if (BufferPos == BufferEnd)
            { 
                // We're done!
                break;
            }
            // Copy a line of text to the line buffer.  Skip spaces, and stop at carriage return or newline.
            
            BufferPos = CopyBufferLine(TextBuffer, BufferPos, BufferEnd, LineBuffer, 0);
            LineNumber += 1;
            PrevLineFilePos = LineFilePos;
            LineFilePos += BufferPos;
            //printf("Line %d starts at %d\n", LineNumber, LineFilePos);
            // Now, move the remaining text to the start of the buffer:
            memmove(TextBuffer, TextBuffer + BufferPos, BufferEnd - BufferPos);
            BufferEnd -= BufferPos;
            BufferPos = 0;
            // Now, process this line of text!
            if (!LineBuffer[0])
            {
                continue;
            }
            BitA = strtok(LineBuffer, "\t\r\n");
            BitB = strtok(NULL, "\t\r\n");
            if (!BitB)
            {
                continue;
            }
            BitC = strtok(NULL, "\t\r\n");
            if (!BitC)
            {
                continue;
            }
            if (!strcmp(BitC, "Probability:") && !strcmp(BitB, "Tag"))
            {
                WithinTagsFlag = 1;
            }
            if (WithinTagsFlag)
            {
                PrefixMass = (float)atof(BitA);
                Probability = (float)atof(BitC);
                if (Probability < (float)0.1)
                {
                    continue;
                }
                NewTag = _TagsOnlyTagList + TotalTagCount;
                memset(NewTag, 0, sizeof(TrieTag));
                // Special code:
                // PepNovo may include MODIFICATIONS in its tags - so, we must parse them.
                // We assume that (a) modifications are written in the form %+d, and (b) we
                // know of the modification type from the inspect input file.
                TempAA = BitB;
                AminoIndex = 0;
                ModBufferPos = 0;
                
                while (*TempAA)
                {
                    if (*TempAA >= 'A' && *TempAA <= 'Z')
                    {
                        // an amino acid - so, finish the modification-in-progress, if there is one.
                        if (ModBufferPos && AminoIndex)
                        {
                            if (NewTag->ModsUsed == MAX_PT_MODS)
                            {
                                printf("** Error tagging scan %d from file %s: Too many PTMs!\n", Node->ScanNumber, Node->InputFile->FileName);
                                break;
                            }
                            ModBuffer[ModBufferPos] = '\0';
                            Delta = FindPTModByName(NewTag->Tag[AminoIndex - 1], ModBuffer);
                            if (Delta)
                            {
                                NewTag->AminoIndex[NewTag->ModsUsed] = AminoIndex - 1;
                                NewTag->ModType[NewTag->ModsUsed] = Delta;
                                NewTag->ModsUsed++;
                            }
                            else
                            {
                                printf("** Error tagging scan %d from file %s: Modification %s not understood!\n", Node->ScanNumber, Node->InputFile->FileName, ModBuffer);
                            }
                        }
                        ModBufferPos = 0;
                        // Add the AA:
                        NewTag->Tag[AminoIndex++] = *TempAA;
                    }// aa
                    else
                    {
                        ModBuffer[ModBufferPos++] = *TempAA;
                    } // not aa
                    TempAA++;
                }
                NewTag->Tag[AminoIndex] = '\0';
                // Finish any pending mod (COPY-PASTA FROM ABOVE)
                if (ModBufferPos && AminoIndex)
                {
                    if (NewTag->ModsUsed == MAX_PT_MODS)
                    {
                        printf("** Error tagging scan %d from file %s: Too many PTMs!\n", Node->ScanNumber, Node->InputFile->FileName);
                        break;
                    }
                    ModBuffer[ModBufferPos] = '\0';
                    Delta = FindPTModByName(NewTag->Tag[AminoIndex - 1], ModBuffer);
                    if (Delta)
                    {
                        NewTag->AminoIndex[NewTag->ModsUsed] = AminoIndex - 1;
                        NewTag->ModType[NewTag->ModsUsed] = Delta;
                        NewTag->ModsUsed++;
                    }
                    else
                    {
                        printf("** Error tagging scan %d from file %s: Modification %s not understood!\n", Node->ScanNumber, Node->InputFile->FileName, ModBuffer);
                    }
                }

                //strncpy(NewTag->Tag, BitB, MAX_TAG_LENGTH);
                NewTag->Charge = Node->Tweaks[TweakIndex].Charge;
                NewTag->ParentMass = Node->Tweaks[TweakIndex].ParentMass;
                NewTag->PSpectrum = Node->Spectrum;
                NewTag->Tweak = Node->Tweaks;
                NewTag->PrefixMass = (int)(PrefixMass * MASS_SCALE + 0.5);
                NewTag->SuffixMass = Node->Spectrum->ParentMass - NewTag->PrefixMass;
                NewTag->Score = Probability;
                for (TempAA = NewTag->Tag; *TempAA; TempAA++)
                {
                    NewTag->SuffixMass -= PeptideMass[*TempAA];
                }
                for (ModIndex = 0; ModIndex < MAX_PT_MODS; ModIndex++)
                {
                    if (NewTag->AminoIndex[ModIndex] >= 0 && NewTag->ModType[ModIndex])
                    {
                        NewTag->SuffixMass -= NewTag->ModType[ModIndex]->RealDelta;
                    }
                }
                TotalTagCount++;
            } // Handle a line AFTER the tag header
        } // Loop over file lines
        fclose(TempTagOutputFile);
    } // Loop over tweaks
    qsort(_TagsOnlyTagList, TotalTagCount, sizeof(TrieTag), (QSortCompare)CompareTagScores);
    TotalTagCount = MergeIdenticalTags(_TagsOnlyTagList, TotalTagCount);
    TotalTagCount = min(TotalTagCount, GlobalOptions->GenerateTagCount);
    for (TagIndex = 0; TagIndex < TotalTagCount; TagIndex++)
    {
        NewTag = _TagsOnlyTagList + TagIndex;
        if (GlobalOptions->RunMode & RUN_MODE_TAGS_ONLY)
        {
            fprintf(GlobalOptions->OutputFile, "%s\t", Node->InputFile->FileName);
            fprintf(GlobalOptions->OutputFile, "%d\t", Node->ScanNumber);
            fprintf(GlobalOptions->OutputFile, "%d\t", Node->FilePosition);
            fprintf(GlobalOptions->OutputFile, "%d\t", TagIndex);
            fprintf(GlobalOptions->OutputFile, "%.2f\t", NewTag->PrefixMass / (float)MASS_SCALE);
            WriteTagToString(NewTag, AnnotationBuffer, 1);
            fprintf(GlobalOptions->OutputFile, "%s\t", AnnotationBuffer);
            fprintf(GlobalOptions->OutputFile, "%.2f\t", NewTag->SuffixMass / (float)MASS_SCALE);
            fprintf(GlobalOptions->OutputFile, "%.2f\t", NewTag->Score);
            fprintf(GlobalOptions->OutputFile, "\n");
        }
        else
        {
            Root = AddTagToTrie(Root, NewTag, &DuplicateFlag);
        }
    }
    return Root;
}

int SearchSpectrumBlockMSAlignment(SearchInfo* Info, SpectrumNode* FirstBlockSpectrum, 
    SpectrumNode* LastBlockSpectrum, DatabaseFile* DB)
{
    SpectrumNode* BlockSpectrum;
    int SpectraSearched = 0;
    int TweakIndex;
    //
    for (BlockSpectrum = FirstBlockSpectrum; BlockSpectrum != LastBlockSpectrum; BlockSpectrum = BlockSpectrum->Next)
    {
        if (!BlockSpectrum->Spectrum)
        {
            continue;
        }
        WindowFilterPeaks(BlockSpectrum->Spectrum, 0, 0);
        IntensityRankPeaks(BlockSpectrum->Spectrum);
        if (!BlockSpectrum->PMCFlag)
        {
            TweakSpectrum(BlockSpectrum);
        }

        printf("Unrestrictive search spectrum '%s' %d\n", BlockSpectrum->InputFile->FileName, BlockSpectrum->ScanNumber);
        fflush(stdout);

        if (!BlockSpectrum->Spectrum)
        {
            continue;
        }
        Info->Spectrum = BlockSpectrum->Spectrum;
        MutationModeSearch(Info);
        MQScoreSpectralMatches(BlockSpectrum);
        //OutputMatchesForSpectrum(BlockSpectrum, GlobalOptions->OutputFile);
        //FreeMatchList(BlockSpectrum);
        // Free PRM scores:
        for (TweakIndex = 0; TweakIndex < TWEAK_COUNT; TweakIndex++)
        {
            if (BlockSpectrum->Tweaks[TweakIndex].PRMScores)
            {
                SafeFree(BlockSpectrum->Tweaks[TweakIndex].PRMScores);
                BlockSpectrum->Tweaks[TweakIndex].PRMScores = NULL;
            }
        }
        if (BlockSpectrum->Spectrum->Graph)
        {
            FreeTagGraph(BlockSpectrum->Spectrum->Graph);
            BlockSpectrum->Spectrum->Graph = NULL;
        }
        SpectraSearched++;
    }
    return SpectraSearched;
}

int SearchSpectrumBlockTrie(SearchInfo* Info, SpectrumNode* FirstBlockSpectrum, SpectrumNode* LastBlockSpectrum, DatabaseFile* DB)
{
    int TagCount;
    SpectrumNode* BlockSpectrum;
    char TagBuffer[256];
    int SpectraSearched = 0;
    //
    // Construct tags for these spectra, and scan with trie:
    TagCount = GlobalOptions->GenerateTagCount;
    for (BlockSpectrum = FirstBlockSpectrum; BlockSpectrum != LastBlockSpectrum; BlockSpectrum = BlockSpectrum->Next)
    {
        if (!BlockSpectrum->Spectrum)
        {
            continue;
        }
        SpectraSearched++;
        // We perform peak filtering AFTER calling the external tagger. 
        if (GlobalOptions->ExternalTagger)
        {
            Info->Root = ConstructTagsExternalTagger(Info->Root, BlockSpectrum, TagCount);
        }
        WindowFilterPeaks(BlockSpectrum->Spectrum, 0, 0);
        IntensityRankPeaks(BlockSpectrum->Spectrum);
        if (!BlockSpectrum->PMCFlag)
        {
            TweakSpectrum(BlockSpectrum);
        }

        if (!GlobalOptions->ExternalTagger)
        {
            Info->Root = ConstructTagsForSpectrum(Info->Root, BlockSpectrum, TagCount);
        }
    }
    // Special case for external tagger, tags only:
    if (GlobalOptions->ExternalTagger && (GlobalOptions->RunMode & RUN_MODE_TAGS_ONLY))
    {
        //
    }
    else
    {
        memset(TagBuffer, 0, sizeof(char)*256);
        InitializeTrieFailureNodes(Info->Root, Info->Root, TagBuffer);
        //printf("Scan file with trie...\n");
        fflush(stdout);
        fseek(Info->DB->DBFile, 0, 0);
        switch (DB->Type)
        {
        case evDBTypeMS2DB:
            SearchMS2DB(Info);
            break;
        case evDBTypeSpliceDB:
            SearchSplicableGenes(Info);
            break;
        case evDBTypeTrie:
            ScanFileWithTrie(Info);
            break;
        default:
            break;
        }
    }
    return SpectraSearched;
}
// Return number of spectra searched
int SearchSpectrumBlockAgainstDB(SpectrumNode* FirstBlockSpectrum, SpectrumNode* LastBlockSpectrum, DatabaseFile* DB)
{
    SearchInfo* Info;
    int SpectraSearched;
    //
    Info = (SearchInfo*)calloc(1, sizeof(SearchInfo));
    Info->DB = DB;

    // MutationMode search is 'unrestricted, but not blind' mode.
    if (GlobalOptions->RunMode & (RUN_MODE_MUTATION | RUN_MODE_BLIND))
    {
        SpectraSearched = SearchSpectrumBlockMSAlignment(Info, FirstBlockSpectrum, LastBlockSpectrum, DB);
    }
    else
    {
        SpectraSearched = SearchSpectrumBlockTrie(Info, FirstBlockSpectrum, LastBlockSpectrum, DB);
        if (GlobalOptions->RunMode & RUN_MODE_BLINDTAG)
        {
            DebugPrintBlindTagExtensions(Info);
        }
    }
    FreeTrieNode(Info->Root);
    Info->Root = NULL;
    return SpectraSearched;
}

// I want to see some basic information about the onesided tag extension
// like the number of DB hits/tag and the number of onesided extends/tag
//SpectrumFileName-ScanCount-prefixMass-prefixExtends-Tag-TagHits-SuffixMass-SuffixExtends-score
void DebugPrintBlindTagExtensions(SearchInfo* Info)
{
    TrieNode* Root = Info->Root;
    TrieNode* L1 = NULL;
    TrieNode* L2 = NULL;
    TrieNode* L3 = NULL;
    TrieTagHanger* Hanger = NULL;
    TrieTag* Tag = NULL;
    int LevelOneKids;
    int LevelTwoKids;
    int LevelThreeKids;
    SpectrumNode* SNode = NULL;
    FILE* OutputFile;
    
    printf("I GOT TO THE DEBUG\n");
    OutputFile = fopen("BlindTaggingInfo.txt", "wb");
    //fprintf(OutputFile, "SpectrumFileName\tScanCount\tPrefixMass\tPrefixExtends\t");
    //fprintf(OutputFile, "Tag\tTagHits\tSuffixMass\tSuffixExtends\tScore\n");
    if (!OutputFile)
    {
        printf("Unable to upen the output file. BlindTaggingInfo\n");
        return;
    }
    for (LevelOneKids = 0; LevelOneKids < TRIE_CHILD_COUNT; LevelOneKids++)
    {
        //every node here is of depth 1, and has a single letter word
        if (LevelOneKids == ('I'-'A') || LevelOneKids == ('Q'-'A'))
        { //don't print out both nodes for I and L, or for Q and K
            continue;
        }
        L1 = Root->Children[LevelOneKids];
        if (L1 != NULL)
        {
            for (LevelTwoKids = 0; LevelTwoKids < TRIE_CHILD_COUNT; LevelTwoKids++)
            {
                if (LevelTwoKids == ('I'-'A') || LevelTwoKids == ('Q'-'A'))
                {
                    continue;
                }
                L2 = L1->Children[LevelTwoKids];
                if(L2 != NULL)
                {
                    for (LevelThreeKids = 0; LevelThreeKids < TRIE_CHILD_COUNT; LevelThreeKids++)
                    {
                        if (LevelThreeKids == ('I'-'A') || LevelThreeKids == ('Q'-'A'))
                        {
                            continue;
                        }
                        L3 = L2->Children[LevelThreeKids];
                        if (L3 != NULL)
                        {
                            //Level three kids should be a tripeptide, with a hanger and tags
                            //printf("My depth is %d\n",L3->Depth);
                            Hanger = L3->FirstTag;
                            while (Hanger) // != NULL; Go through all the hangers on a tag
                            {
                                Tag = Hanger->Tag;
                                SNode = Tag->PSpectrum->Node;
                                fprintf(OutputFile, "%s\t",SNode->InputFile->FileName);
                                fprintf(OutputFile, "%d\t",SNode->ScanNumber);
                                fprintf(OutputFile, "%.2f\t", Tag->PrefixMass / (float)MASS_SCALE);
                                fprintf(OutputFile, "%d\t",Tag->PrefixExtends);
                                fprintf(OutputFile, "%s\t", Tag->Tag);
                                fprintf(OutputFile, "%d\t", Tag->DBTagMatches);
                                fprintf(OutputFile, "%.2f\t", Tag->SuffixMass / (float)MASS_SCALE);
                                fprintf(OutputFile, "%d\t",Tag->SuffixExtends);
                                fprintf(OutputFile, "%.2f\t", Tag->Score);
                                fprintf(OutputFile, "\n");
                                fflush(OutputFile);

                                Hanger = Hanger->Next;
                            }// while
                        }
                    } //Level three kids
                }
            }// level 2 kids
        }
    }// level one kids
    fclose(OutputFile);
}

char* FindExtension(char* FileName)
{
    char* ExtensionString;
    //
    ExtensionString = FileName + strlen(FileName);
    while (ExtensionString > FileName)
    {
        ExtensionString--;
        if (*ExtensionString == '.')
        {
            return ExtensionString;
        }
    }
    return FileName + strlen(FileName);
}

// Set DB->IndexFileName, based upon the FileName and database type.
void FindDatabaseIndexFile(DatabaseFile* DB)
{
    char* Extension;
    
    strcpy(DB->IndexFileName, DB->FileName);
    Extension = FindExtension(DB->IndexFileName);
    sprintf(Extension, ".index\0");
    DB->IndexFile = fopen(DB->IndexFileName, "rb");
    if (DB->IndexFile)
    {
        return;
    }
    sprintf(Extension, ".ms2index\0");
    DB->IndexFile = fopen(DB->IndexFileName, "rb");
    if (DB->IndexFile)
    {
        return;
    }
    // No index file; that's ok.
    return;
}

// Search all our SpectrumNodes, one block at a time.
// Once the search is complete, compute p-values and output search results.
void RunSearch() 
{
    int BlockSize;
    SpectrumNode* FirstBlockSpectrum;
    SpectrumNode* LastBlockSpectrum;
    SpectrumNode* BlockSpectrum;
    FILE* SpectrumFile;
    int Result;
    TrieNode* Root = NULL;
    int SpectraSearched = 0;
    int ThisBlockSpectraSearched;
    DatabaseFile* DB;

    // Find index filenames:
    for (DB = GlobalOptions->FirstDatabase; DB; DB = DB->Next)
    {
        FindDatabaseIndexFile(DB);
    }
    // Open database files:
    for (DB = GlobalOptions->FirstDatabase; DB; DB = DB->Next)
    {
        if (!DB->DBFile)
        {
            DB->DBFile = fopen(DB->FileName, "rb");
        }
        if (!DB->IndexFile)
        {
            DB->IndexFile = fopen(DB->IndexFileName, "rb");
        }
    }

    if (GlobalOptions->RunMode & (RUN_MODE_MUTATION | RUN_MODE_BLIND))
    {
        GlobalOptions->TrieBlockSize = 5;  
        GlobalOptions->StoreMatchCount = 100; 
        GlobalOptions->ReportMatchCount = 10; // in production report at MOST 20 even in blind mode
    }
    BuildDecorations();

    FirstBlockSpectrum = GlobalOptions->FirstSpectrum;
    while (FirstBlockSpectrum)
    {
        fflush(stdout);

        // Load one block of spectrum objects:
        BlockSize = 0;
        LastBlockSpectrum = FirstBlockSpectrum;
        for (BlockSize = 0; BlockSize < GlobalOptions->TrieBlockSize; BlockSize++)
        {
            fflush(stdout);
    
            SpectrumFile = fopen(LastBlockSpectrum->InputFile->FileName, "rb");
            if (SpectrumFile)
            {
                LastBlockSpectrum->Spectrum = (MSSpectrum*)calloc(1, sizeof(MSSpectrum));
                LastBlockSpectrum->Spectrum->Node = LastBlockSpectrum;
                fseek(SpectrumFile, LastBlockSpectrum->FilePosition, 0);
                Result = SpectrumLoadFromFile(LastBlockSpectrum->Spectrum, SpectrumFile);
                //printf("Load from '%s' result %d\n", LastBlockSpectrum->InputFile->FileName, Result);
                fclose(SpectrumFile);
                if (!Result)
                {
                    SafeFree(LastBlockSpectrum->Spectrum);
                    LastBlockSpectrum->Spectrum = NULL;
                }
            }
            LastBlockSpectrum = LastBlockSpectrum->Next;
            if (!LastBlockSpectrum)
            {
                BlockSize++;
                break;
            }
        }
        printf("Search block of %d spectra starting with %s:%d\n", BlockSize, FirstBlockSpectrum->InputFile->FileName, FirstBlockSpectrum->ScanNumber);
        fflush(stdout);
        ThisBlockSpectraSearched = 0;
        for (DB = GlobalOptions->FirstDatabase; DB; DB = DB->Next)
        {
            ThisBlockSpectraSearched = SearchSpectrumBlockAgainstDB(FirstBlockSpectrum, LastBlockSpectrum, DB);
        }
        SpectraSearched += ThisBlockSpectraSearched;
        printf("Search progress: %d / %d (%.2f%%)\n", SpectraSearched, GlobalOptions->SpectrumCount, 100 * SpectraSearched / (float)max(1, GlobalOptions->SpectrumCount));
        fflush(stdout);

        // Clean up this block, and move to the next:
        fflush(stdout);
        for (BlockSpectrum = FirstBlockSpectrum; BlockSpectrum != LastBlockSpectrum; BlockSpectrum = BlockSpectrum->Next)
        {
            if (BlockSpectrum->Spectrum)
            {
                OutputMatchesForSpectrum(BlockSpectrum, GlobalOptions->OutputFile);
                FreeSpectrum(BlockSpectrum->Spectrum);
                FreeMatchList(BlockSpectrum);
                BlockSpectrum->Spectrum = NULL;
            }
        }
        fflush(stdout);

        FreeTrieNode(Root);
        Root = NULL;
        FirstBlockSpectrum = LastBlockSpectrum;
        fflush(stdout);
    }
    ///////////////////////////////////////////////////////////
    // After searching, we compute p-values and output matches:
    fclose(GlobalOptions->OutputFile);
    GlobalOptions->OutputFile = NULL;

    /// Compute p-values, write them to final output file:
    if (GlobalOptions->ExternalTagger && (GlobalOptions->RunMode & RUN_MODE_TAGS_ONLY))
    {
        // do nothing
    }
    else
    {
        CalculatePValues(GlobalOptions->OutputFileName, GlobalOptions->FinalOutputFileName);
    }

    // Close database files:
    for (DB = GlobalOptions->FirstDatabase; DB; DB = DB->Next)
    {
        if (DB->DBFile)
        {
            fclose(DB->DBFile);
            DB->DBFile = NULL;
        }
        if (DB->IndexFile)
        {
            fclose(DB->IndexFile);
            DB->IndexFile = NULL;
        }
    }
}

// Special run mode: Perform parent mass correction on our input spectra.  Output the
// parent masses and charge states.
void PerformSpectrumTweakage()
{
    SpectrumNode* Node;
    FILE* SpectrumFile;
    int TweakIndex;
    int Result;
    //
    for (Node = GlobalOptions->FirstSpectrum; Node; Node = Node->Next)
    {
        SpectrumFile = fopen(Node->InputFile->FileName, "rb");
        fseek(SpectrumFile, Node->FilePosition, 0);
        if (SpectrumFile)
        {
            Node->Spectrum = (MSSpectrum*)calloc(1, sizeof(MSSpectrum));
            Node->Spectrum->Node = Node;
            Result = SpectrumLoadFromFile(Node->Spectrum, SpectrumFile);
            fclose(SpectrumFile);
            if (!Result)
            {
                FreeSpectrum(Node->Spectrum);
                Node->Spectrum = NULL;
            }
            else
            {
                WindowFilterPeaks(Node->Spectrum, 0, 0);
                IntensityRankPeaks(Node->Spectrum);
                PrepareSpectrumForIonScoring(PRMModelCharge2, Node->Spectrum, 0);
                //SpectrumComputeBinnedIntensities(Node);
                Node->Spectrum->Node = Node;
                TweakSpectrum(Node);
                
                fprintf(GlobalOptions->OutputFile, "%s\t", Node->InputFile->FileName);
                fprintf(GlobalOptions->OutputFile, "%d\t", Node->ScanNumber);
                fprintf(GlobalOptions->OutputFile, "%d\t", Node->FilePosition);
                for (TweakIndex = 0; TweakIndex < TWEAK_COUNT; TweakIndex++)
                {
                    if (Node->Tweaks[TweakIndex].Charge)
                    {
                        fprintf(GlobalOptions->OutputFile, "%.2f\t%d\t", Node->Tweaks[TweakIndex].ParentMass / (float)DALTON, Node->Tweaks[TweakIndex].Charge);
                    }
                }
                fprintf(GlobalOptions->OutputFile, "\n");
            }
            FreeSpectrum(Node->Spectrum);
            Node->Spectrum = NULL;
        }
        
    }
}

//For phosphorylated spectra, the superprominent M-p peak can 
//fritz the charge state guessing, and tagging.  So we remove it.
void AttemptParentMassPeakRemoval(MSSpectrum* Spectrum)
{
    int MostIntensePeakIndex;
    int MostIntenseMass;
    int PeakIndex;
    float MostIntense = 0.0;
    float NextMostIntense = 0.0;
    int Diff;
    int ExpectedDiff;
    int ExpectedDiff2;
    int Epsilon = HALF_DALTON;
	int CalculatedMZ;
    //
    for (PeakIndex = 0; PeakIndex < Spectrum->PeakCount; PeakIndex++)
    {
        if (Spectrum->Peaks[PeakIndex].Intensity > MostIntense)
        {
            NextMostIntense = MostIntense;
            MostIntense = Spectrum->Peaks[PeakIndex].Intensity;
            MostIntensePeakIndex = PeakIndex;
            MostIntenseMass = Spectrum->Peaks[PeakIndex].Mass;
        }
        else if(Spectrum->Peaks[PeakIndex].Intensity > NextMostIntense)
        {
            NextMostIntense = Spectrum->Peaks[PeakIndex].Intensity;
        }
    }
    //printf("Most intense %f, next %f\n",MostIntense, NextMostIntense);
    //if more than 2 times great, and in the right place, remove peak.
  //  if (MostIntense < 2 * NextMostIntense)
  //  {
		//Spectrum->RemovedPeakIndex = -1;//dummy holder
  //      return;
  //  }
    //printf ("MZ of %d, charge %d\n", Spectrum->MZ, Spectrum->Charge);
	//Set m/z with the new parentmass and charge that was just assigned in ConstructTags
    CalculatedMZ = (Spectrum->ParentMass + (Spectrum->Charge - 1) * HYDROGEN_MASS) / Spectrum->Charge;
    Diff = abs(CalculatedMZ - MostIntenseMass);
    ExpectedDiff = PHOSPHATE_WATER_MASS / Spectrum->Charge;
    ExpectedDiff2 = (PHOSPHATE_WATER_MASS + WATER_MASS) / Spectrum->Charge;
    if (abs (Diff - ExpectedDiff) < Epsilon)
    { //remove peak
        Spectrum->RemovedPeakIndex = MostIntensePeakIndex;
        Spectrum->RemovedPeakIntensity = Spectrum->Peaks[MostIntensePeakIndex].Intensity;
        Spectrum->Peaks[MostIntensePeakIndex].Intensity = 1.0; //cut to ground
    }
    else if (abs(Diff - ExpectedDiff2) < Epsilon)
    { //remove peak
        Spectrum->RemovedPeakIndex = MostIntensePeakIndex;
        Spectrum->RemovedPeakIntensity = Spectrum->Peaks[MostIntensePeakIndex].Intensity;
        Spectrum->Peaks[MostIntensePeakIndex].Intensity = 1.0; //cut to ground
    }
	else
	{
		Spectrum->RemovedPeakIndex = -1;//dummy holder
	}
}

void RestoreParentMassPeakRemoval(MSSpectrum* Spectrum)
{
	if (Spectrum->RemovedPeakIndex == -1)
	{
		return;
	}
	Spectrum->Peaks[Spectrum->RemovedPeakIndex].Intensity = Spectrum->RemovedPeakIntensity;
}