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


// Mod-tolerant matching of peptides to spectra.  See header file FreeMod.h for overview.
#include "CMemLeak.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Utils.h"
#include "Inspect.h"
#include "Trie.h"
#include "Mods.h"
#include "Tagger.h"
#include "Score.h"
#include "FreeMod.h"
#include "Scorpion.h"
#include "SVM.h"
#include "IonScoring.h"

// SkewPenalty[n] is a score penalty applied to a node that is n/100 daltons
// away from where it should be.  Size 100.  (Derived from a functional
// fit to empirical histogram)
int g_SkewPenalty[] = {0, 0, 0, 0, -1, -2, -2, -3, -4, -4, -4, -4, -4, -5,
        -5, -6, -7, -8, -8, -9, -9, -9, -9, -9, -10, -11, -11, -12, -12,
        -12, -12, -12, -12, -13, -13, -14, -14, -15, -15, -15, -15, -15,
        -15, -15, -16, -16, -16, -17, -17, -17, -17, -17, -17, -17, -18,
        -18, -18, -19, -19, -19, -19, -19, -19, -19, -20, -20, -20, -21,
        -21, -22, -22, -22, -22, -22, -22, -23, -24, -24, -25, -25, -25,
        -25, -25, -26, -26, -28, -29, -30, -31, -32, -32, -32, -32, -34,
        -35, -39, -41, -48, -57, -65};

int g_SkewPenaltySize = sizeof(g_SkewPenalty) / sizeof(int);
int g_SkewPenaltyMax = sizeof(g_SkewPenalty) / sizeof(int) - 1;

// MassDeltas[AminoAcid][n] is the nth modification (normally sorted by size) possible on AminoAcid. 
// e.g. MassDeltas[0][0] is the smallest (or most negative) modification allowed on alanine
MassDelta** MassDeltas = NULL;

// MassDeltaByMass[AminoAcid][PRM] is a linked list of MassDeltaNodes corresponding to MassDeltas matching 
// PRM on amino acid.  
MassDeltaNode*** MassDeltaByMass = NULL;

// For user-supplied (limited) PTMs, MassDeltaByIndex[AA*MAX_PT_MODTYPE + n] points to an instance of the PTM with 
// index n, attached to AA.  Used for decorations!  A special case: Always, for AA of 26 (MDBI_ALL_MODS), store 
// a pointer to a valid PTM.
// The array AllKnownPTMods holds *one* entry for all modifications of the same type (e.g. all phosphorylations),
// but there's a separate MassDelta instance for serine-phos, threonine-phos, and tyrosine-phos.  (And this is probably
// as it should be, since we might attach a different penalty to phosphotyrosine than to phosphoserine, to reflect the
// fact that serines are more commonly phosphorylated)
MassDelta** MassDeltaByIndex = NULL;

/////////////////////////////////////////////////////////
// Forward declarations:
void DebugPrintMultiModTable(TagGraph* Graph, char* Buffer, int MaxX, int MaxY, int MaxZ);
void AddMultiModMatch(MSSpectrum* Spectrum, int CellIndex, int Bonus,
    char* Buffer, int StartPos, int ModBlockSize, int AminoBlockSize, int BonusLength,
    int BufferEnd);
int ExtendMatchRightwardDuo(SearchInfo* Info, char* Buffer, int BufferEnd, int MatchMass, 
    int MaxMods, int ScoreToBeat, int FilePos, SpectrumTweak* Tweak);

int MSAlignmentGeneral(SearchInfo* Info, char* Buffer, int BufferEnd, int MatchMass, 
    int MaxMods, int ScoreToBeat, int FilePos, SpectrumTweak* Tweak);

void AddNewMatchDuo(SearchInfo* Info, SpectrumTweak* Tweak, char* Buffer, int Score, int* PrevCellTable, MassDelta** DeltaTable, 
    int CellIndex, MassDelta* FinalDelta, int AminoBlockSize, int AminoIndex, int EndAminoIndex,
    int FilePos);

void DebugPrintPRMScores(MSSpectrum* Spectrum, SpectrumTweak* Tweak);


// Allocate the array MassDeltaByIndex, if it's not already allocated.
void AllocMassDeltaByIndex()
{
    int MallocSize;
    if (MassDeltaByIndex)
    {
        return;
    }
    MallocSize = (MAX_PT_MODTYPE * (AMINO_ACIDS + 1)) * sizeof(MassDelta*);
    MassDeltaByIndex = (MassDelta**)malloc(MallocSize);
}

// Free the 2-dimensional table MassDeltaByMass.  It's a big table, so don't forget to free it :)
void FreeMassDeltaByMass()
{
    int AA;
    int PRM;
    MassDeltaNode* Node;
    MassDeltaNode* Prev = NULL;
    if (MassDeltaByMass)
    {
        for (AA = 0; AA < AMINO_ACIDS; AA++)
        {
            for (PRM = 0; PRM < GlobalOptions->DeltaBinCount; PRM++)
            {
                // MassDeltaByMass[AA][PRM] is either null, or it points to the head of a 
                // linked list of MassDeltaNode objects.  
                if (MassDeltaByMass[AA][PRM])
                {
                    // Free each node of the list:
                    Node = MassDeltaByMass[AA][PRM];
                    Prev = NULL;
                    while (Node)
                    {
                        SafeFree(Prev);
                        Prev = Node;
                        Node = Node->Next;
                    }
                    SafeFree(Prev);
                }
            }
            SafeFree(MassDeltaByMass[AA]);
        }
        SafeFree(MassDeltaByMass);
        MassDeltaByMass = NULL;
    }
}

// Free all the mods in MassDeltas array.
void FreeMassDeltas()
{
    int AA;
    if (MassDeltas)
    {
        for (AA = 0; AA < AMINO_ACIDS; AA++)
        {
            SafeFree(MassDeltas[AA]);
        }
        SafeFree(MassDeltas);
        MassDeltas = NULL;
    }
    FreeMassDeltaByMass();
    //SafeFree(MassDeltaByIndex);
    //MassDeltaByIndex = NULL;
}

// Initialize the hash MassDeltaByMass.  The table entry MassDeltaByMass[AA][Delta] points to a linked list of
// mass deltas for amino acid AA matching Delta.  
// In some cases, it makes sense to consider two mass deltas of the same size. Example: Mutation
// to Q or to K.  We keep a *list* of mass deltas in all cases. 
void InitMassDeltaByMass()
{
    int MassDeltaIndex;
    int Fudge;
    int Delta;
    int AA;
    // We populate adjacent cells in MassDeltas[AA] as well.  FudgeMax == how many bins away from the "right" bin to consider.
    // FudgeMax is usually 1, to handle roundoff error.  But FudgeMax can be 2-3 if the parent mass epsilon is quite large.
    int FudgeMax = 1 + (GlobalOptions->ParentMassEpsilon / DALTON);
    MassDeltaNode* OldNode;
    MassDeltaNode* NewNode;
    //
    FreeMassDeltaByMass();
    MassDeltaByMass = (MassDeltaNode***)calloc(AMINO_ACIDS, sizeof(MassDeltaNode**));
    for (AA = 0; AA < AMINO_ACIDS; AA++)
    {
        MassDeltaByMass[AA] = (MassDeltaNode**)calloc(GlobalOptions->DeltaBinCount + 1, sizeof(MassDeltaNode**));
        for (MassDeltaIndex = 0; MassDeltaIndex < GlobalOptions->DeltaBinCount; MassDeltaIndex++)
        {
            if (!MassDeltas[AA][MassDeltaIndex].Flags)
            {
                // Null array entry.
                break;
            }
            ROUND_MASS_TO_DELTA_BIN(MassDeltas[AA][MassDeltaIndex].RealDelta, Delta);
            // Add our MassDelta to the bin (and to neighboring bins), either filling the bin
            // or adding the new MassDelta to the end of the bin's linked list of MassDeltaNodes:
            for (Fudge = max(0, Delta - FudgeMax); Fudge < min(GlobalOptions->DeltaBinCount, Delta + FudgeMax + 1); Fudge++)
            {
                NewNode = (MassDeltaNode*)calloc(1, sizeof(MassDeltaNode));
                NewNode->Delta = &MassDeltas[AA][MassDeltaIndex];
                //NewNode->RealDelta = NewNode->Delta->RealDelta;
                OldNode = MassDeltaByMass[AA][Fudge];
                if (!OldNode)
                {
                    MassDeltaByMass[AA][Fudge] = NewNode;
                }
                else
                {
                    while (OldNode->Next)
                    {
                        OldNode = OldNode->Next;
                    }
                    OldNode->Next = NewNode;
                }
            }
        }
    }
}

// Read, from the binary file Mutations.dat, the definitions of all mass modifications we will consider.
// (It's faster to consider a large but LIMITED set of modifications than to consider every feasible value
// of delta.  Also, this limited set lets us assign a SCORE and a NAME to each delta, which is very useful)
// The file Mutations.dat is written out by the scaffold script PrepBlosum.py 
// Any mass delta with Flags == 0 is a dummy record, which is included simply to pad out the 
// array to a uniform size; such deltas should *never* actually be used! 
// If ReadFlag is false, then don't actually read anything from a file - just init the structure.
// After calling this, the caller should also call InitMassDeltaByMass to init the hash.
void LoadMassDeltas(char* FileName, int ReadFlag)
{
    int AA;
    int DeltaIndex;
    FILE* MassDeltaFile;

    FreeMassDeltas(); // Free up any pre-conceived notions

    MassDeltas = (MassDelta**)calloc(AMINO_ACIDS, sizeof(MassDelta*));
    for (AA = 0; AA < AMINO_ACIDS; AA++)
    {
        MassDeltas[AA] = (MassDelta*)calloc(GlobalOptions->DeltasPerAA, sizeof(MassDelta));
    }
    if (!ReadFlag)
    {
        // That was a freebie.
        return;
    }
    if (!FileName || !*FileName)
    {
        // No file to open.
        return;
    }
    MassDeltaFile = fopen(FileName, "rb");
    if (!MassDeltaFile)
    {
        printf("Error: Unable to open mutation data file '%s'", FileName);
        return;
    }
    for (AA = 0; AA < AMINO_ACIDS; AA++)
    {
        MassDeltas[AA] = (MassDelta*)calloc(GlobalOptions->DeltasPerAA, sizeof(MassDelta));
        for (DeltaIndex = 0; DeltaIndex < GlobalOptions->DeltasPerAA; DeltaIndex++)
        {
            ReadBinary(&MassDeltas[AA][DeltaIndex].RealDelta, sizeof(int), 1, MassDeltaFile);
            MassDeltas[AA][DeltaIndex].Delta = MASS_TO_BIN(MassDeltas[AA][DeltaIndex].RealDelta);
            ReadBinary(&MassDeltas[AA][DeltaIndex].Score, sizeof(float), 1, MassDeltaFile);
            ReadBinary(&MassDeltas[AA][DeltaIndex].Name, sizeof(char), 20, MassDeltaFile);
            ReadBinary(&MassDeltas[AA][DeltaIndex].Flags, sizeof(int), 1, MassDeltaFile);
            ReadBinary(&MassDeltas[AA][DeltaIndex].Amino, sizeof(char), 1, MassDeltaFile);
        }
    }
    fclose(MassDeltaFile);
    // The caller should now invoke InitMassDeltaByMass()
}

// Enrich MassDeltas[] to include one modification for any (reasonable) mass change applicable to any
// amino acid.  Also, update InitMassDeltaByMass()
void AddBlindMods()
{
    int AA;
    int DeltaMass;
    int Bin;
    int FoundFlag;
    int Index;
    int MaxDeltaMass;
    //
    for (AA = 0; AA < AMINO_ACIDS; AA++)
    {
        DeltaMass = PeptideMass['A' + AA];
        if (!DeltaMass)
        {
            continue; // bogus amino like B or Z
        }
        DeltaMass = (DeltaMass / 1000) * 1000;
        // The largest *negative* modification permitted is one that takes us down to the mass of glycine:
        DeltaMass = GLYCINE_MASS - DeltaMass;
        MaxDeltaMass = GlobalOptions->MaxPTMDelta * MASS_SCALE;
        while (DeltaMass < MaxDeltaMass)
        {
            // Don't add a mutation for mass delta ~0:
            if (abs(DeltaMass) < MASS_SCALE)
            {
                DeltaMass += MASS_SCALE;
                continue;
            }
            ROUND_MASS_TO_DELTA_BIN(DeltaMass, Bin);
            FoundFlag = 0;
            // If we already know a PTM that matches this mass closely enough, don't add another:
            for (Index = 0; Index < GlobalOptions->DeltasPerAA; Index++)
            {
                if (!MassDeltas[AA][Index].Flags)
                {
                    break;
                }
                if (abs(MassDeltas[AA][Index].RealDelta - DeltaMass) < HALF_DALTON)
                {
                    FoundFlag = 1;
                    break;
                }
            }
            if (!FoundFlag)
            {
                MassDeltas[AA][Index].RealDelta = DeltaMass;
                MassDeltas[AA][Index].Delta = Bin;
                MassDeltas[AA][Index].Flags = 1;
                MassDeltas[AA][Index].Score = -1; // Somewhat magical score!
                sprintf(MassDeltas[AA][Index].Name, "%+d", DeltaMass / MASS_SCALE);
            }
            DeltaMass += DALTON;
        }
    }
    InitMassDeltaByMass();
}

// For development only: Print out the scores of all PRMs, as well as the b and y scores
// and witness scores from which the PRMScores were derived.  Requires some slow
// business to set Spectrum->PRMDebugStrings, and so is not enabled in normal builds.
void DebugPrintPRMScores(MSSpectrum* Spectrum, SpectrumTweak* Tweak)
{
    FILE* PRMFile = NULL;
    int PRM;
    ///
    PRMFile = fopen("PRMScores.xls", "w");
    if (!PRMFile)
    {
        printf("NO DEBUG PRINT OF PRM SCORES DONE.\n");
        return;
    }
    fprintf(PRMFile, "#PRM\tMass\tScore\tBScore\tYScore\tWitnessScore\n");
    for (PRM = 0; PRM < Tweak->PRMScoreMax; PRM++)
    {
        fprintf(PRMFile, "%d\t%.2f\t%d\t", PRM, PRM / 10.0, Tweak->PRMScores[PRM]);
        //fprintf(PRMFile, "%s\n", Spectrum->PRMDebugStrings[PRM]);
        fprintf(PRMFile, "\n");
    }
    fclose(PRMFile);
}

// When doing 2-mutant extension, MAX_RIGHT_EXTENSIONS needs to be large:
#define MAX_RIGHT_EXTENSIONS 512
Peptide LeftExtensions[MAX_RIGHT_EXTENSIONS];
int LeftExtensionCount;
Peptide RightExtensions[MAX_RIGHT_EXTENSIONS];
int RightExtensionCount;

Peptide* Add1ModMatch(SearchInfo* Info, char* Buffer, int BufferLength, int SuffixEndPos, int SuffixStartPos, int PrefixLength, 
    int Score, MassDelta* Delta, SpectrumTweak* Tweak, int FilePos, char ExtraPrefixChar)
{
    Peptide* Match;
    int Length;
    int Pos;
    MSSpectrum* Spectrum = Info->Spectrum;
    //
    Length = SuffixEndPos - SuffixStartPos + PrefixLength + 1;
    Match = NewPeptideNode();
    Match->Tweak = Tweak;
    strncpy(Match->Bases, Buffer + SuffixStartPos - PrefixLength, Length);

    Match->InitialScore = Score;
    Match->RecordNumber = Info->RecordNumber;
    Match->FilePos = FilePos + SuffixStartPos - PrefixLength;
    if (SuffixStartPos - PrefixLength > 0)
    {
        Match->PrefixAmino = Buffer[SuffixStartPos - PrefixLength - 1];
    }
    else
    {
        Match->PrefixAmino = ExtraPrefixChar;
    }
    Pos = SuffixStartPos - PrefixLength + Length;
    if (Pos < BufferLength)
    {
        Match->SuffixAmino = Buffer[Pos];
    }
    if (Delta)
    {
        Match->AminoIndex[0] = PrefixLength;
        Match->ModType[0] = Delta;
    }
    Match->DB = Info->DB;
    GetPeptideParentMass(Match);
    return StoreSpectralMatch(Spectrum, Match, Length, 0);
}

// SeekMatch1PTM performs a blind search with at most one PTM permitted.
// Schematic of SeekMatch1PTM:
//                SuffixStartPos = PrefixEndPos, PTM attaches here
//               /
//              / SuffixEndPos
//             / / 
//            * *  
//  IKKWLSLPGEMTRPLIL
//     *       
//      \--PrefixStartPos

// Kludge: If Buffer points to the middle of a long peptide,
// ExtraPrefixChar is the character that precedes Buffer.
#define MAX_1MOD_PEPTIDE_LEN 64

int SeekMatch1PTM(SearchInfo* Info, char* Buffer, int BufferLen, int MatchMass, int ScoreToBeat,
                  SpectrumTweak* Tweak, int FilePos, char ExtraPrefixChar)
{
    static int* PrefixScores = NULL;
    static int* PrefixMasses = NULL;
    int PrefixStartPos;
    int PrefixEndPos;
    int PRM;
    int PRMBin;
    int Score;
    int MatchScore;
    int MaxPrefix;
    int SkipBases;
    int ArrayIndex;
    int PrefixLength;
    int MaxPrefixLength;
    int MinPossibleDelta = -130000;
    int MaxPossibleDelta = GlobalOptions->MaxPTMDelta * MASS_SCALE;
    int MaxMass = MatchMass + GlobalOptions->ParentMassEpsilon + MaxPossibleDelta;
    int Delta;
    int DeltaBin;
    MassDeltaNode* DeltaNode;
    char AA;
    int AAIndex;
    int AbsSkew;
    int AAMass = 0;
    int Skew;
    int SuffixEndPos;
    int SuffixStartPos;
    int MatchScoreWithDelta;    
    Peptide* Match;
    MSSpectrum* Spectrum = Info->Spectrum;
    //
    if (!BufferLen)
    {
        return 1;
    }
    if (!PrefixScores)
    {
        PrefixScores = (int*)calloc(512 * MAX_1MOD_PEPTIDE_LEN, sizeof(int));
        PrefixMasses = (int*)calloc(512 * MAX_1MOD_PEPTIDE_LEN, sizeof(int));
    }
    // The prefix of our peptide will extend from 
    // PrefixStartPos...PrefixEndPos, NOT including PrefixEndPos
    
    // By default, we cover up to 450 bases in each call.  If we hit the end of a protein, we stop there
    // and handle the next protein in the next call to this function.
    SkipBases = min(BufferLen, 450); 
    for (PrefixStartPos = 0; PrefixStartPos < SkipBases; PrefixStartPos++)
    {
        if (Buffer[PrefixStartPos]=='*')
        {
            SkipBases = PrefixStartPos + 1;
            break;
        }
        PRM = 0;
        Score = 0;
        MaxPrefix = min(SkipBases, PrefixStartPos + MAX_1MOD_PEPTIDE_LEN);
        for (PrefixEndPos = PrefixStartPos; PrefixEndPos < MaxPrefix; PrefixEndPos++)
        {
            ArrayIndex = PrefixEndPos * MAX_1MOD_PEPTIDE_LEN + (PrefixEndPos - PrefixStartPos);
            if (ArrayIndex < 0 || ArrayIndex > 512 * MAX_1MOD_PEPTIDE_LEN)
            {
                printf("** error: Array index for prefix is %d\n", ArrayIndex);
            }
            PrefixScores[ArrayIndex] = Score;
            PrefixMasses[ArrayIndex] = PRM;
            //printf("%d: Prefix %d-%d score %d PRM %d\n", ArrayIndex, PrefixStartPos, PrefixEndPos, Score, PRM); //Verbose1Mod
            if (PRM > MaxMass)
            {
                break;
            }
            AAMass = PeptideMass[Buffer[PrefixEndPos]];
            if (AAMass)
            {
                PRM += AAMass;
            }
            else
            {
                Score = -9999999;
            }
            if (PRM > MaxMass)
            {
                // Modless prefix is too long!
                Score = -9999999;
                break;
            }
            else
            {
                PRMBin = (PRM + 50) / 100;
                if (PRMBin >= 0 && PRMBin < Tweak->PRMScoreMax)
                {
                    Score += Tweak->PRMScores[PRMBin];
                }
            }
        }
    }
    // Now that the prefix table's complete, consider all possible suffixes.
    // The suffix of our peptide will extend from SuffixStartPos...SuffixEndPos, INCLUSIVE.
    for (SuffixEndPos = SkipBases - 1; SuffixEndPos > 0; SuffixEndPos--)
    {
        //printf("Try ending at pos'n %d (%c)\n", SuffixEndPos, Buffer[SuffixEndPos]); //Verbose1Mod
        PRM = MatchMass;
        Score = 0;
        if (Spectrum->Node->LastMatch)
        {
            ScoreToBeat = Spectrum->Node->LastMatch->InitialScore;
        }
        for (SuffixStartPos = SuffixEndPos; SuffixStartPos >= 0; SuffixStartPos--)
        {   
            // Grow the c-terminal suffix by one residue:
            AA = Buffer[SuffixStartPos];
            AAIndex = AA - 'A';
            AAMass = PeptideMass[AA];
            if (AAMass)
            {
                PRM -= AAMass;
            }
            else
            {
                break; // bogus AA encountered
            }
            
            //printf("Suffix %d (%c) to %d (%c), mass remaining %.2f\n", SuffixStartPos, Buffer[SuffixStartPos],
            //    SuffixEndPos, Buffer[SuffixEndPos], PRM / 1000.0); //Verbose1Mod
            if (PRM < -GlobalOptions->ParentMassEpsilon)
            {
                break; // modless suffix is too long!
            }
            // Try to hook up a to prefix:
            ArrayIndex = SuffixStartPos * MAX_1MOD_PEPTIDE_LEN;
            MaxPrefixLength = min(MAX_1MOD_PEPTIDE_LEN, SuffixStartPos + 1);
            for (PrefixLength = 0; PrefixLength < MaxPrefixLength; PrefixLength++)
            {
                Delta = PRM - PrefixMasses[ArrayIndex];
                if (Delta < MinPossibleDelta)
                {
                    break;
                }
                if (Delta < MaxPossibleDelta)
                {
                    MatchScore = Score + PrefixScores[ArrayIndex];
                    //printf("Prefix %d-%d, suffix %d-%d, delta %.2f, score %d\n", SuffixStartPos - PrefixLength, SuffixStartPos, SuffixStartPos, SuffixEndPos, Delta / (float)DALTON, MatchScore);
                    if (MatchScore > ScoreToBeat)
                    {
                        // Look for the delta that can hook these together:
                        if (abs(Delta) < GlobalOptions->ParentMassEpsilon)
                        {
                            Match = Add1ModMatch(Info, Buffer, BufferLen, SuffixEndPos, SuffixStartPos, PrefixLength, MatchScore, NULL, Tweak, FilePos, ExtraPrefixChar);
                            //after every call to add a match, the ScoreToBeat MUST be updated.
                            ScoreToBeat = Spectrum->Node->LastMatch->InitialScore;
                        }
                        else
                        {
                            ROUND_MASS_TO_DELTA_BIN(Delta, DeltaBin);
                            DeltaNode = MassDeltaByMass[AAIndex][DeltaBin];
                            while (DeltaNode)
                            {
                                Skew = Delta - DeltaNode->Delta->RealDelta;
                                AbsSkew = abs(Skew);
                                if (AbsSkew <= GlobalOptions->Epsilon)
                                {                                    
                                    MatchScoreWithDelta = MatchScore + rint(DeltaNode->Delta->Score * DELTA_SCORE_SCALER);
                                    if (MatchScoreWithDelta > ScoreToBeat)
                                    {
                                        Match = Add1ModMatch(Info, Buffer, BufferLen, SuffixEndPos, SuffixStartPos, PrefixLength, MatchScoreWithDelta, DeltaNode->Delta, Tweak, FilePos, ExtraPrefixChar);
                                        ScoreToBeat = Spectrum->Node->LastMatch->InitialScore;
                                    }
                                }
                                DeltaNode = DeltaNode->Next;
                            }
                            // If the modification mass was small, then ALSO try the unmodified peptide:
                            if (abs(Delta) < 5 * DALTON)
                            {
                                Add1ModMatch(Info, Buffer, BufferLen, SuffixEndPos, SuffixStartPos, PrefixLength, MatchScore, NULL, Tweak, FilePos, ExtraPrefixChar);
                                ScoreToBeat = Spectrum->Node->LastMatch->InitialScore;

                            }
                        }
                    }
                }
                ArrayIndex++;
            }
            // If we didn't just link up, then accumulate some score:
            if (PRM >= 0)
            {
                Score += Tweak->PRMScores[MASS_TO_BIN(PRM)];
                //printf("Accumulate score %.2f from PRM %d\n", Tweak->PRMScores[MASS_TO_BIN(PRM)], PRM);
            }
        } // SuffixEndPos loop
    } // SuffixStartPos loop
    return SkipBases;

}

int* PTMScoreTable = NULL;
int* PrevCellIndexTable = NULL;
MassDelta** DeltaTable = NULL;
int* MassDeltaTable = NULL;

#define DB_BUFFER_SIZE 1024000
#define DB_SHUNT_BOUNDARY 900000
#define DB_READ_BOUNDARY 900000

// Search a database, using *no* tag-based filtering at all.  This is much slower than searching
// with tag-based filters, but also more sensitive, particularly since we haven't handled the
// problem of tagging in the presence of mutations.  
void SearchDatabaseTagless(SearchInfo* Info, int MaxMods, int VerboseFlag, SpectrumTweak* Tweak)
{
    static char* Buffer = NULL;  // will be big
    int IsEOF = 0;
    FILE* DBFile;
    int FilePos = 0;
    int BufferEnd = 0;
    int BufferPos = 0;
    int BytesRead;
    MSSpectrum* Spectrum = Info->Spectrum;
    // We require all peptide candidates to be long enough to be meaningful (at least 500 Da, or large
    // enough to equal the parent mass after maximum mod mass, whichever is largest)
    // We also stop considering peptide candidates after they are too long to match
    // the spectrum (even after deducting some mass due to modifications)
    int ParentResidueMass = Spectrum->ParentMass - PARENT_MASS_BOOST;
    int ScoreToBeat = -999999;
    int SkipBases;
    char PrefixChar;
    //
    if (!Buffer)
    {
        Buffer = (char*)malloc(sizeof(char) * DB_BUFFER_SIZE);
        if (!Buffer)
        {
            printf("** ERROR: Unable to allocate buffer in SearchDatabaseTagless()!\n");
            return;
        }
    }

    // Ensure that the PRM scores of this spectrum are set, so that we can score candidates:
    if (!Tweak->PRMScores)
    {
        if (VerboseFlag)
        {
            printf("[V] SetPRMScores()\n");
        }
        SetSpectrumPRMScores(Spectrum, Tweak);
    }
    // Open the database, and start reading:

    DBFile = Info->DB->DBFile;

    //DebugPrintPRMScores(Spectrum, Tweak); 
    Info->RecordNumber = 0;
    while (1)
    {
        if (VerboseFlag)
        {
            printf("[V] Bufferpos %d BufferEnd %d IsEOF %d FilePos %d Record# %d\n", BufferPos, BufferEnd, IsEOF, FilePos, Info->RecordNumber);
        }

        // Shunt bases toward front of buffer:
        if (BufferPos > DB_SHUNT_BOUNDARY)
        {
            memmove(Buffer, Buffer + BufferPos, DB_BUFFER_SIZE - BufferPos);
            BufferEnd = DB_BUFFER_SIZE - BufferPos;
            BufferPos = 0;
        }
        // Read more bases:
        if (!IsEOF && BufferEnd < DB_READ_BOUNDARY)
        {
            BytesRead = ReadBinary(Buffer + BufferEnd, sizeof(char), DB_BUFFER_SIZE - BufferEnd, DBFile);
            BufferEnd += BytesRead;
            if (!BytesRead)
            {
                IsEOF = 1;
            }
        }
        if (BufferPos >= BufferEnd) // hit the end if the database
        {
            break;
        }

        // If this isn't an amino acid, skip onward:
        if (!PeptideMass[Buffer[BufferPos]])
        {
            BufferPos++;
            FilePos++;
            continue;
        }
        if (Buffer[BufferPos]=='*')
        { 
            BufferPos++;
            FilePos++;
            Info->RecordNumber++;
            continue;
        }
        // Try to find peptide matches from a prefix of Buffer[BufferPos:]
        if (MaxMods > 2)
        {
            // The SLOW way!
            SkipBases = MSAlignmentGeneral(Info, Buffer + BufferPos, BufferEnd - BufferPos, ParentResidueMass, 
                MaxMods, ScoreToBeat, FilePos, Tweak);
            if (VerboseFlag)
            {
                printf("[V] General() return.  SkipBases %d\n", SkipBases);
            }

        }
        else if (MaxMods > 1)
        {
            // Extend into a match, possibly up to 2 mods.  
            SkipBases = ExtendMatchRightwardDuo(Info, Buffer + BufferPos, BufferEnd - BufferPos, ParentResidueMass, 
                min(2, MaxMods), ScoreToBeat, FilePos, Tweak);
        }
        else
        {
            // Extend into a match using at most one mod.  
            if (BufferPos)
            {
                PrefixChar = Buffer[BufferPos - 1];
            }
            else
            {
                PrefixChar = '*';
            }
            SkipBases = SeekMatch1PTM(Info, Buffer + BufferPos, BufferEnd - BufferPos, ParentResidueMass, ScoreToBeat, Tweak, FilePos, PrefixChar);
        }
        // RightExtensionCount is set to -1 if there's an error and this spectrum can't be searched.
        if (RightExtensionCount < 0)
        {
            break;
        }
        BufferPos += SkipBases;
        FilePos += SkipBases;
        if (Buffer[BufferPos-1] == '*')
        {
            Info->RecordNumber++;
        }
        if (Spectrum->Node->MatchCount == GlobalOptions->StoreMatchCount)
        {
            ScoreToBeat = Spectrum->Node->LastMatch->InitialScore;
        }
        
    }
    SafeFree(Buffer);
    Buffer = NULL;

    // At this point, we have a list of candidates.  They've been quick-scored, but we can sort them better if
    // we score them more meticulously.  The *caller* will call ScoreSpectralMatches(Spectrum) to re-score them
    // (we could do it here, but that would be wrong in the MultiCharge case)
    //fclose(DBFile);
}

void DebugPrintMatch(Peptide* Match)
{
    int ModIndex;
    char* Amino;
    int Mass = 0;
    printf("Match '%s' ", Match->Bases);
    // Show the mods:
    for (ModIndex = 0; ModIndex < MAX_PT_MODS; ModIndex++)
    {
        if (Match->AminoIndex[ModIndex] < 0)
        {
            break;
        }
        printf(" %c%d:%s(%.2f)", Match->Bases[Match->AminoIndex[ModIndex]], Match->AminoIndex[ModIndex],
            Match->ModType[ModIndex]->Name, Match->ModType[ModIndex]->RealDelta/100.0);
        Mass += Match->ModType[ModIndex]->RealDelta;
    }
    for (Amino = Match->Bases; *Amino; Amino++)
    {
        Mass += PeptideMass[*Amino];
    }
    printf(" mass %.2f score %d:%.3f dcn%.2f dcno%.2f\n", Mass/(float)MASS_SCALE, Match->InitialScore, Match->MatchQualityScore, Match->DeltaCN, Match->DeltaCNOther);
}

// Print out a list of matches for the spectrum node (Spectrum->FirstMatch through Spectrum->LastMatch).
void DebugPrintMatchList(SpectrumNode* Spectrum)
{
    Peptide* Match;
    //
    printf("Spectrum has %d matches:\n", Spectrum->MatchCount);
    for (Match = Spectrum->FirstMatch; Match; Match = Match->Next)
    {
        DebugPrintMatch(Match);
    }
}

// Re-score spectral matches.  The matches in the list Spectrum->FirstMatch have been 
// quick-scored, but we can sort them better if we score them more meticulously.  
// Let's do so, and re-sort the list based on the new scores.
void MQScoreSpectralMatches(SpectrumNode* Node)
{
    Peptide* PrevMatch;
    Peptide* Match;
    int OldScore;
    int VerboseFlag = 0;
    MSSpectrum* Spectrum = Node->Spectrum;
    //
    if (!Node->FirstMatch)
    {
        return; // that was easy - we scored 0 of 0 :)
    }
    PrevMatch = Node->FirstMatch;

    Match = PrevMatch->Next;
    Node->FirstMatch = NULL;
    Node->LastMatch = NULL;
    Node->MatchCount = 0;
    while (PrevMatch)
    {
        PrevMatch->Prev = NULL;
        PrevMatch->Next = NULL;
        OldScore = PrevMatch->InitialScore;
        ComputeMQScoreFeatures(Spectrum, PrevMatch, PrevMatch->ScoreFeatures, VerboseFlag);
#ifdef MQSCORE_USE_SVM
        PrevMatch->MatchQualityScore = SVMComputeMQScore(Spectrum, PrevMatch, PrevMatch->ScoreFeatures);
#else
        PrevMatch->MatchQualityScore = LDAComputeMQScore(Spectrum, PrevMatch, PrevMatch->ScoreFeatures);
#endif
        StoreSpectralMatch(Spectrum, PrevMatch, strlen(PrevMatch->Bases), 1);
        PrevMatch = Match;
        if (!Match)
        {
            break;
        }
        Match = Match->Next;
    }
    //SetMatchDeltaCN(Spectrum);
}

void PrunePoorGraphNodes(TagGraph* Graph)
{
    int NodeIndex = 0;
    float* NodeScores;
    float CutoffNodeScore;
    TagGraphNode* Node;
    TagGraphNode* NextNode = NULL;
    TagGraphEdge* Edge;
    TagGraphEdge* NextEdge = NULL;
    TagGraphEdge* PrevEdge = NULL;
    //
    // Write the node scores to array NodeScores, sort them, and select the cutoff score.
    NodeScores = (float*)malloc(sizeof(float) * Graph->NodeCount);
    
    for (NodeIndex = 0, Node = Graph->FirstNode; Node; Node = Node->Next,NodeIndex++)
    {
        NodeScores[NodeIndex] = Node->Score;
    }
    qsort(NodeScores, NodeIndex, sizeof(float), (QSortCompare)CompareFloats);
    CutoffNodeScore = NodeScores[498]; // Allow two endpoint nodes to survive
    SafeFree(NodeScores);
    // Eliminate every node whose score is <= the cutoff.  Start by eliminating all EDGES to such nodes!
    for (Node = Graph->FirstNode; Node; Node = Node->Next)
    {
        PrevEdge = NULL;
        for (Edge= Node->FirstEdge; Edge; Edge = NextEdge)
        {
            NextEdge = Edge->Next;
            if (Edge->ToNode->Score <= CutoffNodeScore && (Edge->ToNode->NodeType == evGraphNodeB || Edge->ToNode->NodeType == evGraphNodeY))
            {
                // Free this edge:
                if (PrevEdge)
                {
                    PrevEdge->Next = Edge->Next;
                }
                if (Node->FirstEdge == Edge)
                {
                    Node->FirstEdge = Edge->Next;
                }
                if (Node->LastEdge == Edge)
                {
                    Node->LastEdge = PrevEdge;
                }
                SafeFree(Edge);
            }
            else
            {
                PrevEdge = Edge;
            }
        }
    }
    // Now free the nodes themselves:
    for (Node = Graph->FirstNode; Node; Node = NextNode)
    {
        NextNode = Node->Next;
        if (Node->Score <= CutoffNodeScore && (Node->NodeType == evGraphNodeB || Node->NodeType == evGraphNodeY))
        {
            if (Node->Prev)
            {
                Node->Prev->Next = Node->Next;
            }
            if (Node->Next)
            {
                Node->Next->Prev = Node->Prev;
            }
            if (Graph->FirstNode == Node)
            {
                Graph->FirstNode = Node->Next;
            }
            if (Graph->LastNode == Node)
            {
                Graph->LastNode = Node->Prev;
            }
            FreeTagGraphNode(Node);
            Graph->NodeCount--;
        }
    }
    if (Graph->NodeCount > 500)
    {
        printf("* ERROR: Failed to prune excess graph nodes!\n");
    }
    // Fix node numbering:
    for (NodeIndex = 0, Node = Graph->FirstNode; Node; Node = Node->Next, NodeIndex++)
    {
        Node->Index = NodeIndex;
    }
    // And now, rebuild the node index:
    TagGraphBuildNodeIndex(Graph);
}

// Called after populating the tag graph with nodes.
// Now we add edges between any two nodes that can be linked by a JUMP (an amino acid, or 
// an amino acid plus a decoration, or two amino acids plus 0-1 decorations)
void TagGraphPopulateBackEdges(TagGraph* Graph)
{
    TagGraphNode* Node;
    TagGraphNode* OtherNode;
    TagGraphBackEdge* Edge;
    TagGraphBackEdge* OldEdge;
    int AA1;
    int AA2;
    int AA3;
    int Mass;
    int AA1Mass;
    int AA2Mass;
    int AA3Mass;
    int Skew;
    int AbsSkew;
    int NextBackEdgeIndex = 0;
    int BackEdgeBufferSize;
    //

    if (!Graph->NodeIndex)
    {
        TagGraphBuildNodeIndex(Graph);
    }
    SafeFree(Graph->BackEdgeBuffer);
    BackEdgeBufferSize = min(5000000, 8420 * Graph->NodeCount);
    Graph->BackEdgeBuffer = (TagGraphBackEdge*)calloc(BackEdgeBufferSize, sizeof(TagGraphBackEdge));
    if (!Graph->BackEdgeBuffer)
    {
        printf("*** ERROR: Unable to allocate BackEdgeBuffer!\n");
        fflush(stdout);
    }
    // NB: We can't easily realloc the BackEdgeBuffer, because there are many many pointers into it.  If
    // we overflow the buffer, we just complain and then bail out to avoid crashing.

    // Ensure that there aren't too many PRMNodes.  The array in ExtendMatchRightwardDuo assumes that there
    // are at most 500.  That should be *plenty*, since at most 20-30 of them can be true.
    if (Graph->NodeCount > 500)
    {
        PrunePoorGraphNodes(Graph);
    }

    for (Node = Graph->FirstNode; Node; Node = Node->Next)
    {
        Node->BackEdge = (TagGraphBackEdge**)calloc(AMINO_ACIDS, sizeof(TagGraphBackEdge*));
        Node->BackEdgeDouble = (TagGraphBackEdge**)calloc(AMINO_ACIDS*AMINO_ACIDS, sizeof(TagGraphBackEdge*));
        Node->BackEdgeTriple = (TagGraphBackEdge**)calloc(AMINO_ACIDS*AMINO_ACIDS*AMINO_ACIDS, sizeof(TagGraphBackEdge*));
        for (AA1 = 0; AA1 < AMINO_ACIDS; AA1++)
        {
            AA1Mass = PeptideMass[AA1 + 'A'];
            if (!AA1Mass)
            {
                continue;
            }
            // Try to jump back by this amino acid's mass:
            Mass = Node->Mass - AA1Mass;
            if (Mass < -GlobalOptions->Epsilon)
            {
                continue;
            }
            Mass = max(Mass, 0);
            OtherNode = Graph->NodeIndex[Mass / MASS_SCALE];
            while (OtherNode)
            {
                Skew = OtherNode->Mass - Mass;
                if (Skew > GlobalOptions->Epsilon)
                {
                    break;
                }
                if (Skew < -GlobalOptions->Epsilon)
                {
                    OtherNode = OtherNode->Next;
                    continue;
                }
                AbsSkew = abs(Skew) / 10;
                Edge = Graph->BackEdgeBuffer + NextBackEdgeIndex;
                NextBackEdgeIndex++;
                if (NextBackEdgeIndex >= BackEdgeBufferSize)
                {
                    printf("** Too many BackEdges for buffer - bailing out!\n");
                    return;
                }
                //Edge = (TagGraphBackEdge*)calloc(1, sizeof(TagGraphBackEdge));
                Edge->FromNode = Node;
                Edge->ToNode = OtherNode;
                Edge->Skew = Skew;
                Edge->Next = NULL;
                if (AbsSkew > g_SkewPenaltyMax)
                {
                    Edge->Score = g_SkewPenalty[g_SkewPenaltyMax];
                }
                else
                {
                    Edge->Score = g_SkewPenalty[AbsSkew];
                }
                OldEdge = Node->BackEdge[AA1];
                if (!OldEdge)
                {
                    Node->BackEdge[AA1] = Edge;
                }
                else
                {
                    while (OldEdge->Next)
                    {
                        OldEdge = OldEdge->Next;
                    }
                    OldEdge->Next = Edge;
                }
                OtherNode = OtherNode->Next;
            }

            // Try to jump back by a pair of amino acids:
            for (AA2 = 0; AA2 < AMINO_ACIDS; AA2++)
            {
                AA2Mass = PeptideMass[AA2 + 'A'];
                if (!AA2Mass)
                {
                    continue;
                }

                Mass = Node->Mass - AA1Mass - AA2Mass;
                if (Mass < -GlobalOptions->Epsilon)
                {
                    continue;
                }
                Mass = max(Mass, 0);
                OtherNode = Graph->NodeIndex[Mass / MASS_SCALE];
                while (OtherNode)
                {
                    Skew = OtherNode->Mass - Mass;
                    if (Skew > GlobalOptions->Epsilon)
                    {
                        break;
                    }
                    if (Skew < -GlobalOptions->Epsilon)
                    {
                        OtherNode = OtherNode->Next;
                        continue;
                    }
                    AbsSkew = abs(Skew) / 10;

                    Edge = Graph->BackEdgeBuffer + NextBackEdgeIndex;
                    NextBackEdgeIndex++;
                    if (NextBackEdgeIndex >= BackEdgeBufferSize)
                    {
                        printf("** Too many BackEdges for buffer - bailing out!\n");
                        return;
                    }
                
                    //Edge = (TagGraphBackEdge*)calloc(1, sizeof(TagGraphBackEdge));
                    Edge->FromNode = Node;
                    Edge->ToNode = OtherNode;
                    Edge->Skew = Skew;
                    Edge->Next = NULL;
                    Edge->HalfMass = Node->Mass - AA1Mass;
                    if (AbsSkew > g_SkewPenaltyMax)
                    {
                        Edge->Score = g_SkewPenalty[g_SkewPenaltyMax];
                    }
                    else
                    {
                        Edge->Score = g_SkewPenalty[AbsSkew];
                    }
                    OldEdge = Node->BackEdgeDouble[AA1 * AMINO_ACIDS + AA2];
                    if (!OldEdge)
                    {
                        Node->BackEdgeDouble[AA1 * AMINO_ACIDS + AA2] = Edge;
                    }
                    else
                    {
                        while (OldEdge->Next)
                        {
                            OldEdge = OldEdge->Next;
                        }
                        OldEdge->Next = Edge;
                    }
                    OtherNode = OtherNode->Next;
                }

                // Triple-jump (three amino acids) 
                for (AA3 = 0; AA3 < AMINO_ACIDS; AA3++)
                {
                    AA3Mass = PeptideMass[AA3+'A'];
                    if (!AA3Mass)
                    {
                        continue;
                    }

                    Mass = Node->Mass - AA1Mass - AA2Mass - AA3Mass;
                    if (Mass < -GlobalOptions->Epsilon)
                    {
                        continue;
                    }
                    Mass = max(Mass, 0);
                    OtherNode = Graph->NodeIndex[Mass / MASS_SCALE];
                    while (OtherNode)
                    {
                        Skew = OtherNode->Mass - Mass;
                        if (Skew > GlobalOptions->Epsilon)
                        {
                            break;
                        }
                        if (Skew < -GlobalOptions->Epsilon)
                        {
                            OtherNode = OtherNode->Next;
                            continue;
                        }
                        AbsSkew = abs(Skew) / 10;

                        Edge = Graph->BackEdgeBuffer + NextBackEdgeIndex;
                        NextBackEdgeIndex++;
                        if (NextBackEdgeIndex >= BackEdgeBufferSize)
                        {
                            printf("** Too many BackEdges for buffer - bailing out!\n");
                            return;
                        }
                    
                        //Edge = (TagGraphBackEdge*)calloc(1, sizeof(TagGraphBackEdge));
                        Edge->FromNode = Node;
                        Edge->ToNode = OtherNode;
                        Edge->Skew = Skew;
                        Edge->Next = NULL;
                        Edge->HalfMass = Node->Mass - AA1Mass;
                        Edge->HalfMass2 = Node->Mass - AA1Mass - AA2Mass;
                        if (AbsSkew > g_SkewPenaltyMax)
                        {
                            Edge->Score = g_SkewPenalty[g_SkewPenaltyMax];
                        }
                        else
                        {
                            Edge->Score = g_SkewPenalty[AbsSkew];
                        }
                        OldEdge = Node->BackEdgeTriple[AA1*676 + AA2*26 + AA3];
                        if (!OldEdge)
                        {
                            Node->BackEdgeTriple[AA1*676 + AA2*26 + AA3] = Edge;
                        }
                        else
                        {
                            while (OldEdge->Next)
                            {
                                OldEdge = OldEdge->Next;
                            }
                            OldEdge->Next = Edge;
                        }
                        OtherNode = OtherNode->Next;
                    }
                }
            }
        }
    }
    return;
}

// The maximum dimensions of the Duo Table are 512 rows (for amino acids) and 500 columns (for the nodes).
#define MAX_ROWS 512
#define MAX_NODES 500

// Find the end of this peptide block.  Returns 1 if the block is valid, 0 if we
// needn't bother searching this block at all.
int FindPeptideBlockEnd(MSSpectrum* Spectrum, char* Buffer,int BufferEnd, int* pMaxAmino, int* pReturnAmino)
{
    int AccumMass = 0;
    int AminoIndex;
    char Amino;
    //
    *pMaxAmino = MAX_ROWS - 1; // default;
    *pReturnAmino = -1; // uninitialized
    // Iterate over the amino acids, keeping track of the total mass (AccumMass), and watching for
    // record boundaries.
    for (AminoIndex = 0; AminoIndex < MAX_ROWS; AminoIndex++)
    {
        if (AminoIndex >= BufferEnd)
        {
            *pMaxAmino = AminoIndex;
            *pReturnAmino = *pMaxAmino;
            break;
        }
        Amino = Buffer[AminoIndex];
        if (Amino == '*')
        {
            // No peptide can span record boundaries, so we'll stop the block here.
            *pMaxAmino = AminoIndex;// + 1;
            *pReturnAmino = *pMaxAmino;
            break;
        }
        if (!PeptideMass[Amino])
        {
            // A bogus amino acid in the database.  Stop the block here.
            *pMaxAmino = AminoIndex;// + 1;
            *pReturnAmino = *pMaxAmino;
            break;
        }
        AccumMass += PeptideMass[Amino];
        if (AminoIndex >= MAX_ROWS-1)
        {
            // Our block is as large as it can get.
            *pMaxAmino = AminoIndex-1;
            *pReturnAmino = *pMaxAmino - 20;
            break;
        }
    }
    if (*pMaxAmino < 5)
    {
        // Not enough amino acids to make a reasonable peptide candidate, so just exit.
        // (The longest we could get is 4aa, which is too small)
        *pReturnAmino = *pMaxAmino;
        return 0; 
    }
    // Check to see whether the amino acids we've got are large enough - with PTMs included - to match the target:
    AccumMass += GlobalOptions->MaxPTMDelta*100*2 + PARENT_MASS_BOOST + GlobalOptions->ParentMassEpsilon;
    if (AccumMass < Spectrum->ParentMass)
    {
        // There's not enough peptide sequence left to match our target.
        *pReturnAmino = *pMaxAmino;
        return 0; 
    }
    if (*pReturnAmino < 0)
    {
        // Shift forward by most of the length of the block.  (Leave some overlap, because a valid peptide
        // may start near the end of block #1)
        *pReturnAmino = max(1, *pMaxAmino - 20);
    }
    return 1;
}

void DebugPrintPrefixSuffixTable(FILE* TableFile, char* Buffer, int MaxAmino, int* ScoreTable, int* MassTable)
{
    int AminoIndexI;
    int AminoIndexJ;
    int CellIndex;

    // Header line:
    fprintf(TableFile, "<  >\t");
    for (AminoIndexJ = 0; AminoIndexJ <= MaxAmino; AminoIndexJ++)
    {
        if (AminoIndexJ > 0)
        {
            fprintf(TableFile, "[%d %c]\t", AminoIndexJ, Buffer[AminoIndexJ-1]);
        }
        else
        {
            fprintf(TableFile, "[0]\t");
        }
    }
    fprintf(TableFile, "\n");
    // Other lines:
    for (AminoIndexI = 0; AminoIndexI <= MaxAmino; AminoIndexI++)
    {
        if (AminoIndexI > 0)
        {
            fprintf(TableFile, "[%d %c]\t", AminoIndexI, Buffer[AminoIndexI-1]);
        }
        else
        {
            fprintf(TableFile, "[0]\t");
        }

        for (AminoIndexJ = 0; AminoIndexJ <= MaxAmino; AminoIndexJ++)
        {
            if (AminoIndexJ < AminoIndexI)
            {
                fprintf(TableFile, "\t");
                continue;
            }
            CellIndex = AminoIndexI*MAX_ROWS + AminoIndexJ;
            fprintf(TableFile, "%d : %d\t", MassTable[CellIndex], ScoreTable[CellIndex]);
        }
        fprintf(TableFile, "\n");
    }
}

void DebugPrintPrefixSuffixTables(int MaxAmino, char* Buffer, 
                                  int* PrefixTable, int* SuffixTable, int* PrefixMassTable, int* SuffixMassTable)
{
    FILE* TableFile;
    TableFile = fopen("PrefixSuffix.xls", "w");
    if (!TableFile)
    {
        return;
    }
    fprintf(TableFile, "Prefix table:\n");
    DebugPrintPrefixSuffixTable(TableFile, Buffer, MaxAmino, PrefixTable, PrefixMassTable);
    fprintf(TableFile, "\n\nSuffix table:\n");
    DebugPrintPrefixSuffixTable(TableFile, Buffer, MaxAmino, SuffixTable, SuffixMassTable);
    fclose(TableFile);
}

// Fill in the score tables PrefixTable and SuffixTable, plus the mass tables PrefixMassTable and SuffixMassTable.
// The entry PrefixTable[i, j] is the score that one obtains by matching Buffer[i..j] against the spectrum as
// a prefix.  PrefixTable[i,i] uses one PRM score, PrefixTable[i, i + 1] uses two PRM scores, and so on.
// Most candidate peptides will have a PrefixTable entry as part of their final score.
void FillPrefixSuffixTables(MSSpectrum* Spectrum, SpectrumTweak* Tweak, int MatchMass, char* Buffer, int MaxAmino, 
    int* PrefixTable, int* SuffixTable, int* PrefixMassTable, int* SuffixMassTable)
{
    int AminoIndexI;
    int AminoIndexJ;
    int CellIndex;
    int PrevCellIndex = 0;
    int PRM;

    int MaxPRM = Tweak->PRMScoreMax;
    //
    // Brute force initializiation.  (Note: Don't use memset here, because setting every byte
    // to -1 (a) is hacky, and (b) makes scores that can easily wrap around to become HUGE POSITIVE)
    for (CellIndex = 0; CellIndex < MAX_ROWS * MAX_ROWS; CellIndex++)
    {
        PrefixTable[CellIndex] = FORBIDDEN_PATH;
        SuffixTable[CellIndex] = FORBIDDEN_PATH;
        PrefixMassTable[CellIndex] = -999999;
        SuffixMassTable[CellIndex] = -999999;
    }
    
    for (AminoIndexI = 1; AminoIndexI <= MaxAmino; AminoIndexI++)
    {
        for (AminoIndexJ = AminoIndexI; AminoIndexJ <= MaxAmino; AminoIndexJ++)
        {
            CellIndex = AminoIndexI * MAX_ROWS + AminoIndexJ;
            /////////////////////////////
            // Prefix table:
            if (AminoIndexJ == AminoIndexI)
            {
                PrefixMassTable[CellIndex] = PeptideMass[Buffer[AminoIndexI-1]];
            }
            else
            {
                PrevCellIndex = AminoIndexI*MAX_ROWS + (AminoIndexJ-1);
                PrefixMassTable[CellIndex] = PrefixMassTable[PrevCellIndex] + PeptideMass[Buffer[AminoIndexJ-1]];
            }
            PRM = MASS_TO_BIN(PrefixMassTable[CellIndex]);

            // Allow PRMs that are slightly too big or small:
            if (PRM > -PRM_ARRAY_SLACK)
            {
                PRM = max(PRM, 0);
            }
            if (PRM < MaxPRM + 5)
            {
                PRM = min(MaxPRM, PRM);
            }
            if (PRM >= 0 && PRM <= MaxPRM)
            {
                PrefixTable[CellIndex] = Tweak->PRMScores[PRM];
                if (AminoIndexJ != AminoIndexI)
                {
                    PrefixTable[CellIndex] += PrefixTable[PrevCellIndex];
                }
            }
            else
            {
                PrefixTable[CellIndex] = FORBIDDEN_PATH;
                break;
            }
        }
    }
    for (AminoIndexJ = MaxAmino; AminoIndexJ; AminoIndexJ--)
    {
        for (AminoIndexI = AminoIndexJ; AminoIndexI; AminoIndexI--)
        {
            CellIndex = AminoIndexI*MAX_ROWS + AminoIndexJ;

            /////////////////////////////
            // Suffix table:
            if (AminoIndexJ == AminoIndexI)
            {
                SuffixMassTable[CellIndex] = MatchMass - PeptideMass[Buffer[AminoIndexI-1]];
            }
            else
            {
                PrevCellIndex = (AminoIndexI + 1)*MAX_ROWS + AminoIndexJ;
                SuffixMassTable[CellIndex] = SuffixMassTable[PrevCellIndex] - PeptideMass[Buffer[AminoIndexI-1]];
            }
            PRM = MASS_TO_BIN(SuffixMassTable[CellIndex]);
            if (PRM > -PRM_ARRAY_SLACK)
            {
                PRM = max(PRM, 0);
            }
            if (PRM < MaxPRM+5)
            {
                PRM = min(MaxPRM, PRM);
            }
            if (PRM >= 0 && PRM <= MaxPRM)
            {
                SuffixTable[CellIndex] = Tweak->PRMScores[PRM];
                if (AminoIndexI!=AminoIndexJ)
                {
                    SuffixTable[CellIndex] += SuffixTable[PrevCellIndex];
                }
            }
            else
            {
                SuffixTable[CellIndex] = FORBIDDEN_PATH;
                break;
            }
        }
    }

}

// Print the Duo table to a file, for debugging.  (This is most useful when searching a very small database, since then
// the table has a managable height)
void DebugPrintDTable(MSSpectrum* Spectrum, char* Buffer, int* DTable, MassDelta** DeltaTable, int* PrevCellTable, int MaxAmino)
{
    int AminoIndex;
    int NodeIndex;
    TagGraphNode* Node;
    FILE* DTableFile = NULL;
    int CellIndex;
    int PrevCellIndex = 0;
    
    int AminoBlockSize = Spectrum->Graph->NodeCount;
    //
    DTableFile = fopen("DTable.xls", "w");
    if (!DTableFile)
    {
        return;
    }
    // Header:
    fprintf(DTableFile, "\t");
    for (NodeIndex = 0, Node = Spectrum->Graph->FirstNode; Node; Node = Node->Next, NodeIndex++)
    {
        fprintf(DTableFile, "%d (%.2f)\t", NodeIndex, Node->Mass / 100.0);
    }
    fprintf(DTableFile, "\n");
    // Body:
    for (AminoIndex = 0; AminoIndex < MaxAmino; AminoIndex++)
    {
        if (AminoIndex)
        {
            fprintf(DTableFile, "%c %d\t", Buffer[AminoIndex-1], AminoIndex);
        }
        else
        {
            fprintf(DTableFile, "%d\t", AminoIndex);
        }
        for (NodeIndex = 0, Node = Spectrum->Graph->FirstNode; Node; Node = Node->Next, NodeIndex++)
        {
            CellIndex = AminoIndex*AminoBlockSize + NodeIndex;
            fprintf(DTableFile, "%d ", DTable[CellIndex]);
            if (DeltaTable[CellIndex])
            {
                fprintf(DTableFile, "%s", DeltaTable[CellIndex]->Name);
            }
            PrevCellIndex = PrevCellTable[CellIndex];
            if (PrevCellIndex >0)
            {
                fprintf(DTableFile, "-> (%d, %d)", PrevCellIndex/AminoBlockSize, PrevCellIndex%AminoBlockSize);
            }
            fprintf(DTableFile, "\t");
        }
        fprintf(DTableFile, "\n");
    }

    fclose(DTableFile);
}

void DebugPrintGeneralTable(MSSpectrum* Spectrum, char* Buffer, int MaxAmino, int MaxMods,
    int* ScoreTable, int* PrevCellTable)
{
    int AminoIndex;
    int NodeIndex;
    int CellIndex;
    int ModCountIndex;
    FILE* DebugFile;
    TagGraphNode* Node;
    int AminoBlockSize;
    int ZSize = MaxMods + 1;
    //
    DebugFile = fopen("DPTable.txt", "wb");
    if (!DebugFile)
    {
        printf("Unable to open DPTable.txt - not debugprinting.\n");
        return;
    }
    AminoBlockSize = Spectrum->Graph->NodeCount * ZSize;
    for (ModCountIndex = 0; ModCountIndex < ZSize; ModCountIndex++)
    {
        fprintf(DebugFile, "\nZ = %d\n", ModCountIndex);
        /////////////////////////////////
        // Column headers:
        fprintf(DebugFile, "\t");
        for (AminoIndex = 0; AminoIndex < MaxAmino; AminoIndex++)
        {
            fprintf(DebugFile, "%d\t", AminoIndex);
        }
        fprintf(DebugFile, "\n");
        fprintf(DebugFile, "\t\t");
        for (AminoIndex = 1; AminoIndex < MaxAmino; AminoIndex++)
        {
            fprintf(DebugFile, "%c\t", Buffer[AminoIndex - 1]);
        }
        fprintf(DebugFile, "\n");
        /////////////////////////////////
        // Body:
        for (NodeIndex = 0, Node = Spectrum->Graph->FirstNode; Node; Node = Node->Next, NodeIndex++)
        {
            // Print a ROW for this node:
            fprintf(DebugFile, "%d (%.2f)\t", NodeIndex, Node->Mass / 100.0);
            for (AminoIndex = 0; AminoIndex < MaxAmino; AminoIndex++)
            {
                CellIndex = AminoIndex * AminoBlockSize + NodeIndex * ZSize + ModCountIndex;
                fprintf(DebugFile, "%d (c%d)\t", ScoreTable[CellIndex], CellIndex);
            }
            fprintf(DebugFile, "\n");
        }
    }
    fclose(DebugFile);
}


static int* PrefixTable = NULL;
static int* SuffixTable = NULL;
static int* PrefixMassTable = NULL;
static int* SuffixMassTable = NULL;

// MS-Alignment algorithm, general case (handles k>2).  For most purposes, this code is
// unacceptably slow and non-selective.  But, for completeness, it is implemented.
// In practice, one should use "mods,1" or "mods,2" and find the corpus of 
// available PTMs that way.
int MSAlignmentGeneral(SearchInfo* Info, char* Buffer, int BufferEnd, int MatchMass, 
    int MaxMods, int ScoreToBeat, int FilePos, SpectrumTweak* Tweak)
{
    static int* PrevCellTable = NULL;
    static int* ScoreTable = NULL;
    // StartPointPenalty and EndPointPenalty provide a simple protease specificity.
    int StartPointPenalty[MAX_ROWS];
    int EndPointPenalty[MAX_ROWS];
    int Result;
    int ReturnAmino = -1;
    int AminoIndex;
    int AA = 0;
    int AminoBlock = 0;
    int NodeIndex;
    int CellIndex;
    int SliceIndex;
    int AA2;
    int AA3;
    TagGraphNode* Node;
    TagGraphBackEdge* Edge;
    int ModCountIndex;
    int AminoBlockSize;
    int ZSize = MaxMods + 1;
    int MaxAmino;
    int AAMass = 0;
    int AA2Mass;
    int BackEdgeDoubleIndex = 0;
    int BackEdgeTripleIndex = 0;
    int PrevCellIndex;
    int Score;
    int PRM;
    TagGraphNode* BackNode;
    //char MatchBuffer[256];
    char MatchBufferPos;
    Peptide* Match;
    int Mass;
    int BackNodeIndex;
    int BackNodeDirection;
    int Delta = 0;
    int MinPossibleDelta = -130000; // W->G mutation
    int MaxPossibleDelta = GlobalOptions->MaxPTMDelta * MASS_SCALE;
    int DeltaBin;
    MassDeltaNode* DeltaNode;
    int Skew;
    int AbsSkew;
    int MaxPRM = Tweak->PRMScoreMax - 1;
    int X;
    int Y;
    int NextY;
    int Z;
    int MatchPTMCount;
    int ModIndex;
    int TokenDropped;
    MSSpectrum* Spectrum = Info->Spectrum;
    //

    // Allocate tables, if necessary:
    if (!PrevCellTable)
    {
        PrevCellTable = (int*)malloc(sizeof(int) * MAX_ROWS * MAX_ROWS * (MaxMods + 1));
        ScoreTable = (int*)malloc(sizeof(int) * MAX_ROWS * MAX_NODES * (MaxMods + 1));
        DeltaTable = (MassDelta**)malloc(sizeof(MassDelta*) * MAX_ROWS * MAX_NODES * (MaxMods+1));
        MassDeltaTable = (int*)malloc(sizeof(int) * MAX_ROWS * MAX_NODES * (MaxMods + 1));
    }

    /////////////////////////////////////////
    // Find MaxAmino:
    //VerboseFlag = 1;
    if (Info->VerboseFlag)
    {
        printf("[V] FindPeptideBlockEnd:\n");
    }
    Result = FindPeptideBlockEnd(Spectrum, Buffer, BufferEnd, &MaxAmino, &ReturnAmino);
    if (!Result)
    {
        // No extension necessary.  Advance the database pointer:
        return ReturnAmino;
    }
    if (Info->VerboseFlag)
    {
        printf("[V] FindPeptideBlockEnd: MaxAmino %d returnamino %d\n", MaxAmino, ReturnAmino);
    }

    // Apply a slap-on-the-wrist for using non-tryptic peptides:
    for (AminoIndex = 0; AminoIndex < MaxAmino; AminoIndex++)
    {
        StartPointPenalty[AminoIndex] = 0;
        if (AminoIndex)
        {
            AA = Buffer[AminoIndex - 1];
            if (AA != 'R' && AA != 'K' && AA != '*')
            {
                StartPointPenalty[AminoIndex] = -500;
            }
        }

        EndPointPenalty[AminoIndex] = 0;
        if (AminoIndex)
        {
            AA = Buffer[AminoIndex - 1];
            if ((AA != 'R' && AA != 'K') && (AminoIndex <= MaxAmino-1 && Buffer[AminoIndex + 1]!='*'))
            {
                EndPointPenalty[AminoIndex] = -500;
            }
        }
    }
    AminoBlockSize = Spectrum->Graph->NodeCount * ZSize;
    // Loop over the d.p. table to populate scores and path.
    // Iterate by amino acid index (row), by node (column), then by PTMCount (z).
    for (AminoIndex = 0; AminoIndex <= MaxAmino; AminoBlock += AminoBlockSize, AminoIndex++)
    {
        CellIndex = AminoBlock;
        AA2 = 0;
        AA3 = 0;
        if (AminoIndex)
        {
            AA = Buffer[AminoIndex-1] - 'A';
            AAMass = PeptideMass[Buffer[AminoIndex - 1]];
        }
        if (AminoIndex > 1)
        {
            AA2 = Buffer[AminoIndex-2] - 'A';
            AA2Mass = PeptideMass[Buffer[AminoIndex - 2]];
            BackEdgeDoubleIndex = AA*AMINO_ACIDS + AA2;
        }
        if (AminoIndex > 2)
        {
            AA3 = Buffer[AminoIndex - 3] - 'A';
            BackEdgeTripleIndex = AA*676 + AA2*26 + AA3;
        }
        for (NodeIndex = 0, Node = Spectrum->Graph->FirstNode; Node; Node = Node->Next, NodeIndex++)
        {
            SliceIndex = CellIndex;
            for (ModCountIndex = 0; ModCountIndex < ZSize; ModCountIndex++)
            {
                // Check our cell index:
                if (CellIndex != AminoIndex * AminoBlockSize + NodeIndex * ZSize + ModCountIndex)
                {
                    printf("Bad cell index %d, %d, %d -> %d (%d)\n", NodeIndex, AminoIndex, ModCountIndex, CellIndex,
                        AminoIndex * AminoBlockSize + NodeIndex * ZSize + ModCountIndex);
                }
                ScoreTable[CellIndex] = FORBIDDEN_PATH; // default
                DeltaTable[CellIndex] = NULL;
                PrevCellTable[CellIndex] = -1;
                
                ///////////////
                // Free rides:
                if (ModCountIndex == 0 && Node->Mass < GlobalOptions->ParentMassEpsilon)
                {
                    ScoreTable[CellIndex] = StartPointPenalty[AminoIndex];
                    PrevCellTable[CellIndex] = -1;
                    DeltaTable[CellIndex] = NULL;
                    MassDeltaTable[CellIndex] = Node->Mass;
                    CellIndex++;
                    continue;
                }
                ///////////////
                // Drop a token:
                TokenDropped = 0;
                if (ModCountIndex)
                {
                    PrevCellIndex = CellIndex - 1; 
                    ScoreTable[CellIndex] = ScoreTable[PrevCellIndex];
                    PrevCellTable[CellIndex] = PrevCellIndex;
                    TokenDropped = 1;
                }

                // And that's all we do on the top row:
                if (AminoIndex == 0)
                {
                    CellIndex++;
                    continue; 
                }
               
                ///////////////
                // One unmodified amino acid:
                Edge = Node->BackEdge[AA];
                while (Edge)
                {
                    PrevCellIndex = AminoBlock - AminoBlockSize + (Edge->ToNode->Index * ZSize) + ModCountIndex;
                    Score = ScoreTable[PrevCellIndex] + Edge->Score; 
                    if (Score > ScoreTable[CellIndex])
                    {
                        ScoreTable[CellIndex] = Score;
                        PrevCellTable[CellIndex] = PrevCellIndex;
                        DeltaTable[CellIndex] = NULL;
                        MassDeltaTable[CellIndex] = Edge->Skew + MassDeltaTable[PrevCellIndex];
                        TokenDropped = 0;
                    }
                    Edge = Edge->Next;
                }

                ///////////////
                // Two unmodified amino acids:
                if (AminoIndex > 1)
                {
                    Edge = Node->BackEdgeDouble[BackEdgeDoubleIndex];
                    while (Edge)
                    {
                        PrevCellIndex = AminoBlock - AminoBlockSize - AminoBlockSize + (Edge->ToNode->Index * ZSize) + ModCountIndex;
                        // Accumulate points for the middle of the jump:
                        PRM = MASS_TO_BIN(Edge->HalfMass);
                        Score = Tweak->PRMScores[PRM] + ScoreTable[PrevCellIndex] + Edge->Score;
                        if (Score > ScoreTable[CellIndex])
                        {
                            ScoreTable[CellIndex] = Score;
                            PrevCellTable[CellIndex] = PrevCellIndex;
                            DeltaTable[CellIndex] = NULL;
                            MassDeltaTable[CellIndex] = Edge->Skew + MassDeltaTable[PrevCellIndex];
                            TokenDropped = 0;
                        }
                        Edge = Edge->Next;
                    }
                }
                ///////////////
                // Three unmodified amino acids:
                if (AminoIndex > 2)
                {
                    Edge = Node->BackEdgeTriple[BackEdgeTripleIndex];
                    while (Edge)
                    {
                        PrevCellIndex = AminoBlock - AminoBlockSize - AminoBlockSize - AminoBlockSize + (Edge->ToNode->Index * ZSize) + ModCountIndex;
                        // Accumulate points for the middle of the jump:
                        PRM = MASS_TO_BIN(Edge->HalfMass);
                        Score = Tweak->PRMScores[PRM] + ScoreTable[PrevCellIndex] + Edge->Score;
                        PRM = MASS_TO_BIN(Edge->HalfMass2);
                        Score += Tweak->PRMScores[PRM];
                        if (Score > ScoreTable[CellIndex])
                        {
                            ScoreTable[CellIndex] = Score;
                            PrevCellTable[CellIndex] = PrevCellIndex;
                            DeltaTable[CellIndex] = NULL;
                            MassDeltaTable[CellIndex] = Edge->Skew + MassDeltaTable[PrevCellIndex];
                            TokenDropped = 0;
                        }
                        Edge = Edge->Next;
                    }
                }

                if (ModCountIndex)
                {
                    ///////////////
                    // Modification!  
                    // Remember: there may be no node corresponding to
                    // the peptide with a PTM removed.  Example: Assume the peptide is 
                    // AFKDEDTQAM+16PFR and we're at the node at 1152 for AFKDEDTQAM+16.
                    // We cannot place the M+16 PTM and jump to a node.  We must place
                    // the M+16 PTM while placing the M amino acid in order to jump.  
                    // If (due to poor fragmentation) there's no node available, then 
                    // we cannot place the PTM at all, but hopefully we will place the correct
                    // PTM mass at another node.
                    Mass = Node->Mass - AAMass;
                    BackNodeIndex = NodeIndex;
                    BackNodeDirection = -1;
                    BackNode = Node;
                    while (1)
                    {
                        // Bouncing iteration: Iterate back along the node list until
                        // you hit the start of the list (or mass becomes too small).
                        // Then iterate forward along the list until you hit the end of the
                        // list (or mass becomes too large).  We iterate over a "neighborhood"
                        // to save time.
                        if (BackNodeDirection < 0)
                        {
                            BackNode = BackNode->Prev;
                            BackNodeIndex--;
                            if (!BackNode)
                            {
                                BackNodeDirection = 1;
                                BackNodeIndex = NodeIndex;
                                BackNode = Node;
                            }
                            else
                            {
                                Delta = Mass - BackNode->Mass;
                                if (Delta > MaxPossibleDelta)
                                {
                                    BackNodeDirection = 1;
                                    BackNodeIndex = NodeIndex;
                                    BackNode = Node;
                                }
                            }
                        }
                        if (BackNodeDirection > 0)
                        {
                            BackNode = BackNode->Next;
                            BackNodeIndex++;
                            if (!BackNode)
                            {
                                break;
                            }
                            else
                            {
                                Delta = Mass - BackNode->Mass;
                                if (Delta < MinPossibleDelta)
                                {
                                    break;
                                }
                            }
                        }
                        
                        ROUND_MASS_TO_DELTA_BIN(Delta, DeltaBin);
                        DeltaNode = MassDeltaByMass[AA][DeltaBin];
                        while (DeltaNode)
                        {
                            Skew = Delta - DeltaNode->Delta->RealDelta;
                            AbsSkew = abs(Skew) / 10;
                            if (abs(Skew) <= GlobalOptions->Epsilon)
                            {
                                PrevCellIndex = AminoBlock - AminoBlockSize + (BackNodeIndex * ZSize) + ModCountIndex - 1;
                                Score = g_SkewPenalty[AbsSkew] + rint(DeltaNode->Delta->Score * DELTA_SCORE_SCALER) + ScoreTable[PrevCellIndex];
                                if (Score > ScoreTable[CellIndex])
                                {
                                    ScoreTable[CellIndex] = Score;
                                    PrevCellTable[CellIndex] = PrevCellIndex;
                                    DeltaTable[CellIndex] = DeltaNode->Delta;
                                    MassDeltaTable[CellIndex] = MassDeltaTable[PrevCellIndex] + Skew;
                                    TokenDropped = 0;
                                }
                            }
                            DeltaNode = DeltaNode->Next;
                        }
                    } // loop over back-nodes
                    
                } // if ModCount
                //////////////////////////////////////////////////////////////////////////////
                // We now have our move backwards (or our FORBIDDEN_PATH).  Get points for this node's PRM:
                if (!TokenDropped)
                {
                    Mass = MASS_TO_BIN(Node->Mass + MassDeltaTable[CellIndex]);
                    Mass = min(MaxPRM, max(0, Mass));
                    ScoreTable[CellIndex] += Tweak->PRMScores[Mass];
                }
                CellIndex++;
            } // ModCountIndex loop
        } // NodeIndex loop
    } // AminoIndex loop

    if (Info->VerboseFlag)
    {
        DebugPrintGeneralTable(Spectrum, Buffer, MaxAmino, MaxMods,
            ScoreTable, PrevCellTable);
    }

    ////////////////////////////////////////////////////////////////////////////////////
    // The d.p. table has been populated.  Now we must read off the candidate(s).
    for (Node = Spectrum->Graph->LastNode; Node; Node = Node->Prev)
    {
        NodeIndex = Node->Index;
        // We want Node->Mass to be about equal to Tweak->ParentMass - 1900
        // If it's too small, STOP.  If it's too large, CONTINUE.
        if (Node->Mass > Tweak->ParentMass - 1900 + GlobalOptions->ParentMassEpsilon)
        {
            continue;
        }
        if (Node->Mass < Tweak->ParentMass - 1900 - GlobalOptions->ParentMassEpsilon)
        {
            break;
        }
        for (AminoIndex = 0; AminoIndex < MaxAmino; AminoIndex++)
        {
            CellIndex = AminoIndex * AminoBlockSize + NodeIndex * ZSize + MaxMods;
            Score = ScoreTable[CellIndex] + EndPointPenalty[AminoIndex];
            if (Score > ScoreToBeat)
            {
                ///////////////////////////////////////
                // We have a match.  Read off the sequence and the PTMs.
                X = NodeIndex;
                Y = AminoIndex;
                Z = MaxMods;
                if (Info->VerboseFlag)
                {
                    printf("\nMatch found at (%d, %d, %d)\n", X, Y, Z);
                }
                Match = NewPeptideNode();
                Match->Tweak = Tweak;
                MatchBufferPos = 0;
                MatchPTMCount = 0;
                while (1)
                {
                    if (Info->VerboseFlag)
                    {
                        printf("Move to (%d, %d, %d), match is '%s'\n", X, Y, Z, Match->Bases);
                    }
                    //CellIndex = Y * AminoBlock + Y * ZSize + Z;
                    if (Y)
                    {
                        if (DeltaTable[CellIndex])
                        {
                            Match->AminoIndex[MatchPTMCount] = MatchBufferPos;
                            Match->ModType[MatchPTMCount] = DeltaTable[CellIndex];
                            MatchPTMCount++;
                            if (Info->VerboseFlag)
                            {
                                printf("Apply PTM '%s' (%d)\n", DeltaTable[CellIndex]->Name, DeltaTable[CellIndex]->RealDelta);
                            }
                        }
                    }
                    CellIndex = PrevCellTable[CellIndex];
                    if (CellIndex < 0)
                    {
                        break;
                    }
                    NextY = CellIndex / AminoBlockSize;
                    X = (CellIndex - NextY * AminoBlockSize) / ZSize;
                    Z = (CellIndex - NextY * AminoBlockSize) % ZSize;
                    while (Y > NextY)
                    {
                        Match->Bases[MatchBufferPos] = Buffer[Y - 1];
                        MatchBufferPos++;
                        Y--;
                    }
                }
                Match->FilePos = FilePos + Y;
                Match->RecordNumber = Info->RecordNumber;
                Match->Bases[MatchBufferPos] = '\0';
                Match->SuffixAmino = Buffer[AminoIndex];
                if (Y)
                {
                    Match->PrefixAmino = Buffer[Y - 1];
                }
                ReverseString(Match->Bases);
                for (ModIndex = 0; ModIndex < MatchPTMCount; ModIndex++)
                {
                    Match->AminoIndex[ModIndex] = MatchBufferPos - 1 - Match->AminoIndex[ModIndex];
                }
                Match->DB = Info->DB;
                //Match->Score = Score;
                Match->InitialScore = Score;
                Match = StoreSpectralMatch(Spectrum, Match, MatchBufferPos, 0);
                if (Info->VerboseFlag && Match)
                {
                    printf("Store match '%c.%s.%c' score %d endpointpenalty %d\n", Match->PrefixAmino, Match->Bases, Match->SuffixAmino, Match->InitialScore, EndPointPenalty[AminoIndex]);
                }
            } // final AminoIndex loop
        } // final node loop
    }
    return ReturnAmino;
}

// The MS-Alignment algorithm.
// New and improved version of the d.p. algorithm for generating candidates with 0-2 PTMs
int ExtendMatchRightwardDuo(SearchInfo* Info, char* Buffer, int BufferEnd, int MatchMass, 
    int MaxMods, int ScoreToBeat, int FilePos, SpectrumTweak* Tweak)
{
    static int* PrevCellTable = NULL;
    static int* DTable = NULL;
    int MaxAmino;
    int ReturnAmino = -1;
    int Result;
    int NodeIndex;
    TagGraphNode* Node;
    int AminoBlockSize;
    int AminoBlock;
    int CellIndex;
    int PrevCellIndex = 0;
    int Score;
    int CellMass;
    TagGraphBackEdge* Edge;
    int StartAminoIndex;
    int DeltaBin;
    int Delta;
    int MinPossibleDelta = -13000; // W->G mutation
    int MaxPossibleDelta = GlobalOptions->MaxPTMDelta * 100;
    MassDeltaNode* DeltaNode;
    int AA;
    int AAMass = 0;
    int AA2;
    int AA2Mass;
    int AA3;
    int EndAminoIndex;
    int ComplementMass;
    int CellScore;
    int AminoIndex;
    int Skew;
    int AbsSkew;
    int Mass;
    int PRM;
    int BackEdgeDoubleIndex = 0;
    int BackEdgeTripleIndex = 0;
    int MaxPRM = Tweak->PRMScoreMax - 1;
    int StartPointPenalty[MAX_ROWS];
    int EndPointPenalty[MAX_ROWS];
    MSSpectrum* Spectrum = Info->Spectrum;
    //
    // Allocate tables, if necessary:
    if (!PrefixTable)
    {
        PrefixTable = (int*)malloc(sizeof(int) * MAX_ROWS * MAX_ROWS);
        SuffixTable = (int*)malloc(sizeof(int) * MAX_ROWS * MAX_ROWS);
        PrefixMassTable = (int*)malloc(sizeof(int) * MAX_ROWS * MAX_ROWS);
        SuffixMassTable = (int*)malloc(sizeof(int) * MAX_ROWS * MAX_ROWS);
        DTable = (int*)malloc(sizeof(int) * MAX_ROWS * MAX_NODES);
        // MassDeltaTable stores mass delta used in reaching a cell of DTable:
        MassDeltaTable = (int*)malloc(sizeof(int) * MAX_ROWS * MAX_NODES);
        PrevCellTable = (int*)malloc(sizeof(int) * MAX_ROWS * MAX_NODES);
        DeltaTable = (MassDelta**)malloc(sizeof(MassDelta*) * MAX_ROWS * MAX_NODES);
    }
    /////////////////////////////////////////
    // Find MaxAmino:
    Result = FindPeptideBlockEnd(Spectrum, Buffer, BufferEnd, &MaxAmino, &ReturnAmino);
    if (!Result)
    {
        // No extension necessary.  Advance the database pointer:
        return ReturnAmino;
    }

    for (AminoIndex = 0; AminoIndex < MaxAmino; AminoIndex++)
    {
        StartPointPenalty[AminoIndex] = 0;
        if (AminoIndex)
        {
            AA = Buffer[AminoIndex - 1];
            if (AA != 'R' && AA != 'K' && AA != '*')
            {
                StartPointPenalty[AminoIndex] = -500;
            }
        }
        EndPointPenalty[AminoIndex] = 0;
        AA = Buffer[AminoIndex];
        if ((AA != 'R' && AA != 'K') && (AminoIndex <= MaxAmino-1 && Buffer[AminoIndex + 1]!='*'))
        {
            EndPointPenalty[AminoIndex] = -500;
        }
    }
    /////////////////////////////////////////
    // Fill the Forward and Suffix tables:

    FillPrefixSuffixTables(Spectrum, Tweak, MatchMass, Buffer, MaxAmino, PrefixTable, SuffixTable,
        PrefixMassTable, SuffixMassTable);
#ifdef VERBOSE_DEBUGGING
    DebugPrintPrefixSuffixTables(MaxAmino, Buffer, PrefixTable, SuffixTable, PrefixMassTable, SuffixMassTable); 
#endif
    /////////////////////////////////////////
    // Fill table D[]
    AminoBlockSize = Spectrum->Graph->NodeCount;
    AminoBlock = 0;
    for (AminoIndex = 0; AminoIndex <= MaxAmino; AminoBlock += AminoBlockSize, AminoIndex++)
    {

        CellIndex = AminoBlock;
        AA2 = 0;
        AA3 = 0;
        if (AminoIndex)
        {
            AA = Buffer[AminoIndex-1] - 'A';
            AAMass = PeptideMass[Buffer[AminoIndex-1]];
        }
        if (AminoIndex>1)
        {
            AA2 = Buffer[AminoIndex-2] - 'A';
            AA2Mass = PeptideMass[Buffer[AminoIndex-2]];
            BackEdgeDoubleIndex = AA*AMINO_ACIDS + AA2;
        }
        if (AminoIndex>2)
        {
            AA3 = Buffer[AminoIndex-3] - 'A';
            BackEdgeTripleIndex = AA*676 + AA2*26 + AA3;
        }
        for (NodeIndex = 0, Node = Spectrum->Graph->FirstNode; Node; Node = Node->Next, NodeIndex++, CellIndex++)
        {
            DTable[CellIndex] = FORBIDDEN_PATH; // default
            DeltaTable[CellIndex] = NULL;
            PrevCellTable[CellIndex] = -1;
            MassDeltaTable[CellIndex] = 0;

            ///////////////
            // Free rides:
            if (Node->Mass < GlobalOptions->ParentMassEpsilon)
            {
                DTable[CellIndex] = StartPointPenalty[AminoIndex];
                PrevCellTable[CellIndex] = -1;
                DeltaTable[CellIndex] = NULL;
                MassDeltaTable[CellIndex] = Node->Mass;
                continue;
            }
            if (AminoIndex == 0)
            {
                continue; // And that's all we do on the top row.
            }
           
            ///////////////
            // One unmodified amino acid:
            Edge = Node->BackEdge[AA];
            while (Edge)
            {
                PrevCellIndex = AminoBlock - AminoBlockSize + Edge->ToNode->Index;
                Score = DTable[PrevCellIndex] + Edge->Score; 
                if (Score > DTable[CellIndex])
                {
                    DTable[CellIndex] = Score;
                    PrevCellTable[CellIndex] = PrevCellIndex;
                    DeltaTable[CellIndex] = NULL;
                    MassDeltaTable[CellIndex] = Edge->Skew + MassDeltaTable[PrevCellIndex];
                }
                Edge = Edge->Next;
            }

            ///////////////
            // Two unmodified amino acids:
            if (AminoIndex > 1)
            {
                Edge = Node->BackEdgeDouble[BackEdgeDoubleIndex];
                while (Edge)
                {
                    PrevCellIndex = AminoBlock - AminoBlockSize - AminoBlockSize + Edge->ToNode->Index;
                    // Accumulate points for the middle of the jump:
                    PRM = MASS_TO_BIN(Edge->HalfMass);
                    Score = Tweak->PRMScores[PRM] + DTable[PrevCellIndex] + Edge->Score;
                    if (Score > DTable[CellIndex])
                    {
                        DTable[CellIndex] = Score;
                        PrevCellTable[CellIndex] = PrevCellIndex;
                        DeltaTable[CellIndex] = NULL;
                        MassDeltaTable[CellIndex] = Edge->Skew + MassDeltaTable[PrevCellIndex];
                    }
                    Edge = Edge->Next;
                }
            }
            ///////////////
            // Three unmodified amino acids:
            if (AminoIndex > 2)
            {
                Edge = Node->BackEdgeTriple[BackEdgeTripleIndex];
                while (Edge)
                {
                    PrevCellIndex = AminoBlock - AminoBlockSize - AminoBlockSize - AminoBlockSize + Edge->ToNode->Index;
                    // Accumulate points for the middle of the jump:
                    PRM = MASS_TO_BIN(Edge->HalfMass);
                    Score = Tweak->PRMScores[PRM] + DTable[PrevCellIndex] + Edge->Score;
                    PRM = MASS_TO_BIN(Edge->HalfMass2);
                    Score += Tweak->PRMScores[PRM];
                    if (Score > DTable[CellIndex])
                    {
                        DTable[CellIndex] = Score;
                        PrevCellTable[CellIndex] = PrevCellIndex;
                        DeltaTable[CellIndex] = NULL;
                        MassDeltaTable[CellIndex] = Edge->Skew + MassDeltaTable[PrevCellIndex];
                    }
                    Edge = Edge->Next;
                }
            }

            ///////////////
            // No amino acid at all, or modified amino acid.  Try using the prefix StartAminoIndex...EndAminoIndex.
            // Also, try using an EMPTY prefix (the case where StartAminoIndex == AminoIndex > EndAminoIndex)
            EndAminoIndex = AminoIndex - 1;
            Mass = Node->Mass - AAMass;
            for (StartAminoIndex = AminoIndex; StartAminoIndex>0; StartAminoIndex--)
            {
                if (StartAminoIndex == AminoIndex)
                {
                    Delta = Mass; // Modification on the first amino acid of the peptide
                }
                else
                {
                    Delta = Mass - PrefixMassTable[StartAminoIndex*MAX_ROWS + EndAminoIndex];
                }
                if (Delta > MaxPossibleDelta)
                {
                    continue;
                }
                if (Delta < MinPossibleDelta)
                {
                    break;
                }
                ROUND_MASS_TO_DELTA_BIN(Delta, DeltaBin);
                DeltaNode = MassDeltaByMass[AA][DeltaBin];
                while (DeltaNode)
                {
                    Skew = Delta - DeltaNode->Delta->RealDelta;
                    //Skew = Delta - DeltaNode->RealDelta;
                    AbsSkew = abs(Skew) / 10;
                    if (AbsSkew <= GlobalOptions->Epsilon)
                    {
                        if (StartAminoIndex == AminoIndex)
                        {
                            Score = g_SkewPenalty[AbsSkew] + rint(DeltaNode->Delta->Score * DELTA_SCORE_SCALER);
                        }
                        else
                        {
                            Score = g_SkewPenalty[AbsSkew] + rint(DeltaNode->Delta->Score * DELTA_SCORE_SCALER + PrefixTable[StartAminoIndex*MAX_ROWS + EndAminoIndex]);
                        }
                        Score += StartPointPenalty[StartAminoIndex - 1];
                        //Score += Spectrum->PRMScores[PRM];
                        if (Score > DTable[CellIndex])
                        {
                            DTable[CellIndex] = Score;
                            PrevCellTable[CellIndex] = (StartAminoIndex-1) * AminoBlockSize;
                            DeltaTable[CellIndex] = DeltaNode->Delta;
                            MassDeltaTable[CellIndex] = Skew;
                        }
                    }
                    DeltaNode = DeltaNode->Next;
                }
                Skew = abs(Delta) / 10;
                if (Skew < GlobalOptions->Epsilon)
                {
                    if (StartAminoIndex > EndAminoIndex)
                    {
                        Score = g_SkewPenalty[Skew];
                    }
                    else
                    {
                        Score = g_SkewPenalty[Skew] + PrefixTable[StartAminoIndex*MAX_ROWS + EndAminoIndex];
                    }
                    Score += StartPointPenalty[StartAminoIndex - 1];
                    if (Score > DTable[CellIndex])
                    {
                        DTable[CellIndex] = Score;
                        PrevCellTable[CellIndex] = (StartAminoIndex-1) * AminoBlockSize;
                        DeltaTable[CellIndex] = NULL;
                    }
                }
            }
            //////////////////////////////////////////////////////////////////////////////
            // We now have our move backwards (or our FORBIDDEN_PATH).  Get points for this node's PRM:
            Mass = MASS_TO_BIN(Node->Mass + MassDeltaTable[CellIndex]);
            Mass = min(MaxPRM, max(0, Mass));
            DTable[CellIndex] += Tweak->PRMScores[Mass];
        }
    }
#ifdef VERBOSE_DEBUGGING
    DebugPrintDTable(Spectrum, Buffer, DTable, DeltaTable, PrevCellTable, MaxAmino); 
#endif

    /////////////////////////////////////////
    // Find candidate peptides, using tables PrefixTable, SuffixTable, and D:
    AminoBlock = 0;
    for (AminoIndex = 0; AminoIndex <= MaxAmino; AminoIndex++, AminoBlock += AminoBlockSize)
    {
        AA = Buffer[AminoIndex] - 'A'; // amimoindex + 1 - 1
        CellIndex = AminoBlock;
        for (NodeIndex = 0, Node = Spectrum->Graph->FirstNode; Node; Node = Node->Next, NodeIndex++, CellIndex++)
        {
            CellMass = Node->Mass + MassDeltaTable[CellIndex];
            CellScore = DTable[CellIndex];
            ComplementMass = MatchMass - CellMass;
            // We can end right here:
            if (abs(ComplementMass) < GlobalOptions->FlankingMassEpsilon)
            {
                Spectrum->CandidatesScored++;
                GlobalStats->CandidatesScored++;
                Score = CellScore + EndPointPenalty[AminoIndex-1];
                if (CellScore > ScoreToBeat)
                {
                    AddNewMatchDuo(Info, Tweak, Buffer, CellScore, PrevCellTable, DeltaTable, CellIndex, NULL,
                        AminoBlockSize, AminoIndex, AminoIndex, FilePos);
                    if (Spectrum->Node->MatchCount == GlobalOptions->StoreMatchCount)
                    {
                        ScoreToBeat = Spectrum->Node->LastMatch->InitialScore;
                    }
                }
            }
            for (EndAminoIndex = AminoIndex + 1; EndAminoIndex <= MaxAmino; EndAminoIndex++)
            {
                Delta = SuffixMassTable[(AminoIndex + 1)*MAX_ROWS + EndAminoIndex] - CellMass;
                
                if (Delta > MaxPossibleDelta)
                {
                    continue;
                }
                if (Delta < MinPossibleDelta)
                {
                    break;
                }
                //EndAA = Buffer[EndAminoIndex] - 'A'; // amimoindex+1-1
                // Maybe we match a suffix mass:
                if (abs(Delta) < GlobalOptions->Epsilon)
                {
                    Skew = abs(Delta) / 10;
                    Score = CellScore + SuffixTable[(AminoIndex + 1)*MAX_ROWS + EndAminoIndex] + g_SkewPenalty[Skew];
                    Score += EndPointPenalty[EndAminoIndex - 1];
                    Spectrum->CandidatesScored++;
                    GlobalStats->CandidatesScored++;
                    if (Score > ScoreToBeat)
                    {
                        AddNewMatchDuo(Info, Tweak, Buffer, Score, PrevCellTable, DeltaTable, CellIndex, NULL,
                            AminoBlockSize, AminoIndex, EndAminoIndex, FilePos);
                        if (Spectrum->Node->MatchCount == GlobalOptions->StoreMatchCount)
                        {
                            ScoreToBeat = Spectrum->Node->LastMatch->InitialScore;
                        }

                    }
                }
                ROUND_MASS_TO_DELTA_BIN(Delta, DeltaBin);
                DeltaNode = MassDeltaByMass[AA][DeltaBin];
                while (DeltaNode)
                {
                    Skew = abs(DeltaNode->Delta->RealDelta - Delta);
                    //Skew = abs(DeltaNode->RealDelta - Delta);
                    if (Skew < GlobalOptions->Epsilon)
                    {
                        Score = CellScore + rint(DeltaNode->Delta->Score * DELTA_SCORE_SCALER) + SuffixTable[(AminoIndex + 1)*MAX_ROWS + EndAminoIndex] + g_SkewPenalty[Skew / 10];
                        Score += EndPointPenalty[EndAminoIndex - 1];
                        Spectrum->CandidatesScored++;
                        GlobalStats->CandidatesScored++;
                        if (Score > ScoreToBeat)
                        {
                            AddNewMatchDuo(Info, Tweak, Buffer, Score, PrevCellTable, DeltaTable, CellIndex, DeltaNode->Delta,
                                AminoBlockSize, AminoIndex, EndAminoIndex, FilePos);
                            if (Spectrum->Node->MatchCount == GlobalOptions->StoreMatchCount)
                            {
                                ScoreToBeat = Spectrum->Node->LastMatch->InitialScore;
                            }

                        }
                    }
                    DeltaNode = DeltaNode->Next;
                }
            }
        }
    }
    return ReturnAmino;
}
Peptide* ClonePeptide(Peptide* Match)
{
    Peptide* NewMatch = NewPeptideNode();
    memcpy(NewMatch, Match, sizeof(Peptide));
    return NewMatch;
}

// AddNewMatchDuo considers an unmodified peptide whenever it considers a modification of size +1..+n or -1..-n.
// We only want to add one "undecorated" peptide, not re-add the same thing multiple times,
// Store the FilePos, StartAminoIndex, EndAminoIndex of the last match 
// We've got a new match!  It ends at EndAminoIndex, and its D-path ends at the cell CellIndex.  
void AddNewMatchDuo(SearchInfo* Info, SpectrumTweak* Tweak, char* Buffer, int Score, int* PrevCellTable, MassDelta** DeltaTable, 
    int CellIndex, MassDelta* FinalDelta, int AminoBlockSize, int AminoIndex, int EndAminoIndex,
    int FilePos)
{
    int StartAminoIndex;
    int PeptideLength;
    int OldCellIndex;
    Peptide* Match;
    Peptide* VariantMatch;
    int SlideLeftIndex = -1;
    int SlideRightIndex = -1;
    int ModCount = 0;
    int PlainMass;
    int ModdedMass;
    float RunningScore;
    float PlainScore;
    float ModdedScore;
    int Diff;
    MassDeltaNode* Node;
    int PRM;
    int BestDiff;
    MSSpectrum* Spectrum = Info->Spectrum;
    //
    Match = NewPeptideNode();
    Match->Tweak = Tweak;
    // Trace back through the d.p. table to find our starting amino index:

    StartAminoIndex = AminoIndex;
    OldCellIndex = CellIndex;
    while (OldCellIndex >= 0)
    {
        StartAminoIndex = OldCellIndex / AminoBlockSize;
        if (DeltaTable[OldCellIndex])
        {
            Match->ModType[0] = DeltaTable[OldCellIndex];
            Match->AminoIndex[0] = StartAminoIndex;
            ModCount++;
        }
        OldCellIndex = PrevCellTable[OldCellIndex];
    }
    if (ModCount)
    {
        Match->AminoIndex[0] -= (StartAminoIndex + 1);
        SlideLeftIndex = Match->AminoIndex[0];
    }
    PeptideLength = EndAminoIndex - StartAminoIndex;
    strncpy(Match->Bases, Buffer + StartAminoIndex, PeptideLength);
    Match->Bases[PeptideLength] = '\0';
    Match->InitialScore = Score;
    Match->FilePos = FilePos + StartAminoIndex;
    Match->RecordNumber = Info->RecordNumber;
    if (StartAminoIndex)
    {
        Match->PrefixAmino = Buffer[StartAminoIndex - 1];
    }
    if (FinalDelta)
    {
        Match->AminoIndex[ModCount] = AminoIndex - StartAminoIndex;
        SlideRightIndex = Match->AminoIndex[ModCount];
        Match->ModType[ModCount] = FinalDelta;
        ModCount++;
    }
    GetPeptideParentMass(Match);

    Match->SuffixAmino = Buffer[EndAminoIndex];
#ifdef VERBOSE_DEBUGGING
    DebugPrintMatch(Match);
#endif
    // STRIP DECORATION:
    // If we placed a small PTM (mass -3...4), then be sure to consider a match with no modification.
    // If we placed the PTM only to make the parent mass match up, then the modless peptide will get a 
    // better score, and we'll filter the spurious +1 modification.  (There are a few real +1 modifications,
    // such as deamidation of N, but spurious +1 modifications are much more common.
    if (FinalDelta && FinalDelta->RealDelta >= -300 && FinalDelta->RealDelta < 500)
    {
        VariantMatch = ClonePeptide(Match);
        VariantMatch->InitialScore = Score;
        VariantMatch->AminoIndex[0] = -1;
        VariantMatch->ModType[0] = NULL;
        VariantMatch->DB = Info->DB;
        StoreSpectralMatch(Spectrum, VariantMatch, PeptideLength, 0);
    }
    // SLIDE LEFT:
    // If we placed a PTM at the edge of our prefix, but the PTM could just as easily have been placed earlier,
    // then do so:
    if (SlideLeftIndex > 0)
    {
        PlainMass = 0;
        //ModdedMass = Match->ModType[0]->RealDelta;
        for (AminoIndex = 0; AminoIndex <= Match->AminoIndex[0]; AminoIndex++)
        {
            PlainMass += PeptideMass[Match->Bases[AminoIndex]];
        }
        ModdedMass = PlainMass + Match->ModType[0]->RealDelta;
        RunningScore = (float)Match->InitialScore;
        for (AminoIndex = Match->AminoIndex[0]; AminoIndex > 0; AminoIndex--)
        {
            PlainMass -= PeptideMass[Match->Bases[AminoIndex]];
            ModdedMass -= PeptideMass[Match->Bases[AminoIndex]];
            PRM = MASS_TO_BIN(PlainMass);
            if (PRM < -PRM_ARRAY_SLACK || PRM >= Tweak->PRMScoreMax)
            {
                break;
            }
            PRM = max(0, PRM);
            PlainScore = (float)Tweak->PRMScores[PRM];
            PRM = MASS_TO_BIN(ModdedMass);
            if (PRM < -PRM_ARRAY_SLACK || PRM >= Tweak->PRMScoreMax)
            {
                break;
            }
            PRM = max(0, PRM);
            ModdedScore = (float)Tweak->PRMScores[PRM];
            if (ModdedScore > 0)
            {
                // We've already had the chance to attach this ptm here.
                break;
            }
            RunningScore += (ModdedScore - PlainScore);
            if (RunningScore < Match->InitialScore - 100)
            {
                // We've hurt our score quite a bit; let's stop.
                break;
            }
            // Make a variant-match:
            VariantMatch = ClonePeptide(Match);
            VariantMatch->InitialScore = rint(RunningScore);
            BestDiff = -1;
            VariantMatch->ModType[0] = NULL;
            for (Node = MassDeltaByMass[Match->Bases[AminoIndex-1]-'A'][Match->ModType[0]->Delta]; Node; Node = Node->Next)
            {
                Diff = abs(Node->Delta->RealDelta - Match->ModType[0]->RealDelta);
                if (BestDiff < 0 || Diff < BestDiff)
                {
                    BestDiff = Diff;
                    VariantMatch->ModType[0] = Node->Delta;
                }
            }
            VariantMatch->AminoIndex[0] = AminoIndex - 1;
            if (VariantMatch->ModType[0])
            {
#ifdef VERBOSE_DEBUGGING
                printf("Variant:\n");
                DebugPrintMatch(VariantMatch);
#endif
                VariantMatch->DB = Info->DB;
                StoreSpectralMatch(Spectrum, VariantMatch, PeptideLength, 0);
            }
            else
            {
                FreePeptideNode(VariantMatch);
            }
        }
    }
    // SLIDE RIGHT:
    // If we placed a PTM at the edge of our prefix, but the PTM could just as easily have been placed earlier,
    // then do so:
    if (SlideRightIndex > 0)
    {
        PlainMass = 0;
        if (ModCount>1)
        {
            PlainMass += Match->ModType[0]->RealDelta;
        }
        //ModdedMass = Match->ModType[0]->RealDelta;
        for (AminoIndex = 0; AminoIndex < Match->AminoIndex[ModCount-1]; AminoIndex++)
        {
            PlainMass += PeptideMass[Match->Bases[AminoIndex]];
        }
        ModdedMass = PlainMass + Match->ModType[ModCount-1]->RealDelta;
        RunningScore = (float)Match->InitialScore;

        for (AminoIndex = Match->AminoIndex[ModCount-1]; Match->Bases[AminoIndex]; AminoIndex++)
        {
            PlainMass += PeptideMass[Match->Bases[AminoIndex]];
            ModdedMass += PeptideMass[Match->Bases[AminoIndex]];
            PRM = MASS_TO_BIN(PlainMass);
            if (PRM < -PRM_ARRAY_SLACK || PRM >= Tweak->PRMScoreMax)
            {
                break;
            }
            PlainScore = (float)Tweak->PRMScores[PRM];
            PRM = MASS_TO_BIN(ModdedMass);
            if (PRM < -PRM_ARRAY_SLACK || PRM >= Tweak->PRMScoreMax)
            {
                break;
            }
            ModdedScore = (float)Tweak->PRMScores[PRM];
            RunningScore += (PlainScore - ModdedScore);
            if (RunningScore < Match->InitialScore - 100)
            {
                // We've hurt our score quite a bit; let's stop.
                break;
            }
            if (AminoIndex > Match->AminoIndex[ModCount-1])
            {
                if (ModdedScore > 0)
                {
                    // We've already had the chance to attach this ptm here.
                    break;
                }

                // Make a variant-match:
                VariantMatch = ClonePeptide(Match);
                VariantMatch->InitialScore = rint(RunningScore);
                BestDiff = -1;
                VariantMatch->ModType[ModCount-1] = NULL;
                for (Node = MassDeltaByMass[Match->Bases[AminoIndex]-'A'][Match->ModType[ModCount-1]->Delta]; Node; Node = Node->Next)
                {
                    Diff = abs(Node->Delta->RealDelta - Match->ModType[ModCount-1]->RealDelta);
                    if (BestDiff < 0 || Diff < BestDiff)
                    {
                        BestDiff = Diff;
                        VariantMatch->ModType[ModCount-1] = Node->Delta;
                    }
                }
                if (VariantMatch->ModType[ModCount-1])
                {
                    VariantMatch->AminoIndex[ModCount-1] = AminoIndex;
#ifdef VERBOSE_DEBUGGING
                    printf("Variant:\n");
                    DebugPrintMatch(VariantMatch);
#endif
                    VariantMatch->DB = Info->DB;
                    StoreSpectralMatch(Spectrum, VariantMatch, PeptideLength, 0);
                }
                else
                {
                    FreePeptideNode(VariantMatch);
                }
            }
        }

    }
    Match->DB = Info->DB;
    StoreSpectralMatch(Spectrum, Match, PeptideLength, 0);
}
