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
#include <string.h>
#include <math.h> 
#include <stdlib.h>
#include "Mods.h"
#include "Score.h"
#include "Spectrum.h"
#include "Inspect.h"
#include "Tagger.h"
#include "SVM.h"
#include "IonScoring.h"
#include "ParentMass.h"

////////////////////////////////////////////////////////////////////////////////////////
// Score.c and Score.h support an 'alignment-based' scoring method which has been 
// discarded in favor of SVM-based scoring.  The two methods have similar performance, but
// the SVM was slightly better and slightly faster.  This code is kept around for reference,
// but isn't executed in practice.

// The max length permitted for any peptide match:
#define MAX_PEPTIDE_LENGTH 256

// Set SHOW_DP_TABLE to enable verbose printout of the d.p. table from scoring.
//#define SHOW_DP_TABLE 1

// TheoPeak has a mass, an ion type, and a score.  It's a peak
// in the *theoretical fragmentation spectrum*.  Some peaks (e.g.
// an a-ion peak late in the peptide) have bad scores; others get very 
// good scores.  The scores are log-odds scores.  Scores are based on the
// scores in the ScoringModel.
typedef struct TheoPeak
{
    int Mass;
    int IonType; 
    int LossType;
    int Score;
    int AntiScore;
    // CutIndex is the number of amino acids in this fragment.
    // For instance, in EAMAPK, b1 and y5 are the same cut point.
    // CutIndex is always at least 1, since a fragment has at least one amino acid.
    int CutIndex; 
    // TrueCutIndex is the index of the cut.  (So, b and y fragments with same TrueCutIndex represent
    // breakage of the same peptide bond)
    int TrueCutIndex; 
    // AssignedPeak is built during spectral scoring, when backtracking along the DP table
    SpectralPeak* AssignedPeak;
} TheoPeak;

// Our scoring model uses several type of cuts, and computes probabilities for each.
// Here's a diagram of the cut points for peptide "SPECTRUM":
// LeftEdge   L1   L2   Mid   Mid   Mid   R2   R1   RightEdge 
//          S    P    E     C     T     R    U    M
typedef enum CutPointType
{
    CutPointTypeL1 = 0,
    CutPointTypeL2,
    CutPointTypeR1,
    CutPointTypeR2,
    CutPointTypeMid,
    CutPointTypeLeftEdge,
    CutPointTypeRightEdge,
    CutPointTypeCount
} CutPointType;

// Skew odds are modeled by bins (0 to 0.05 Da, 0.05 to 0.1, up to a bin for 0.5+).
#define SKEW_BIN_COUNT 10

#define SECTOR_COUNT 3

// Scoring model, built by DPTrainer.py
typedef struct ScoringModel
{
    int BScore[CutPointTypeCount];
    int YScore[CutPointTypeCount];
    int BH2OBoostScore[CutPointTypeCount];
    int BH2OScore[CutPointTypeCount];
    int BNH3BoostScore[CutPointTypeCount];
    int BNH3Score[CutPointTypeCount];
    int ABoostScore[CutPointTypeCount];
    int AScore[CutPointTypeCount];
    int YH2OBoostScore[CutPointTypeCount];
    int YH2OScore[CutPointTypeCount];
    int YNH3Score[CutPointTypeCount];
    int B2BoostScore[CutPointTypeCount];
    int B2Score[CutPointTypeCount];
    int Y2Score[CutPointTypeCount];
    int AH2OScore[CutPointTypeCount];
    int ANH3Score[CutPointTypeCount];
    int NoisePenalty[21];
    int PresenceScore[SECTOR_COUNT];
    int AbsenceScore[SECTOR_COUNT];
    int InexplicableScore[SECTOR_COUNT];
    int SkewScore[SKEW_BIN_COUNT];
} ScoringModel;

// We have an array of models - one for each charge state (1, 2, 3+)
ScoringModel* Models;

// P-values are assigned based on a histogram of Match Quality Scores for false matches.
// Histogram bin-count (and edges) are hard-coded.
#define PVALUE_BIN_COUNT 300
#define PVALUE_BIN_BOOST 100
float g_MatchPValueShort[PVALUE_BIN_COUNT];
float g_MatchPValueMedium[PVALUE_BIN_COUNT];
float g_MatchPValueLong[PVALUE_BIN_COUNT];
float g_MatchPValueLongLong[PVALUE_BIN_COUNT];

int InitPValue(char* FileName)
{
    FILE* File;
    //
    File = fopen(FileName, "rb");
    if (!File)
    {
        return 0;
    }
    ReadBinary(g_MatchPValueShort, sizeof(float), PVALUE_BIN_COUNT, File);
    ReadBinary(g_MatchPValueMedium, sizeof(float), PVALUE_BIN_COUNT, File);
    ReadBinary(g_MatchPValueLong, sizeof(float), PVALUE_BIN_COUNT, File);
    ReadBinary(g_MatchPValueLongLong, sizeof(float), PVALUE_BIN_COUNT, File);
    fclose(File);
    return 1;
}


// For debug output: Return a description of an ion type code.
char* GetIonTypeName(int IonType)
{
    switch (IonType)
    {
    case evIonTypeNone:
        return "None";
    case evIonTypeB:
        return "B";
    case evIonTypeY:
        return "Y";
    case evIonTypeA:
        return "A";
    case evIonTypeBH2O:
        return "B-H2O";
    case evIonTypeAH2O:
        return "A-H2O";
    case evIonTypeBNH3:
        return "B-NH3";
    case evIonTypeANH3:
        return "A-NH3";
    case evIonTypeYH2O:
        return "Y-H2O";
    case evIonTypeYNH3:
        return "Y-NH3";
    case evIonTypeB2:
        return "B2";
    case evIonTypeY2:
        return "Y2";
    case evIonTypeNoise:
        return "<noise>";
    case evIonTypeBPhos:
        return "b-p";
    case evIonTypeYPhos:
        return "y-p";
    default:
        return "BROKEN****";
    }
}

// For debug output: Return a description of an ion type code.
char* GetShortIonTypeName(int IonType)
{
    switch (IonType)
    {
    case evIonTypeNone:
        return "-";
    case evIonTypeB:
        return "b";
    case evIonTypeY:
        return "y";
    case evIonTypeA:
        return "a";
    case evIonTypeBH2O:
        return "b-H2O";
    case evIonTypeAH2O:
        return "a-H2O";
    case evIonTypeBNH3:
        return "b-NH3";
    case evIonTypeANH3:
        return "a-NH3";
    case evIonTypeYH2O:
        return "y-H2O";
    case evIonTypeYNH3:
        return "y-NH3";
    case evIonTypeB2:
        return "b2";
    case evIonTypeY2:
        return "y2";
    case evIonTypeNoise:
        return "<noise>";
    case evIonTypeBPhos:
        return "b-p";
    case evIonTypeYPhos:
        return "y-p";
    default:
        return "*";
    }
}

// For debugging: Print a theoretical fragmentation spectrum.  
void DebugPrintPeaks(TheoPeak* Peaks, int PeakCount)
{
    int PeakIndex;
    //
    printf("\n-----Peak list:-----\n");
    for (PeakIndex = 0; PeakIndex < PeakCount; PeakIndex++)
    {
        printf("%d: m/z %.2f %s %d %d\n", PeakIndex, Peaks[PeakIndex].Mass / (float)MASS_SCALE, GetIonTypeName(Peaks[PeakIndex].IonType),
            Peaks[PeakIndex].CutIndex, Peaks[PeakIndex].Score);
    }
}

// Theoretical spectrum builder - simple struct for remembering sector edges
// and the current proline bonus.
typedef struct TheoPeakBuilder
{
    int SectorEdgeA;
    int SectorEdgeB;
    int BProlineBonus;
    int YProlineBonus;
    int MatchLength;
} TheoPeakBuilder;

// Add a new theoretical peak to Peaks.
void AddTheoPeak(TheoPeakBuilder* Theo, TheoPeak* Peaks, int PeakCount, int IonType, int LossType, int Mass,
                 int CutIndex, int PrefixFlag, ScoringModel* Model, int Score)
{
    int SectorNumber;
    //
    Peaks[PeakCount].IonType = IonType;
    Peaks[PeakCount].LossType = LossType;
    Peaks[PeakCount].Mass = Mass;
    Peaks[PeakCount].CutIndex = CutIndex;
    Peaks[PeakCount].Score = Score;
    if (PrefixFlag)
    {
        Peaks[PeakCount].TrueCutIndex = CutIndex;
        Peaks[PeakCount].Score = min(0, Peaks[PeakCount].Score + Theo->BProlineBonus);
    }
    else
    {
        Peaks[PeakCount].TrueCutIndex = Theo->MatchLength - CutIndex;
        Peaks[PeakCount].Score = min(0, Peaks[PeakCount].Score + Theo->YProlineBonus);
    }
    if (Mass > Theo->SectorEdgeB)
    {
        SectorNumber = 2;
    }
    else if (Mass > Theo->SectorEdgeA)
    {
        SectorNumber = 1;
    }
    else
    {
        SectorNumber = 0;
    }
    // Compare against the null model right here:
    Peaks[PeakCount].Score -= Model->PresenceScore[SectorNumber];
}

// Give a bonus for peptides that seem plausible given the cleavage type.
// For instance: If our sample was subjected to trypsin digest, then *most* fragments will
// end in K or R (and be preceded by a K or R), so give a bonus to such peptides.
// The DEFAULT behavior is to assume no digest and give no points.
// This code is NO LONGER USED in production; instead, number of tryptic termini (NTT) is used
// as one feature for LDA.
int ApplyDigestBonus(Peptide* Match)
{
    int Score = 0;
    int AminoCount;
    int AminoIndex;
    // SWT 3/9/5: Use a somewhat HEAVIER penalty for broken endpoints.  And, penalize two bad endpoints
    // superadditively
    int MissedCleavagePenalty = 100;
    int BrokenSuffixPenalty = 550;
    int BrokenPrefixPenalty = 550;
    int BrokenBothPenalty = 400; // extra penalty if both endpoints broken
    char MutantBases[256];
    int BrokenTermini = 0;
    int ModIndex;
    // Write MutantBases string.  This contains the *real* amino acids of the match, with any mutations
    // applied.  For instance, if Match->Bases is "EAMAPK" but match has an M->Q mutation in position 2,
    // then MutantBases will be "EAQAPK".  
    strcpy(MutantBases, Match->Bases);
    if (!GlobalOptions->TaglessSearchFlag)
    {
        for (ModIndex = 0; ModIndex < MAX_PT_MODS; ModIndex++)
        {
            if (!Match->ModType[ModIndex])
            {
                break;
            }
            if (Match->ModType[ModIndex]->Amino)
            {
                MutantBases[Match->AminoIndex[ModIndex]] = Match->ModType[ModIndex]->Amino;
            }
        }
    }
    switch (GlobalOptions->DigestType)
    {
    case 0:
        // No digest (or unknown digest), so no points
        return 0;
    case 1:
        // A tryptic peptide gets a minor bonus, and 
        // missed cleavages get a minor penalty.
        AminoCount = strlen(MutantBases);
        for (AminoIndex = 1; AminoIndex < AminoCount - 1; AminoIndex++)
        {
            if ((MutantBases[AminoIndex] == 'K' || MutantBases[AminoIndex] == 'R') && (MutantBases[AminoIndex + 1]!='P'))
                Score -= MissedCleavagePenalty;
        }
        if (MutantBases[AminoCount - 1] != 'K' && MutantBases[AminoCount - 1] != 'R')
        {
            Score -= BrokenSuffixPenalty;
            BrokenTermini++;
        }
        if (Match->PrefixAmino && (Match->PrefixAmino!='K' && Match->PrefixAmino!='R'))
        {
            Score -= BrokenPrefixPenalty;
            BrokenTermini++;
        }
        if (BrokenTermini==2)
        {
            Score -= BrokenBothPenalty;
        }
        return Score;
    case 2:
        // Chymotrypsin: Cleaves C-terminal side of FYWL (if not followed by P)
        AminoCount = strlen(MutantBases);
        for (AminoIndex = 1; AminoIndex < AminoCount - 1; AminoIndex++)
        {
            if ((MutantBases[AminoIndex] == 'F' || MutantBases[AminoIndex] == 'Y' ||
                MutantBases[AminoIndex] == 'W' || MutantBases[AminoIndex] == 'L') && (MutantBases[AminoIndex + 1]!='P'))
                Score -= MissedCleavagePenalty;
        }
        if (MutantBases[AminoCount - 1] != 'F' && MutantBases[AminoCount - 1] != 'Y' &&
            MutantBases[AminoCount - 1] != 'W' && MutantBases[AminoCount - 1] != 'L')
        {
            BrokenTermini++;
            Score -= BrokenSuffixPenalty;
        }
        if (Match->PrefixAmino && (Match->PrefixAmino != 'F' && Match->PrefixAmino != 'Y' &&
            Match->PrefixAmino != 'W' && Match->PrefixAmino != 'L'))
        {
            BrokenTermini++;
            Score -= BrokenPrefixPenalty;
        }
        if (BrokenTermini==2)
        {
            Score -= BrokenBothPenalty;
        }

        return Score;
    case 3:
        // Lys-C - similar to trypsin.  Cleaves after K if not before P.
        // missed cleavages get a minor penalty.
        AminoCount = strlen(MutantBases);
        for (AminoIndex = 1; AminoIndex < AminoCount - 1; AminoIndex++)
        {
            if ((MutantBases[AminoIndex] == 'K') && (MutantBases[AminoIndex+1]!='P'))
                Score -= MissedCleavagePenalty;
        }
        if (MutantBases[AminoCount - 1] != 'K')
        {
            Score -= BrokenSuffixPenalty;
            BrokenTermini++;
        }
        if (Match->PrefixAmino && (Match->PrefixAmino!='K'))
        {
            Score -= BrokenPrefixPenalty;
            BrokenTermini++;
        }
        if (BrokenTermini==2)
        {
            Score -= BrokenBothPenalty;
        }

        return Score;
    case 4:
        // Asp-N - Cleaves before (on N-terminal side of) DE
        AminoCount = strlen(MutantBases);
        // Penalty for missed cleavages:
        for (AminoIndex = 1; AminoIndex < AminoCount - 1; AminoIndex++)
        {
            if ((MutantBases[AminoIndex] == 'D') || (MutantBases[AminoIndex]=='E'))
                Score -= MissedCleavagePenalty;
        }
        if (Match->SuffixAmino && (Match->SuffixAmino!='K' && Match->SuffixAmino!='E'))
        {
            Score -= BrokenSuffixPenalty;
            BrokenTermini++;
        }
        if (MutantBases[0]!='D' && MutantBases[0]!='E')
        {
            Score -= BrokenPrefixPenalty;
            BrokenTermini++;
        }
        if (BrokenTermini==2)
        {
            Score -= BrokenBothPenalty;
        }

        return Score;
    case 5:
        // GluC cleaves c-terminal of E
        AminoCount = strlen(MutantBases);
        for (AminoIndex = 1; AminoIndex < AminoCount - 1; AminoIndex++)
        {
            if (MutantBases[AminoIndex] == 'E')
            {
                Score -= MissedCleavagePenalty;
            }
        }
        if (MutantBases[AminoCount - 1] != 'E')
        {
            Score -= BrokenSuffixPenalty;
            BrokenTermini++;
        }
        if (Match->PrefixAmino && (Match->PrefixAmino!='E'))
        {
            Score -= BrokenPrefixPenalty;
            BrokenTermini++;
        }
        if (BrokenTermini==2)
        {
            Score -= BrokenBothPenalty;
        }
        return Score;
    default:
        printf("Unknown digest type '%d' encountered, no scoring adjustment applied.\n", GlobalOptions->DigestType);
        return 0;
    }
}

int DiffPeptides(char* AA1, char* AA2)
{
    int DiffCount = 0;
    while (*AA1 && *AA2)
    {
        if (*AA1 != *AA2)
        {
            DiffCount++;
        }
        AA1++;
        AA2++;
    }
    return DiffCount;
}

void SetMatchDeltaCN(SpectrumNode* Spectrum)
{
    Peptide* Match;
    Peptide* OtherMatch;
    int MatchNumber = 0;
    int IsSame;

    // Init DeltaCN and DeltaCNOther:
    for (Match = Spectrum->FirstMatch; Match; Match = Match->Next)
    {
        Match->DeltaCN = (float)FORBIDDEN_PATH;
        Match->DeltaCNOther = (float)FORBIDDEN_PATH;
    }

    // Properly set DeltaCN and DeltaCNOther:
    for (Match = Spectrum->FirstMatch; Match; Match = Match->Next)
    {
        MatchNumber++;
        if (Match != Spectrum->FirstMatch)
        {
            Match->DeltaCN = Match->MatchQualityScore - Spectrum->FirstMatch->MatchQualityScore;
        }
        else
        {
            if (Match->Next)
            {
                Match->DeltaCN = Match->MatchQualityScore - Match->Next->MatchQualityScore;
            }
            else
            {
                Match->DeltaCN = max(0, Match->MatchQualityScore);
            }
        }
        // If this match is already dissimilar to a higher-scoring one, stop now:
        if (Match->DeltaCNOther != FORBIDDEN_PATH)
        {
            continue;
        }
        if (Match->FamilyLeader)
        {
            Match->DeltaCNOther = Match->FamilyLeader->DeltaCNOther + (Match->MatchQualityScore - Match->FamilyLeader->MatchQualityScore);
            continue;
        }
        if (MatchNumber > GlobalOptions->ReportMatchCount)
        {
            // We won't bother computing DeltaCNOther for any poorer matches, because we'll drop them anyway.
            break;
        }
        for (OtherMatch = Match->Next; OtherMatch; OtherMatch = OtherMatch->Next)
        {
            IsSame = 0;
            if (abs(Match->FilePos - OtherMatch->FilePos) < 3)
            {
                IsSame = 1;
            }
            if (DiffPeptides(Match->Bases, OtherMatch->Bases) < 2)
            {
                IsSame = 1;
            }
            if (DiffPeptides(Match->Bases, OtherMatch->Bases + 1) < 2)
            {
                IsSame = 1;
            }
            if (DiffPeptides(Match->Bases + 1, OtherMatch->Bases) < 2)
            {
                IsSame = 1;
            }
            if (IsSame)
            {
                OtherMatch->FamilyLeader = Match;
            }
            else
            {
                OtherMatch->DeltaCNOther = OtherMatch->MatchQualityScore - Match->MatchQualityScore;
                if (Match->DeltaCNOther == FORBIDDEN_PATH)
                {
                    Match->DeltaCNOther = Match->MatchQualityScore - OtherMatch->MatchQualityScore;
                }
            }
        }
        if (Match->DeltaCNOther == FORBIDDEN_PATH)
        {
            if (Match == Spectrum->LastMatch)
            {
                Match->DeltaCNOther = max(Match->MatchQualityScore, 0);
            }
            else
            {
                Match->DeltaCNOther = Match->MatchQualityScore - Spectrum->LastMatch->MatchQualityScore;
            }
        }
    }
}


// Get PeptideMatchFeatures having to do with cut scores (mean, median...)
int GetCutScorePeptideMatchFeatures(MSSpectrum* Spectrum, Peptide* Match, float* FeatureArray, PRMBayesianModel* Model)
{
    int FeatureIndex = 0;
    float CutScores[256];
    int PRMCount;
    int AminoIndex;
    float ScoreTotal;
    int PeptideLength;
    //
    PeptideLength = strlen(Match->Bases);
    //for (NodeIndex = 0, Node = Model->Head; Node; NodeIndex++, Node = Node->Next)
    //{
    //    PRM = 0;
    //    for (AminoIndex = 0; AminoIndex <= PeptideLength; AminoIndex++)
    //    {
    //        ///////////////////////////////////////////////////////////////////////////////////////
    //        // Set values, and accumulate table entries:
    //        Node->Values[AminoIndex] = IonScoringGetNodeValue(Model, Node, Spectrum, PRM, Match, AminoIndex);
    //        ///////////////////////////////////////////////////////////////////////////////////////
    //        // Add to PRM:
    //        if (AminoIndex == PeptideLength)
    //        {
    //            break;
    //        }
    //        PRM += PeptideMass[Match->Bases[AminoIndex]];
    //        for (ModIndex = 0; ModIndex < MAX_PT_MODS; ModIndex++)
    //        {
    //            if (Match->AminoIndex[ModIndex] == AminoIndex)
    //            {
    //                PRM += Match->ModType[ModIndex]->RealDelta;
    //            }
    //        }
    //    } // Amino loop
    //} // NodeIndex loop

    //// Populate the CutScores array:
    //for (AminoIndex = 0; AminoIndex <= PeptideLength; AminoIndex++)
    //{
    //    CutScores[AminoIndex] = PRMBNGetCutScore(Spectrum, Model, AminoIndex);
    //}
    PopulateCutScores(Model, Spectrum, Match, CutScores);

    // Compute features based upon cut scores:
    // Total/mean for ALL cut scores:
    ScoreTotal = 0;
    PRMCount = 0;
    for (AminoIndex = 0; AminoIndex <= PeptideLength; AminoIndex++)
    {
        ScoreTotal += CutScores[AminoIndex];
        PRMCount++;
    }
    FeatureArray[FeatureIndex++] = ScoreTotal;
    FeatureArray[FeatureIndex++] = ScoreTotal / (float)PRMCount;

    // Total/mean for CENTRAL cut scores:
    ScoreTotal = 0;
    PRMCount = 0;
    for (AminoIndex = 1; AminoIndex < PeptideLength; AminoIndex++)
    {
        ScoreTotal += CutScores[AminoIndex];
        PRMCount++;
    }
    FeatureArray[FeatureIndex++] = ScoreTotal;
    FeatureArray[FeatureIndex++] = ScoreTotal / (float)max(1, PRMCount);

    // Total/mean for CENTRAL cut scores:
    ScoreTotal = 0;
    PRMCount = 0;
    for (AminoIndex = 2; AminoIndex < (PeptideLength - 1); AminoIndex++)
    {
        ScoreTotal += CutScores[AminoIndex];
        PRMCount++;
    }
    FeatureArray[FeatureIndex++] = ScoreTotal;
    FeatureArray[FeatureIndex++] = ScoreTotal / (float)max(1, PRMCount);

    // Median cut score:
    PRMCount = PeptideLength + 1;
    FeatureArray[FeatureIndex++] = GetMedian(CutScores + 2, PRMCount - 4);
    FeatureArray[FeatureIndex++] = GetMedian(CutScores + 1, PRMCount - 2);
    FeatureArray[FeatureIndex++] = GetMedian(CutScores, PRMCount);

    return FeatureIndex;
}

// Helper for GetPeptideMatchFeaturesFull: Compute features having to do with the percentage of peaks
// and peak intensity explained by the match.
int GetExplainedPeakPeptideMatchFeatures(MSSpectrum* Spectrum, Peptide* Match, float* FeatureArray)
{
    int PeakIndex;
    float IntensityB = 0;
    float IntensityY = 0;
    float IntensityBSeries = 0;
    float IntensityYSeries = 0;
    float TotalIntensity = 0;
    int PeakCountB = 0;
    int PeakCountY = 0;
    int StrongPeakCountB = 0;
    int StrongPeakCountY = 0;
    float WeightedPeakCountTotal = 0;
    float WeightedPeakCountB = 0;
    float WeightedPeakCountY = 0;
    int StrongPeakCount;
    int FeatureIndex = 0;
    int FragmentType;
    float PeakIntensity;
    float WeightedPeakIndex;
    int BFlag[256];
    int YFlag[256];
    int PeptideLength;
    int PresentCount;
    int AminoIndex;
    //
    PeptideLength = strlen(Match->Bases);
    memset(BFlag, 0, sizeof(int) * (PeptideLength + 1));
    memset(YFlag, 0, sizeof(int) * (PeptideLength + 1));
    StrongPeakCount = PeptideLength * 2;
    for (PeakIndex = 0; PeakIndex < Spectrum->PeakCount; PeakIndex++)
    {
        FragmentType = Spectrum->Peaks[PeakIndex].IonType;
        PeakIntensity = Spectrum->Peaks[PeakIndex].Intensity;
        TotalIntensity += PeakIntensity;
        WeightedPeakIndex = (float)1.0 / (Spectrum->Peaks[PeakIndex].IntensityRank + 1);
        WeightedPeakCountTotal += WeightedPeakIndex;
        switch (FragmentType)
        {
            case evFragmentY:
                PeakCountY++;
                IntensityY += PeakIntensity;
                IntensityYSeries += PeakIntensity;
                WeightedPeakCountY += WeightedPeakIndex;
                YFlag[Spectrum->Peaks[PeakIndex].AminoIndex] = 1;
                if (Spectrum->Peaks[PeakIndex].IntensityRank < StrongPeakCount)
                {
                    StrongPeakCountY++;
                }
                break;
            case evFragmentYLoss:
                IntensityYSeries += PeakIntensity;
                break;
            case evFragmentB:
                PeakCountB++;
                IntensityB += PeakIntensity;
                IntensityBSeries += PeakIntensity;
                WeightedPeakCountB += WeightedPeakIndex;
                BFlag[Spectrum->Peaks[PeakIndex].AminoIndex] = 1;
                if (Spectrum->Peaks[PeakIndex].IntensityRank < StrongPeakCount)
                {
                    StrongPeakCountB++;
                }
                break;
            case evFragmentBLoss:
                IntensityBSeries += PeakIntensity;
                break;
        }
    }
    // Fraction of B, Y present:
    PresentCount = 0;
    for (AminoIndex = 0; AminoIndex <= PeptideLength; AminoIndex++)
    {
        PresentCount += YFlag[AminoIndex];
    }
    FeatureArray[FeatureIndex++] = PresentCount / (float)(PeptideLength + 1);
    PresentCount = 0;
    for (AminoIndex = 0; AminoIndex <= PeptideLength; AminoIndex++)
    {
        PresentCount += BFlag[AminoIndex];
    }
    FeatureArray[FeatureIndex++] = PresentCount / (float)(PeptideLength + 1);
    PresentCount = 0;
    for (AminoIndex = 1; AminoIndex < PeptideLength; AminoIndex++)
    {
        PresentCount += YFlag[AminoIndex];
    }
    FeatureArray[FeatureIndex++] = PresentCount / (float)(PeptideLength - 1);
    PresentCount = 0;
    for (AminoIndex = 1; AminoIndex < PeptideLength; AminoIndex++)
    {
        PresentCount += BFlag[AminoIndex];
    }
    FeatureArray[FeatureIndex++] = PresentCount / (float)(PeptideLength - 1);

    // Fraction of top peaks:
    FeatureArray[FeatureIndex++] = (StrongPeakCountY + StrongPeakCountB) / (float)StrongPeakCount;
    FeatureArray[FeatureIndex++] = StrongPeakCountY / (float)StrongPeakCount;
    FeatureArray[FeatureIndex++] = StrongPeakCountB / (float)StrongPeakCount;

    FeatureArray[FeatureIndex++] = (WeightedPeakCountY + WeightedPeakCountB) / WeightedPeakCountTotal;
    FeatureArray[FeatureIndex++] = WeightedPeakCountY / WeightedPeakCountTotal;
    FeatureArray[FeatureIndex++] = WeightedPeakCountB / WeightedPeakCountTotal;
    
    // Fraction of intensity:
    FeatureArray[FeatureIndex++] = (IntensityY + IntensityB) / TotalIntensity;
    FeatureArray[FeatureIndex++] = IntensityY / TotalIntensity;
    FeatureArray[FeatureIndex++] = IntensityB / TotalIntensity;

    // Fraction of intensity:
    FeatureArray[FeatureIndex++] = (IntensityYSeries + IntensityBSeries) / TotalIntensity;
    FeatureArray[FeatureIndex++] = IntensityYSeries / TotalIntensity;
    FeatureArray[FeatureIndex++] = IntensityBSeries / TotalIntensity;

    return FeatureIndex;
}

// Compute features rating the quality of this annotation for the spectrum.
// Set feature values in FeatureArray, return the number of features set.
int GetPeptideMatchFeaturesFull(MSSpectrum* Spectrum, Peptide* Match, float* FeatureArray)
{
    int FeatureIndex = 0;
    int PeptideLength;
    PRMBayesianModel* Model;
    PMCSpectrumInfo* SpectrumInfo;
    PMCInfo* Info;
    //

    // Length:
    PeptideLength = strlen(Match->Bases);
    FeatureArray[FeatureIndex++] = (float)PeptideLength;
    
    if (Spectrum->Charge < 3)
    {
        Model = TAGModelCharge2;
    }
    else
    {
        Model = TAGModelCharge3;
    }

    Spectrum->ParentMass = GetPeptideParentMass(Match);

    // Compute cut scores:
    FeatureIndex += GetCutScorePeptideMatchFeatures(Spectrum, Match, FeatureArray + FeatureIndex, Model);

    // Compute features based on the fraction of top peaks / intensity explained:
    FeatureIndex += GetExplainedPeakPeptideMatchFeatures(Spectrum, Match, FeatureArray + FeatureIndex);

    ///////////////////////////////
    // Spectral convolution:
    SpectrumInfo = GetPMCSpectrumInfo(Spectrum);
    Info = (PMCInfo*)calloc(1, sizeof(PMCInfo));
    Info->Charge = SpectrumInfo->Charge;
    Info->ParentMass = Spectrum->ParentMass;
    SpectrumInfo->Head = Info;
    SpectrumInfo->Tail = Info;
    ConvolveMassCorrectedSpectrum(Info, SpectrumInfo);
    FeatureArray[FeatureIndex++] = Info->Convolve[2];
    FeatureArray[FeatureIndex++] = Info->Convolve2[0];
    FreePMCSpectrumInfo(SpectrumInfo);

    /////////////////////////////////
    // Number of tryptic termini:
    FeatureArray[FeatureIndex++] = (float)CountTrypticTermini(Match);

    ////////////////////////////////
    // Fancy length feature:
    FeatureArray[FeatureIndex++] = (float)log(max(1, PeptideLength - 5));
    FeatureArray[FeatureIndex++] = (float)log(max(1, PeptideLength - 4));
    FeatureArray[FeatureIndex++] = (float)log(max(1, PeptideLength - 3));

    return FeatureIndex;
}

float GetPenalizedScore(MSSpectrum* Spectrum, Peptide* Match, float Score)
{
    int ModIndex;
    if (strlen(Match->Bases) < MIN_VALID_PEPTIDE_LENGTH)
    {
        Score -= 1.0;
    }
    for (ModIndex = 0; ModIndex < MAX_PT_MODS; ModIndex++)
    {
        if (Match->ModType[ModIndex])
        {
            Score -= 0.25;
        }
    }
    return Score;
}
