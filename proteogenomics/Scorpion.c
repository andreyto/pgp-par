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
#include <stdlib.h>
#include <limits.h>
#include <assert.h>
#include <math.h> 
#include "Mods.h"
#include "Scorpion.h"
#include "Spectrum.h"
#include "Inspect.h"
#include "Tagger.h"
#include "SVM.h"
#include "BN.h"
#include "FreeMod.h"
#include "IonScoring.h"

#define PRM_FEATURE_COUNT 32
#define MAX_PEPTIDE_LENGTH 256
#define INTENSITY_LEVEL_COUNT 4

int g_CutFeatures[MAX_PEPTIDE_LENGTH * CUT_FEATURE_COUNT];
float g_VerboseCutFeatures[MAX_PEPTIDE_LENGTH * CUT_FEATURE_COUNT];
int g_PRMFeatures[PRM_FEATURE_COUNT];
//float g_PRMBScore; // hax
//float g_PRMYScore; // hax
float g_CutScores[MAX_PEPTIDE_LENGTH];
extern PRMBayesianModel* PRMModelCharge2;

int SeizePeaks(MSSpectrum* Spectrum, int TargetMass, int IonType, int AminoIndex, float* pIntensity, float* pSkew, float* pAbsSkew);

FILE* g_ScorpionScoringFile = NULL;

float GetExplainedPeakPercent(MSSpectrum* Spectrum, int PeakCount, int BYOnly)
{
    int PeaksSeen = 0;
    int AnnotatedCount = 0;
    int PeakIndex;
    SpectralPeak* Peak;
    int VerboseFlag = 0;
    //
    for (PeakIndex = 0; PeakIndex < Spectrum->PeakCount; PeakIndex++)
    {
        Peak = Spectrum->Peaks + PeakIndex;
        if (PeakCount>0 && Peak->IntensityRank >= PeakCount)
        {
            continue;
        }
        PeaksSeen++;
        switch (Peak->IonType)
        {
        case IonB:
        case IonB2:
        case IonBI:
        case IonY:
        case IonY2:
        case IonYI:
            AnnotatedCount++;
            if (VerboseFlag)
            {
                printf("Peak index %d at %.2f: %d\n", Peak->IntensityRank, Peak->Mass / (float)MASS_SCALE, Peak->IonType);
            }

            break;
        case 0:
            if (VerboseFlag)
            {
                printf("* Peak index %d at %.2f NOT annotated\n", Peak->IntensityRank, Peak->Mass / (float)MASS_SCALE);
            }
            break; // No annotated intensity for you!
        default:
            if (VerboseFlag)
            {
                printf("Peak index %d at %.2f: %d\n", Peak->IntensityRank, Peak->Mass / (float)MASS_SCALE, Peak->IonType);
            }

            if (!BYOnly)
            {
                AnnotatedCount++;
            }
        }
    }
    if (!PeaksSeen)
    {
        return 0;
    }
    if (PeakCount > 0)
    {
        return AnnotatedCount / (float)PeakCount;
    }
    return AnnotatedCount / (float)PeaksSeen;
}


float GetExplainedIntensityPercent(MSSpectrum* Spectrum, int PeakCount, int BYOnly)
{
    float PeakIntensity = 0;
    float AnnotatedIntensity = 0;
    int PeakIndex;
    SpectralPeak* Peak;
    //
    for (PeakIndex = 0; PeakIndex < Spectrum->PeakCount; PeakIndex++)
    {
        Peak = Spectrum->Peaks + PeakIndex;
        if (PeakCount>0 && Peak->IntensityRank >= PeakCount)
        {
            continue;
        }
        
        PeakIntensity += Peak->Intensity;
        //printf("%.2f\t%.2f\t%d\n", Peak->Mass / (float)MASS_SCALE, Peak->Intensity, Peak->IonType);
        switch (Peak->IonType)
        {
        case IonB:
        case IonB2:
        case IonBI:
        case IonY:
        case IonY2:
        case IonYI:
            AnnotatedIntensity += Peak->Intensity;
            break;
        case 0:
            break; // No annotated intensity for you!
        default:
            if (!BYOnly)
            {
                AnnotatedIntensity += Peak->Intensity;
            }
        }
    }
    if (PeakIntensity == 0)
    {
        return 0;
    }
    return AnnotatedIntensity / PeakIntensity;
}

int GetFlankBFeature(char Left, char Right)
{
    // H on right: Strong suppression
    if (Right == 'H')
    {
        return 0;
    }
    // G or P on left: Strong suppression
    if (Left == 'G' || Left == 'P')
    {
        return 1;
    }
    // K or R on left: Augmentation
    if (Left == 'K' || Left == 'R')
    {
        return 2;
    }
    // P on right: Augmentation
    if (Right == 'P')
    {
        return 3;
    }

    return 4;
}

int GetFlankYFeature(char Left, char Right)
{
    // K or R on right: Strong suppression
    if (Right == 'K' || Right == 'R')
    {
        return 0;
    }
    // G or P on left: Strong suppression
    if (Left == 'G' || Left == 'P')
    {
        return 1;
    }
    // K or R on left: Augmentation
    if (Left == 'K' || Left == 'R')
    {
        return 2;
    }
    // P on right: Augmentation
    if (Right == 'P')
    {
        return 3;
    }

    return 4;
}

// SECTOR_COUNT
#define GET_SECTOR(Mass) \
if (Mass > GlobalOptions->DynamicRangeMax || Mass < GlobalOptions->DynamicRangeMin) \
{\
    Sector = -1;\
}\
else if (Mass > SectorCutoffA) \
{\
    Sector = 1;\
}\
else\
{\
    Sector = 0;\
}

// SECTOR_COUNT
#define GET_CUT_SECTOR(Mass) \
if (Mass > GlobalOptions->DynamicRangeMax || Mass < GlobalOptions->DynamicRangeMin) \
{\
    Sector = -1;\
}\
else if (Mass > SectorCutoffA) \
{\
    Sector = 1;\
}\
else\
{\
    Sector = 0;\
}

// Helper macro for GetPRMFeatures.  NOT applicable in peptide context, only in tagging context
#define GET_BIN_INTENSITY(Mass) \
    Bin = (Mass + 50) / 100;\
    if (Bin >= 0 && Bin < Spectrum->IntensityBinCount) \
    { \
    IntensityLevel = Spectrum->BinnedIntensityLevels[Bin]; \
    } \
    else \
    { \
        IntensityLevel = 0; \
    } 

// Given a putative prefix residue mass for an unannotated spectrum, compute features for scoring its quality.
// Used in building a score array for blind search, and in quick-scoring tagged search.
// This code has SIGNIFICANT OVERLAP with the code in GetCutFeatures()
// Here, PRM is a mass (in thousandths-of-a-dalton).
float GetPRMFeatures(MSSpectrum* Spectrum, SpectrumTweak* Tweak, BayesianModel* Model, int PRM, int VerboseFlag)
{
    int ParentMass;
    int MassB;
    int MassY;
    int Mass;
    int IntensityLevel;
    int Sector;
    int SectorCutoffA;
    float Score = 0;
    int Bin;
    //Spectrum->Charge = Tweak->Charge;
    //Spectrum->ParentMass = Tweak->ParentMass;
    ParentMass = Tweak->ParentMass;
    SectorCutoffA = rint(ParentMass * 0.5);
    // SECTOR_COUNT
    if (PRM > SectorCutoffA)
    {
        g_PRMFeatures[SISector] = 1;
    }
    else
    {
        g_PRMFeatures[SISector] = 0;
    }
    MassB = PRM + DALTON;
    MassY = ParentMass - PRM;
    // Compute the vector of features.  Compute parent features BEFORE computing children.
    g_PRMFeatures[SICharge] = Spectrum->Charge;

    // Alterations to PRM scoring in context of a phosphopeptide search:
    // - Don't try to use a phosphate-loss peak as a b or y peak
    // - Give a bonus for phosphate loss peaks, maybe

    // Find the intensity level for the y peak, and store it in the feature-vector:
    GET_BIN_INTENSITY(MassY);
    //IntensityLevel = SeizePeaks(Spectrum, MassY, 0);
    g_PRMFeatures[IonY] = IntensityLevel;

    // If the y peak is outside dynamic range, then don't adjust the score.
    // If it's in range: Add the y node's log-probability, and subtract the null model's log-probability.
    GET_SECTOR(MassY);
    if (Model && Sector >= 0)
    {
        Score += ComputeBNProbability(Model->Nodes + IonY, IonY, g_PRMFeatures);
        Score -= Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel];
        if (VerboseFlag)
        {
            printf("Y\t%.1f\t\t%d\t%.2f\t%.2f\n", MassY / (float)MASS_SCALE, IntensityLevel, 
                ComputeBNProbability(Model->Nodes + IonY, IonY, g_PRMFeatures),
                Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel]);
        }
    }

    // FOR PHOSOPHOPEPTIDE TAGGING:
    if (GlobalOptions->PhosphorylationFlag)
    {
        GET_BIN_INTENSITY(MassY - PHOSPHATE_WATER_MASS);
        if (IntensityLevel)
        {
            Score += 0.5;
        }
    }

    // b peak:
    //IntensityLevel = SeizePeaks(Spectrum, MassB, 0);
    GET_BIN_INTENSITY(MassB);
    g_PRMFeatures[IonB] = IntensityLevel;

    GET_SECTOR(MassB);
    if (Model && Sector >= 0)
    {
        Score += ComputeBNProbability(Model->Nodes + IonB, IonB, g_PRMFeatures);
        Score -= Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel];
        if (VerboseFlag)
        {
            printf("B\t%.1f\t\t%d\t%.2f\t%.2f\n", MassB / (float)MASS_SCALE, IntensityLevel, 
                ComputeBNProbability(Model->Nodes + IonB, IonB, g_PRMFeatures),
                Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel]);
        }
    }

    // FOR PHOSOPHOPEPTIDE TAGGING:
    if (GlobalOptions->PhosphorylationFlag)
    {
        GET_BIN_INTENSITY(MassB - PHOSPHATE_WATER_MASS);
        if (IntensityLevel)
        {
            Score += 0.5;
        }
    }

    // y isotopic peak:
    Mass = MassY + DALTON;
    //IntensityLevel = SeizePeaks(Spectrum, Mass, 0);
    GET_BIN_INTENSITY(Mass);
    g_PRMFeatures[IonYI] = IntensityLevel;
    GET_SECTOR(Mass);
    if (Model && Sector >= 0)
    {
        Score += ComputeBNProbability(Model->Nodes + IonYI, IonYI, g_PRMFeatures);
        Score -= Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel];
        if (VerboseFlag)
        {
            printf("YI\t%.1f\t\t%d\t%.2f\t%.2f\n", Mass / (float)MASS_SCALE, IntensityLevel, 
                ComputeBNProbability(Model->Nodes + IonYI, IonYI, g_PRMFeatures),
                Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel]);
        }
    }

    // b isotopic peak:
    Mass = MassB + DALTON;
    //IntensityLevel = SeizePeaks(Spectrum, Mass, 0);
    GET_BIN_INTENSITY(Mass);
    g_PRMFeatures[IonBI] = IntensityLevel;
    GET_SECTOR(Mass);
    if (Model && Sector >= 0)
    {
        Score += ComputeBNProbability(Model->Nodes + IonBI, IonBI, g_PRMFeatures);
        Score -= Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel];
        if (VerboseFlag)
        {
            printf("BI\t%.1f\t\t%d\t%.2f\t%.2f\n", Mass / (float)MASS_SCALE, IntensityLevel, 
                ComputeBNProbability(Model->Nodes + IonBI, IonBI, g_PRMFeatures),
                Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel]);
        }
    }

    // doubly-charged y:
    Mass = rint((MassY + HYDROGEN_MASS)/2);
    //IntensityLevel = SeizePeaks(Spectrum, Mass, 0);
    GET_BIN_INTENSITY(Mass);
    g_PRMFeatures[IonY2] = IntensityLevel;
    GET_SECTOR(Mass);
    if (Model && Sector >= 0)
    {
        Score += ComputeBNProbability(Model->Nodes + IonY2, IonY2, g_PRMFeatures);
        Score -= Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel];
        if (VerboseFlag)
        {
            printf("Y2\t%.1f\t\t%d\t%.2f\t%.2f\n", Mass / (float)MASS_SCALE, IntensityLevel, 
                ComputeBNProbability(Model->Nodes + IonY2, IonY2, g_PRMFeatures),
                Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel]);
        }
    }

    // doubly-charged b:
    Mass = rint((MassB + HYDROGEN_MASS)/2);
    //IntensityLevel = SeizePeaks(Spectrum, Mass, 0);
    GET_BIN_INTENSITY(Mass);
    g_PRMFeatures[IonB2] = IntensityLevel; 
    if (Model && Sector >= 0)
    {
        Score += ComputeBNProbability(Model->Nodes + IonB2, IonB2, g_PRMFeatures);
        Score -= Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel];
        if (VerboseFlag)
        {
            printf("B2\t%.1f\t\t%d\t%.2f\t%.2f\n", Mass / (float)MASS_SCALE, IntensityLevel, 
                ComputeBNProbability(Model->Nodes + IonB2, IonB2, g_PRMFeatures),
                Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel]);
        }
    }

    // Y-H2O:
    Mass = MassY - WATER_MASS;
    //IntensityLevel = SeizePeaks(Spectrum, Mass, 0);
    GET_BIN_INTENSITY(Mass);
    g_PRMFeatures[IonYH2O] = IntensityLevel;
    GET_SECTOR(Mass);
    if (Model && Sector >= 0)
    {
        Score += ComputeBNProbability(Model->Nodes + IonYH2O, IonYH2O, g_PRMFeatures);
        Score -= Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel];
        if (VerboseFlag)
        {
            printf("Y-h2o\t%.1f\t\t%d\t%.2f\t%.2f\n", Mass / (float)MASS_SCALE, IntensityLevel, 
                ComputeBNProbability(Model->Nodes + IonYH2O, IonYH2O, g_PRMFeatures),
                Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel]);
        }
    }

    // a:
    Mass = MassB - 27000;
    //IntensityLevel = SeizePeaks(Spectrum, Mass, IonA);
    GET_BIN_INTENSITY(Mass);
    g_PRMFeatures[IonA] = IntensityLevel;
    GET_SECTOR(Mass);
    if (Model && Sector >= 0)
    {
        Score += ComputeBNProbability(Model->Nodes + IonA, IonA, g_PRMFeatures);
        Score -= Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel];
        if (VerboseFlag)
        {
            printf("a\t%.1f\t\t%d\t%.2f\t%.2f\n", Mass / (float)MASS_SCALE, IntensityLevel, 
                ComputeBNProbability(Model->Nodes + IonA, IonA, g_PRMFeatures),
                Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel]);
        }
    }

    // b-H2O:
    Mass = MassB - WATER_MASS;
    //IntensityLevel = SeizePeaks(Spectrum, Mass, IonBH2O);
    GET_BIN_INTENSITY(Mass);
    g_PRMFeatures[IonBH2O] = IntensityLevel;
    GET_SECTOR(Mass);
    if (Model && Sector >= 0)
    {
        Score += ComputeBNProbability(Model->Nodes + IonBH2O, IonBH2O, g_PRMFeatures);
        Score -= Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel];
        if (VerboseFlag)
        {
            printf("b-h2o\t%.1f\t\t%d\t%.2f\t%.2f\n", Mass / (float)MASS_SCALE, IntensityLevel, 
                ComputeBNProbability(Model->Nodes + IonBH2O, IonBH2O, g_PRMFeatures),
                Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel]);
        }
    }

    // y-NH3:
    Mass = MassY - 17000;
    //IntensityLevel = SeizePeaks(Spectrum, Mass, IonYNH3);
    GET_BIN_INTENSITY(Mass);
    g_PRMFeatures[IonYNH3] = IntensityLevel;
    GET_SECTOR(Mass);
    if (Model && Sector >= 0)
    {
        Score += ComputeBNProbability(Model->Nodes + IonYNH3, IonYNH3, g_PRMFeatures);
        Score -= Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel];
        if (VerboseFlag)
        {
            printf("y-nh3\t%.1f\t\t%d\t%.2f\t%.2f\n", Mass / (float)MASS_SCALE, IntensityLevel, 
                ComputeBNProbability(Model->Nodes + IonYNH3, IonYNH3, g_PRMFeatures),
                Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel]);
        }
    }

    // b-NH3:
    Mass = MassB - 17000;
    //IntensityLevel = SeizePeaks(Spectrum, Mass, IonBNH3);
    GET_BIN_INTENSITY(Mass);
    g_PRMFeatures[IonBNH3] = IntensityLevel;
    GET_SECTOR(Mass);
    if (Model && Sector >= 0)
    {
        Score += ComputeBNProbability(Model->Nodes + IonBNH3, IonBNH3, g_PRMFeatures);
        Score -= Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel];
        if (VerboseFlag)
        {
            printf("b-nh3\t%.1f\t\t%d\t%.2f\t%.2f\n", Mass / (float)MASS_SCALE, IntensityLevel, 
                ComputeBNProbability(Model->Nodes + IonBNH3, IonBNH3, g_PRMFeatures),
                Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel]);
        }
    }

    if (GlobalOptions->PhosphorylationFlag)
    {
        Mass = (Spectrum->ParentMass + HYDROGEN_MASS * (Tweak->Charge - 1) - PHOSPHATE_WATER_MASS) / 2;
        if (abs(PRM - Mass) < 1000)
        {
            Score = min(Score, (float)2.0);
        }
    }

    return Score;
}

float g_BIntensity[MAX_PEPTIDE_LENGTH];
float g_YIntensity[MAX_PEPTIDE_LENGTH];
float g_BSkew[MAX_PEPTIDE_LENGTH];
float g_YSkew[MAX_PEPTIDE_LENGTH];
float g_BAbsSkew[MAX_PEPTIDE_LENGTH];
float g_YAbsSkew[MAX_PEPTIDE_LENGTH];

#define FRAGMENTATION_NORMAL 0
#define FRAGMENTATION_PHOSPHO 1

// Given a spectrum and a peptide, generate the feature-vector for each cut-point along the backbone.
// This is separate from GetPRMFeatures, which is done in the context of a spectrum WITHOUT a peptide.
// (The difference: Here we know the flanking peptide, and we have 5 possible sectors rather than 3)
//
// Either FeatureFile is non-null (in which case we should write our feature-vectors to the file), or 
// ScoringNetwork is non-null (in which case we should save the array of cut-point probabilities)
void GetCutFeatures(MSSpectrum* Spectrum, SpectrumTweak* Tweak, Peptide* Match, 
    BayesianModel* Model)
{
    int Mass;
    int PRM;
    int AminoIndex;
    int Length;
    int CutMasses[MAX_PEPTIDE_LENGTH];
    //int BaseFlags[MAX_PEPTIDE_LENGTH];
    //int AcidFlags[MAX_PEPTIDE_LENGTH];
    int ModIndex;
    int MassB;
    int MassY;
    int SectorCutoffA;
    //int SectorCutoffB;
    //int SectorCutoffC;
    int ParentMass;
    int FeatureValue;
    int Sector;
    int IntensityLevel;
    int VerboseFlag = 0;
    int AminoIndexY;
    int AminoIndexB;
    int CutFeaturesBaseIndex;
    //
    Spectrum->Charge = Tweak->Charge;
    Spectrum->ParentMass = Tweak->ParentMass;

    // Check whether we're using special fragmentation models.
    // Use phosphopeptide fragmentation rules if this is Sphos or Tphos (but not for Yphos)
    for (ModIndex = 0; ModIndex < MAX_PT_MODS; ModIndex++)
    {
        if (!Match->ModType[ModIndex])
        {
            break;
        }
        if (Match->ModType[ModIndex]->Flags & DELTA_FLAG_PHOSPHORYLATION && Match->Bases[Match->AminoIndex[ModIndex]]!='Y')
        {
            Match->SpecialFragmentation = FRAGMENTATION_PHOSPHO;
            Match->SpecialModPosition = Match->AminoIndex[ModIndex];
        }
    }
    Length = strlen(Match->Bases);
    Mass = 0;
    //memset(BaseFlags, 0, sizeof(int)*MAX_PEPTIDE_LENGTH);
    //memset(AcidFlags, 0, sizeof(int)*MAX_PEPTIDE_LENGTH);
    // Get the array of masses:
    for (AminoIndex = 0; AminoIndex < Length; AminoIndex++)
    {
        Mass += PeptideMass[Match->Bases[AminoIndex]];
        
        for (ModIndex = 0; ModIndex < MAX_PT_MODS; ModIndex++)
        {
            if (Match->AminoIndex[ModIndex] != AminoIndex)
            {
                continue;
            }
            if (!Match->ModType[ModIndex])
            {
                break;
            }
            Mass += Match->ModType[ModIndex]->RealDelta;
        }
        CutMasses[AminoIndex] = Mass;
    }
    ParentMass = Mass + PARENT_MASS_BOOST;
    // Set dynamic-range feature:
    for (AminoIndex = 0; AminoIndex < Length-1; AminoIndex++)
    {
        Mass = CutMasses[AminoIndex];
        MassB = Mass + DALTON;
        MassY = ParentMass - Mass;
        if (MassB < GlobalOptions->DynamicRangeMin || MassB > GlobalOptions->DynamicRangeMax)
        {
            // B is out.
            if (MassY < GlobalOptions->DynamicRangeMin || MassY > GlobalOptions->DynamicRangeMax)
            {
                // Both out:
                g_CutFeatures[AminoIndex * CUT_FEATURE_COUNT] = 0;
            }
            else
            {
                // Y, no B:
                g_CutFeatures[AminoIndex * CUT_FEATURE_COUNT] = 2;
            }

        }
        else
        {
            if (MassY < GlobalOptions->DynamicRangeMin || MassY > GlobalOptions->DynamicRangeMax)
            {
                // B, no Y:
                g_CutFeatures[AminoIndex*CUT_FEATURE_COUNT] = 1;
            }
            else
            {
                // Both lie inside dynamic range
                g_CutFeatures[AminoIndex*CUT_FEATURE_COUNT] = 3;
            }
        }
    }

    SectorCutoffA = rint(ParentMass * 0.5);
    //SectorCutoffB = (int)(ParentMass * 0.667 + 0.5);
    //SectorCutoffC = (int)(ParentMass * 0.667 + 0.5); // SECTOR_COUNT
    memset(g_VerboseCutFeatures, 0, sizeof(int) * MAX_PEPTIDE_LENGTH * CUT_FEATURE_COUNT);
    memset(g_BSkew, 0, sizeof(float) * MAX_PEPTIDE_LENGTH);
    memset(g_YSkew, 0, sizeof(float) * MAX_PEPTIDE_LENGTH);
    memset(g_BAbsSkew, 0, sizeof(float) * MAX_PEPTIDE_LENGTH);
    memset(g_YAbsSkew, 0, sizeof(float) * MAX_PEPTIDE_LENGTH);

    // Annotate the M-P peak, if it's a phosphopeptide:
    if (Match->SpecialFragmentation == FRAGMENTATION_PHOSPHO)
    {
        Mass = (ParentMass + (Tweak->Charge - 1) * HYDROGEN_MASS - PHOSPHATE_WATER_MASS) / Tweak->Charge;
        SeizePeaks(Spectrum, Mass, IonParentLoss, -1, NULL, NULL, NULL);
    }

    // Capture ions:
    for (AminoIndex = 0; AminoIndex < Length; AminoIndex++)
    {
        // We number the CUTS from 0 up through length-1.  We set AminoIndexY and AminoIndexB
        // to the length (in amino acids) of the b and y fragments, for easy reading by humans.
        AminoIndexY = Length - AminoIndex - 1;
        AminoIndexB = AminoIndex + 1;
        CutFeaturesBaseIndex = AminoIndex * CUT_FEATURE_COUNT;
        PRM = CutMasses[AminoIndex];
        g_CutScores[AminoIndex] = 0;
        MassB = PRM + DALTON;
        MassY = ParentMass - PRM;
        // Compute the vector of features.  Compute parent features BEFORE computing children.
        g_CutFeatures[CutFeaturesBaseIndex + SICharge] = Spectrum->Charge;
        g_CutFeatures[CutFeaturesBaseIndex + SIFlankB] = GetFlankBFeature(Match->Bases[AminoIndex], Match->Bases[AminoIndex + 1]);
        g_CutFeatures[CutFeaturesBaseIndex + SIFlankY] = GetFlankYFeature(Match->Bases[AminoIndex], Match->Bases[AminoIndex + 1]);
        // 2+3 value sector (first cut, last cut, and three regions)
        
        if (PRM > ParentMass * 0.5)
        {
            FeatureValue = 1;
        }
        else
        {
            FeatureValue = 0;
        }
        g_CutFeatures[CutFeaturesBaseIndex + SISector] = FeatureValue;

        // Find the intensity level for the y peak, and store it in the feature-vector:
        IntensityLevel = SeizePeaks(Spectrum, MassY, IonY, AminoIndexY, g_YIntensity + AminoIndex, g_YSkew + AminoIndex, g_YAbsSkew + AminoIndex);
        if (Match->SpecialFragmentation && AminoIndex < Match->SpecialModPosition)
        {
            IntensityLevel = max(IntensityLevel, SeizePeaks(Spectrum, MassY - PHOSPHATE_WATER_MASS, IonY, AminoIndexY, g_YIntensity + AminoIndex, g_YSkew + AminoIndex, g_YAbsSkew + AminoIndex));
        }
        g_CutFeatures[CutFeaturesBaseIndex + IonY] = IntensityLevel;
        // If the y peak is outside dynamic range, then don't adjust the score.
        // If it's in range: Add the y node's log-probability, and subtract the null model's log-probability.
        GET_SECTOR(MassY);
        if (Model && Sector >= 0)
        {
            g_CutScores[AminoIndex] += ComputeBNProbability(Model->Nodes + IonY, IonY, g_CutFeatures + CutFeaturesBaseIndex);
            g_CutScores[AminoIndex] -= Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel];
            g_VerboseCutFeatures[AminoIndex*CUT_FEATURE_COUNT + IonY] = ComputeBNProbability(Model->Nodes + IonY, IonY, g_CutFeatures + CutFeaturesBaseIndex) - Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel];
            if (VerboseFlag)
            {
                printf("%s\t%s\t%d\t%.1f\tY\t%d\t%.2f\t%.2f\t\n", Spectrum->Node->InputFile->FileName,
                    Match->Bases, AminoIndex, MassY / (float)MASS_SCALE, IntensityLevel, 
                    ComputeBNProbability(Model->Nodes + IonY, IonY, g_CutFeatures + CutFeaturesBaseIndex),
                    Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel]);
            }
        }
        // b peak:
        IntensityLevel = SeizePeaks(Spectrum, MassB, IonB, AminoIndexB, g_BIntensity + AminoIndex, g_BSkew + AminoIndex, g_BAbsSkew + AminoIndex);
        if (Match->SpecialFragmentation && AminoIndex >= Match->SpecialModPosition)
        {
            IntensityLevel = max(IntensityLevel, SeizePeaks(Spectrum, MassB - PHOSPHATE_WATER_MASS, IonB, AminoIndexB, g_BIntensity + AminoIndex, g_BSkew + AminoIndex, g_BAbsSkew + AminoIndex));
        }

        g_CutFeatures[CutFeaturesBaseIndex + IonB] = IntensityLevel;
        GET_SECTOR(MassB);
        if (Model && Sector >= 0)
        {
            g_CutScores[AminoIndex] += ComputeBNProbability(Model->Nodes + IonB, IonB, g_CutFeatures + CutFeaturesBaseIndex);
            g_CutScores[AminoIndex] -= Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel];
            g_VerboseCutFeatures[AminoIndex*CUT_FEATURE_COUNT + IonB] = ComputeBNProbability(Model->Nodes + IonB, IonB, g_CutFeatures + CutFeaturesBaseIndex) - Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel];
            if (VerboseFlag)
            {
                printf("%s\t%s\t%d\t%.1f\tB\t%d\t%.2f\t%.2f\t\n", Spectrum->Node->InputFile->FileName,
                    Match->Bases, AminoIndex, MassB / (float)MASS_SCALE, IntensityLevel, 
                    ComputeBNProbability(Model->Nodes + IonB, IonB, g_CutFeatures + CutFeaturesBaseIndex),
                    Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel]);
            }
        }
        // y isotopic peak:
        Mass = MassY + DALTON;
        IntensityLevel = SeizePeaks(Spectrum, Mass, IonYI, AminoIndexY, 0, 0, 0);
        g_CutFeatures[CutFeaturesBaseIndex + IonYI] = IntensityLevel;
        GET_SECTOR(Mass);
        if (Model && Sector >= 0)
        {
            g_CutScores[AminoIndex] += ComputeBNProbability(Model->Nodes + IonYI, IonYI, g_CutFeatures + CutFeaturesBaseIndex);
            g_CutScores[AminoIndex] -= Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel];
            g_VerboseCutFeatures[AminoIndex*CUT_FEATURE_COUNT + IonYI] = ComputeBNProbability(Model->Nodes + IonYI, IonYI, g_CutFeatures + CutFeaturesBaseIndex) - Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel];
            if (VerboseFlag)
            {
                printf("%s\t%s\t%d\t%.1f\tYI\t%d\t%.2f\t%.2f\t\n", Spectrum->Node->InputFile->FileName,
                    Match->Bases, AminoIndex, Mass / (float)MASS_SCALE, IntensityLevel, 
                    ComputeBNProbability(Model->Nodes + IonYI, IonYI, g_CutFeatures + CutFeaturesBaseIndex),
                    Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel]);
            }
        }

        // b isotopic peak:
        Mass = MassB + DALTON;
        IntensityLevel = SeizePeaks(Spectrum, Mass, IonBI, AminoIndexB, 0, 0, 0);
        g_CutFeatures[CutFeaturesBaseIndex + IonBI] = IntensityLevel;
        GET_SECTOR(Mass);
        if (Model && Sector >= 0)
        {
            g_CutScores[AminoIndex] += ComputeBNProbability(Model->Nodes + IonBI, IonBI, g_CutFeatures + CutFeaturesBaseIndex);
            g_CutScores[AminoIndex] -= Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel];
            g_VerboseCutFeatures[AminoIndex*CUT_FEATURE_COUNT + IonBI] = ComputeBNProbability(Model->Nodes + IonBI, IonBI, g_CutFeatures + CutFeaturesBaseIndex) - Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel];
            if (VerboseFlag)
            {
                printf("%s\t%s\t%d\t%.1f\tBI\t%d\t%.2f\t%.2f\t\n", Spectrum->Node->InputFile->FileName,
                    Match->Bases, AminoIndex, Mass / (float)MASS_SCALE, IntensityLevel, 
                    ComputeBNProbability(Model->Nodes + IonBI, IonBI, g_CutFeatures + CutFeaturesBaseIndex),
                    Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel]);
            }
        }

        // doubly-charged y:
        Mass = rint((MassY + HYDROGEN_MASS)/2);
        IntensityLevel = SeizePeaks(Spectrum, Mass, IonY2, AminoIndexY, 0, 0, 0);
        g_CutFeatures[CutFeaturesBaseIndex + IonY2] = IntensityLevel;
        GET_SECTOR(Mass);
        if (Model && Sector >= 0)
        {
            g_CutScores[AminoIndex] += ComputeBNProbability(Model->Nodes + IonY2, IonY2, g_CutFeatures + CutFeaturesBaseIndex);
            g_CutScores[AminoIndex] -= Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel];
            g_VerboseCutFeatures[CutFeaturesBaseIndex + IonY2] = ComputeBNProbability(Model->Nodes + IonY2, IonY2, g_CutFeatures + CutFeaturesBaseIndex) - Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel];
            if (VerboseFlag)
            {
                printf("%s\t%s\t%d\t%.1f\tY2\t%d\t%.2f\t%.2f\t\n", Spectrum->Node->InputFile->FileName,
                    Match->Bases, AminoIndex, Mass / (float)MASS_SCALE, IntensityLevel, 
                    ComputeBNProbability(Model->Nodes + IonY2, IonY2, g_CutFeatures + CutFeaturesBaseIndex),
                    Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel]);
            }
        }

        Mass = rint((MassB + HYDROGEN_MASS)/2);
        IntensityLevel = SeizePeaks(Spectrum, Mass, IonB2, AminoIndexB, 0, 0, 0);
        g_CutFeatures[CutFeaturesBaseIndex + IonB2] = IntensityLevel;
        GET_SECTOR(Mass);
        if (Model && Sector >= 0)
        {
            g_CutScores[AminoIndex] += ComputeBNProbability(Model->Nodes + IonB2, IonB2, g_CutFeatures + CutFeaturesBaseIndex);
            g_CutScores[AminoIndex] -= Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel];
            g_VerboseCutFeatures[AminoIndex*CUT_FEATURE_COUNT + IonB2] = ComputeBNProbability(Model->Nodes + IonB2, IonB2, g_CutFeatures) - Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel];
        }

        // Y-H2O:
        Mass = MassY - WATER_MASS;
        IntensityLevel = SeizePeaks(Spectrum, Mass, IonYH2O, AminoIndexY, 0, 0, 0);
        g_CutFeatures[CutFeaturesBaseIndex + IonYH2O] = IntensityLevel;
        GET_SECTOR(Mass);
        if (Model && Sector >= 0)
        {
            g_CutScores[AminoIndex] += ComputeBNProbability(Model->Nodes + IonYH2O, IonYH2O, g_CutFeatures + CutFeaturesBaseIndex);
            g_CutScores[AminoIndex] -= Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel];
            g_VerboseCutFeatures[AminoIndex*CUT_FEATURE_COUNT + IonYH2O] = ComputeBNProbability(Model->Nodes + IonYH2O, IonYH2O, g_CutFeatures) - Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel];
        }

        // a:
        Mass = MassB - 27000;
        IntensityLevel = SeizePeaks(Spectrum, Mass, IonA, AminoIndexB, 0, 0, 0);
        g_CutFeatures[CutFeaturesBaseIndex + IonA] = IntensityLevel;
        GET_SECTOR(Mass);
        if (Model && Sector >= 0)
        {
            g_CutScores[AminoIndex] += ComputeBNProbability(Model->Nodes + IonA, IonA, g_CutFeatures + CutFeaturesBaseIndex);
            g_CutScores[AminoIndex] -= Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel];
            g_VerboseCutFeatures[AminoIndex*CUT_FEATURE_COUNT + IonA] = ComputeBNProbability(Model->Nodes + IonA, IonA, g_CutFeatures) - Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel];
        }

        // b-H2O:
        Mass = MassB - WATER_MASS;
        IntensityLevel = SeizePeaks(Spectrum, Mass, IonBH2O, AminoIndexB, 0, 0, 0);
        g_CutFeatures[CutFeaturesBaseIndex + IonBH2O] = IntensityLevel;
        GET_SECTOR(Mass);
        if (Model && Sector >= 0)
        {
            g_CutScores[AminoIndex] += ComputeBNProbability(Model->Nodes + IonBH2O, IonBH2O, g_CutFeatures + CutFeaturesBaseIndex);
            g_CutScores[AminoIndex] -= Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel];
            g_VerboseCutFeatures[AminoIndex*CUT_FEATURE_COUNT + IonBH2O] = ComputeBNProbability(Model->Nodes + IonBH2O, IonBH2O, g_CutFeatures) - Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel];
        }

        // y-NH3:
        Mass = MassY - 17000;
        IntensityLevel = SeizePeaks(Spectrum, Mass, IonYNH3, AminoIndexY, 0, 0, 0);
        g_CutFeatures[CutFeaturesBaseIndex + IonYNH3] = IntensityLevel;
        GET_SECTOR(Mass);
        if (Model && Sector >= 0)
        {
            g_CutScores[AminoIndex] += ComputeBNProbability(Model->Nodes + IonYNH3, IonYNH3, g_CutFeatures + CutFeaturesBaseIndex);
            g_CutScores[AminoIndex] -= Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel];
            g_VerboseCutFeatures[AminoIndex*CUT_FEATURE_COUNT + IonYNH3] = ComputeBNProbability(Model->Nodes + IonYNH3, IonYNH3, g_CutFeatures) - Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel];
        }

        // b-NH3:
        Mass = MassB - 17000;
        IntensityLevel = SeizePeaks(Spectrum, Mass, IonBNH3, AminoIndexB, 0, 0, 0);
        g_CutFeatures[CutFeaturesBaseIndex + IonBNH3] = IntensityLevel;
        GET_SECTOR(Mass);
        if (Model && Sector >= 0)
        {
            g_CutScores[AminoIndex] += ComputeBNProbability(Model->Nodes + IonBNH3, IonBNH3, g_CutFeatures + CutFeaturesBaseIndex);
            g_CutScores[AminoIndex] -= Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel];
            g_VerboseCutFeatures[AminoIndex*CUT_FEATURE_COUNT + IonBNH3] = ComputeBNProbability(Model->Nodes + IonBNH3, IonBNH3, g_CutFeatures) - Tweak->Intensities[Sector * INTENSITY_LEVEL_COUNT + IntensityLevel];
        }

    }
}

// Take all unlabeled peaks in a radius of the target m/z.  Annotate them with this ion type, 
// and return the cumulative intensity level.
int SeizePeaks(MSSpectrum* Spectrum, int TargetMass, int IonType, int AminoIndex, float* pIntensity, float* pSkew, float *pAbsSkew)
{
    int PeakIndex;
    int MaxMass;
    float Intensity = 0;
    int Bin;
    float WeightedSkew = 0;
    float AbsWeightedSkew = 0;
    //
    Bin = (TargetMass + 50) / 100;
    MaxMass = TargetMass + INTENSITY_BIN_RADIUS;

    // If the mass is off the scale, then you get no peaks:
    if (Bin >= Spectrum->IntensityBinCount || Bin < 0)
    {
        return 0;
    }
    PeakIndex = Spectrum->BinPeakIndex[Bin];
    if (PeakIndex >= 0)
    {
        for ( ; PeakIndex < Spectrum->PeakCount; PeakIndex++)
        {
            if (Spectrum->Peaks[PeakIndex].Mass > MaxMass)
            {
                break;
            }
            Intensity += Spectrum->Peaks[PeakIndex].Intensity;
            Spectrum->Peaks[PeakIndex].IonType = IonType;
            Spectrum->Peaks[PeakIndex].AminoIndex = AminoIndex;
            WeightedSkew += Spectrum->Peaks[PeakIndex].Intensity * (Spectrum->Peaks[PeakIndex].Mass - TargetMass);
            AbsWeightedSkew += Spectrum->Peaks[PeakIndex].Intensity * abs(Spectrum->Peaks[PeakIndex].Mass - TargetMass);
        }
    }
    if (pIntensity)
    {
        *pIntensity = Intensity;
        *pSkew = WeightedSkew;
        *pAbsSkew = AbsWeightedSkew;
    }

    if (Intensity < Spectrum->IntensityCutoffLow)
    {
        return 0;
    }
    if (Intensity < Spectrum->IntensityCutoffMedium)
    {
        return 1;
    }
    if (Intensity < Spectrum->IntensityCutoffHigh)
    {
        return 2;
    }
    return 3;
}

FILE* g_TrainFile2;
FILE* g_TrainFile3;

// Callback for trining PRM scorer in peptide context.
void TrainPepPRMCallback(SpectrumNode* Node, int Charge, int ParentMass, Peptide* Annotation)
{
    FILE* OutputFile;
    int AminoIndex;
    int Length;
    int FeatureIndex;
    MSSpectrum* Spectrum;

    Spectrum = Node->Spectrum;
    Length = strlen(Annotation->Bases);
    WindowFilterPeaks(Spectrum, 0, 0);
    IntensityRankPeaks(Spectrum);
    // Use the charge+PM oracle:
    Node->Tweaks[0].Charge = Charge;
    Node->Tweaks[0].ParentMass = ParentMass;
    
    PrepareSpectrumForIonScoring(PRMModelCharge2, Node->Spectrum, 0);
    //SpectrumComputeBinnedIntensities(Node);
    if (Charge < 3)
    {
        OutputFile = g_TrainFile2;
    }
    else
    {
        OutputFile = g_TrainFile3;
    }
    GetCutFeatures(Node->Spectrum, Node->Tweaks, Annotation, NULL);

    for (AminoIndex = 0; AminoIndex < Length; AminoIndex++)
    {
        for (FeatureIndex = 0; FeatureIndex < SIMax; FeatureIndex++)
        {
            fprintf(OutputFile, "%d\t", g_CutFeatures[AminoIndex*CUT_FEATURE_COUNT + FeatureIndex]);
        }
        fprintf(OutputFile, "\n");
    }
    fflush(OutputFile);
}

void TrainPepPRM(char* OracleFile, char* OracleDir)
{
    g_TrainFile2 = fopen("TrainingFiles\\PEPPRM2.txt", "w");
    g_TrainFile3 = fopen("TrainingFiles\\PEPPRM3.txt", "w");
    TrainOnOracleFile(OracleFile, OracleDir, TrainPepPRMCallback);
}

void ScorpionSetPRMScores(MSSpectrum* Spectrum, SpectrumTweak* Tweak)
{
    BayesianModel* Model;
    int PRM;
    float fScore;
    //
    // Ensure models are loaded:
    if (!BNCharge2TaggingBN)
    {
        InitBayesianModels();
    }
    Tweak->PRMScoreMax = Tweak->ParentMass;
    if (Spectrum->Graph)
    {
        Tweak->PRMScoreMax = max(Tweak->PRMScoreMax, Spectrum->Graph->LastNode->Mass);
    }
    Tweak->PRMScoreMax = PRM_ARRAY_SLACK + (Tweak->PRMScoreMax / PRM_BIN_SIZE);
    SafeFree(Tweak->PRMScores);
    Tweak->PRMScores = (int*)calloc(Tweak->PRMScoreMax + 5, sizeof(int)); // extra slack in alloc
    if (Tweak->Charge > 2)
    {
        Model = BNCharge3TaggingBN;
    }
    else
    {
        Model = BNCharge2TaggingBN;
    }
    for (PRM = 0; PRM < Tweak->PRMScoreMax; PRM++)
    {
        fScore = GetPRMFeatures(Spectrum, Tweak, Model, PRM * PRM_BIN_SIZE, 0);
        Tweak->PRMScores[PRM] = rint(fScore * 1000);
    }
}

void FinishPRMTestRecord(char* RememberFileName, int* Scores, int MatchCount, int* RankHistogram, char* CandidateAnnotations)
{
    int TrueScore;
    int ScoreIndex;
    int BestScore = -9999;
    int BestScoreIndex = 0;
    int HistogramPoint = 0;
    //
    // Find the best score:
    for (ScoreIndex = 0; ScoreIndex < MatchCount; ScoreIndex++)
    {
        if (Scores[ScoreIndex] > BestScore)
        {
            BestScore = Scores[ScoreIndex];
            BestScoreIndex = ScoreIndex;
        }
    }
    TrueScore = Scores[0];

    qsort(Scores, MatchCount, sizeof(int), (QSortCompare)CompareInts);
    for (ScoreIndex = 0; ScoreIndex < MatchCount; ScoreIndex++)
    {
        if (Scores[ScoreIndex] <= TrueScore)
        {
            // Found it!
            RankHistogram[ScoreIndex] += 1;
            HistogramPoint = ScoreIndex;
            break;
        }
    }
    // Verbose output:
    printf("%s\t%s\t%s\t%d\t%d\t%d\n", RememberFileName, CandidateAnnotations, CandidateAnnotations + 128*BestScoreIndex,
        BestScore, TrueScore, HistogramPoint);

}


// TestPepPRMCallback:
// Print the minimum, maximum, and average score of cut-point scores for this annotation.
// The primary goal here is to evaluate whether minor changes to the PepPRM scoring
// model, such as changing the intensity cutoffs or adding new nodes and edges,
// improve the model's effectiveness.
// This function also serves as a 'sanity check' that true matches are getting
// reasonably good PepPRM scores.
void TestPepPRMCallback(SpectrumNode* Node, int Charge, int ParentMass, Peptide* Annotation)
{
    float MinPRMScore = 9999;
    float MaxPRMScore = -9999;
    int AminoIndex;
    int PRMCount;
    int PRM;
    int Len;
    float Score;
    float TotalScore;
    BayesianModel* Model;
    ////////////////////////////////////////////////////////////////////////
    // Main
    Node->Tweaks[0].ParentMass = Annotation->ParentMass;
    Node->Tweaks[0].Charge = Charge;
    Node->Spectrum->ParentMass = Annotation->ParentMass;
    Node->Spectrum->Charge = Charge;
    WindowFilterPeaks(Node->Spectrum, 0, 0);
    IntensityRankPeaks(Node->Spectrum);
    PrepareSpectrumForIonScoring(PRMModelCharge2, Node->Spectrum, 0);
    //SpectrumComputeBinnedIntensities(Node);
    //SpectrumComputeNoiseDistributions(Node);
    if (Charge > 2)
    {
        Model = BNCharge3ScoringBN;
    }
    else
    {
        Model = BNCharge2ScoringBN;
    }
    PRMCount = 0;
    PRM = 0;
    TotalScore = 0;
    Len = strlen(Annotation->Bases);
    GetCutFeatures(Node->Spectrum, Node->Tweaks, Annotation, Model);
    for (AminoIndex = 0; AminoIndex < Len; AminoIndex++)
    {
        Score = g_CutScores[AminoIndex];
        TotalScore += Score;
        MinPRMScore = min(MinPRMScore, Score);
        MaxPRMScore = max(MaxPRMScore, Score);
        PRMCount++;
    }
    Score = TotalScore / PRMCount;
    printf("%s\t%s\t%.2f\t%.2f\t%.2f\t\n", Node->InputFile->FileName, Annotation->Bases, Score, MinPRMScore, MaxPRMScore);
}

void TestPRMQuickScoringCallback(SpectrumNode* Node, int Charge, int ParentMass, Peptide* Annotation)
{
    static int* Scores;
    static int MatchCount;
    static char* CurrentFile;
    static int* RankHistogram;
    static int RowsProcessed;
    int Score;
    int PRM;
    int AminoIndex;
    int ModIndex;
    int Len;
    BayesianModel* Model;
    BayesianModel* PepPRMModel;
    int Cumulative;
    int TotalPeptides;
    int RankIndex;
    int PRMCount;
    static char CandidateAnnotations[512*128];
    static char RememberFileName[1024];

    // If Node is null, then we've been called in initialize / cleanup mode:
    if (!Node)
    {
        if (!Charge)
        {
            CurrentFile = (char*)calloc(256, sizeof(char));
            Scores = (int*)calloc(512, sizeof(int));
            RankHistogram = (int*)calloc(512, sizeof(int));
            RowsProcessed = 0;
            return;
        }
        // Finish the current peptide, if any:
        if (*CurrentFile)
        {
            FinishPRMTestRecord(RememberFileName, Scores, MatchCount, RankHistogram, CandidateAnnotations);
        }
        // Now report:
        printf("Histogram of PRM quick score pack positions:\n");
        TotalPeptides = 0;
        Cumulative = 0;
        for (RankIndex = 0; RankIndex< 512; RankIndex++)
        {
            TotalPeptides += RankHistogram[RankIndex];
        }
        for (RankIndex = 0; RankIndex< 512; RankIndex++)
        {
            Cumulative += RankHistogram[RankIndex];
            printf("%d\t%d\t%d\t%.2f\t%.2f\t\n", RankIndex, RankHistogram[RankIndex], Cumulative,
                RankHistogram[RankIndex] / (float)TotalPeptides, Cumulative / (float)TotalPeptides);
        }
        // Lastly, free memory:
        SafeFree(Scores);
        SafeFree(RankHistogram);
        SafeFree(CurrentFile);
        return;
    }
    ////////////////////////////////////////////////////////////////////////
    // Main
    Node->Tweaks[0].ParentMass = Annotation->ParentMass;
    Node->Tweaks[0].Charge = Charge;
    Node->Spectrum->ParentMass = Annotation->ParentMass;
    Node->Spectrum->Charge = Charge;
    WindowFilterPeaks(Node->Spectrum, 0, 0);
    IntensityRankPeaks(Node->Spectrum);
    PrepareSpectrumForIonScoring(PRMModelCharge2, Node->Spectrum, 0);
    //SpectrumComputeBinnedIntensities(Node);
    //SpectrumComputeNoiseDistributions(Node);
    // OLD PRM SCORING: 
    //SetPRMScores(Node->Spectrum); 
    RowsProcessed++;
    if (strcmp(CurrentFile, Node->InputFile->FileName))
    {
        if (MatchCount)
        {
            FinishPRMTestRecord(RememberFileName, Scores, MatchCount, RankHistogram, CandidateAnnotations);
        }
        MatchCount = 0;
        strcpy(CurrentFile, Node->InputFile->FileName);
    }
    sprintf(RememberFileName, Node->InputFile->FileName);
    // Compute score for these PRM values:
    if (Charge > 2)
    {
        PepPRMModel = BNCharge3ScoringBN;
        Model = BNCharge3TaggingBN;
    }
    else
    {
        PepPRMModel = BNCharge2ScoringBN;
        Model = BNCharge2TaggingBN;
    }
    Len = strlen(Annotation->Bases);
    PRM = 0;
    Score = 0;
    PRMCount = 0;
    
    // Verify that using flanking aminos and such really does improve things:
    GetCutFeatures(Node->Spectrum, Node->Tweaks, Annotation, PepPRMModel);

    for (AminoIndex = 0; AminoIndex < Len; AminoIndex++)
    {
        PRM += PeptideMass[Annotation->Bases[AminoIndex]];
        for (ModIndex = 0; ModIndex < MAX_PT_MODS; ModIndex++)
        {
            if (Annotation->AminoIndex[ModIndex] == AminoIndex && Annotation->ModType[ModIndex])
            {
                PRM += Annotation->ModType[ModIndex]->RealDelta;
            }
        }
        PRMCount++;


        Score += rint(1000 * g_CutScores[AminoIndex]); 
    }
    Score = Score / PRMCount;
    // Cheese to prevent running off the edge of the array:
    if (MatchCount < 512)
    {
        Scores[MatchCount] = Score;
        strcpy(CandidateAnnotations + 128 * MatchCount, Annotation->Bases);
        MatchCount++;
    }
}

void TestPRMQuickScoring(char* OracleFile, char* OracleDir)
{
    InitBayesianModels();
    TestPRMQuickScoringCallback(NULL, 0, 0, NULL); // initialization
    TrainOnOracleFile(OracleFile, OracleDir, TestPRMQuickScoringCallback);
    TestPRMQuickScoringCallback(NULL, 1, 0, NULL); // completion
}

void TestPepPRM(char* OracleFile, char* OracleDir)
{
    InitBayesianModels();
    TrainOnOracleFile(OracleFile, OracleDir, TestPepPRMCallback);
}

// TestLDACallback:
void TestLDACallback(SpectrumNode* Node, int Charge, int ParentMass, Peptide* Annotation)
{
    BayesianModel* Model;
    // 
    Node->Tweaks[0].ParentMass = Annotation->ParentMass;
    Node->Tweaks[0].Charge = Charge;
    Annotation->Tweak = Node->Tweaks;
    Node->Spectrum->ParentMass = Annotation->ParentMass;
    Node->Spectrum->Charge = Charge;
    WindowFilterPeaks(Node->Spectrum, 0, 0);
    IntensityRankPeaks(Node->Spectrum);
    PrepareSpectrumForIonScoring(PRMModelCharge2, Node->Spectrum, 0);
    //SpectrumComputeBinnedIntensities(Node);
    //SpectrumComputeNoiseDistributions(Node);
    GlobalOptions->DigestType = DIGEST_TYPE_TRYPSIN;
    if (Charge > 2)
    {
        Model = BNCharge3ScoringBN;
    }
    else
    {
        Model = BNCharge2ScoringBN;
    }
    ScorpionSetPRMScores(Node->Spectrum, Node->Tweaks); 

}

void TestLDA(char* OracleFile, char* OracleDir)
{
    InitBayesianModels();
    TrainOnOracleFile(OracleFile, OracleDir, TestLDACallback);
}

// For debug output: Return a description of an ion type code.
char* GetScorpIonName(int IonType)
{
    switch (IonType)
    {
    case IonY:
        return "y";
    case IonB:
        return "b";
    case IonYI:
        return "yi";
    case IonBI:
        return "bi";
    case IonY2:
        return "y2";
    case IonB2:
        return "b2";
    case IonYH2O:
        return "y-h2o";
    case IonA:
        return "a";
    case IonBH2O:
        return "b-h2o";
    case IonYNH3:
        return "y-nh3";
    case IonBNH3:
        return "b-nh3";
    default:
        return "";
    }
}
