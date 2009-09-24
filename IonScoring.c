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
#include <math.h>
#include "Inspect.h"
#include "Utils.h"
#include "Errors.h"
#include "IonScoring.h"
#include "Spectrum.h"
#include "Tagger.h"

// Global variables: Bayesian networks for PRM scoring (for MS-Alignment) and
// for cut scoring (for tagging and match-scoring)
PRMBayesianModel* PRMModelCharge2 = NULL;
PRMBayesianModel* PRMModelCharge3 = NULL;
PRMBayesianModel* TAGModelCharge2 = NULL;
PRMBayesianModel* TAGModelCharge3 = NULL;
PRMBayesianModel* PhosCutModelCharge2 = NULL;
PRMBayesianModel* PhosCutModelCharge3 = NULL;

// Forward declarations:
int IonScoringGetPrefixContainPhos(PRMBayesianNode* Node, Peptide* Match, int AminoIndex);
int IonScoringGetSuffixContainPhos(PRMBayesianNode* Node, Peptide* Match, int BreakIndex);
PRMBayesianModel* GetScoringModel(Peptide* Match, int Charge);
void AnnotateParentPeaks(MSSpectrum* Spectrum, Peptide* Match, PRMBayesianModel* Model);
void ClaimParentPeak(MSSpectrum* Spectrum, Peptide* Match, int Mass, PRMBayesianModel* Model);

// Free a node from a Bayesian network; helper for FreePRMBayesianModel
void FreePRMBayesianNode(PRMBayesianNode* Node)
{
    if (!Node)
    {
        return;
    }
    SafeFree(Node->Parents);
    SafeFree(Node->ParentBlocks);
    SafeFree(Node->CountTable);
    SafeFree(Node->ProbTable);
    SafeFree(Node);
}

// Free a Bayesian network model.
void FreePRMBayesianModel(PRMBayesianModel* Model)
{
    PRMBayesianNode* Node;
    PRMBayesianNodeHolder* Holder;
    PRMBayesianNodeHolder* PrevHolder = NULL;
    PRMBayesianNode* Prev = NULL;
    //
    if (!Model)
    {
        return;
    }
    // Free the linked list of node-holders that require flanking amino acid info:
    for (Holder = Model->FirstFlank; Holder; Holder = Holder->Next)
    {
        SafeFree(PrevHolder);
        PrevHolder = Holder;
    }
    SafeFree(PrevHolder);

    // Free the linked list of all nodes:
    for (Node = Model->Head; Node; Node = Node->Next)
    {
        FreePRMBayesianNode(Prev);
        Prev = Node;
    }
    FreePRMBayesianNode(Prev);
    SafeFree(Model->Nodes);
    SafeFree(Model);
}

// Add a node to a Bayesian network.  Called from PyInspect when building up 
// a network (semi)interactively, not used in production.
void AddPRMBayesianNode(PRMBayesianModel* Model, char* Name, int NodeType, int NodeFlag, float NodeMassOffset, 
    int FragmentType)
{
    PRMBayesianNode* Node;
    //
    // Create the node:
    Node = (PRMBayesianNode*)calloc(1, sizeof(PRMBayesianNode));
    Node->Type = NodeType;
    strncpy(Node->Name, Name, 256);
    Node->MassOffset = (int)(NodeMassOffset * DALTON);
    Node->Flag = NodeFlag;
    Node->Index = Model->NodeCount;
    Node->FragmentType = FragmentType;
    Model->NodeCount++;
    // Insert the node into the list:
    if (Model->Tail)
    {
        Model->Tail->Next = Node;
    }
    else
    {
        Model->Head = Node;
    }
    Model->Tail = Node;
    
    // Insert the node into the array:
    if (Model->Nodes)
    {
        Model->Nodes = (PRMBayesianNode**)realloc(Model->Nodes, Model->NodeCount * sizeof(PRMBayesianNode*));
    }
    else
    {
        Model->Nodes = (PRMBayesianNode**)calloc(sizeof(PRMBayesianNode*), 1);
    }
    Model->Nodes[Model->NodeCount - 1] = Node;

    // Now set the value count:
    switch (Node->Type)
    {
    case evPRMBPrefix:
    case evPRMBPrefix2:
    case evPRMBSuffix:
    case evPRMBSuffix2:
        // The number of values is determined by the intensity scheme:
        switch (Model->IntensityScheme)
        {
        case 0:
        case 1:
        case 4:
            Node->ValueCount = 4;
            break;
        case 2:
        case 3:
            Node->ValueCount = 3;
            break;
        default:
            REPORT_ERROR(0);
            break;
        }
        break;
    case evSector:
        // The number of values is 2, 3, 4, or 5, depending on our sector count:
        switch (Node->Flag)
        {
        case 0:
            Node->ValueCount = 2;
            break;
        case 1:
            Node->ValueCount = 3;
            break;
        case 2:
            Node->ValueCount = 4;
            break;
        case 3:
            Node->ValueCount = 5;
            break;
        case 4:
            Node->ValueCount = 5;
            break;
        default:
            REPORT_ERROR(0);
            break;
        }
        break;
    case evFlank:
        // The number of values depends on the flank scheme flag:
        switch (Node->Flag)
        {
        case 0:
            Node->ValueCount = 4;
            break;
        case 1:
            Node->ValueCount = 4;
            break;
        case 2:
            Node->ValueCount = 3;
            break;
        case 3:
            Node->ValueCount = 3;
            break;
        default:
            REPORT_ERROR(0);
            break;
        }
        break;
    case evPrefixAA:
    case evSuffixAA:
        // PrefixAA and SuffixAA nodes are simple binary nodes.
        Node->ValueCount = 2;
        break;
    case evPrefixContain:
        switch (Node->Flag)
        {
        case 0:
            // Acid residue (flag)
            Node->ValueCount = 2;
            break;
        case 1:
            // Acid residue (0, 1, many)
            Node->ValueCount = 3;
            break;
        case 2:
            // Basic residue (flag)
            Node->ValueCount = 2;
            break;
        case 3:
            // Basic residue (0, 1, many)
            Node->ValueCount = 2;
            break;
        default:
            REPORT_ERROR(0);
            break;
        }
        break;
    case evSuffixContain:
        switch (Node->Flag)
        {
        case 0:
            // Acid residue (flag)
            Node->ValueCount = 2;
            break;
        case 1:
            // Acid residue (0, 1, many)
            Node->ValueCount = 3;
            break;
        case 2:
            // Basic residue (flag)
            Node->ValueCount = 2;
            break;
        case 3:
            // Basic residue (0, 1, many)
            Node->ValueCount = 2;
            break;
        default:
            REPORT_ERROR(0);
            break;
        }
        break;
    case evPrefixContainPhos:
    case evSuffixContainPhos:
        //has phosphate on the fragment (flag)
        Node->ValueCount = 2;
        break;
    default:
        printf("* Error: Unknown Node->Type in AddPRMBayesianNode\n");
        break;
    }
    // Allocate initial count/probability tables, assuming no parents for the node:
    Node->CountTable = (int*)calloc(Node->ValueCount, sizeof(int));
    Node->ProbTable = (float*)calloc(Node->ValueCount, sizeof(float));
    Node->TableSize = Node->ValueCount;
}

// Given a spectrum, compute the intensity-thresholds for level 0 (strongest)
// through level n (absent).
int ComputeSpectrumIntensityThresholds(PRMBayesianModel* Model, MSSpectrum* Spectrum)
{
    int ThresholdCount;
    int CutoffRank;
    int WeakRank;
    int PeakIndex;
    float SortedIntensity[200];
    int WeakPeakCount = 0;
    float TotalIntensity = 0;
    float GrassIntensity;
    float StrongPeakIntensity;
    //

    switch (Model->IntensityScheme)
    {
    case 0:
    case 1:
        // Scheme 1: Top N peaks, high, low, absent
        ThresholdCount = 4;
        Spectrum->IntensityThresholds = (float*)calloc(5, sizeof(float));
        StrongPeakIntensity = -1;
        CutoffRank = (int)(Spectrum->ParentMass / (50 * DALTON));
        WeakRank = max(CutoffRank, Spectrum->PeakCount - 200);
        for (PeakIndex = 0; PeakIndex < Spectrum->PeakCount; PeakIndex++)
        {
            if (Spectrum->Peaks[PeakIndex].IntensityRank >= WeakRank)
            {
                SortedIntensity[WeakPeakCount] = Spectrum->Peaks[PeakIndex].Intensity;
                WeakPeakCount++;
            }
            else
            {
                if (StrongPeakIntensity < 0 || StrongPeakIntensity > Spectrum->Peaks[PeakIndex].Intensity)
                {
                    StrongPeakIntensity = Spectrum->Peaks[PeakIndex].Intensity;
                }
            }
            TotalIntensity += Spectrum->Peaks[PeakIndex].Intensity;
            if (WeakPeakCount == 200)
            {
                break;
            }
        }
        if (!WeakPeakCount)
        {
            GrassIntensity = TotalIntensity / (2 * Spectrum->PeakCount);
        }
        else
        {
            qsort(SortedIntensity, WeakPeakCount, sizeof(float), (QSortCompare)CompareFloats);
            GrassIntensity = SortedIntensity[WeakPeakCount / 2];
        }
        Spectrum->IntensityThresholds[0] = StrongPeakIntensity;
        Spectrum->IntensityThresholds[1] = (float)min(StrongPeakIntensity * 0.5, GrassIntensity * 2);
        //Spectrum->IntensityThresholds[2] = (float)0.5 * GrassIntensity;
        Spectrum->IntensityThresholds[2] = 0;
        Spectrum->IntensityThresholds[3] = -1;
        break;
    case 2:
    case 3:
        // Scheme 1: Top N peaks, present, absent
        ThresholdCount = 3;
        Spectrum->IntensityThresholds = (float*)calloc(5, sizeof(float));
        StrongPeakIntensity = -1;
        CutoffRank = (int)(Spectrum->ParentMass / (50 * DALTON));
        WeakRank = max(CutoffRank, Spectrum->PeakCount - 200);
        for (PeakIndex = 0; PeakIndex < Spectrum->PeakCount; PeakIndex++)
        {
            if (Spectrum->Peaks[PeakIndex].IntensityRank >= WeakRank)
            {
                SortedIntensity[WeakPeakCount] = Spectrum->Peaks[PeakIndex].Intensity;
                WeakPeakCount++;
            }
            else
            {
                if (StrongPeakIntensity < 0 || StrongPeakIntensity > Spectrum->Peaks[PeakIndex].Intensity)
                {
                    StrongPeakIntensity = Spectrum->Peaks[PeakIndex].Intensity;
                }
            }
            TotalIntensity += Spectrum->Peaks[PeakIndex].Intensity;
            if (WeakPeakCount == 200)
            {
                break;
            }
        }
        if (!WeakPeakCount)
        {
            GrassIntensity = TotalIntensity / (2 * Spectrum->PeakCount);
        }
        else
        {
            qsort(SortedIntensity, WeakPeakCount, sizeof(float), (QSortCompare)CompareFloats);
            GrassIntensity = SortedIntensity[WeakPeakCount / 2];
        }
        Spectrum->IntensityThresholds[0] = StrongPeakIntensity;
        Spectrum->IntensityThresholds[1] = 0; //GrassIntensity * 0.5;
        //Spectrum->IntensityThresholds[2] = (float)0.5 * GrassIntensity;
        Spectrum->IntensityThresholds[2] = -1;
        //Spectrum->IntensityThresholds[3] = -1;
        break;
    case 4:
        //Scheme 4: partitioned by ratio to grass
        ThresholdCount = 4;
        Spectrum->IntensityThresholds = (float*)calloc(5, sizeof(float));
        WeakRank = (Spectrum->PeakCount / 3 ); //AverageGrass = median of bottom 1/3 of peaks
        WeakRank = min(200, WeakRank); //at most 200, limited by array size
        for (PeakIndex = 0; PeakIndex < Spectrum->PeakCount; PeakIndex++)
        {
            if (Spectrum->Peaks[PeakIndex].IntensityRank >= WeakRank)
            {
                SortedIntensity[WeakPeakCount] = Spectrum->Peaks[PeakIndex].Intensity;
                WeakPeakCount++;
            }
        }
        if (!WeakPeakCount)
        {
            GrassIntensity = TotalIntensity / (2 * Spectrum->PeakCount);
        }
        else
        {
            qsort(SortedIntensity, WeakPeakCount, sizeof(float), (QSortCompare)CompareFloats);
            GrassIntensity = SortedIntensity[WeakPeakCount / 2];
        }
        Spectrum->IntensityThresholds[0] = GrassIntensity * (float)10.0;
        Spectrum->IntensityThresholds[1] = GrassIntensity * 2;
        Spectrum->IntensityThresholds[2] = GrassIntensity * (float)0.1;
        Spectrum->IntensityThresholds[3] = -1;
        break;
    default:
        REPORT_ERROR(0);
        return 0;
    }
    return ThresholdCount;
}

// Prepare a spectrum for PRM and cut scoring.  Compute intensity cutoffs, compute binned
// intensity, comput binned intensity levels.  
void PrepareSpectrumForIonScoring(PRMBayesianModel* Model, MSSpectrum* Spectrum, int ForceRefresh)
{
    int WeakPeakCount = 0;
    float TotalIntensity = 0;
    int ThresholdCount;
    int PeakIndex;
    int IntensityLevel;
    int BinScalingFactor = 100; // One bin per 0.1Da
    int CountByIntensityLevel[16];
    int Bin;
    int NearBin;
    SpectralPeak* Peak;
    int MaxParentMass;
    float Intensity;
    float Probability;
    float Multiplier;
    int Skew;
    //
    if (Spectrum->IntensityThresholds && !ForceRefresh)
    {
        return; // Already set!
    }
    if (!Spectrum->PeakCount)
    {
        return;
    }
    if (!Model)
    {
        return;
    }

    ///////////////////////////////
    // Free any old info:
    SafeFree(Spectrum->BinnedIntensities);
    Spectrum->BinnedIntensities = NULL;
    SafeFree(Spectrum->BinnedIntensitiesTight);
    Spectrum->BinnedIntensitiesTight = NULL;
    SafeFree(Spectrum->BinnedIntensityLevels);
    Spectrum->BinnedIntensityLevels = NULL;
    SafeFree(Spectrum->BinPeakIndex);
    Spectrum->BinPeakIndex = NULL;
    SafeFree(Spectrum->IonScoringNoiseProbabilities);
    Spectrum->IonScoringNoiseProbabilities = NULL;
    SafeFree(Spectrum->IntensityThresholds);
    Spectrum->IntensityThresholds = NULL;
    ///////////////////////////////
    ThresholdCount = ComputeSpectrumIntensityThresholds(Model, Spectrum);

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // We know our intensity thresholds; now let's compute the binned intensities:
    MaxParentMass = Spectrum->MZ * 3 + (2 * HYDROGEN_MASS);
    
    Spectrum->IntensityBinCount = (MaxParentMass + DALTON) / BinScalingFactor; 
    Spectrum->BinnedIntensities = (float*)calloc(Spectrum->IntensityBinCount, sizeof(float));
    Spectrum->BinnedIntensitiesTight = (float*)calloc(Spectrum->IntensityBinCount, sizeof(float));
    Spectrum->BinnedIntensityLevels = (int*)calloc(Spectrum->IntensityBinCount, sizeof(int));
    
    Spectrum->BinPeakIndex = (int*)calloc(Spectrum->IntensityBinCount, sizeof(int));
    for (Bin = 0; Bin < Spectrum->IntensityBinCount; Bin++)
    {
        Spectrum->BinPeakIndex[Bin] = -1;
    }
    // Iterate over spectral peaks, putting intensity into bins:
    for (PeakIndex = 0; PeakIndex < Spectrum->PeakCount; PeakIndex++)
    {
        Peak = Spectrum->Peaks + PeakIndex;
        Bin = (Peak->Mass + 50) / BinScalingFactor;
        for (NearBin = Bin - 6; NearBin < Bin + 7; NearBin++)
        {
            if (NearBin < 0 || NearBin >= Spectrum->IntensityBinCount)
            {
                continue;
            }
            Skew = abs(Peak->Mass - (NearBin * BinScalingFactor));
            if (Skew > Model->IntensityRadius)
            {
                continue;
            }
            Multiplier = 1.0; // default
            if (Model->IntensityScheme == 1 || Model->IntensityScheme == 3)
            {
                if (Skew >= Model->HalfIntensityRadius)
                {
                    Multiplier = 0.5;
                }
            }
            Spectrum->BinnedIntensities[NearBin] += Peak->Intensity * Multiplier;
            if (Skew < INTENSITY_BIN_RADIUS_TIGHT)
            {
                Spectrum->BinnedIntensitiesTight[NearBin] += Peak->Intensity;
            }
            if (Spectrum->BinPeakIndex[NearBin] < 0)
            {
                Spectrum->BinPeakIndex[NearBin] = PeakIndex;
            }
        }
    }
    // Compute the intensity level (absent, lo, med, hi) for each bin:
    //ComputeSpectrumIntensityCutoffs(Spectrum);
    memset(CountByIntensityLevel, 0, sizeof(int) * 16);
    for (Bin = 0; Bin < Spectrum->IntensityBinCount; Bin++)
    {
        Intensity = Spectrum->BinnedIntensities[Bin];
        for (IntensityLevel = 0; IntensityLevel < 99; IntensityLevel++)
        {
            if (Intensity > Spectrum->IntensityThresholds[IntensityLevel])
            {
                Spectrum->BinnedIntensityLevels[Bin] = IntensityLevel;
                CountByIntensityLevel[IntensityLevel]++;
                break;
            }
        }
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Now let's compute the fraction of mass bins which attain these intensity thresholds 'by chance'.
    // This fraction is used for scoring PRMs; the bonus for having a y peak is smaller for a very
    // thick spectrum than for a very sparse spectrum.
    Spectrum->IonScoringNoiseProbabilities = (float*)calloc(ThresholdCount + 1, sizeof(float));
    for (IntensityLevel = 0; IntensityLevel < ThresholdCount; IntensityLevel++)
    {
        Probability = (CountByIntensityLevel[IntensityLevel] + 1) / (float)Spectrum->IntensityBinCount;
        Spectrum->IonScoringNoiseProbabilities[IntensityLevel] = (float)log(Probability);
    }
}

// Return the intensity level for this mass.  If this is a cut, claim the peaks; otherwise,
// just return the intensity level.
int IonScoringGetPeakIntensity(PRMBayesianModel* Model, MSSpectrum* Spectrum, int Mass, int FragmentType, int SeizePeakAminoIndex)
{
    int Bin;
    int MinMass;
    int MaxMass;
    float Intensity = 0;
    int IntensityLevelIndex;
    int PeakIndex;
    int Skew;
    float Multiplier;
    //
    Bin = (Mass + 50) / 100; // Bin width 0.1Da
    MinMass = Mass - Model->IntensityRadius;
    MaxMass = Mass + Model->IntensityRadius;

    // If the mass is off the scale, then you get no peaks:
    if (Bin >= Spectrum->IntensityBinCount || Bin < 0)
    {
        return Model->MinIntensityLevel;
    }
    
    // If this is a PRM (not a cut), then look up the intensity level
    // in the spectrum's array:
    if (SeizePeakAminoIndex < 0)
    {
        return Spectrum->BinnedIntensityLevels[Bin];
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
            if (Spectrum->Peaks[PeakIndex].Mass < MinMass)
            {
                continue;
            }

            Multiplier = 1.0; // default
            Skew = abs(Mass - Spectrum->Peaks[PeakIndex].Mass);
            if (Model->IntensityScheme == 1 || Model->IntensityScheme == 3)
            {
                if (Skew >= Model->HalfIntensityRadius)
                {
                    Multiplier = 0.5;
                }
            }
            if (Spectrum->Peaks[PeakIndex].IonType)
            {
                // This peak has already been CLAIMED by another ion type:
                continue;
            }
            Intensity += Spectrum->Peaks[PeakIndex].Intensity * Multiplier;
            // CLAIM this spectrum:
            Spectrum->Peaks[PeakIndex].IonType = FragmentType;
            Spectrum->Peaks[PeakIndex].AminoIndex = SeizePeakAminoIndex;
        }
    }
    for (IntensityLevelIndex = 0; IntensityLevelIndex < 99; IntensityLevelIndex++)
    {
        if (Intensity > Spectrum->IntensityThresholds[IntensityLevelIndex])
        {
            return IntensityLevelIndex;
        }
    }
    return 0;
}

// Compute the sector for a given mass.  The sector is a simple partition of
// the mass range (low/high, or low/medium/high, etc).  
int IonScoringGetSector(PRMBayesianNode* Node, int ParentMass, int Mass)
{
    switch (Node->Flag)
    {
    case 0:
        // Two sectors, LOW and HIGH:
        if (Mass < ParentMass / 2)
        {
            return 0;
        }
        else
        {
            return 1;
        }
        break;
    case 1:
        // Three sectors, LOW and MEDIUM and HIGH:
        if (Mass < ParentMass * 0.33)
        {
            return 0;
        }
        if (Mass < ParentMass * 0.66)
        {
            return 1;
        }
        return 2;
        break;
    case 2:
        // Four sectors:
        if (Mass < ParentMass * 0.25)
        {
            return 0;
        }
        if (Mass < ParentMass * 0.5)
        {
            return 1;
        }
        if (Mass < ParentMass * 0.75)
        {
            return 2;
        }
        return 3;
        break;
    case 3:
        // Five sectors:
        if (Mass < ParentMass * 0.2)
        {
            return 0;
        }
        if (Mass < ParentMass * 0.4)
        {
            return 1;
        }
        if (Mass < ParentMass * 0.6)
        {
            return 2;
        }
        if (Mass < ParentMass * 0.8)
        {
            return 3;
        }
        return 4;
        break;
    default:
        REPORT_ERROR(0);
        break;
    }
    return 0;
}

// Compute the value of the Flank feature.  These features reflect flanking amino acids
// which have effects on fragment intensities.
int IonScoringGetFlank(PRMBayesianNode* Node, char Left, char Right)
{
    //
    switch (Node->Flag)
    {
    case 0:
        // Default B flank:
        // G or P on left: Strong suppression
        if (Left == 'G' || Left == 'P')
        {
            return 0;
        }
        // P on right: Augmentation
        if (Right == 'P')
        {
            return 1;
        }
        // H or R on right: Suppression
        if (Right == 'H' || Right == 'R')
        {
            return 2;
        }
        return 3;
        break;
    case 1:
        // Default Y flank:
        // P on right: Strong augmentation
        if (Right == 'P')
        {
            return 0;
        }
        // K or R on right: Strong suppression
        if (Right == 'R' || Right == 'K')
        {
            return 1;
        }
        // H on right or P on left: suppression
        if (Left == 'P' || Right == 'H')
        {
            return 2;
        }
        return 3;
    default:
        REPORT_ERROR(0);
        break;
    }
    return 0;
}

// Compute a feature for whether the N- or C-terminal portion of a peptide contains acidic
// or basic residues.  (Not used in production)
int IonScoringGetFragmentContain(PRMBayesianNode* Node, Peptide* Match, int AminoIndex, int SuffixFlag)
{
    int MinIndex;
    int MaxIndex;
    int CheckIndex;
    int Count = 0;
    //
    if (SuffixFlag)
    {
        MinIndex = AminoIndex;
        MaxIndex = strlen(Match->Bases);
    }
    else
    {
        MinIndex = 0;
        MaxIndex = AminoIndex;
    }
    for (CheckIndex = MinIndex; CheckIndex < MaxIndex; CheckIndex++)
    {
        switch (Match->Bases[CheckIndex])
        {
        case 'D':
        case 'E':
            if (Node->Flag == 0 || Node->Flag == 1)
            {
                Count++;
            }
            break;
        case 'R':
        case 'K':
        case 'H':
            if (Node->Flag == 2 || Node->Flag == 3)
            {
                Count++;
            }
            break;
        default:
            break;
        }
    }
    switch (Node->Flag)
    {
    case 0:
    case 2:
        if (Count)
        {
            return 1;
        }
        else
        {
            return 0;
        }
        break;
    case 1:
    case 3:
        if (Count > 1)
        {
            return 2;
        }
        else if (Count)
        {
            return 1;
        }
        else
        {
            return 0;
        }
        break;
    default:
        REPORT_ERROR(0);
        break;
    }
    return 0; // unreachable
}

int IonScoringGetPrefixContainPhos(PRMBayesianNode* Node, Peptide* Match, int AminoIndex)
{
    int ModIndex;
    int ModifiedResidueIndex = -1;
    //
    for (ModIndex = 0; ModIndex < MAX_PT_MODS; ModIndex++)
    {
        if (!Match->ModType[ModIndex])
        {
            break;
        }
        if (Match->ModType[ModIndex]->Flags & DELTA_FLAG_PHOSPHORYLATION)
        {
            ModifiedResidueIndex = Match->AminoIndex[ModIndex];
            if (ModifiedResidueIndex < AminoIndex)
            {
                return 1;
            }
        }
    }
    //got all the way here without returning anything.
    return 0;
}

int IonScoringGetSuffixContainPhos(PRMBayesianNode* Node, Peptide* Match, int BreakIndex)
{
    int ModIndex;
    int ModifiedResidueIndex = -1;
    //
    for (ModIndex = 0; ModIndex < MAX_PT_MODS; ModIndex++)
    {
        if (!Match->ModType[ModIndex])
        {
            break;
        }
        if (Match->ModType[ModIndex]->Flags & DELTA_FLAG_PHOSPHORYLATION)
        {
            ModifiedResidueIndex = Match->AminoIndex[ModIndex];
            if (BreakIndex <= ModifiedResidueIndex)
            {
                return 1;
            }
        }
    }
    //got all the way here without returning anything.
    return 0;
}

// Return the value for a particular PRM or cut.  This function calls the appropriate setter based 
// on the node type.  
// Important special note: AminoIndex should be -1 if we're getting PRM scores!
int IonScoringGetNodeValue(PRMBayesianModel* Model, PRMBayesianNode* Node, MSSpectrum* Spectrum, int PRM,
    Peptide* Match, int AminoIndex)
{
    int Suffix;
    int PeptideLen;
    char PrefixAA;
    char SuffixAA;
    ///////////////////////////////////////////////////////////////////////////////////////
    // Set values for the current PRM:
    switch (Node->Type)
    {
    case evPRMBPrefix:
        // Handle b peak, or other N-terminal fragment:
        return IonScoringGetPeakIntensity(Model, Spectrum, PRM + Node->MassOffset, Node->FragmentType, AminoIndex);
        //Node->Value = IonScoringSetIntensityLevel(Spectrum, Intensity);
        break;
    case evPRMBPrefix2:
        // Handle doubly-charged N-terminal fragment:
        return IonScoringGetPeakIntensity(Model, Spectrum, (PRM + Node->MassOffset + HYDROGEN_MASS) / 2, Node->FragmentType, AminoIndex);
        //Node->Value = IonScoringSetIntensityLevel(Spectrum, Intensity);
        break;
    case evPRMBSuffix:
        // Handle C-terminal fragment:
        Suffix = Spectrum->ParentMass - PRM;
        return IonScoringGetPeakIntensity(Model, Spectrum, Suffix + Node->MassOffset, Node->FragmentType, AminoIndex);
        //Node->Value = IonScoringSetIntensityLevel(Spectrum, Intensity);
        break;
    case evPRMBSuffix2:
        // Handle doubly-charged C-terminal fragment:
        Suffix = Spectrum->ParentMass - PRM;
        return IonScoringGetPeakIntensity(Model, Spectrum, (Suffix + Node->MassOffset + HYDROGEN_MASS) / 2, Node->FragmentType, AminoIndex);
        //Node->Value = IonScoringSetIntensityLevel(Spectrum, Intensity);
        break;
    case evSector:
        // Handle "sector" (which part of the mass range this mass lies in)
        return IonScoringGetSector(Node, Spectrum->ParentMass, PRM);
        break;
    case evFlank:
        // Handle "flank" (for cuts only: based on prefix and suffix amino acids)
        // If no peptide, return 0 (always the "default" intensity)
        if (!Match)
        {
            return 0;
        }
        PeptideLen = strlen(Match->Bases);
        if (AminoIndex > 0)
        {
            PrefixAA = Match->Bases[AminoIndex - 1];
        }
        else
        {
            PrefixAA = '\0';
        }
        if (AminoIndex < PeptideLen)
        {
            SuffixAA = Match->Bases[AminoIndex];
        }
        else
        {
            SuffixAA = '\0';
        }
        return IonScoringGetFlank(Node, PrefixAA, SuffixAA);
        break;
    case evPrefixAA:
        if (AminoIndex > 0)
        {
            PrefixAA = Match->Bases[AminoIndex - 1];
        }
        else
        {
            PrefixAA = '\0';
        }
        if ((PrefixAA - 'A') == Node->Flag)
        {
            return 1;
        }
        else
        {
            return 0;
        }
        break;
    case evSuffixAA:
        PeptideLen = strlen(Match->Bases);
        if (AminoIndex < PeptideLen)
        {
            SuffixAA = Match->Bases[AminoIndex];
        }
        else
        {
            SuffixAA = '\0';
        }
        if ((SuffixAA - 'A') == Node->Flag)
        {
            return 1;
        }
        else
        {
            return 0;
        }
        break;
    case evPrefixContain:
        return IonScoringGetFragmentContain(Node, Match, AminoIndex, 0);
    case evSuffixContain:
        return IonScoringGetFragmentContain(Node, Match, AminoIndex, 1);
    case evPrefixContainPhos:
        return IonScoringGetPrefixContainPhos(Node, Match, AminoIndex);
    case evSuffixContainPhos:
        return IonScoringGetSuffixContainPhos(Node, Match, AminoIndex);
    default:
        REPORT_ERROR(0);
        break;
    }
    return 0;
}

// For debugging purposes, print out the definition of a PRMBayesianModel.
void DebugPrintPRMBayesianModel(PRMBayesianModel* Model)
{
    PRMBayesianNode* Node;
    PRMBayesianNode* Parent;
    int NodeIndex;
    int ParentIndex;
    printf(">>>DebugPrintPRMBayesianModel\n");
    printf("CutFlag %d IntensityRadius %.2f IntensityScheme %d\n", Model->CutFlag, Model->IntensityRadius / (float)DALTON, Model->IntensityScheme);
    for (Node = Model->Head, NodeIndex = 0; Node; Node = Node->Next, NodeIndex++)
    {
        printf(">>Node %d of %d %s:\n", NodeIndex, Model->NodeCount, Node->Name);
        printf("  Type %d flag %d mass offset %.2f\n", Node->Type, Node->Flag, Node->MassOffset / (float)DALTON);
        printf("  Valuecount %d\n", Node->ValueCount);
        for (ParentIndex = 0; ParentIndex < Node->ParentCount; ParentIndex++)
        {
            Parent = Node->Parents[ParentIndex];
            printf("  Parent %d of %d: %s\n", ParentIndex, Node->ParentCount, Parent->Name);
        }
    }
    printf(">>> End of model <<<\n");
}

// The tag-scoring Bayesian network includes some features which rely upon flanking amino acids.
// These nodes must be visited during tag generation, when the flanking amino acids are finally 
// learned.  To save time, we build up a singly-linked list (Model->FirstFlank...Model->LastFlank)
// to keep track of such nodes.
void BuildModelFlankList(PRMBayesianModel* Model)
{
    int NodeIndex;
    PRMBayesianNode* Node;
    PRMBayesianNode* Parent;
    PRMBayesianNodeHolder* Holder;
    int ParentIndex;
    //
    // Set flank flags of all nodes:
    for (NodeIndex = 0; NodeIndex < Model->NodeCount; NodeIndex++)
    {
        Node = Model->Nodes[NodeIndex];
        Node->FlankFlag = 0; // default
        if (Node->Type == evFlank || Node->Type == evPrefixAA || Node->Type == evSuffixAA)
        {
            Node->FlankFlag = 1;
        }
        for (ParentIndex = 0; ParentIndex < Node->ParentCount; ParentIndex++)
        {
            Parent = Node->Parents[ParentIndex];
            if (Parent->Type == evFlank || Parent->Type == evPrefixAA || Parent->Type == evSuffixAA)
            {
                Node->FlankFlag = 1;
            }
        }
    }
    // Build linked list of nodes which rely upon flanking amino acid info:
    for (NodeIndex = 0; NodeIndex < Model->NodeCount; NodeIndex++)
    {
        Node = Model->Nodes[NodeIndex];
        if (Node->FlankFlag)
        {
            // Add a NodeHolder for this node:
            Holder = (PRMBayesianNodeHolder*)calloc(1, sizeof(PRMBayesianNodeHolder));
            Holder->Node = Node;
            if (Model->FirstFlank)
            {
                Model->LastFlank->Next = Holder;
            }
            else
            {
                Model->FirstFlank = Holder;
            }
            Model->LastFlank = Holder;
        }
    }
}

// Save a PRMBayesianNode to a binary file.  Helper function for SavePRMBayesianModel.
void SavePRMBayesianNode(PRMBayesianNode* Node, FILE* ModelFile)
{
    int ParentIndex;
    WriteBinary(&Node->Name, sizeof(char), 256, ModelFile);
    WriteBinary(&Node->Type, sizeof(int), 1, ModelFile);
    WriteBinary(&Node->Flag, sizeof(int), 1, ModelFile);
    WriteBinary(&Node->FragmentType, sizeof(int), 1, ModelFile);
    WriteBinary(&Node->MassOffset, sizeof(int), 1, ModelFile);
    WriteBinary(&Node->ValueCount, sizeof(int), 1, ModelFile);
    WriteBinary(&Node->ParentCount, sizeof(int), 1, ModelFile);
    // Write parent indices:
    for (ParentIndex = 0; ParentIndex < Node->ParentCount; ParentIndex++)
    {
        WriteBinary(&Node->Parents[ParentIndex]->Index, sizeof(int), 1, ModelFile);
    }
    WriteBinary(Node->ParentBlocks, sizeof(int), Node->ParentCount, ModelFile);
    WriteBinary(&Node->TableSize, sizeof(int), 1, ModelFile);
    WriteBinary(Node->CountTable, sizeof(int), Node->TableSize, ModelFile);
    WriteBinary(Node->ProbTable, sizeof(float), Node->TableSize, ModelFile);
}

// Load a PRMBayesianNode from a binary file.  Helper function for LoadPRMBayesianModel.
PRMBayesianNode* LoadPRMBayesianNode(PRMBayesianModel* Model, FILE* ModelFile)
{
    PRMBayesianNode* Node;
    int ParentIndex;
    int ParentNodeIndex;
    //
    Node = (PRMBayesianNode*)calloc(1, sizeof(PRMBayesianNode));
    ReadBinary(&Node->Name, sizeof(char), 256, ModelFile);
    ReadBinary(&Node->Type, sizeof(int), 1, ModelFile);
    ReadBinary(&Node->Flag, sizeof(int), 1, ModelFile);
    ReadBinary(&Node->FragmentType, sizeof(int), 1, ModelFile);
    ReadBinary(&Node->MassOffset, sizeof(int), 1, ModelFile);
    ReadBinary(&Node->ValueCount, sizeof(int), 1, ModelFile);
    ReadBinary(&Node->ParentCount, sizeof(int), 1, ModelFile);
    if (Node->ParentCount < 0 || Node->ParentCount > 100)
    {
        REPORT_ERROR(0);
        return NULL;
    }
    if (Node->ParentCount)
    {
        Node->Parents = (PRMBayesianNode**)calloc(Node->ParentCount, sizeof(PRMBayesianNode*));
        for (ParentIndex = 0; ParentIndex < Node->ParentCount; ParentIndex++)
        {
            ReadBinary(&ParentNodeIndex, sizeof(int), 1, ModelFile);
            if (ParentNodeIndex < 0 || ParentNodeIndex >= Model->NodeCount)
            {
                REPORT_ERROR(0);
                return NULL;
            }
            Node->Parents[ParentIndex] = Model->Nodes[ParentNodeIndex];
        }
        Node->ParentBlocks = (int*)calloc(Node->ParentCount, sizeof(int));
        ReadBinary(Node->ParentBlocks, sizeof(int), Node->ParentCount, ModelFile);
    }
    ReadBinary(&Node->TableSize, sizeof(int), 1, ModelFile);
    if (Node->TableSize <= 0 || Node->TableSize > 10000)
    {
        REPORT_ERROR(0);
        return NULL;
    }
    Node->CountTable = (int*)calloc(Node->TableSize, sizeof(int));
    ReadBinary(Node->CountTable, sizeof(int), Node->TableSize, ModelFile);
    Node->ProbTable = (float*)calloc(Node->TableSize, sizeof(float));
    ReadBinary(Node->ProbTable, sizeof(float), Node->TableSize, ModelFile);
    return Node;
}

// Save a PRMBayesian model to a binary file.  In production, the model
// is loaded (using LoadPRMBayesianModel) and then used.
void SavePRMBayesianModel(PRMBayesianModel* Model, char* FileName)
{
    FILE* ModelFile;
    PRMBayesianNode* Node;
    //
    if (!Model)
    {
        REPORT_ERROR(0);
        return;
    }
    ModelFile = fopen(FileName, "wb");
    if (!ModelFile)
    {
        REPORT_ERROR_S(8, FileName);
        return;
    }
    WriteBinary(&Model->CutFlag, sizeof(int), 1, ModelFile);
    WriteBinary(&Model->IntensityScheme, sizeof(int), 1, ModelFile);
    WriteBinary(&Model->MinIntensityLevel, sizeof(int), 1, ModelFile);
    WriteBinary(&Model->IntensityRadius, sizeof(int), 1, ModelFile);
    fwrite(&Model->NoiseModel, sizeof(int), 1, ModelFile);
    WriteBinary(Model->RandomIntensityCounts, sizeof(int), 10, ModelFile);
    WriteBinary(Model->RandomIntensityScores, sizeof(float), 10, ModelFile);
    WriteBinary(&Model->NodeCount, sizeof(int), 1, ModelFile);
    for (Node = Model->Head; Node; Node = Node->Next)
    {
        SavePRMBayesianNode(Node, ModelFile);
    }
    fclose(ModelFile);
}

// Load a PRMBayesianModel from a binary file.
PRMBayesianModel* LoadPRMBayesianModel(char* FileName)
{
    PRMBayesianModel* Model;
    FILE* ModelFile;
    int NodeIndex;
    PRMBayesianNode* Node;

    //
    ModelFile = fopen(FileName, "rb");
    if (!ModelFile)
    {
        REPORT_ERROR_S(8, FileName);
        return NULL;
    }
    Model = (PRMBayesianModel*)calloc(1, sizeof(PRMBayesianModel));
    ReadBinary(&Model->CutFlag, sizeof(int), 1, ModelFile);
    ReadBinary(&Model->IntensityScheme, sizeof(int), 1, ModelFile);
    ReadBinary(&Model->MinIntensityLevel, sizeof(int), 1, ModelFile);
    ReadBinary(&Model->IntensityRadius, sizeof(int), 1, ModelFile);
    ReadBinary(&Model->NoiseModel, sizeof(int), 1, ModelFile);
    ReadBinary(Model->RandomIntensityCounts, sizeof(int), 10, ModelFile);
    ReadBinary(Model->RandomIntensityScores, sizeof(float), 10, ModelFile);
    ReadBinary(&Model->NodeCount, sizeof(int), 1, ModelFile);
    Model->Nodes = (PRMBayesianNode**)calloc(Model->NodeCount, sizeof(PRMBayesianNode*));
    for (NodeIndex = 0; NodeIndex < Model->NodeCount; NodeIndex++)
    {
        Node = LoadPRMBayesianNode(Model, ModelFile);
        Node->Index = NodeIndex;
        Model->Nodes[NodeIndex] = Node;
        if (Model->Tail)
        {
            Model->Tail->Next = Node;
        }
        else
        {
            Model->Head = Node;
        }
        Model->Tail = Node;
    }
    BuildModelFlankList(Model);
    fclose(ModelFile);
    return Model;
}

// Translate the CountTables for this model's nodes into probability tables.
// We use a "buffer" count for each node to pad out the probabilities; if our training
// set was small, it may have left ZERO entries in some nodes, and we don't want 
// probabilities of zero (since then we can't take their natural logarithm).
void ComputePRMBayesianModelProbabilityTables(PRMBayesianModel* Model, int PaddingCount)
{
    PRMBayesianNode* Node;
    int TotalEntries;
    int TableIndex;
    float Probability;
    int Count;
    int BlockStartIndex;
    int IntensityLevel;
    //

    // Set global noise probabilities:
    Count = 0;
    for (IntensityLevel = 0; IntensityLevel <= Model->MinIntensityLevel; IntensityLevel++)
    {
        Count += (1 + Model->RandomIntensityCounts[IntensityLevel]);
    }
    for (IntensityLevel = 0; IntensityLevel <= Model->MinIntensityLevel; IntensityLevel++)
    {
        Probability = (1 + Model->RandomIntensityCounts[IntensityLevel]) / (float)Count;
        Model->RandomIntensityScores[IntensityLevel] = (float)log(Probability);
    }

    // Set probabilities for each node:
    for (Node = Model->Head; Node; Node = Node->Next)
    {
        // Compute the probability that this node will have a value,
        // GIVEN the values of any parent nodes:
        for (BlockStartIndex = 0; BlockStartIndex < Node->TableSize; BlockStartIndex += Node->ValueCount)
        {
            TotalEntries = 0;
            for (TableIndex = BlockStartIndex; TableIndex < BlockStartIndex + Node->ValueCount; TableIndex++)
            {
                TotalEntries += Node->CountTable[TableIndex] + PaddingCount;
            }
            for (TableIndex = BlockStartIndex; TableIndex < BlockStartIndex + Node->ValueCount; TableIndex++)
            {
                if (TableIndex >= Node->TableSize)
                {
                    REPORT_ERROR(0);
                }
                Count = Node->CountTable[TableIndex] + PaddingCount;
                Probability = Count / (float)TotalEntries;
                Node->ProbTable[TableIndex] = (float)log(Probability);
            }
        }
    }
}

// PRMBNGetCutScore returns the score for a cut-point.  
// It's called AFTER setting the Values array for each node, with calls to IonScoringGetNodeValue
float PRMBNGetCutScore(MSSpectrum* Spectrum, PRMBayesianModel* Model, int AminoIndex)
{
    float Score = 0;
    float NodeScore = 0;
    int TableIndex;
    int ParentIndex;
    PRMBayesianNode* Node;
    int VerboseFlag = 0;
    //
    for (Node = Model->Head; Node; Node = Node->Next)
    {
        switch (Node->Type)
        {
        case evPRMBPrefix:
        case evPRMBPrefix2:
        case evPRMBSuffix:
        case evPRMBSuffix2:
            TableIndex = Node->Values[AminoIndex];
            for (ParentIndex = 0; ParentIndex < Node->ParentCount; ParentIndex++)
            {
                TableIndex += Node->Parents[ParentIndex]->Values[AminoIndex] * Node->ParentBlocks[ParentIndex];
            }
            if (TableIndex >= Node->TableSize)
            {
                REPORT_ERROR(0);
            }
            NodeScore = Node->ProbTable[TableIndex];
            if (Model->NoiseModel)
            {
                // GLOBAL noise model, based on all spectra
                NodeScore -= Model->RandomIntensityScores[Node->Values[AminoIndex]];
            }
            else
            {
                // SPECTRUM noise model:
                NodeScore -= Spectrum->IonScoringNoiseProbabilities[Node->Values[AminoIndex]];
            }
            if (VerboseFlag)
            {
                printf("  AA %d: Node %d (%s) contributes %.3f - %.3f = %.5f\n", AminoIndex, Node->Index, Node->Name, 
                    Node->ProbTable[TableIndex], Spectrum->IonScoringNoiseProbabilities[Node->Values[AminoIndex]], NodeScore);
            }
            Score += NodeScore;
            break;
        default:
            // Other node-types don't contribute to the score.
            break;
        }
    }
    return Score;
}

// Compute the score for a PRM, using a bayesian network.  Sum the log-probabilities over
// all ion fragment nodes.  
float GetIonPRMFeatures(MSSpectrum* Spectrum, SpectrumTweak* Tweak, PRMBayesianModel* Model, int PRM, int VerboseFlag)
{
    PRMBayesianNode* Node;
    int ParentIndex;
    int TableIndex;
    float Score = 0;
    float NodeScore;
    //
    // Compute each node value:
    for (Node = Model->Head; Node; Node = Node->Next)
    {
        Node->Value = IonScoringGetNodeValue(Model, Node, Spectrum, PRM, NULL, -1);
        if (VerboseFlag)
        {
            printf("Score(%.2f): Node %d (%s) has value %d\n", PRM / (float)DALTON, Node->Index, Node->Name, Node->Value);
        }
    }
    // Compute a SCORE for this collection of values:
    for (Node = Model->Head; Node; Node = Node->Next)
    {
        switch (Node->Type)
        {
        case evPRMBPrefix:
        case evPRMBPrefix2:
        case evPRMBSuffix:
        case evPRMBSuffix2:
            TableIndex = Node->Value;
            for (ParentIndex = 0; ParentIndex < Node->ParentCount; ParentIndex++)
            {
                TableIndex += Node->Parents[ParentIndex]->Value * Node->ParentBlocks[ParentIndex];
            }
            if (TableIndex >= Node->TableSize)
            {
                REPORT_ERROR(0);
            }
            NodeScore = Node->ProbTable[TableIndex];
            if (Model->NoiseModel)
            {
                // GLOBAL noise model, based on all spectra
                NodeScore -= Model->RandomIntensityScores[Node->Value];
            }
            else
            {
                // SPECTRUM noise model:
                NodeScore -= Spectrum->IonScoringNoiseProbabilities[Node->Value];
            }
            if (VerboseFlag)
            {
                printf("  Node %d (%s) contributes %.3f - %.3f = %.5f\n", Node->Index, Node->Name, 
                    Node->ProbTable[TableIndex], Spectrum->IonScoringNoiseProbabilities[Node->Value], NodeScore);
            }
            Score += NodeScore;
            break;
        default:
            // Other node-types don't contribute to the score.
            break;
        }
    }
    return Score;
}

// Iterate over all the nodes in our TagGraph, and assign a score to each.
void TagGraphScorePRMNodes(PRMBayesianModel* Model, TagGraph* Graph, MSSpectrum* Spectrum, SpectrumTweak* Tweak)
{
    TagGraphNode* Node;

    if (!Model)
    {
        if (Tweak->Charge < 3)
        {
            Model = TAGModelCharge2;
        }
        else
        {
            Model = TAGModelCharge3;
        }
    }

    for (Node = Graph->FirstNode; Node; Node = Node->Next)
    {
        if (Node->NodeType != evGraphNodeB && Node->NodeType != evGraphNodeY)
        {
            Node->Score = 0;
            continue;
        }
        Node->Score = GetIonPRMFeatures(Spectrum, Tweak, Model, Node->Mass, 0);
        continue; 
    }
}

void FreeBayesianModels()
{
    FreePRMBayesianModel(PRMModelCharge2);
    PRMModelCharge2 = NULL;
    FreePRMBayesianModel(PRMModelCharge3);
    PRMModelCharge3 = NULL;
    FreePRMBayesianModel(TAGModelCharge2);
    TAGModelCharge2 = NULL;
    FreePRMBayesianModel(TAGModelCharge3);
    TAGModelCharge3 = NULL;
    FreePRMBayesianModel(PhosCutModelCharge2);
    PhosCutModelCharge2 = NULL;
    FreePRMBayesianModel(PhosCutModelCharge3);
    PhosCutModelCharge3 = NULL;
}

// Load PRMBayesianModel objects for scoring PRMs and for scoring tags.
void InitBayesianModels()
{
    char FilePath[2048];
    // Return immediately, if models are loaded already:
    if (PRMModelCharge2)
    {
        return;
    }
    sprintf(FilePath, "%s%s", GlobalOptions->ResourceDir, "PRM2.bn");
    PRMModelCharge2 = LoadPRMBayesianModel(FilePath);
    sprintf(FilePath, "%s%s", GlobalOptions->ResourceDir, "PRM3.bn");
    PRMModelCharge3 = LoadPRMBayesianModel(FilePath);
    sprintf(FilePath, "%s%s", GlobalOptions->ResourceDir, "TAG2.bn");
    TAGModelCharge2 = LoadPRMBayesianModel(FilePath);
    sprintf(FilePath, "%s%s", GlobalOptions->ResourceDir, "TAG3.bn");
    TAGModelCharge3 = LoadPRMBayesianModel(FilePath);
    sprintf(FilePath, "%s%s", GlobalOptions->ResourceDir, "PhosCut2.bn");
    PhosCutModelCharge2 = LoadPRMBayesianModel(FilePath);
    sprintf(FilePath, "%s%s", GlobalOptions->ResourceDir, "PhosCut3.bn");
    PhosCutModelCharge3 = LoadPRMBayesianModel(FilePath);
}

// Replace a PRMScoring model with one specified in the input file (the "PRMModel" option).
// Useful for handling new instrument types, etc.
int ReplacePRMScoringModel(int Charge, char* FileName)
{
    PRMBayesianModel* Model;
    //
    Model = LoadPRMBayesianModel(FileName);
    if (!Model)
    {
        return 0;
    }
    if (Charge == 2)
    {
        FreePRMBayesianModel(PRMModelCharge2);
        PRMModelCharge2 = Model;
    }
    else if (Charge == 3)
    {
        FreePRMBayesianModel(PRMModelCharge3);
        PRMModelCharge3 = Model;
    }
    else
    {
        REPORT_ERROR(0);
    }

    return 1;
}

// Replace a tag scoring model with one specified in the input file (the "TAGModel" option).
// Useful for handling new instrument types, etc.
int ReplaceTAGScoringModel(int Charge, char* FileName)
{
    PRMBayesianModel* Model;
    //
    Model = LoadPRMBayesianModel(FileName);
    if (!Model)
    {
        return 0;
    }
    if (Charge == 2)
    {
        FreePRMBayesianModel(TAGModelCharge2);
        TAGModelCharge2 = Model;
    }
    else if (Charge == 3)
    {
        FreePRMBayesianModel(TAGModelCharge3);
        TAGModelCharge3 = Model;
    }
    else
    {
        REPORT_ERROR(0);
    }
    return 1;
}

// Set the array Tweak->PRMScores.  This is used in unrestrictive ("blind") searches.
void SetSpectrumPRMScores(MSSpectrum* Spectrum, SpectrumTweak* Tweak)
{
    PRMBayesianModel* Model;
    int PRM;
    float fScore;
    //
    // Ensure models are loaded:
    if (!PRMModelCharge2)
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
        Model = PRMModelCharge3;
    }
    else
    {
        Model = PRMModelCharge2;
    }
    for (PRM = 0; PRM < Tweak->PRMScoreMax; PRM++)
    {
        fScore = GetIonPRMFeatures(Spectrum, Tweak, Model, PRM * PRM_BIN_SIZE, 0);
        //GetPRMFeatures(Spectrum, Tweak, Model, PRM * PRM_BIN_SIZE, 0);
        Tweak->PRMScores[PRM] = (int)(fScore * 1000);
    }
    //DebugPrintPRMScores(Spectrum, Tweak);
}

int CountTrypticTermini(Peptide* Match)
{
    int NTT = 0;
    int PeptideLength = strlen(Match->Bases);
    switch (GlobalOptions->DigestType)
    {
    case DIGEST_TYPE_TRYPSIN:
        /////////////////////////////////
        // Number of tryptic termini:
        NTT = 0;
        if (Match->PrefixAmino == '\0' || Match->PrefixAmino == '-' || Match->PrefixAmino == '*')
        {
            NTT++;
        }
        else if ((Match->PrefixAmino == 'K' || Match->PrefixAmino == 'R') && Match->Bases[0] != 'P')
        {
            NTT++;
        }
        if (Match->SuffixAmino == '\0' || Match->SuffixAmino == '-' || Match->SuffixAmino == '*')
        {
            NTT++;
        }
        else if ((Match->Bases[PeptideLength - 1] == 'K' || Match->Bases[PeptideLength - 1] == 'R') && Match->SuffixAmino != 'P')
        {
            NTT++;
        }
        break;
    case DIGEST_TYPE_UNKNOWN:
    default:
        NTT = 2;
        break;
    }
    return NTT;
}

void PopulateCutScores(PRMBayesianModel* Model, MSSpectrum* Spectrum, Peptide* Match, float* CutScores)
{
    int PRM = 0;  
    int NodeIndex;
    PRMBayesianNode* Node;
    int AminoIndex;
    int ModIndex;
    int PeptideLength = strlen(Match->Bases);
    int PeakIndex;

    // Reset all peak annotations:
    for (PeakIndex = 0; PeakIndex < Spectrum->PeakCount; PeakIndex++)
    {
        Spectrum->Peaks[PeakIndex].IonType = evFragmentTypeNone;
    }
    if (Match->SpecialFragmentation)
    { // phosphorylated spectra
        AnnotateParentPeaks(Spectrum, Match, Model);
    }

    for (NodeIndex = 0, Node = Model->Head; Node; NodeIndex++, Node = Node->Next)
    {
        PRM = 0;
        for (AminoIndex = 0; AminoIndex <= PeptideLength; AminoIndex++)
        {
            ///////////////////////////////////////////////////////////////////////////////////////
            // Set values, and accumulate table entries:
            Node->Values[AminoIndex] = IonScoringGetNodeValue(Model, Node, Spectrum, PRM, Match, AminoIndex);
            ///////////////////////////////////////////////////////////////////////////////////////
            // Add to PRM:
            if (AminoIndex == PeptideLength)
            {
                break;
            }
            PRM += PeptideMass[Match->Bases[AminoIndex]];
            for (ModIndex = 0; ModIndex < MAX_PT_MODS; ModIndex++)
            {
                if (Match->AminoIndex[ModIndex] == AminoIndex)
                {
                    PRM += Match->ModType[ModIndex]->RealDelta;
                }
            }
        } // Amino loop
    } // NodeIndex loop

    // Populate the CutScores array:
    for (AminoIndex = 0; AminoIndex <= PeptideLength; AminoIndex++)
    {
        CutScores[AminoIndex] = PRMBNGetCutScore(Spectrum, Model, AminoIndex);
    }
}

// Compute MQScore features, in preparation for MQScore calculation
int ComputeMQScoreFeatures(MSSpectrum* Spectrum, Peptide* Match, float* MQFeatures, int VerboseFlag)
{
    int FeatureIndex = 0;
    PRMBayesianModel* Model;
    int PeptideLength;
    float CutScores[256];
    int PRM = 0;  
    int AminoIndex;
    int PRMCount;
    float ScoreTotal;
    int YFlag[256];
    int BFlag[256];
    int PeakIndex;
    int PresentCount;
    int FragmentType;
    float PeakIntensity;
    float TotalIntensity;
    float IntensityY = 0;
    float IntensityYSeries = 0;
    float IntensityB = 0;
    float IntensityBSeries = 0;
    //
    Spectrum->ParentMass = Match->ParentMass;
    Model = GetScoringModel(Match, Spectrum->Charge);
    PeptideLength = strlen(Match->Bases);
    // If the peptide is very short (length 5 or less), wey may not even want to bother
    // computing features.  Peptides that short are not informative!

    MQFeatures[FeatureIndex++] = (float)PeptideLength; // #2
    
    ///////////////////////////////////////
    // Cut score features (5, 11):
    PopulateCutScores(Model, Spectrum, Match, CutScores);

    // Total/mean for CENTRAL cut scores:
    ScoreTotal = 0;
    PRMCount = 0;
    for (AminoIndex = 1; AminoIndex < PeptideLength; AminoIndex++)
    {
        ScoreTotal += CutScores[AminoIndex];
        PRMCount++;
        if (VerboseFlag)
        {
            printf("  Cut score %d: %.2f\n", AminoIndex, CutScores[AminoIndex]);
        }
    }
    MQFeatures[FeatureIndex++] = ScoreTotal; // #5

    // Median cut score: 
    PRMCount = PeptideLength + 1;
    MQFeatures[FeatureIndex++] = GetMedian(CutScores, PRMCount); // #11

    // Count b and y peak presence:
    memset(BFlag, 0, sizeof(int) * (PeptideLength + 1));
    memset(YFlag, 0, sizeof(int) * (PeptideLength + 1));
    TotalIntensity = 0;
    for (PeakIndex = 0; PeakIndex < Spectrum->PeakCount; PeakIndex++)
    {
        FragmentType = Spectrum->Peaks[PeakIndex].IonType;
        PeakIntensity = Spectrum->Peaks[PeakIndex].Intensity;
        if (FragmentType == evParentLoss)
        { 
            // I don't want the parent loss peaks to count against
            // phosphorylation, because it is typically very strong.
            PeakIntensity = 0;
        }
        TotalIntensity += PeakIntensity;
        switch (FragmentType)
        {
            case evFragmentY:
                IntensityY += PeakIntensity;
                IntensityYSeries += PeakIntensity;
                YFlag[Spectrum->Peaks[PeakIndex].AminoIndex] = 1;
                break;
            case evFragmentYLoss:
                IntensityYSeries += PeakIntensity;
                break;
            case evFragmentB:
                IntensityB += PeakIntensity;
                IntensityBSeries += PeakIntensity;
                BFlag[Spectrum->Peaks[PeakIndex].AminoIndex] = 1;
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
    MQFeatures[FeatureIndex++] = PresentCount / (float)(PeptideLength + 1); // #12
    PresentCount = 0;
    for (AminoIndex = 0; AminoIndex <= PeptideLength; AminoIndex++)
    {
        PresentCount += BFlag[AminoIndex];
    }
    MQFeatures[FeatureIndex++] = PresentCount / (float)(PeptideLength + 1); // #13
    
    // Fraction of total intensity in B and Y series:
    MQFeatures[FeatureIndex++] = (IntensityY + IntensityB) / TotalIntensity; // #25

    MQFeatures[FeatureIndex++] = (float)CountTrypticTermini(Match); // #30
    return FeatureIndex;
}

// This is currently only called for Phosphorylated spectra, for that reason we
// claim the peaks corresponding to Parent-phosphate, and parent-phosphate-water
// this bears resemblance to IonScoringGetPeakIntensity, when we claim the peak
// as this does not correspond to an AminoIndex, yet requires peak claiming.
// we have to rewrite it here.
void AnnotateParentPeaks(MSSpectrum* Spectrum, Peptide* Match, PRMBayesianModel* Model)
{
    int Loss;
    //the loss of phosphate from the spectrum is actually 98 daltons.(80 + 18)
    int PMMinusPhosphate;
    //actually phosphate and twoWaters (80+18+18)
    int PMMinusPhosphateAndWater;
    //
    //set mz according to the current parent mass of the spectrum
    Spectrum->MZ = (Spectrum->ParentMass + (Spectrum->Charge - 1) * HYDROGEN_MASS) / Spectrum->Charge;
    Loss = PHOSPHATE_WATER_MASS / Spectrum->Charge;
    PMMinusPhosphate = Spectrum->MZ - Loss;
    Loss = (PHOSPHATE_WATER_MASS + WATER_MASS)/Spectrum->Charge;
    PMMinusPhosphateAndWater = Spectrum->MZ - Loss;
    ClaimParentPeak(Spectrum, Match, PMMinusPhosphate, Model);
    ClaimParentPeak(Spectrum, Match, PMMinusPhosphateAndWater, Model);
    // Now look for +1 isotopes
    Loss = PHOSPHATE_WATER_MASS / Spectrum->Charge;
    PMMinusPhosphate = Spectrum->MZ - Loss + (HYDROGEN_MASS/Spectrum->Charge);
    Loss = (PHOSPHATE_WATER_MASS + WATER_MASS) / Spectrum->Charge;
    PMMinusPhosphateAndWater = Spectrum->MZ - Loss + (HYDROGEN_MASS / Spectrum->Charge);
    ClaimParentPeak(Spectrum, Match, PMMinusPhosphate, Model);
    ClaimParentPeak(Spectrum, Match, PMMinusPhosphateAndWater, Model);
}

void ClaimParentPeak(MSSpectrum* Spectrum, Peptide* Match, int Mass, PRMBayesianModel* Model)
{
    int Bin;
    int MinMass;
    int MaxMass;
    float Intensity = 0;
    int PeakIndex;
    int Skew;
    float Multiplier;
    //
    Bin = (Mass + 50) / 100; // Bin width 0.1Da
    MinMass = Mass - Model->IntensityRadius;
    MaxMass = Mass + Model->IntensityRadius;
    
    // If the mass is off the scale, then you get no peaks:
    if (Bin >= Spectrum->IntensityBinCount || Bin < 0)
    {
        return;
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
            if (Spectrum->Peaks[PeakIndex].Mass < MinMass)
            {
                continue;
            }

            Multiplier = 1.0; // default
            Skew = abs(Mass - Spectrum->Peaks[PeakIndex].Mass);
            if (Model->IntensityScheme == 1 || Model->IntensityScheme == 3)
            {
                if (Skew >= Model->HalfIntensityRadius)
                {
                    Multiplier = 0.5;
                }
            }
            if (Spectrum->Peaks[PeakIndex].IonType)
            {
                // This peak has already been CLAIMED by another ion type:
                continue;
            }
            Intensity += Spectrum->Peaks[PeakIndex].Intensity * Multiplier;
            // CLAIM this spectrum:
            Spectrum->Peaks[PeakIndex].IonType = evParentLoss;
            Spectrum->Peaks[PeakIndex].AminoIndex = -1; //not an amino index.  is this a problem?
        }
    }
}

PRMBayesianModel* GetScoringModel(Peptide* Match, int Charge)
{
    int ModIndex;
    //
    for (ModIndex = 0; ModIndex < MAX_PT_MODS; ModIndex++)
    {
        if (!Match->ModType[ModIndex])
        {
            break;
        }
        if (Match->ModType[ModIndex]->Flags & DELTA_FLAG_PHOSPHORYLATION)
        {
            Match->SpecialFragmentation = FRAGMENTATION_PHOSPHO;
            Match->SpecialModPosition = Match->AminoIndex[ModIndex];
            break;
        }
    }
    if (Match->SpecialFragmentation)
    {
        if (Charge > 2)
        {
            return PhosCutModelCharge3;
        }
        return PhosCutModelCharge2;
    }
    if (Charge > 2)
    {
        return TAGModelCharge3;
    }
    return TAGModelCharge2; //default
}

char* GetFragmentTypeName(int FragmentType)
{
    switch (FragmentType)
    {
    case evFragmentY:
        return "Y";
    case evFragmentB:
        return "B";
    case evFragmentYLoss:
        return "Y loss";
    case evFragmentBLoss:
        return "B loss";
    case evParentLoss:
        return "Parent loss";
    case evFragmentTypeNone:
    default:
        return "";
    }
}

