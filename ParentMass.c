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

// ParentMass.c: Routines for parent mass correction.  The precursor mass, as supplied,
// may be off by up to 1 Da (or more, depending on the experiment).  Here we determine
// which parent mass is correct by considering the spectrum's self-convolution: The
// overlap between b and y peaks should be highest when the parent mass is exactly right.
//
// Our implementation: We construct a PMCSpectrumInfo object for the spectrum, which keeps
// track of PMCInfo nodes.  We build one PMCInfo node for each mass we're testing.  We
// compute self-convolutions for each PMCInfo node, and we also compare convolutions across
// PMCInfo nodes and across mass offsets.  Finally, we feed these features into a model
// which assigns each PMCInfo a score, and we keep the best PMCInfo.

#include "CMemLeak.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include "Utils.h"
#include "ChargeState.h"
#include "Spectrum.h"
#include "Inspect.h"
#include "SVM.h"
#include "Errors.h"
#include "LDA.h"

#ifdef _WIN32
#include <windows.h>
#include <sys/types.h>
#include <sys/stat.h>
#else
#include <dirent.h>
#include <sys/stat.h>
#endif

#define EPSILON (float)0.000001

// Models, for parent mass correction in various charge states:
extern LDAModel* PMCCharge1LDA;
extern LDAModel* PMCCharge2LDA;
extern LDAModel* PMCCharge3LDA;
SVMModel* PMCCharge1SVM;
SVMModel* PMCCharge2SVM;
SVMModel* PMCCharge3SVM;

// For converting parts-per-million:
#define ONE_MILLION 1000000

///////////////////////////////////////////////////
// Forward declarations:

///////////////////////////////////////////////////
// Functions:
void CharacterizePhosphatePeaks(PMCInfo* Info, PMCSpectrumInfo* SpectrumInfo, int Offset, int FeatureIndex);

// Free PMCSpectrumInfo, which is only kept around during
// parent mass and charge state correction.
void FreePMCSpectrumInfo(PMCSpectrumInfo* SpectrumInfo)
{
    PMCInfo* Info;
    PMCInfo* Prev;
    SelfConvolutionNode* Node;
    SelfConvolutionNode* PrevNode;
    int HashIndex;
    //
    if (!SpectrumInfo)
    {
        return;
    }
    // Free PMCInfo list:
    Prev = NULL;
    for (Info = SpectrumInfo->Head; Info; Info = Info->Next)
    {
        SafeFree(Prev);
        Prev = Info;
    }
    SafeFree(Prev);
    // Free SelfConvolution list:
    for (HashIndex = 0; HashIndex < SC_HASH_SIZE; HashIndex++)
    {
        PrevNode = NULL;
        for (Node = SpectrumInfo->SCHash[HashIndex]; Node; Node = Node->Next)
        {
            SafeFree(PrevNode);
            PrevNode = Node;
        }
        SafeFree(PrevNode);
    }
    // Free SelfConvolution2 list:
    for (HashIndex = 0; HashIndex < SC_HASH_SIZE; HashIndex++)
    {
        PrevNode = NULL;
        for (Node = SpectrumInfo->SC2Hash[HashIndex]; Node; Node = Node->Next)
        {
            SafeFree(PrevNode);
            PrevNode = Node;
        }
        SafeFree(PrevNode);
    }
    // Free the parent:
    SafeFree(SpectrumInfo);
}

// Build PMCInfo nodes for the masses we'll consider adjusting to.  We'll add one node
// for the core mass, and we'll add some more nodes in the neighborhood (MinMass, MaxMass).
// The PMCInfo nodes are children of SpectrumInfo.
void AddPMCNodes(PMCSpectrumInfo* SpectrumInfo, int CoreMass, int MinMass, int MaxMass)
{
    PMCInfo* Info;
    int MassChange;
    int Mass;
    //
    // Iterate from the core mass downward.  When you reach the end, iterate
    // from the core mass (+0.1Da) upward.  (Use the 'two-way iteration' instead
    // of two loops)
    Mass = CoreMass;
    MassChange = -DECI_DALTON;
    while (1)
    {
        if (Mass < MinMass)
        {
            MassChange = DECI_DALTON;
            Mass = CoreMass + MassChange;
        }
        if (Mass > MaxMass)
        {
            break;
        }
        Info = (PMCInfo*)calloc(1, sizeof(PMCInfo));
        Info->Charge = SpectrumInfo->Charge;
        Info->ParentMass = Mass;
        if (!SpectrumInfo->Head)
        {
            SpectrumInfo->Head = Info;
        }
        else
        {
            SpectrumInfo->Tail->Next = Info;
        }
        SpectrumInfo->Tail = Info;
        Mass += MassChange;
    }
}

// Compute features for performing parent mass correction.
// Assumes that the charge state is set!
void ComputePMCFeatures(PMCSpectrumInfo* SpectrumInfo)
{
    int OffsetIndex;
    int FeatureIndex;
    int Charge;
    float PMRadius;
    float AverageConvolution;
    PMCInfo* Info;
    MSSpectrum* Spectrum;
    float MaxConvolution;
    int InfoCount;
    float Diff;
    //
    Spectrum = SpectrumInfo->Spectrum;

    // Set the spectrum's mass:
    Spectrum->ParentMass = Spectrum->MZ * SpectrumInfo->Charge - HYDROGEN_MASS * (SpectrumInfo->Charge - 1); // base mass
    Charge = min(3, SpectrumInfo->Charge);

    ////////////////////////////////////////////////////////////
    // Build PMCInfo structs for allowed masses.  We're always allowed a +1 or -1 isotope.
    // And we're allowed to move around by 0.1Da until our mass error (in PPM) becomes too large.
    PMRadius = (float)Spectrum->ParentMass;
    PMRadius *= GlobalOptions->ParentMassPPM / (float)ONE_MILLION;

    AddPMCNodes(SpectrumInfo, Spectrum->ParentMass,
        rint(Spectrum->ParentMass - PMRadius), rint(Spectrum->ParentMass + PMRadius));

    // We're always allowed a +1 and -1 shift:
    if (PMRadius < DALTON)
    {
        AddPMCNodes(SpectrumInfo, Spectrum->ParentMass - DALTON,
            rint(Spectrum->ParentMass - DALTON - PMRadius),
            rint(min(Spectrum->ParentMass - DALTON + PMRadius, Spectrum->ParentMass - PMRadius)));
        AddPMCNodes(SpectrumInfo, Spectrum->ParentMass + DALTON,
            rint(max(Spectrum->ParentMass + DALTON - PMRadius, Spectrum->ParentMass + PMRadius)),
            rint(Spectrum->ParentMass + DALTON + PMRadius));
    }
    // Ok, PMCInfo nodes have now been created.
    // Perform self-conovolution, at various parent masses.  This populates Info->Convolve, Info->Convolve2.
    // Along the way, track the *average* and *maximum* self-convolutions.
    InfoCount = 0;
    for (Info = SpectrumInfo->Head; Info; Info = Info->Next)
    {
        ConvolveMassCorrectedSpectrum(Info, SpectrumInfo);
        InfoCount++;
    }

    // Use the self-convolution info to populate the feature-vector for each PMCInfo:
    for (Info = SpectrumInfo->Head; Info; Info = Info->Next)
    {
        FeatureIndex = 0;

        // First feature is derived from the mass offset:
        if (SpectrumInfo->Charge == 1)
        {
            // Absolute Mass offset
            Info->Features[FeatureIndex++] = (float)fabs((Spectrum->ParentMass - Info->ParentMass) / (float)MASS_SCALE);
        }
        else
        {
            // Mass offset:
            Diff = (Spectrum->ParentMass - Info->ParentMass) / (float)MASS_SCALE;
            Info->Features[FeatureIndex++] = Diff;
            // Absolute mass offset:
            Info->Features[FeatureIndex++] = Diff * Diff;
        }

        ////////////////////////////////////////////////////////////////
        // Convolution features:
        // Find the average convolution for several masses:
        AverageConvolution = 0;
        MaxConvolution = 0;
        for (OffsetIndex = 0; OffsetIndex < SELF_CONVOLVE_OFFSETS; OffsetIndex++)
        {
            AverageConvolution += Info->Convolve[OffsetIndex];
            MaxConvolution = max(MaxConvolution, Info->Convolve[OffsetIndex]);
        }
//        printf("ParentMass.c::ComputePMCFeatures. Iteration on PM %d\n", Info->ParentMass);
//        printf("Average convolution (summation) %f\n", AverageConvolution);
        AverageConvolution /= (float)SELF_CONVOLVE_OFFSETS;
        AverageConvolution = max(EPSILON, AverageConvolution);
//        printf("Average convolution (final) %f\n", AverageConvolution);
        MaxConvolution = max(EPSILON, MaxConvolution);
        for (OffsetIndex = 0; OffsetIndex < SELF_CONVOLVE_OFFSETS; OffsetIndex++)
        {
            if (OffsetIndex < 4)
            {
                Info->Features[FeatureIndex++] = Info->Convolve[OffsetIndex];
                Info->Features[FeatureIndex++] = Info->Convolve[OffsetIndex] / AverageConvolution;
            }
        }
        // Convolutions of singly- and doubly-charged peaks.
        // (These features aren't computed for charge 1!)
        if (SpectrumInfo->Charge > 1)
        {
            ////////////////////////////////////////////////////////////////
            // Convolution2 features:
            // Find the average convolution for several masses:
            AverageConvolution = 0;
            MaxConvolution = 0;
            for (OffsetIndex = 0; OffsetIndex < SELF_CONVOLVE2_OFFSETS; OffsetIndex++)
            {
                AverageConvolution += Info->Convolve2[OffsetIndex];
                MaxConvolution = max(MaxConvolution, Info->Convolve2[OffsetIndex]);
            }
            AverageConvolution /= (float)SELF_CONVOLVE2_OFFSETS;
            AverageConvolution = max(EPSILON, AverageConvolution);
            MaxConvolution = max(EPSILON, MaxConvolution);
            for (OffsetIndex = 0; OffsetIndex < SELF_CONVOLVE2_OFFSETS; OffsetIndex++)
            {
                if (OffsetIndex < 3)
                {
                    Info->Features[FeatureIndex++] = Info->Convolve2[OffsetIndex];
                    Info->Features[FeatureIndex++] = Info->Convolve2[OffsetIndex] / AverageConvolution;
                }
            }
        }
		if (GlobalOptions->PhosphorylationFlag) // sam's new insertion
		{ //feature is simple sum of M-p and M-p-h20 intensity and skew
			CharacterizePhosphatePeaks(Info, SpectrumInfo, PHOSPHATE_WATER_MASS / Info->Charge, 0);
			CharacterizePhosphatePeaks(Info, SpectrumInfo, (PHOSPHATE_WATER_MASS + WATER_MASS) / Info->Charge, 1);
			Info->Features[FeatureIndex++] = Info->IntensePeakIntensity[0] + Info->IntensePeakIntensity[1];
			Info->Features[FeatureIndex++] = (float)(Info->IntensePeakSkew[0] + Info->IntensePeakSkew[1]);
			//save this information for the charge state model.
			Info->IntensePeakIntensity[2] = Info->IntensePeakIntensity[0] + Info->IntensePeakIntensity[1];
		}
    }
}

// Get Features of possible phosphorylated peak
// looking for the most intense peak within a given range
void CharacterizePhosphatePeaks(PMCInfo* Info, PMCSpectrumInfo* SpectrumInfo, int Offset, int FeatureIndex)
{
    MSSpectrum* Spectrum;
    int PeakIndex = -1;
    int MZ; // Of this PMCInfo guess at parent mass, not what is listed in the file.
    int Epsilon = rint(0.5 * DALTON);
    int SavedPeakIndex = -1;
    int Difference;
    int Skew = 0;
    float Intensity = 0;
    float TotalIntensity = 0;
    int ExpectedPeakMass;
    //
    Spectrum = SpectrumInfo->Spectrum;
    MZ = (Info->ParentMass + (Info->Charge - 1) * HYDROGEN_MASS) / Info->Charge;
    ExpectedPeakMass = MZ - Offset;
    for (PeakIndex = 0; PeakIndex < Spectrum->PeakCount; PeakIndex++)
    {
        TotalIntensity += Spectrum->Peaks[PeakIndex].Intensity;
        Difference = abs(Spectrum->Peaks[PeakIndex].Mass - ExpectedPeakMass);
        if (Difference > Epsilon)
        {
            continue;
        }
        if (Spectrum->Peaks[PeakIndex].Intensity > Intensity)
        {
            Intensity = Spectrum->Peaks[PeakIndex].Intensity;
            Skew = Difference;
            SavedPeakIndex = PeakIndex;
        }
    }
    if (SavedPeakIndex > 0)
    {
        Info->IntensePeakIndex[FeatureIndex] = SavedPeakIndex;
        Info->IntensePeakIntensity[FeatureIndex] = Intensity / TotalIntensity; //percent total
		Info->IntensePeakSkew[FeatureIndex] = Skew;
    }
	else
	{//nothing found, save zeros
        Info->IntensePeakIndex[FeatureIndex] = 0;
        Info->IntensePeakIntensity[FeatureIndex] = 0;
		Info->IntensePeakSkew[FeatureIndex] = 0;
	}
}


typedef enum {
    DebugLevelLOW,
    DebugLevelHIGH
} DebugLevel;

void DebugPrintPMCSpectrumInfo(PMCSpectrumInfo* SpectrumInfo, DebugLevel debugLevel) {
    PMCInfo* Info;
    for (Info = SpectrumInfo->Head; Info; Info = Info->Next)
    {
        printf("Charge %d ParentMass %d SVMScore %f\n",
                Info->Charge, Info->ParentMass, Info->SVMScore );

        if ( debugLevel == DebugLevelHIGH ) {
            int i;
            for(i = 0; i < 64; i++) {
                printf("%f ",Info->Features[i]);
                if ((i+1) % 8 == 0)
                    printf("\n");
            }
            printf("\n");
        }
    }
}

// Carry out parent mass correction on this spectrum.
void PerformPMC(PMCSpectrumInfo* SpectrumInfo)
{
    PMCInfo* Info;

    //
#ifdef PMC_USE_SVM
    LoadPMCSVM(0);
#else
    LoadPMCLDA(0);
#endif
    ComputePMCFeatures(SpectrumInfo);
    // If we don't have a model (yet), then give the FIRST mass the best score:
    switch (SpectrumInfo->Head->Charge)
    {
    case 1:
        if (!PMCCharge1LDA && !PMCCharge1SVM)
        {
            SpectrumInfo->BestInfo = SpectrumInfo->Head;
            SpectrumInfo->Head->SVMScore = 1.0;
            return;
        }
        break;
    case 3:
        if (!PMCCharge3LDA && !PMCCharge3SVM)
        {
            SpectrumInfo->BestInfo = SpectrumInfo->Head;
            SpectrumInfo->Head->SVMScore = 1.0;
            return;
        }
        break;
    default:
        if (!PMCCharge2LDA && !PMCCharge2SVM)
        {
            SpectrumInfo->BestInfo = SpectrumInfo->Head;
            SpectrumInfo->Head->SVMScore = 1.0;
            return;
        }
        break;
    }

    // Apply the machine learning model to each one:
    for (Info = SpectrumInfo->Head; Info; Info = Info->Next)
    {
        if (Info->Charge == 1)
        {
#ifdef PMC_USE_SVM
            Info->SVMScore = SVMClassify(PMCCharge1SVM, Info->Features, 0);
#else
            Info->SVMScore = ApplyLDAModel(PMCCharge1LDA, Info->Features);
#endif
        }
        else if (Info->Charge == 2)
        {
#ifdef PMC_USE_SVM
            Info->SVMScore = SVMClassify(PMCCharge2SVM, Info->Features, 0);
#else
            Info->SVMScore = ApplyLDAModel(PMCCharge2LDA, Info->Features);
#endif
        }
        else
        {
#ifdef PMC_USE_SVM
            Info->SVMScore = SVMClassify(PMCCharge3SVM, Info->Features, 0);
#else
            Info->SVMScore = ApplyLDAModel(PMCCharge3LDA, Info->Features);
#endif
        }
    }
//    DebugPrintPMCSpectrumInfo(SpectrumInfo,DebugLevelHIGH);
    // Remember the best one:
    for (Info = SpectrumInfo->Head; Info; Info = Info->Next)
    {
        if (!SpectrumInfo->BestInfo || Info->SVMScore > SpectrumInfo->BestInfo->SVMScore)
        {
            SpectrumInfo->BestInfo = Info;
        }
    }
    // Remember the second-best one:
    for (Info = SpectrumInfo->Head; Info; Info = Info->Next)
    {
        //if (Info == SpectrumInfo->BestInfo)
        //{
        //    continue;
        //}
        if (abs(Info->ParentMass - SpectrumInfo->BestInfo->ParentMass) <= 400)
        {
            continue;
        }

        if (!SpectrumInfo->RunnerUpInfo || (Info->SVMScore > SpectrumInfo->RunnerUpInfo->SVMScore))
        {
            SpectrumInfo->RunnerUpInfo = Info;
        }
    }
}

// Load parent mass correction SVM models.
int LoadPMCSVM()
{
    char FilePath[1024];
    if (PMCCharge1SVM)
    {
        return 1;
    }
    sprintf(FilePath, "%s%s", GlobalOptions->ResourceDir, "PMCSVM1.model");
    PMCCharge1SVM = ReadSVMModel(FilePath);
    sprintf(FilePath, "%s%s", GlobalOptions->ResourceDir, "PMCSVM1.range");
    ReadSVMScaling(PMCCharge1SVM, FilePath);
    sprintf(FilePath, "%s%s", GlobalOptions->ResourceDir, "PMCSVM2.model");
    PMCCharge2SVM = ReadSVMModel(FilePath);
    sprintf(FilePath, "%s%s", GlobalOptions->ResourceDir, "PMCSVM2.range");
    ReadSVMScaling(PMCCharge2SVM, FilePath);
    sprintf(FilePath, "%s%s", GlobalOptions->ResourceDir, "PMCSVM3.model");
    PMCCharge3SVM = ReadSVMModel(FilePath);
    sprintf(FilePath, "%s%s", GlobalOptions->ResourceDir, "PMCSVM3.range");
    ReadSVMScaling(PMCCharge3SVM, FilePath);
    return 1;
}

// Build a PMCSpectrumInfo instance for the spectrum.  Assumes the charge state is set.
PMCSpectrumInfo* GetPMCSpectrumInfo(MSSpectrum* Spectrum)
{
    PMCSpectrumInfo* SpectrumInfo;
    float SelfConvolve;
    int PeakIndex;
    int Bin;
    float Intensity;
    //
    SpectrumInfo = (PMCSpectrumInfo*)calloc(1, sizeof(PMCSpectrumInfo));
    SpectrumInfo->Spectrum = Spectrum;
    SpectrumInfo->Charge = Spectrum->Charge;
    SpectrumInfo->Mass = (Spectrum->MZ * Spectrum->Charge) - ((Spectrum->Charge - 1) * HYDROGEN_MASS);

    // Scale spectrum peaks to a TOTAL intensity of 100:
    Intensity = 0;
    for (PeakIndex = 0; PeakIndex < Spectrum->PeakCount; PeakIndex++)
    {
        Intensity += Spectrum->Peaks[PeakIndex].Intensity;
    }
    SpectrumInfo->PeakScalingFactor = 100 / Intensity;
    SpectrumInfo->PeakScalingFactor *= SpectrumInfo->PeakScalingFactor;

    // Compute self-convolution:
    SelfConvolve = EPSILON;
    for (PeakIndex = 0; PeakIndex < Spectrum->PeakCount; PeakIndex++)
    {
        Bin = (Spectrum->Peaks[PeakIndex].Mass + 50) / 100;
        if (Bin >= 0 && Bin < Spectrum->IntensityBinCount)
        {
            Intensity = Spectrum->BinnedIntensitiesTight[Bin];
            SelfConvolve += Spectrum->Peaks[PeakIndex].Intensity * Intensity * SpectrumInfo->PeakScalingFactor;
        }
    }
    SpectrumInfo->SelfConvolution = SelfConvolve;
    return SpectrumInfo;
}

float SpectrumGetSelfConvolution(MSSpectrum* Spectrum, PMCSpectrumInfo* SpectrumInfo, int Offset, int DoublyChargedFlag)
{
    SelfConvolutionNode* Node;
    SelfConvolutionNode* OldNode;

    int PeakIndex;
    int OtherMass;
    int Bin;
    float Product;
    int VerboseFlag = 0;
    float Convolution;
    int HashIndex;
    //

    HashIndex = abs(Offset / 100) % SC_HASH_SIZE;
    // If the self-convolution has already been computed, then we simply look it up:
    if (DoublyChargedFlag)
    {
        Node = SpectrumInfo->SC2Hash[HashIndex];
    }
    else
    {
        Node = SpectrumInfo->SCHash[HashIndex];
    }
    for (; Node; Node = Node->Next)
    {
        if (Node->MassOffset == Offset)
        {
            //printf("SGSC%d: Return already-computed for offset %d\n", DoublyChargedFlag, Offset);
            return Node->Value;
        }
    }

    //printf("SGSC%d: Compute for offset %d\n", DoublyChargedFlag, Offset);
    // Compute convolution value for these parameters:
    Convolution = 0;
    for (PeakIndex = 0; PeakIndex < Spectrum->PeakCount; PeakIndex++)
    {
        if (DoublyChargedFlag)
        {
            OtherMass = SpectrumInfo->Mass + 2 * HYDROGEN_MASS - (2 * Spectrum->Peaks[PeakIndex].Mass) + Offset;
        }
        else
        {
            OtherMass = SpectrumInfo->Mass + HYDROGEN_MASS - Spectrum->Peaks[PeakIndex].Mass + Offset;
        }
        Bin = ((OtherMass + 50) / 100);
        if (Bin < 0 || Bin >= Spectrum->IntensityBinCount)
        {
            continue;
        }
        Product = Spectrum->Peaks[PeakIndex].Intensity * Spectrum->BinnedIntensitiesTight[Bin] * SpectrumInfo->PeakScalingFactor;
        if (VerboseFlag && Product)
        {
            printf("Peak@%.2f and binned intensity %d (%.2f) -> %.5f\n", Spectrum->Peaks[PeakIndex].Mass / (float)DALTON,
                Bin, OtherMass / (float)DALTON, Product);
        }
        Convolution += Product;
    }
    Node = (SelfConvolutionNode*)calloc(1, sizeof(SelfConvolutionNode));
    Node->MassOffset = Offset;
    Node->Value = Convolution / SpectrumInfo->SelfConvolution;
    if (DoublyChargedFlag)
    {
        if (SpectrumInfo->SC2Hash[HashIndex])
        {
            OldNode = SpectrumInfo->SC2Hash[HashIndex];
            while (OldNode->Next)
            {
                OldNode = OldNode->Next;
            }
            OldNode->Next = Node;
        }
        else
        {
            SpectrumInfo->SC2Hash[HashIndex] = Node;
        }
    }
    else
    {
        if (SpectrumInfo->SCHash[HashIndex])
        {
            OldNode = SpectrumInfo->SCHash[HashIndex];
            while (OldNode->Next)
            {
                OldNode = OldNode->Next;
            }
            OldNode->Next = Node;
        }
        else
        {
            SpectrumInfo->SCHash[HashIndex] = Node;
        }
    }
    return Node->Value;
}


// ConvolveMassCorrectedSpectrum computes self-convolution for the given charge
// and parent mass.
void ConvolveMassCorrectedSpectrum(PMCInfo* Info, PMCSpectrumInfo* SpectrumInfo)
{
    MSSpectrum* Spectrum;
    int OffsetIndex;
    int VerboseFlag = 0;
    int OverallOffset;

    // Offsets consists of some masses where we expect LARGE self-convolution, followed by others where we expect SMALL
    // self-convolution:
    int Offsets[SELF_CONVOLVE_OFFSETS] = {-18 * DALTON, -17 * DALTON, 0 * DALTON, 1 * DALTON,
        -1 * DALTON, rint(0.5 * DALTON), rint(-16.5 * DALTON)};
    int Offsets2[SELF_CONVOLVE2_OFFSETS] = {rint(0.4 * DALTON), rint(1.2 * DALTON), rint(-17.5 * DALTON),
        -1 * DALTON, 4 * DALTON};

	if(GlobalOptions->PhosphorylationFlag)
	{//for phos searches, these offsets produce much better results.
		Offsets2[0] = rint(0.2 * DALTON);
		Offsets2[2] = rint(-18.0 * DALTON);
	}

    //
    Spectrum = SpectrumInfo->Spectrum;
    if (!Spectrum->BinnedIntensities) // move to caller!
    {
        REPORT_ERROR_S(4, "Error in ConvolveMassCorrectedSpectrum(): Spectrum binned intensities not set.\n");
        return;
    }
    for (OffsetIndex = 0; OffsetIndex < SELF_CONVOLVE_OFFSETS; OffsetIndex++)
    {
        if (VerboseFlag)
        {
            printf("\n>>>Offset %d: %.2f\n", OffsetIndex, Offsets[OffsetIndex] / (float)DALTON);
        }
        OverallOffset = Offsets[OffsetIndex] + (Info->ParentMass - SpectrumInfo->Mass);
        Info->Convolve[OffsetIndex] = SpectrumGetSelfConvolution(Spectrum, SpectrumInfo, OverallOffset, 0);

        //// Compute convolution value for these parameters:
        //Convolution = 0;
        //for (PeakIndex = 0; PeakIndex < Spectrum->PeakCount; PeakIndex++)
        //{
        //    OtherMass = Info->ParentMass + HYDROGEN_MASS - Spectrum->Peaks[PeakIndex].Mass + Offsets[OffsetIndex];
        //    Bin = ((OtherMass + 50) / 100);
        //    if (Bin < 0 || Bin >= Spectrum->IntensityBinCount)
        //    {
        //        continue;
        //    }
        //    Product = Spectrum->Peaks[PeakIndex].Intensity * Spectrum->BinnedIntensitiesTight[Bin] * SpectrumInfo->PeakScalingFactor;
        //    if (VerboseFlag && Product)
        //    {
        //        printf("Peak@%.2f and binned intensity %d (%.2f) -> %.5f\n", Spectrum->Peaks[PeakIndex].Mass / (float)DALTON,
        //            Bin, OtherMass / (float)DALTON, Product);
        //    }
        //    Convolution += Product;
        //}
        //Info->Convolve[OffsetIndex] = Convolution / SpectrumInfo->SelfConvolution;
        //if (VerboseFlag)
        //{
        //    printf(">>Convolve[%d] = %.4f\n", OffsetIndex, Convolution);
        //}
    }

    if (Spectrum->Charge > 1)//I compute these values for phos charge 2, but don't use them for the PMC model.
    {                        //they do go into the ChargeCorrection model, so they are still calculated.
        // Compute convolution of charge-1 and charge-2 peaks:
        for (OffsetIndex = 0; OffsetIndex < SELF_CONVOLVE2_OFFSETS; OffsetIndex++)
        {
            OverallOffset = Offsets2[OffsetIndex] + (Info->ParentMass - SpectrumInfo->Mass);
            Info->Convolve2[OffsetIndex] = SpectrumGetSelfConvolution(Spectrum, SpectrumInfo, OverallOffset, 1);
            //Convolution = 0;
            //for (PeakIndex = 0; PeakIndex < Spectrum->PeakCount; PeakIndex++)
            //{
            //    OtherMass = Info->ParentMass + 2 * HYDROGEN_MASS - (2 * Spectrum->Peaks[PeakIndex].Mass) + Offsets2[OffsetIndex];
            //    Bin = ((OtherMass + 50) / 100);
            //    if (Bin < 0 || Bin >= Spectrum->IntensityBinCount)
            //    {
            //        continue;
            //    }
            //    Convolution += Spectrum->Peaks[PeakIndex].Intensity * Spectrum->BinnedIntensitiesTight[Bin] * SpectrumInfo->PeakScalingFactor;
            //}
            //Info->Convolve2[OffsetIndex] = Convolution / SpectrumInfo->SelfConvolution;
        }
    }
}
