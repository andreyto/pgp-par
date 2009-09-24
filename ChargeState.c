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
#include "IonScoring.h"

#ifdef _WIN32
#include <Windows.h>
#include <sys/types.h>
#include <sys/stat.h>
#else
#include <dirent.h>
#include <sys/stat.h>
#endif

#define CC_USE_SVM

#define EPSILON (float)0.00001

SVMModel** PMCModel = NULL;

extern LDAModel* PMCCharge1LDA;
extern LDAModel* PMCCharge2LDA;
extern LDAModel* PMCCharge3LDA;

extern SVMModel* PMCCharge1SVM;
extern SVMModel* PMCCharge2SVM;
extern SVMModel* PMCCharge3SVM;

extern LDAModel* CCModel1LDA;
extern LDAModel* CCModel2LDA;

extern SVMModel* CCModel1SVM;
extern SVMModel* CCModel2SVM;

extern PRMBayesianModel* PRMModelCharge2;

// For converting parts-per-million:
#define ONE_MILLION 1000000

///////////////////////////////////////////////////
// Forward declarations:
void ConvolveMassCorrectedSpectrum(PMCInfo* Info, PMCSpectrumInfo* SpectrumInfo);

///////////////////////////////////////////////////
// Functions:

// Get charge correction features.  Most of the charge correction features are set
// during parent mass correction - if BY convolution assuming charge 2 is very high,
// then it's most probable that the true spectrum charge is 2.
void GetChargeCorrectionFeatures1(PMCSpectrumInfo* SpectrumInfo1, PMCSpectrumInfo* SpectrumInfo2,
    PMCSpectrumInfo* SpectrumInfo3, float* Features)
{
    float TotalIntensity = 0;
    float LowIntensity = 0; // Below m/z
    float MediumIntensity = 0; // Between m/z and 2*m/z
    float HighIntensity = 0; // Above 2*m/z
    int LowPeakCount = 0;
    int MediumPeakCount = 0;
    int HighPeakCount = 0;
    int PeakIndex;
    int FeatureIndex = 0;
    float Competitor;
    int MZ;

    MSSpectrum* Spectrum = SpectrumInfo1->Spectrum;
    MZ = SpectrumInfo1->BestInfo->ParentMass;
    //
    for (PeakIndex = 0; PeakIndex < Spectrum->PeakCount; PeakIndex++)
    {
        TotalIntensity += Spectrum->Peaks[PeakIndex].Intensity;
        if (Spectrum->Peaks[PeakIndex].Mass <= MZ)
        {
            LowPeakCount++;
            LowIntensity += Spectrum->Peaks[PeakIndex].Intensity;
        }
        else if (Spectrum->Peaks[PeakIndex].Mass <= 2 * MZ)
        {
            MediumPeakCount++;
            MediumIntensity += Spectrum->Peaks[PeakIndex].Intensity;
        }
        else 
        {
            HighPeakCount++;
            HighIntensity += Spectrum->Peaks[PeakIndex].Intensity;
        }
    }

    // Feature: How much of the spectral intensity is above M/z?  
    Features[FeatureIndex++] = (MediumIntensity + HighIntensity) / (float)max(0.001, TotalIntensity);
    Features[FeatureIndex++] = (MediumPeakCount + HighPeakCount) / (float)Spectrum->PeakCount;
        
    // Features: How do the B/Y convolution values compare between charges 1 and 2?
    Competitor = max(SpectrumInfo2->BestInfo->Convolve[0], SpectrumInfo3->BestInfo->Convolve[0]);
    Features[FeatureIndex++] = SpectrumInfo1->BestInfo->Convolve[0] / max(EPSILON, SpectrumInfo1->BestInfo->Convolve[0] + Competitor);
    Features[FeatureIndex++] = SpectrumInfo1->BestInfo->Convolve[0] - Competitor;
    Competitor = max(SpectrumInfo2->BestInfo->Convolve[1], SpectrumInfo3->BestInfo->Convolve[1]);
    Features[FeatureIndex++] = SpectrumInfo1->BestInfo->Convolve[1] / max(EPSILON, SpectrumInfo1->BestInfo->Convolve[1] + Competitor);
    Features[FeatureIndex++] = SpectrumInfo1->BestInfo->Convolve[1] - Competitor;
    Competitor = max(SpectrumInfo2->BestInfo->Convolve[2], SpectrumInfo3->BestInfo->Convolve[2]);
    Features[FeatureIndex++] = SpectrumInfo1->BestInfo->Convolve[2] / max(EPSILON, SpectrumInfo1->BestInfo->Convolve[2] + Competitor);
    Features[FeatureIndex++] = SpectrumInfo1->BestInfo->Convolve[2] - Competitor;
    Competitor = max(SpectrumInfo2->BestInfo->Convolve[3], SpectrumInfo3->BestInfo->Convolve[3]);
    Features[FeatureIndex++] = SpectrumInfo1->BestInfo->Convolve[3] / max(EPSILON, SpectrumInfo1->BestInfo->Convolve[3] + Competitor);
    Features[FeatureIndex++] = SpectrumInfo1->BestInfo->Convolve[3] - Competitor;
}

// Get charge correction features.  Most of the charge correction features are set
// during parent mass correction - if BY convolution assuming charge 2 is very high,
// then it's most probable that the true spectrum charge is 2.
void GetChargeCorrectionFeatures2(PMCSpectrumInfo* SpectrumInfo2, PMCSpectrumInfo* SpectrumInfo3,
    float* Features)
{
    float TotalIntensity = 0;
    float MediumIntensity = 0;
    float HighIntensity = 0;
    float LowIntensity = 0;
    int LowPeakCount = 0;
    int MediumPeakCount = 0;
    int HighPeakCount = 0;
    int PeakIndex;
    int FeatureIndex = 0;
    float MZ;
    float Balance2;
    MSSpectrum* Spectrum = SpectrumInfo2->Spectrum;
    //
    MZ = SpectrumInfo2->BestInfo->ParentMass / (float)2.0;
    //
    for (PeakIndex = 0; PeakIndex < Spectrum->PeakCount; PeakIndex++)
    {
        TotalIntensity += Spectrum->Peaks[PeakIndex].Intensity;
        if (Spectrum->Peaks[PeakIndex].Mass <= MZ)
        {
            LowPeakCount++;
            LowIntensity += Spectrum->Peaks[PeakIndex].Intensity;
        }
        else if (Spectrum->Peaks[PeakIndex].Mass <= 2 * MZ)
        {
            MediumPeakCount++;
            MediumIntensity += Spectrum->Peaks[PeakIndex].Intensity;
        }
        else 
        {
            HighPeakCount++;
            HighIntensity += Spectrum->Peaks[PeakIndex].Intensity;
        }
    }

    // Feature: How much of the spectral intensity is above M/z?  
    Features[FeatureIndex++] = (MediumIntensity + HighIntensity) / TotalIntensity;
    Features[FeatureIndex++] = (MediumPeakCount + HighPeakCount) / (float)Spectrum->PeakCount;
    Features[FeatureIndex++] = (MediumIntensity) / TotalIntensity;
    Features[FeatureIndex++] = (MediumPeakCount) / (float)Spectrum->PeakCount;
    Features[FeatureIndex++] = (LowIntensity) / TotalIntensity;
    Features[FeatureIndex++] = (LowPeakCount) / (float)Spectrum->PeakCount;
    //Features[FeatureIndex++] = (HighIntensity) / TotalIntensity;
    //Features[FeatureIndex++] = (HighPeakCount) / (float)Spectrum->PeakCount;

    // Features: Balance between low and med-to-high:
    Balance2 = (float)fabs((MediumIntensity + HighIntensity) - LowIntensity) / TotalIntensity;
    Features[FeatureIndex++] = Balance2;
    
    // Features: How do the B/Y convolution values compare between charges 2 and 3?
    Features[FeatureIndex++] = SpectrumInfo2->BestInfo->Convolve[0] / max(EPSILON, SpectrumInfo2->BestInfo->Convolve[0] + SpectrumInfo3->BestInfo->Convolve[0]);
    Features[FeatureIndex++] = SpectrumInfo2->BestInfo->Convolve[0] - SpectrumInfo3->BestInfo->Convolve[0];
    Features[FeatureIndex++] = SpectrumInfo2->BestInfo->Convolve[1] / max(EPSILON, SpectrumInfo2->BestInfo->Convolve[1] + SpectrumInfo3->BestInfo->Convolve[1]);
    Features[FeatureIndex++] = SpectrumInfo2->BestInfo->Convolve[1] - SpectrumInfo3->BestInfo->Convolve[1];
    Features[FeatureIndex++] = SpectrumInfo2->BestInfo->Convolve[2] / max(EPSILON, SpectrumInfo2->BestInfo->Convolve[2] + SpectrumInfo3->BestInfo->Convolve[2]);
    Features[FeatureIndex++] = SpectrumInfo2->BestInfo->Convolve[2] - SpectrumInfo3->BestInfo->Convolve[2];
    Features[FeatureIndex++] = SpectrumInfo2->BestInfo->Convolve[3] / max(EPSILON, SpectrumInfo2->BestInfo->Convolve[3] + SpectrumInfo3->BestInfo->Convolve[3]);
    Features[FeatureIndex++] = SpectrumInfo2->BestInfo->Convolve[3] - SpectrumInfo3->BestInfo->Convolve[3];
    Features[FeatureIndex++] = SpectrumInfo2->BestInfo->Convolve2[0] / max(EPSILON, SpectrumInfo2->BestInfo->Convolve2[0] + SpectrumInfo3->BestInfo->Convolve2[0]);
    Features[FeatureIndex++] = SpectrumInfo2->BestInfo->Convolve2[0] - SpectrumInfo3->BestInfo->Convolve2[0];
    Features[FeatureIndex++] = SpectrumInfo2->BestInfo->Convolve2[1] / max(EPSILON, SpectrumInfo2->BestInfo->Convolve2[1] + SpectrumInfo3->BestInfo->Convolve2[1]);
    Features[FeatureIndex++] = SpectrumInfo2->BestInfo->Convolve2[1] - SpectrumInfo3->BestInfo->Convolve2[1];
    Features[FeatureIndex++] = SpectrumInfo2->BestInfo->Convolve2[2] / max(EPSILON, SpectrumInfo2->BestInfo->Convolve2[2] + SpectrumInfo3->BestInfo->Convolve2[2]);
    Features[FeatureIndex++] = SpectrumInfo2->BestInfo->Convolve2[2] - SpectrumInfo3->BestInfo->Convolve2[2];
    //Features[FeatureIndex++] = Spectrum->PeakCount;
    Features[FeatureIndex++] = (SpectrumInfo2->BestInfo->ParentMass / (float)(1000 * DALTON));
}


//Phosphorylation uses a distinct PMC model, so that means it needs a distinct CC model
//most notably, we are going to use the IntensePeakIntensity and skew
void GetChargeCorrectionFeatures2Phos(PMCSpectrumInfo* SpectrumInfo2, PMCSpectrumInfo* SpectrumInfo3,
    float* Features)
{
    float TotalIntensity = 0;
    float MediumIntensity = 0;
    float HighIntensity = 0;
    float LowIntensity = 0;
    int LowPeakCount = 0;
    int MediumPeakCount = 0;
    int HighPeakCount = 0;
    int PeakIndex;
    int FeatureIndex = 0;
    float MZ;
    float Balance2;
	float PhosPeak2;
	float PhosPeak3;
    MSSpectrum* Spectrum = SpectrumInfo2->Spectrum;
    //
    MZ = SpectrumInfo2->BestInfo->ParentMass / (float)2.0;
    //
	for (PeakIndex = 0; PeakIndex < Spectrum->PeakCount; PeakIndex++)
    {
        TotalIntensity += Spectrum->Peaks[PeakIndex].Intensity;
        if (Spectrum->Peaks[PeakIndex].Mass <= MZ)
        {
            LowPeakCount++;
            LowIntensity += Spectrum->Peaks[PeakIndex].Intensity;
        }
        else if (Spectrum->Peaks[PeakIndex].Mass <= 2 * MZ)
        {
            MediumPeakCount++;
            MediumIntensity += Spectrum->Peaks[PeakIndex].Intensity;
        }
        else 
        {
            HighPeakCount++;
            HighIntensity += Spectrum->Peaks[PeakIndex].Intensity;
        }
    }

    // Feature: How much of the spectral intensity is above M/z?  
    Features[FeatureIndex++] = (MediumIntensity + HighIntensity) / TotalIntensity;
    Features[FeatureIndex++] = (MediumPeakCount + HighPeakCount) / (float)Spectrum->PeakCount;
    Features[FeatureIndex++] = (MediumIntensity) / TotalIntensity;
    Features[FeatureIndex++] = (MediumPeakCount) / (float)Spectrum->PeakCount;
    Features[FeatureIndex++] = (LowIntensity) / TotalIntensity;
    Features[FeatureIndex++] = (LowPeakCount) / (float)Spectrum->PeakCount;
    //Features[FeatureIndex++] = (HighIntensity) / TotalIntensity;
    //Features[FeatureIndex++] = (HighPeakCount) / (float)Spectrum->PeakCount;

    // Features: Balance between low and med-to-high:
    Balance2 = (float)fabs((MediumIntensity + HighIntensity) - LowIntensity) / TotalIntensity;
    Features[FeatureIndex++] = Balance2;
    
    // Features: How do the B/Y convolution values compare between charges 2 and 3?
    Features[FeatureIndex++] = SpectrumInfo2->BestInfo->Convolve[0] / max(EPSILON, SpectrumInfo2->BestInfo->Convolve[0] + SpectrumInfo3->BestInfo->Convolve[0]);
    Features[FeatureIndex++] = SpectrumInfo2->BestInfo->Convolve[0] - SpectrumInfo3->BestInfo->Convolve[0];
    Features[FeatureIndex++] = SpectrumInfo2->BestInfo->Convolve[1] / max(EPSILON, SpectrumInfo2->BestInfo->Convolve[1] + SpectrumInfo3->BestInfo->Convolve[1]);
    Features[FeatureIndex++] = SpectrumInfo2->BestInfo->Convolve[1] - SpectrumInfo3->BestInfo->Convolve[1];
    Features[FeatureIndex++] = SpectrumInfo2->BestInfo->Convolve[2] / max(EPSILON, SpectrumInfo2->BestInfo->Convolve[2] + SpectrumInfo3->BestInfo->Convolve[2]);
    Features[FeatureIndex++] = SpectrumInfo2->BestInfo->Convolve[2] - SpectrumInfo3->BestInfo->Convolve[2];
    Features[FeatureIndex++] = SpectrumInfo2->BestInfo->Convolve[3] / max(EPSILON, SpectrumInfo2->BestInfo->Convolve[3] + SpectrumInfo3->BestInfo->Convolve[3]);
    Features[FeatureIndex++] = SpectrumInfo2->BestInfo->Convolve[3] - SpectrumInfo3->BestInfo->Convolve[3];
    Features[FeatureIndex++] = SpectrumInfo2->BestInfo->Convolve2[0] / max(EPSILON, SpectrumInfo2->BestInfo->Convolve2[0] + SpectrumInfo3->BestInfo->Convolve2[0]);
    Features[FeatureIndex++] = SpectrumInfo2->BestInfo->Convolve2[0] - SpectrumInfo3->BestInfo->Convolve2[0];
    Features[FeatureIndex++] = SpectrumInfo2->BestInfo->Convolve2[1] / max(EPSILON, SpectrumInfo2->BestInfo->Convolve2[1] + SpectrumInfo3->BestInfo->Convolve2[1]);
    Features[FeatureIndex++] = SpectrumInfo2->BestInfo->Convolve2[1] - SpectrumInfo3->BestInfo->Convolve2[1];
    Features[FeatureIndex++] = SpectrumInfo2->BestInfo->Convolve2[2] / max(EPSILON, SpectrumInfo2->BestInfo->Convolve2[2] + SpectrumInfo3->BestInfo->Convolve2[2]);
    Features[FeatureIndex++] = SpectrumInfo2->BestInfo->Convolve2[2] - SpectrumInfo3->BestInfo->Convolve2[2];
    //Features[FeatureIndex++] = Spectrum->PeakCount;
    Features[FeatureIndex++] = (SpectrumInfo2->BestInfo->ParentMass / (float)(1000 * DALTON));
	//M-p peak related stuff
	PhosPeak2 = (float) (max(0.1, SpectrumInfo2->BestInfo->IntensePeakIntensity[2]));
	PhosPeak3 = (float) (max(0.1, SpectrumInfo3->BestInfo->IntensePeakIntensity[2]));
	Features[FeatureIndex++] = PhosPeak2 / (PhosPeak2 + PhosPeak3);
}


int ChargeCorrectSpectrum(SpectrumNode* Node, float* Model1Score, float* Model2Score)
{
    PMCSpectrumInfo* SpectrumInfo1;
    PMCSpectrumInfo* SpectrumInfo2;
    PMCSpectrumInfo* SpectrumInfo3;
    float CCFeatures1[64];
    float CCFeatures2[64];
    float Score1;
    float Score2;
    //
    Score1 = 0;
#ifdef CC_USE_SVM
    LoadCCModelSVM(0);
#else
    LoadCCModelLDA(0);
#endif
    /////////////////////////////////
    // Charge 1 PMC:
    Node->Spectrum->Charge = 1;
    Node->Spectrum->ParentMass = (Node->Spectrum->MZ * 1);
    SpectrumInfo1 = GetPMCSpectrumInfo(Node->Spectrum);
    PerformPMC(SpectrumInfo1);
    /////////////////////////////////
    // Charge 2 PMC:
    Node->Spectrum->Charge = 2;
    Node->Spectrum->ParentMass = (Node->Spectrum->MZ * 2) - HYDROGEN_MASS;
    SpectrumInfo2 = GetPMCSpectrumInfo(Node->Spectrum);
    PerformPMC(SpectrumInfo2);
    /////////////////////////////////
    // Charge 3 PMC:
    Node->Spectrum->Charge = 3;
    Node->Spectrum->ParentMass = (Node->Spectrum->MZ * 3) - 2 * HYDROGEN_MASS;
    SpectrumInfo3 = GetPMCSpectrumInfo(Node->Spectrum);
    PerformPMC(SpectrumInfo3);
    // Get features:
    memset(CCFeatures1, 0, sizeof(float) * 64);
    memset(CCFeatures2, 0, sizeof(float) * 64);
    GetChargeCorrectionFeatures1(SpectrumInfo1, SpectrumInfo2, SpectrumInfo3, CCFeatures1);
    GetChargeCorrectionFeatures2(SpectrumInfo2, SpectrumInfo3, CCFeatures2); //change to Phos function if you need
#ifdef CC_USE_SVM
    Score1 = SVMClassify(CCModel1SVM, CCFeatures1, 0);
    Score2 = SVMClassify(CCModel2SVM, CCFeatures2, 0);
#else
    Score1 = ApplyLDAModel(CCModel1LDA, CCFeatures1);
    Score2 = ApplyLDAModel(CCModel2LDA, CCFeatures2);
#endif
    // If the caller asked for them, return the scores from the two models:
    if (Model1Score)
    {
        *Model1Score = Score1;
    }
    if (Model2Score)
    {
        *Model2Score = Score2;
    }
    // Free temporary structs:
    FreePMCSpectrumInfo(SpectrumInfo1);
    FreePMCSpectrumInfo(SpectrumInfo2);
    FreePMCSpectrumInfo(SpectrumInfo3);
    // Use cutoffs to determine the favorite charge state:
    if (Score1 > 1.0)
    {
        return 1;
    }
    if (Score2 > 0.0)
    {
        return 2;
    }
    return 3;
}

// We've loaded a spectrum.  Now let's adjust its parent mass and its charge to the 
// best possible. 
void TweakSpectrum(SpectrumNode* Node)
{
    MSSpectrum* Spectrum;
    PMCSpectrumInfo* SpectrumInfo;
    PMCSpectrumInfo* SpectrumInfo1;
    PMCSpectrumInfo* SpectrumInfo2;
    PMCSpectrumInfo* SpectrumInfo3;
    float CCFeatures[64];
    float CCScore;
    int TweakIndex;
	int Charge;
    //
    if (!Node->Spectrum || !Node->Spectrum->PeakCount)
    {
        return;
    }
    Spectrum = Node->Spectrum;
    // If our models aren't loaded - which should NEVER happen in production - then we'll
    // trust the input mass and charge.
    if (!PRMModelCharge2)
    {
        if (!Spectrum->Charge)
        {
            Spectrum->Charge = 2;
            Spectrum->ParentMass = (Spectrum->MZ * 2) - HYDROGEN_MASS;
        }
        TweakIndex = (Spectrum->Charge - 1) * 2;
        Node->Tweaks[TweakIndex].Charge = Spectrum->Charge;
        Node->Tweaks[TweakIndex].ParentMass = Spectrum->ParentMass;
        return;
    }
    Node->Spectrum->ParentMass = (Spectrum->MZ * 2) - HYDROGEN_MASS;
    PrepareSpectrumForIonScoring(PRMModelCharge2, Node->Spectrum, 0);
    //SpectrumComputeBinnedIntensities(Node);
    
    if (!GlobalOptions->MultiChargeMode && Spectrum->FileChargeFlag)
    {
        // The spectrum has charge(s) assigned, and we're trusting the charge(s).
		for (Charge = 1; Charge < 5; Charge++)
		{
			if (Spectrum->FileCharge[Charge])
			{
				Spectrum->Charge = Charge;
				SpectrumInfo = GetPMCSpectrumInfo(Spectrum);
				PerformPMC(SpectrumInfo);
				TweakIndex = min(3, Spectrum->Charge - 1) * 2;
				Node->Tweaks[TweakIndex].Charge = Spectrum->Charge;
				Node->Tweaks[TweakIndex].ParentMass = SpectrumInfo->BestInfo->ParentMass;
				if (SpectrumInfo->RunnerUpInfo)
				{
					Node->Tweaks[TweakIndex + 1].Charge = Spectrum->Charge;
					Node->Tweaks[TweakIndex + 1].ParentMass = SpectrumInfo->RunnerUpInfo->ParentMass;
				}
				//SpectrumComputeNoiseDistributions(Node);
				FreePMCSpectrumInfo(SpectrumInfo);
			}
		}
        return;
    }
#ifdef CC_USE_SVM
    LoadCCModelSVM(0);
#else
    LoadCCModelLDA(0);
#endif
    
    // Either the spectrum has no charge set, or we're overriding the file guess 
    // with our charge correction guess.

    // Find the best parent mass if the charge is 1:
    Node->Spectrum->Charge = 1;
    SpectrumInfo1 = GetPMCSpectrumInfo(Spectrum);
    PerformPMC(SpectrumInfo1);
    Node->Tweaks[0].Charge = 1;
    Node->Tweaks[0].ParentMass = SpectrumInfo1->BestInfo->ParentMass;
    if (SpectrumInfo1->RunnerUpInfo)
    {
        Node->Tweaks[1].Charge = 1;
        Node->Tweaks[1].ParentMass = SpectrumInfo1->RunnerUpInfo->ParentMass;
    }

    // Find the best parent mass if the charge is 2:
    Node->Spectrum->Charge = 2;
    SpectrumInfo2 = GetPMCSpectrumInfo(Spectrum);
    PerformPMC(SpectrumInfo2);
    Node->Tweaks[2].Charge = 2;
    Node->Tweaks[2].ParentMass = SpectrumInfo2->BestInfo->ParentMass;
    if (SpectrumInfo2->RunnerUpInfo)
    {
        Node->Tweaks[3].Charge = 2;
        Node->Tweaks[3].ParentMass = SpectrumInfo2->RunnerUpInfo->ParentMass;
    }

    // Find the best parent mass if the charge is 3:
    Node->Spectrum->Charge = 3;
    SpectrumInfo3 = GetPMCSpectrumInfo(Spectrum);
    PerformPMC(SpectrumInfo3);
    Node->Tweaks[4].Charge = 3;
    Node->Tweaks[4].ParentMass = SpectrumInfo3->BestInfo->ParentMass;
    if (SpectrumInfo3->RunnerUpInfo)
    {
        Node->Tweaks[5].Charge = 3;
        Node->Tweaks[5].ParentMass = SpectrumInfo3->RunnerUpInfo->ParentMass;
    }

    GetChargeCorrectionFeatures1(SpectrumInfo1, SpectrumInfo2, SpectrumInfo3, CCFeatures);
    CCScore = SVMClassify(CCModel1SVM, CCFeatures, 0);
    if (CCScore > 0)
    {
        // It's a singly-charged spectrum:
        Node->Tweaks[2].Charge = 0;
        Node->Tweaks[3].Charge = 0;
        Node->Tweaks[4].Charge = 0;
        Node->Tweaks[5].Charge = 0;
    }
    else
    {
        // It's a multiply-charged spectrum:
        Node->Tweaks[0].Charge = 0;
        Node->Tweaks[1].Charge = 0;
		if (GlobalOptions->PhosphorylationFlag)
		{
			GetChargeCorrectionFeatures2Phos(SpectrumInfo2, SpectrumInfo3, CCFeatures);
		}
		else
		{
			GetChargeCorrectionFeatures2(SpectrumInfo2, SpectrumInfo3, CCFeatures);
		}
        CCScore = SVMClassify(CCModel2SVM, CCFeatures, 0);
        if (CCScore >= 0.5)
        {
            // It's clearly not charge-3:
            Node->Tweaks[4].Charge = 0;
            Node->Tweaks[5].Charge = 0;
        }
        if (CCScore <= -0.5)
        {
            // It's clearly not charge-2:
            Node->Tweaks[2].Charge = 0;
            Node->Tweaks[3].Charge = 0;
        }
    }
    
    // cleanup:
    FreePMCSpectrumInfo(SpectrumInfo1);
    FreePMCSpectrumInfo(SpectrumInfo2);
    FreePMCSpectrumInfo(SpectrumInfo3);
    //SpectrumComputeNoiseDistributions(Node);
    return;
}

// Iterate over lines of a training/testing oracle file, and invoke the callback function once for each.
// Line format: Tab-delimited.  Pieces are:
// Spectrum file name (not full path), charge, parent mass, annotation
void TrainOnOracleFile(char* OracleFileName, char* SpectrumDir, TrainingCallback Callback)
{
    int BytesToRead;
    char LineBuffer[MAX_LINE_LENGTH];
    int BufferPos = 0;
    int BytesRead;
    int BufferEnd = 0;
    int LineNumber = 0;
    char TextBuffer[BUFFER_SIZE * 2];
    FILE* OracleFile;
    char* SpectrumFileName;
    char FilePath[2048];
    int Charge;
    int ParentMass;
    FILE* DTAFile;
    char* Field;
    Peptide* Match;
    SpectrumNode* Node;
    InputFileNode* FNode;
    char* ColonPos;
    int SpectrumFilePos;
    char* Extension;
    //
    OracleFile = fopen(OracleFileName, "rb");
    if (!OracleFile)
    {
        printf("** Error: Unable to open training oracle '%s'.\n", OracleFileName);
        return;
    }
    Node = (SpectrumNode*)calloc(1, sizeof(SpectrumNode));
    FNode = (InputFileNode*)calloc(1, sizeof(InputFileNode));
    LineNumber = 0;
    while (1)
    {
        BytesToRead = BUFFER_SIZE - BufferEnd;
        BytesRead = ReadBinary(TextBuffer + BufferEnd, sizeof(char), BytesToRead, OracleFile);
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

        // Now, move the remaining text to the start of the buffer:
        memmove(TextBuffer, TextBuffer + BufferPos, BufferEnd - BufferPos);
        BufferEnd -= BufferPos;
        BufferPos = 0;

        // Now, process this line of text!
        // Skip empty lines:
        if (!LineBuffer[0])
        {
            continue;
        }
        if (LineBuffer[0] == '#')
        {
            continue;
        }
        SpectrumFileName = strtok(LineBuffer, "\t");
        if (!SpectrumFileName)
        {
            continue;
        }

        SpectrumFilePos = 0;
        ColonPos = SpectrumFileName;
        if (SpectrumFileName[1] == ':')
        {
            ColonPos = SpectrumFileName + 2;
        }
        while (*ColonPos)
        {
            if (*ColonPos == ':')
            {
                *ColonPos = '\0';
                SpectrumFilePos = atoi(ColonPos + 1);
                break;
            }
            ColonPos++;
        }
        
        Extension = SpectrumFileName + strlen(SpectrumFileName) - 4;
        if (!CompareStrings(Extension, ".mgf"))
        {
            FNode->Format = SPECTRUM_FORMAT_MGF;
        }
        else if (!CompareStrings(Extension, ".ms2"))
        {
            //FNode->Format = SPECTRUM_FORMAT_MS2;
            FNode->Format = SPECTRUM_FORMAT_MS2_COLONS;
        }
        else if (!CompareStrings(Extension, ".mzxml"))
        {
            FNode->Format = SPECTRUM_FORMAT_MZXML;
        }
        else if (!CompareStrings(Extension, ".mzdata"))
        {
            FNode->Format = SPECTRUM_FORMAT_MZDATA;
        }
        else
        {
            FNode->Format = SPECTRUM_FORMAT_DTA;
        }
        if (SpectrumFileName[1] == ':')
        {
            sprintf(FilePath, "%s", SpectrumFileName);
        }
        else
        {
            sprintf(FilePath, "%s%s", SpectrumDir, SpectrumFileName);
        }
        
        
        DTAFile = fopen(FilePath, "rb");
        if (!DTAFile)
        {
            printf("**Error: Couldn't open training/testing spectrum '%s'\n", FilePath);
            continue;
        }
        fseek(DTAFile, SpectrumFilePos, 0);
        Field = strtok(NULL, "\t");
        if (!Field)
        {
            printf("** Syntax error: Line %d of %s\n", LineNumber, OracleFileName);
            continue;
        }
        Charge = atoi(Field);
        Field = strtok(NULL, "\t");
        if (!Field)
        {
            printf("** Syntax error: Line %d of %s\n", LineNumber, OracleFileName);
            continue;
        }

        ParentMass = (int)(atof(Field) * MASS_SCALE + 0.5);
        Field = strtok(NULL, "\t");
        if (!Field)
        {
            printf("** Syntax error: Line %d of %s\n", LineNumber, OracleFileName);
            continue;
        }

        Match = GetPeptideFromAnnotation(Field);
        Node->Spectrum = (MSSpectrum*)calloc(1, sizeof(MSSpectrum));
        Node->Spectrum->Node = Node;
        strcpy(FNode->FileName, FilePath);
        Node->InputFile = FNode;
        //strcpy(Node->FileName, FilePath);
        SpectrumLoadFromFile(Node->Spectrum, DTAFile);
        fclose(DTAFile);
        (*Callback)(Node, Charge, ParentMass, Match);
        FreeSpectrum(Node->Spectrum);
        FreePeptideNode(Match);
    }
}

