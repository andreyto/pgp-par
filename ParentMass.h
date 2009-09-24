#ifndef PARENT_MASS_H
#define PARENT_MASS_H

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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "Utils.h"
#include "Inspect.h"
#include "Spectrum.h"

#define SELF_CONVOLVE_OFFSETS 7
#define SELF_CONVOLVE2_OFFSETS 5

// A linked list of self-convolutions for a spectrum.  
// We keep this list in the PMCSpectrumInfo, because many
// of the PMCInfo objects will re-use the same self-convolutions;
// it's expensive to re-compute them.
typedef struct SelfConvolutionNode
{
    int MassOffset;
    float Value;
    struct SelfConvolutionNode* Next;
} SelfConvolutionNode;

#define SC_HASH_SIZE 64
// PMCSpectrumInfo is data used during parent mass correction; here we store
// intermediate values which are general across the whole spectrum, so that
// we needn't re-compute them for each PMCInfo
typedef struct PMCSpectrumInfo
{ 
    MSSpectrum* Spectrum;
    int Charge;
    float PeakScalingFactor;
    float SelfConvolution;
    struct PMCInfo* Head;
    struct PMCInfo* Tail;
    struct PMCInfo* BestInfo;
    struct PMCInfo* RunnerUpInfo;
    int Mass; // base mass, from the file
    SelfConvolutionNode* SCHash[SC_HASH_SIZE];
    //SelfConvolutionNode* SCTail;
    //SelfConvolutionNode* SC2Head;
    SelfConvolutionNode* SC2Hash[SC_HASH_SIZE];
} PMCSpectrumInfo;

// We allocate one PMCInfo struct for each candidate parent mass.  We store 
// the SVM features here, along with the mass and other bookkeeping info.
// The PMCInfo structs are kept in a list.  The final tweaks that we keep are 
// the best 1..3 PMCInfos
typedef struct PMCInfo
{
    int Charge;
    int ParentMass;
    
    float Features[64];
    float Convolve[SELF_CONVOLVE_OFFSETS];
    float Convolve2[SELF_CONVOLVE_OFFSETS];
    float SVMScore;
    struct PMCInfo* Next;
    int IntensePeakIndex[6];//these are for keeping track of possible M-p related peaks, which are superintense
    float IntensePeakIntensity[6]; //ratio to mean
    int IntensePeakSkew[6]; //0 = M-p, 1 = m-p-h2o, 2 = feature used
} PMCInfo;

void PerformPMC(PMCSpectrumInfo* SpectrumInfo);
void FreePMCSpectrumInfo(PMCSpectrumInfo* SpectrumInfo);
void ComputePMCFeatures(PMCSpectrumInfo* SpectrumInfo);
PMCSpectrumInfo* GetPMCSpectrumInfo(MSSpectrum* Spectrum);
void ConvolveMassCorrectedSpectrum(PMCInfo* Info, PMCSpectrumInfo* SpectrumInfo);

#endif // PARENT_MASS_H
