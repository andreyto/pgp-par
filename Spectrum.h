
//Copyright 2007, The Regents of the University of California
//All Rights Reserved
//
//Permission to use, copy, modify and distribute any part of this 
//program for educational, research and non-profit purposes, without fee, 
//and without a written agreement is hereby granted, provided that the 
//above copyright notice, this paragraph and the following three paragraphs 
//appear in all copies.
//
//Those desiring to incorporate this <name of work> into commercial 
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

#ifndef SPECTRUM_H
#define SPECTRUM_H

// The basic spectrum object, with array of peaks.
// Structs and functions to support loading of spectra from several
// text-based file formats (.dta files, .mgf files, .ms2 files).

#include <stdio.h>
#include "Inspect.h"

#define DEFAULT_WINDOW_WIDTH 50000 //50Da
#define DEFAULT_WINDOW_KEEP_COUNT 6 // 6

// After filtering, there will probably be just 1 possible isotope neighbor
#define MAX_ISOTOPE_NEIGHBORS 8

#define MAX_NOISY_NEIGHBORS 8

// Intensity bin radii, in thousandths of a dalton
#define INTENSITY_BIN_RADIUS_TIGHT 150 
#define INTENSITY_BIN_RADIUS 500 

// Set VERBOSE_DEBUGGING to true if scoring (in mutation-tolerant mode) is broken.  
// Slows things down a bit, because we write out spdeadsheets:
// PRMScores.xls (verbose annotations for every last PRM bin).
// DTable.xls (in 2-mod mode) The DScore[] d.p. table
// PrefixSuffix.xls (in 2-mod mode) The PrefixTable and SuffixTable
//#define VERBOSE_DEBUGGING

typedef struct SpectralPeak
{
    int Mass;
    float Intensity;
    int IntensityRank;
    int Rank; // binned version of IntensityRank
    int IonType; // for PRM peaks only
    int FilterScore;
    int NoisePenalty;
    float PercentIntensity;
    int HasNeutralLosses; // 0, 1, or 2
    int TheoPeak; // for (greedy) interpretation 
    int Score; // for (greedy) interpretation
    // The IsotopeNeighbors array holds the indices of peaks that are potential 
    // isotopes of this peak.  If a peak was assigned a Noise ion type, but it has 
    // a neighbor peak at -1Da, then we give the peak the Isotope ion type.
    // (The +1 peak gets an IsotopeNeighbors entry)
    int IsotopeNeighbors[MAX_ISOTOPE_NEIGHBORS]; 
    // Sometimes two high-intensity peaks are separated by only 0.1 amu.  That *probably*
    // means there's one big peak that was split by the machine.  
    int NoiseNeighbors[MAX_NOISY_NEIGHBORS];
    int Index;
    int RescuedFlag;
    int AminoIndex; // For labeling purposes only!
} SpectralPeak;


typedef struct ListNode
{
    struct ListNode* Prev;
    struct ListNode* Next;
    int Entry;
} ListNode;

typedef struct MSSpectrum
{
    int MZ;
    int ParentMass;
    float SignalToNoise;
	// Parent mass based on the file (BEFORE correction)
    int FileMass; 
	// The input file may indicate no charge at all (in which case we guess),
	// a single charge (in which case we accept it, OR guess if MultiCharge is set),
	// or multiple charges (in which case we accept it, OR guess if MultiCharge is set).
	char FileCharge[6];
	int FileChargeFlag;
    int Charge;
    int PeakCount;
    // PeakAllocation is the size of the allocated Peaks array; >= PeakCount
    // When we run out of space in the array, we reallocate to double size.
    int PeakAllocation; 
    SpectralPeak* Peaks;
    int UnfilteredPeakCount;
    SpectralPeak* UnfilteredPeaks;
    int PRMPeakCount;
    float MaxIntensity; // max over all peaks
    int PMCorrectedFlag;
    struct TagGraph* Graph;
    int CandidatesScored; 
    int IntensityBinCount;
    float* BinnedIntensitiesTight; // size IntensityBinCount; used for PMC.  Tighter radius
    float* BinnedIntensities; // size IntensityBinCount
    int* BinnedIntensityLevels; // size IntensityBinCount
    int* BinPeakIndex; // size IntensityBinCount
    float IntensityCutoffLow;
    float IntensityCutoffMedium;
    float IntensityCutoffHigh;
    struct SpectrumNode* Node;
    // For use by IonScoring:
    float* IntensityThresholds;
    float* IonScoringNoiseProbabilities;
    //For Phosphorylation trickery.  we remove superdominant peaks for M-p
    int RemovedPeakIndex;
    float RemovedPeakIntensity;

#ifdef VERBOSE_DEBUGGING
    char** PRMDebugStrings;
#endif
} MSSpectrum;

int GuessSpectralCharge(MSSpectrum* Spectrum);
void UnitTestSpectrum();
void WindowFilterPeaks(MSSpectrum* Spectrum, float WindowWidth, int KeepCount);
void SpectrumAssignIsotopeNeighbors(MSSpectrum* Spectrum);
void IntensityRankPeaks();
MSSpectrum* NewSpectrum();
void FreeSpectrum(MSSpectrum* Spectrum);
int SpectrumLoadFromFile(MSSpectrum* Spectrum, FILE* DTAFile);
void SpectrumCorrectParentMass(MSSpectrum* Spectrum);
void SpectrumSetCharge(MSSpectrum* Spectrum, int Charge);
void FreeMatchList(SpectrumNode* Spectrum);
//void SpectrumComputeBinnedIntensities(SpectrumNode* Node);
void SpectrumComputeNoiseDistributions(SpectrumNode* Node);
int GuessSpectrumFormatFromExtension(char* FileName);
void FreeSpectrumNode(SpectrumNode* Node);
int GuessSpectrumFormatFromHeader(char* FilePath, MSSpectrum* Spectrum);
#endif // SPECTRUM_H
