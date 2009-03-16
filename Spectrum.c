
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
#include "Spectrum.h"
#include "Utils.h"
#include "Inspect.h"
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <string.h>
#include <math.h>
#include "Tagger.h"
#include "Errors.h"
#include "ParseXML.h"

#define INITIAL_PEAK_COUNT 100
#define INITIAL_PRM_PEAK_COUNT 500

#define MINIMUM_ALLOWED_PARENT_MASS GLYCINE_MASS
#define MAXIMUM_ALLOWED_PARENT_MASS 6000*DALTON

// This should be MORE than enough peaks for any realistic spectrum.
// If there are more than this, we refuse to parse them all, so there.
#define MAX_PEAKS_PER_SPECTRUM 10000

/////////////////////////////////////////////////////////////////////////////////
// Forward declarations:
int SpectrumLoadHeaderLine(MSSpectrum* Spectrum, char* LineBuffer);
void AttemptParentMassPeakRemoval(MSSpectrum* Spectrum);

/////////////////////////////////////////////////////////////////////////////////
// Functions:

void SpectrumComputeSignalToNoise(MSSpectrum* Spectrum)
{
    int IntensePeakIndex;
    int MedianPeakIndex;
    int PeakIndex;
    float Signal = 0;
    float Noise = 0;
    //
    IntensePeakIndex = min(5, Spectrum->PeakCount) / 2;
    MedianPeakIndex = Spectrum->PeakCount / 2;
    for (PeakIndex = 0; PeakIndex < Spectrum->PeakCount; PeakIndex++)
    {
        if (Spectrum->Peaks[PeakIndex].IntensityRank == IntensePeakIndex)
        {
            Signal = Spectrum->Peaks[PeakIndex].Intensity;
        }
        if (Spectrum->Peaks[PeakIndex].IntensityRank == MedianPeakIndex)
        {
            Noise = Spectrum->Peaks[PeakIndex].Intensity;
        }
    }
    Spectrum->SignalToNoise = Signal / (float)max(1.0, Noise);
}

// Remove peaks that are not reasonably high for their mass window.
// If WindowWidth and KeepCount are <= 0, use reasonable defaults.
void WindowFilterPeaks(MSSpectrum* Spectrum, float WindowWidth, int KeepCount)
{
    int FilterPeakIndex;
    int NewIndex;
    int OtherPeakIndex;
    float* Intensities;
    int Neighbors;
    float WindowStart;
    float WindowEnd;
    int FilteredCount = 0;
    //
    if (Spectrum->UnfilteredPeaks)
    {
        // We've already performed window filtering; don't do it again!
        return;
    }
    if (WindowWidth <= 0)
    {
        WindowWidth = DEFAULT_WINDOW_WIDTH;
    }
    if (KeepCount <= 0)
    {
        KeepCount = DEFAULT_WINDOW_KEEP_COUNT;
    }

    //
    Intensities = (float*)calloc(Spectrum->PeakCount, sizeof(float));
    for (FilterPeakIndex = 0; FilterPeakIndex < Spectrum->PeakCount; FilterPeakIndex++)
    {
        WindowStart = Spectrum->Peaks[FilterPeakIndex].Mass - (WindowWidth / (float)2.0);
        WindowEnd = Spectrum->Peaks[FilterPeakIndex].Mass + (WindowWidth / (float)2.0);
        Neighbors = 0;
        for (OtherPeakIndex = 0; OtherPeakIndex < Spectrum->PeakCount; OtherPeakIndex++)
        {
            if (Spectrum->Peaks[OtherPeakIndex].Mass > WindowEnd)
            {
                break;
            }
            if (Spectrum->Peaks[OtherPeakIndex].Mass > WindowStart)
            {
                Intensities[Neighbors] = Spectrum->Peaks[OtherPeakIndex].Intensity;
                Neighbors++;
            }
        }
        qsort(Intensities, Neighbors, sizeof(float), (QSortCompare)CompareFloats);
        if (Neighbors < KeepCount || Spectrum->Peaks[FilterPeakIndex].Intensity >= Intensities[KeepCount - 1])
        {
            Spectrum->Peaks[FilterPeakIndex].FilterScore = 1;
            FilteredCount++;
        }
    }
    SafeFree(Intensities);
    // New array:
    Spectrum->UnfilteredPeakCount = Spectrum->PeakCount;
    Spectrum->UnfilteredPeaks = Spectrum->Peaks;
    Spectrum->PeakCount = FilteredCount;
    Spectrum->Peaks = (SpectralPeak*)calloc(FilteredCount, sizeof(SpectralPeak));
    NewIndex = 0;
    for (FilterPeakIndex = 0; FilterPeakIndex < Spectrum->UnfilteredPeakCount; FilterPeakIndex++)
    {
        if (Spectrum->UnfilteredPeaks[FilterPeakIndex].FilterScore)
        {
            memcpy(Spectrum->Peaks + NewIndex, Spectrum->UnfilteredPeaks + FilterPeakIndex, sizeof(SpectralPeak));
            Spectrum->Peaks[NewIndex].Index = NewIndex;
            NewIndex++;
        }
    }

}

// Sort from MOST to LEAST intense:
int ComparePeaksIntensity(const SpectralPeak* PeakA, const SpectralPeak* PeakB)
{
    if (PeakA->Intensity < PeakB->Intensity)
    {
        return 1;
    }
    if (PeakA->Intensity > PeakB->Intensity)
    {
        return -1;
    }
    return 0;

}

int ComparePeaksByMass(const SpectralPeak* PeakA, const SpectralPeak* PeakB)
{
    if (PeakA->Mass < PeakB->Mass)
    {
        return -1;
    }
    if (PeakA->Mass > PeakB->Mass)
    {
        return 1;
    }
    return 0;
}

void IntensityRankPeaks(MSSpectrum* Spectrum)
{
    int PeakIndex;
    //
    qsort(Spectrum->Peaks, Spectrum->PeakCount, sizeof(SpectralPeak), (QSortCompare)ComparePeaksIntensity);
    for (PeakIndex = 0; PeakIndex < Spectrum->PeakCount; PeakIndex++)
    {
        Spectrum->Peaks[PeakIndex].IntensityRank = PeakIndex;
    }
    qsort(Spectrum->Peaks, Spectrum->PeakCount, sizeof(SpectralPeak), (QSortCompare)ComparePeaksByMass);
    SpectrumComputeSignalToNoise(Spectrum);
}


void FreeMatchList(SpectrumNode* Spectrum)
{
    Peptide* MatchNode;
    Peptide* MatchPrev = NULL;
    for (MatchNode = Spectrum->FirstMatch; MatchNode; MatchNode = MatchNode->Next)
    {
        if (MatchPrev)
        {
            FreePeptideNode(MatchPrev);
        }
        MatchPrev = MatchNode;
    }
    if (MatchPrev)
    {
        FreePeptideNode(MatchPrev);
    }
    Spectrum->MatchCount = 0;
    Spectrum->FirstMatch = NULL;
    Spectrum->LastMatch = NULL;
}

void FreeSpectrum(MSSpectrum* Spectrum)
{
    if (!Spectrum)
    {
        return;
    }
    SafeFree(Spectrum->UnfilteredPeaks);
    SafeFree(Spectrum->Peaks);
    if (Spectrum->Graph)
    {
        FreeTagGraph(Spectrum->Graph);
        Spectrum->Graph = NULL;
    }
    SafeFree(Spectrum->BinnedIntensities);
    SafeFree(Spectrum->BinnedIntensitiesTight);
    SafeFree(Spectrum->BinnedIntensityLevels);
    SafeFree(Spectrum->BinPeakIndex);
    SafeFree(Spectrum->IntensityThresholds);
    SafeFree(Spectrum->IonScoringNoiseProbabilities);
    SafeFree(Spectrum);
}

// Constructor: Allocate and return a spectrum
MSSpectrum* NewSpectrum()
{
    MSSpectrum* Spectrum;
    // Allocate a spectrum, with a reasonable amount of space to store peaks
    Spectrum = (MSSpectrum*)calloc(1, sizeof(MSSpectrum));
    Spectrum->Peaks = (SpectralPeak*)calloc(INITIAL_PEAK_COUNT, sizeof(SpectralPeak));
    Spectrum->PeakAllocation = INITIAL_PEAK_COUNT;
    return Spectrum;
}

int SpectrumAddPeak(MSSpectrum* Spectrum, float Mass, float Intensity)
{
    int OldAllocation;
    
    // If necessary, reallocate:
    if (Spectrum->PeakCount > MAX_PEAKS_PER_SPECTRUM)
    {
        if (Spectrum->Node->InputFile)
        {
            REPORT_ERROR_IS(31, Spectrum->Node->ScanNumber, Spectrum->Node->InputFile->FileName);
        }
        else
        {
            REPORT_ERROR_IS(31, Spectrum->Node->ScanNumber, "??");
        }
        return 0;
    }
    if (Spectrum->PeakCount == Spectrum->PeakAllocation)
    {
        OldAllocation = Spectrum->PeakAllocation;
        Spectrum->PeakAllocation = max(200, Spectrum->PeakAllocation * 2);
        if (OldAllocation)
        {
            Spectrum->Peaks = (SpectralPeak*)realloc(Spectrum->Peaks, sizeof(SpectralPeak) * Spectrum->PeakAllocation);
            memset(Spectrum->Peaks + OldAllocation, 0, sizeof(SpectralPeak) * (Spectrum->PeakAllocation - OldAllocation));
        }
        else
        {
            Spectrum->Peaks = (SpectralPeak*)calloc(Spectrum->PeakAllocation, sizeof(SpectralPeak));
        }
    }
    ROUND_MASS(Mass, Spectrum->Peaks[Spectrum->PeakCount].Mass);
    //Spectrum->Peaks[Spectrum->PeakCount].Mass = Mass;
    Spectrum->Peaks[Spectrum->PeakCount].Intensity = Intensity;
    Spectrum->Peaks[Spectrum->PeakCount].FilterScore = 0; // init
    Spectrum->Peaks[Spectrum->PeakCount].NoisePenalty = 0; // init
    Spectrum->Peaks[Spectrum->PeakCount].PercentIntensity = 0; // init
    memset(Spectrum->Peaks[Spectrum->PeakCount].IsotopeNeighbors, -1, sizeof(int)*MAX_ISOTOPE_NEIGHBORS);
    memset(Spectrum->Peaks[Spectrum->PeakCount].NoiseNeighbors, -1, sizeof(int)*MAX_NOISY_NEIGHBORS);
    //Log("Added peak %d to spectrum.  (Alloc size is now %d)\n", Spectrum->PeakCount, Spectrum->PeakAllocation);
    Spectrum->PeakCount++;
    Spectrum->MaxIntensity = max(Spectrum->MaxIntensity, Intensity);
    return 1;
}

// Handle the header line for .dta, .pkl, .ms2 formats.
int SpectrumLoadHeaderLine(MSSpectrum* Spectrum, char* LineBuffer)
{
    char* StrA;
    char* StrB;
    char* StrC;
    float Mass;
	int Charge;
    // Default case: The first line should be the parent mass and the charge.
    StrA = strtok(LineBuffer, WHITESPACE);
    if (!StrA || !*StrA)
    {
        return 0;
    }
    StrB = strtok(NULL, WHITESPACE);
    if (!StrB)
    {
        return 0;
    }
    StrC = strtok(NULL, WHITESPACE);

    // MS2 file: Z, charge, parent-mass.
    if (!strcmp(StrA, "Z"))
    {
        
        Mass = (float)atof(StrC);
        ROUND_MASS(Mass, Spectrum->ParentMass);
        if (Spectrum->ParentMass < MINIMUM_ALLOWED_PARENT_MASS || 
            Spectrum->ParentMass > MAXIMUM_ALLOWED_PARENT_MASS)
        {
            if (Spectrum->Node && Spectrum->Node->InputFile)
            {
                REPORT_ERROR_IIS(42, Spectrum->ParentMass / MASS_SCALE, Spectrum->Node->FilePosition, Spectrum->Node->InputFile->FileName);
            }
            else
            {
                REPORT_ERROR_I(43, Spectrum->ParentMass / MASS_SCALE);
            }
            return 0;
        }
        Charge = atoi(StrB);
        if (Charge > 6)
        {
            printf("** Invalid charge '%s' - maximum is 6\n", StrB);
            return 0;
        }
		Spectrum->Charge = Charge;
		Spectrum->FileChargeFlag = 1;
		Spectrum->FileCharge[Charge] = 1;
		Spectrum->MZ = (Spectrum->ParentMass + (Charge - 1) * HYDROGEN_MASS) / Charge;
        return 1;
    }

    // Header line of a PKL file: precursor mz, precursor intensity, and charge.
    if (StrC)
    {
        Mass = (float)atof(StrA);
        ROUND_MASS(Mass, Spectrum->ParentMass);
        if (Spectrum->ParentMass < MINIMUM_ALLOWED_PARENT_MASS || 
            Spectrum->ParentMass > MAXIMUM_ALLOWED_PARENT_MASS)
        {
            printf("** Error in SpectrumLoadFromFile: Mass %.2f not legal.\n", Mass);
            return 0;
        }
        Charge = atoi(StrC);
        if (Charge > 6)
        {
            printf("** Invalid charge '%s' - maximum is 6\n", StrC);
            return 0;
        }
		Spectrum->FileCharge[Charge] = 1;
		Spectrum->FileChargeFlag = 1;
		Spectrum->Charge = Charge;
        if (Charge)
        {
            Spectrum->ParentMass = (Spectrum->ParentMass * Charge) - (Charge - 1)*HYDROGEN_MASS;
        }
    }
    else
    {
        // DTA file:
        Mass = (float)atof(StrA);
        if (Mass < 1)
        {
            // Invalid header line - the mass can't be zero!
            return 0;
        }
        ROUND_MASS(Mass, Spectrum->ParentMass);
        Charge = atoi(StrB);
        if (!Charge)
        {
            Spectrum->MZ = Spectrum->ParentMass;
            Spectrum->ParentMass = 0;
        }
        else
        {
            // The file's mass is the residue mass + 19, which is the parent mass.
			Spectrum->FileCharge[Charge] = 1;
			Spectrum->FileChargeFlag = 1;
			Spectrum->Charge = Charge;
            //Spectrum->ParentMass -= HYDROGEN_MASS; // remove one H+
            //Spectrum->ParentMass = (float)atof(StrA) - HYDROGEN_MASS; // remove one H+
        }        
    }
    return 1;
}

int GuessSpectralCharge(MSSpectrum* Spectrum)
{
    int PeakIndex;
    float MeanMass = 0;
    int Charge;
    int BestDiff = 9999999;
    int BestCharge = 2;
    int ParentMass;
    int Diff;
    // Compute the MEDIAN peak mass:
    for (PeakIndex = 0; PeakIndex < Spectrum->PeakCount; PeakIndex++)
    {
        if (PeakIndex >= Spectrum->PeakCount / 2)
        {
            MeanMass = (float)Spectrum->Peaks[PeakIndex].Mass;
            break;
        }
    }
    //MeanMass /= Spectrum->PeakCount;
    // Use a charge that will bring the parent mass as close as possible to the mean peak mass x 2
    for (Charge = 1; Charge < 10; Charge++)
    {
        ParentMass = (Spectrum->MZ * Charge) - (HYDROGEN_MASS * (Charge - 1));
        Diff = abs(ParentMass - (int)(MeanMass*2));
        if (Diff < BestDiff)
        {
            BestDiff = Diff;
            BestCharge = Charge;
        }
    }
    return BestCharge;
}

// Initial parent mass computation.  We do it LATER if our charge is 0 (unknown).
void SpectrumComputeParentMass(MSSpectrum* Spectrum)
{
    //
    if (!Spectrum->Charge)
    {
        return; // We'll tweak later!
    }    
    else
    {
        Spectrum->MZ = (Spectrum->ParentMass + (Spectrum->Charge - 1)*HYDROGEN_MASS) / Spectrum->Charge;
        Spectrum->FileMass = Spectrum->ParentMass;
    }
}

// Return FALSE when we're done loading
int SpectrumHandleMS2ColonLine(int LineNumber, int FilePos, char* LineBuffer, void* UserData)
{
    char* ValueStr;
    float ParentMass;
    float Mass;
    float Intensity;
    MSSpectrum* Spectrum;
	int Charge;
    int Result;
    //
    Spectrum = (MSSpectrum*)UserData;

    if (LineBuffer[0] == ':')
    {
        if (Spectrum->PeakCount)
        {
            // We've loaded some peaks already, so this is the scan number of 
            // the NEXT scan.
            return 0;
        }
        else
        {
            // We know our scan number already, so we do nothing.
            return 1;
        }
    }

    // If we don't know our MZ yet, then we load it now.  Otherwise we add a peak.
    if (!Spectrum->MZ)
    {
        ValueStr = strtok(LineBuffer, WHITESPACE);
        if (!ValueStr)
        {
            return 1; // INVALID LINE, stop now
        }
        ParentMass = (float)(atof(ValueStr) * MASS_SCALE);
        ValueStr = strtok(NULL, WHITESPACE);
        if (!ValueStr)
        {
            return 0; // INVALID LINE, stop now
        }
        Charge = atoi(ValueStr);
        if (Charge)
        {
			Spectrum->Charge = Charge;
			Spectrum->FileCharge[Charge] = 1;
			Spectrum->FileChargeFlag = 1;
            Spectrum->ParentMass = (int)(ParentMass - HYDROGEN_MASS + 0.5); 
            Spectrum->MZ = (int)((ParentMass + (Spectrum->Charge - 1)*HYDROGEN_MASS) / (float)Spectrum->Charge + 0.5);
        }
        else
        {
            Spectrum->ParentMass = (int)(ParentMass + 0.5);
            Spectrum->MZ = (int)(ParentMass + 0.5);
        }
        return 1;
    }
    // Ordinary peak
    ValueStr = strtok(LineBuffer,  WHITESPACE);
    if (!ValueStr)
    {
        return 0; // INVALID LINE, stop now
    }
    Mass = (float)atof(ValueStr);
    ValueStr = strtok(NULL,  WHITESPACE);
    if (!ValueStr)
    {
        return 0; // INVALID LINE, stop now
    }
    Intensity = (float)atof(ValueStr);
    Result = SpectrumAddPeak(Spectrum, Mass, Intensity);
    return Result;
}

int SpectrumLoadMGFCallback(int LineNumber, int FilePos, char* LineBuffer, void* UserData)
{
    MSSpectrum* Spectrum;
    float Mass;
    float Intensity;
    char* WordA;
    char* WordB;
    char* EQWordA;
    int Result;
	int Charge;
	char* AndWord;
    //
    Spectrum = (MSSpectrum*)UserData;

    // If we see a command we recognize, then handle it:
    WordA = strtok(LineBuffer, WHITESPACE);
    WordB = strtok(NULL, WHITESPACE);
    EQWordA = strtok(WordA, "=");

    if (!CompareStrings(WordA, "END"))
    {
        if (WordB && !CompareStrings(WordB, "IONS"))
        {
            // Stop parsing lines now!
            return 0; 
        }
    }
    else if (!CompareStrings(EQWordA, "PEPMASS"))
    {
        Mass = (float)atof(LineBuffer + 8);
        ROUND_MASS(Mass, Spectrum->MZ);
        if (Spectrum->MZ < MINIMUM_ALLOWED_PARENT_MASS || Spectrum->MZ > MAXIMUM_ALLOWED_PARENT_MASS)
        {
            // Illegal mass!
            if (Spectrum->Node->InputFile)
            {
                REPORT_ERROR_IS(32, Spectrum->Node->ScanNumber, Spectrum->Node->InputFile->FileName);
            }
            else
            {
                REPORT_ERROR_IS(32, Spectrum->Node->ScanNumber, "???");
            }
            return 0;
        }
    }
    else if (!CompareStrings(EQWordA, "CHARGE"))
    {
		
        Charge = atoi(LineBuffer + 7);
		if (Charge)
		{
			Spectrum->Charge = Charge;
			Spectrum->FileCharge[Charge] = 1;
			Spectrum->FileChargeFlag = 1;
		}
		// the CHARGE line may have the form "2+ and 3+"
		if (WordB && !CompareStrings(WordB, "and"))
		{
			Charge = atoi(WordB + 4);
			if (Charge)
			{
				Spectrum->Charge = Charge;
				Spectrum->FileCharge[Charge] = 1;
				Spectrum->FileChargeFlag = 1;
			}
		}
    }
    else
    {
        // Default: Try to read an m/z and intensity
        Mass = (float)atof(WordA);
        if (Mass && WordB)
        {
            Intensity = (float)atof(WordB);
            Result = SpectrumAddPeak(Spectrum, Mass, Intensity);
            return Result;
        }
    }
    return 1;
}

// Load spectrum from an MGF file.  See one or more header lines, then some
// peaks, then an "END IONS" line.
int SpectrumLoadMGF(MSSpectrum* Spectrum, FILE* DTAFile)
{

    ParseFileByLines(DTAFile, SpectrumLoadMGFCallback, Spectrum, 0);
    if (Spectrum->Charge)
    {
        Spectrum->ParentMass = Spectrum->MZ * Spectrum->Charge - (HYDROGEN_MASS * (Spectrum->Charge - 1));
    }

    return 1;
}

int GuessMS2FormatFromLine(int LineNumber, int FilePos, char* LineBuffer, void* UserData)
{
    MSSpectrum* Spectrum;
    Spectrum = (MSSpectrum*)UserData;
    if (LineBuffer[0] == ':')
    {
        Spectrum->Node->InputFile->Format = SPECTRUM_FORMAT_MS2_COLONS;
        return 0;
    }
    if (LineBuffer[0] == 'Z')
    {
        Spectrum->Node->InputFile->Format = SPECTRUM_FORMAT_MS2;
        return 0;
    }
    return 1;

}

// Return 1 if we succeeded.
int SpectrumLoadDTAFileLine(int LineNumber, int FilePos, char* LineBuffer, void* UserData)
{
    MSSpectrum* Spectrum;
    int Result;
    float Mass;
    float Intensity;
    char* Str;
    // 
    Spectrum = (MSSpectrum*)UserData;
    // Special case: MS2 format handles one or more "Z" lines
    if (LineBuffer[0] == 'Z' && (LineBuffer[1] == ' ' || LineBuffer[1] == '\t'))
    {
        Result = SpectrumLoadHeaderLine(Spectrum, LineBuffer);
        return Result;
    }
    // Special case: MS2 format skips the first "S" line, and knows the 
	// second "S" line it sees marks the end of the record
    if (LineBuffer[0] == 'S' && (LineBuffer[1] == ' ' || LineBuffer[1] == '\t'))
    {
		if (LineNumber > 1)
		{
			return 0;
		}
		else
		{
			return 1;
		}
    }
    if (LineNumber == 1)
    {
        Result = SpectrumLoadHeaderLine(Spectrum, LineBuffer);
        return Result;
    }

    // After the first line, we expect to see lines of the form "Mass Intensity"
    Str = strtok(LineBuffer, WHITESPACE);
    if (!Str)
    {
        return 1;
    }
    Mass = (float)atof(Str);
    if (!Mass) 
    {   
        return 1;
    }
    Str = strtok(NULL, WHITESPACE);
    if (!Str)
    {
        // This line had only two pieces on it.  Invalid syntax!
        //printf("**Error in file '%s': peak lines must contain mass AND intensity\n", Spectrum->Node->InputFile->FileName);
        REPORT_ERROR_IS(33, LineNumber, Spectrum->Node->InputFile->FileName);
        return 0;
    }
    Intensity = (float)atof(Str);
    if (!Intensity) 
    {   
        // invalid intensity?  Assume that a string staring with "0" really means intensity zero,
        // god help us.
        if (Str[0] != '0')
        {
            REPORT_ERROR_IS(33, LineNumber, Spectrum->Node->InputFile->FileName);
            return 0;
        }
    }
    // If there's a third piece on the line, then stop parsing now.  (That happens if we run
    // off the end of a record in a pkl file, into the start of the next record):
    Str = strtok(NULL, WHITESPACE);
    if (Str)
    {
        return 0;
    }
    Result = SpectrumAddPeak(Spectrum, Mass, Intensity);
    return Result;
}

int GuessSpectrumFormatFromHeader(char* FilePath, MSSpectrum* Spectrum)
{
    FILE* MS2File;
    //
    MS2File = fopen(FilePath, "rb");
    ParseFileByLines(MS2File, GuessMS2FormatFromLine, Spectrum, 0);
    fclose(MS2File);
    return Spectrum->Node->InputFile->Format;
}

// SpectrumLoadFromFile: Return True if the spectrum is valid, False if it's not.
// Example of an invalid spectrum file: Sequest .out files contaminating the .dta directory.
// Iterate over lines, handling the header specially.
int SpectrumLoadFromFile(MSSpectrum* Spectrum, FILE* DTAFile)
{    
    int ReturnCode = 1;
    int MS2ChargeLineSeen = 0;
    //

    // handle XML formats separately from line-based foramts:
    switch (Spectrum->Node->InputFile->Format)
    {
    case SPECTRUM_FORMAT_MZXML:
        ReturnCode = SpectrumLoadMZXML(Spectrum, DTAFile);
        break;
    case SPECTRUM_FORMAT_MZDATA:
        SpectrumLoadMZData(Spectrum, DTAFile);
        break;
    case SPECTRUM_FORMAT_MGF:
        SpectrumLoadMGF(Spectrum, DTAFile);
        break;
    case SPECTRUM_FORMAT_MS2_COLONS:
        ParseFileByLines(DTAFile, SpectrumHandleMS2ColonLine, Spectrum, 0);
        break;
    case SPECTRUM_FORMAT_PKL:
    case SPECTRUM_FORMAT_DTA:
    case SPECTRUM_FORMAT_MS2:
    default:
        ParseFileByLines(DTAFile, SpectrumLoadDTAFileLine, Spectrum, 0);
        break;
    }
    if (ReturnCode)
    {
        SpectrumComputeParentMass(Spectrum);
    }
    //if (GlobalOptions->PhosphorylationFlag)
    //{
    //    AttemptParentMassPeakRemoval(Spectrum);
    //}
    return ReturnCode;
}

////For phosphorylated spectra, the superprominent M-p peak can 
////fritz the charge state guessing, and tagging.  So we remove it.
//void AttemptParentMassPeakRemoval(MSSpectrum* Spectrum)
//{
//    int MostIntensePeakIndex;
//    int MostIntenseMass;
//    int PeakIndex;
//    float MostIntense = 0.0;
//    float NextMostIntense = 0.0;
//    int Diff;
//    int ExpectedDiff;
//    int ExpectedDiff2;
//    int Epsilon = 2 * DALTON;
//    int Charge;
//    //
//    for (PeakIndex = 0; PeakIndex < Spectrum->PeakCount; PeakIndex++)
//    {
//        if (Spectrum->Peaks[PeakIndex].Intensity > MostIntense)
//        {
//            NextMostIntense = MostIntense;
//            MostIntense = Spectrum->Peaks[PeakIndex].Intensity;
//            MostIntensePeakIndex = PeakIndex;
//            MostIntenseMass = Spectrum->Peaks[PeakIndex].Mass;
//        }
//        else if(Spectrum->Peaks[PeakIndex].Intensity > NextMostIntense)
//        {
//            NextMostIntense = Spectrum->Peaks[PeakIndex].Intensity;
//        }
//    }
//    //printf("Most intense %f, next %f\n",MostIntense, NextMostIntense);
//    //if more than 3 times great, and in the right place, remove peak.
//    //if (MostIntense < 2 * NextMostIntense)
//    //{
//    //    return;
//    //}
//    //printf ("MZ of %d, charge %d\n", Spectrum->MZ, Spectrum->Charge);
//    // If the spectrum has a charge, then trust, otherwise try charge 2, 3
//	//Set m/z with the new parentmass and charge that was just assigned in ConstructTags
//	printf("Old MZ %f\n",Spectrum->MZ);
//	Spectrum->MZ = (Spectrum->ParentMass + (Spectrum->Charge - 1) * HYDROGEN_MASS) / Spectrum->Charge;
//	printf("New MZ %f\n",Spectrum->MZ);
//	return;
//    if (Spectrum->Charge)
//    {
//        Diff = abs(Spectrum->MZ - MostIntenseMass);
//        ExpectedDiff = PHOSPHATE_WATER_MASS / Spectrum->Charge;
//        ExpectedDiff2 = (PHOSPHATE_WATER_MASS + WATER_MASS) / Spectrum->Charge;
//        if (abs (Diff - ExpectedDiff) < Epsilon)
//        { //remove peak
//            Spectrum->RemovedPeakIndex = MostIntensePeakIndex;
//            Spectrum->RemovedPeakIntensity = Spectrum->Peaks[MostIntensePeakIndex].Intensity;
//            Spectrum->Peaks[MostIntensePeakIndex].Intensity = 1.0; //cut to ground
//        }
//        else if (abs(Diff - ExpectedDiff2) < Epsilon)
//        { //remove peak
//            Spectrum->RemovedPeakIndex = MostIntensePeakIndex;
//            Spectrum->RemovedPeakIntensity = Spectrum->Peaks[MostIntensePeakIndex].Intensity;
//            Spectrum->Peaks[MostIntensePeakIndex].Intensity = 1.0; //cut to ground
//        }
//    }
//    else
//    {
//        for (Charge = 1; Charge <= 3; Charge++)
//        {
//            Diff = abs(Spectrum->MZ - MostIntenseMass);
//            ExpectedDiff = PHOSPHATE_WATER_MASS/ Charge;
//            ExpectedDiff2 = (PHOSPHATE_WATER_MASS + WATER_MASS)/ Charge;
//            // printf("Charge %d, Diff %d, ExpectedDiff %d\n", Charge, Diff, ExpectedDiff);
//            if (abs (Diff - ExpectedDiff) < Epsilon)
//            { // remove peak
//                Spectrum->RemovedPeakIndex = MostIntensePeakIndex;
//                Spectrum->RemovedPeakIntensity = Spectrum->Peaks[MostIntensePeakIndex].Intensity;
//                Spectrum->Peaks[MostIntensePeakIndex].Intensity = 1.0; //cut to ground
//                Spectrum->Charge = Charge; // This is a big enough clue, that we are going to guess charge
//                Spectrum->MZ = MostIntenseMass + ExpectedDiff; //testing this feature
//                break;
//            }
//            else if (abs(Diff - ExpectedDiff2) < Epsilon)
//            { // remove peak
//                Spectrum->RemovedPeakIndex = MostIntensePeakIndex;
//                Spectrum->RemovedPeakIntensity = Spectrum->Peaks[MostIntensePeakIndex].Intensity;
//                Spectrum->Peaks[MostIntensePeakIndex].Intensity = 1.0; //cut to ground
//                Spectrum->Charge = Charge;
//                Spectrum->MZ = MostIntenseMass + ExpectedDiff2;
//                break;
//            }
//        } // end for
//    } // end else
//
//}

// Called AFTER filtering.  Looks 1Da to the left of peaks for potential isotope neighbors.
void SpectrumAssignIsotopeNeighbors(MSSpectrum* Spectrum)
{
    // Don't worry *too* much about efficiency, as this happens only once during scoring
    int PeakIndex;
    int OldPeakIndex;
    int IsotopeCount;
    int NoiseCount;
    int MaxMass;
    int MinMass;
    int OtherPeakIndex;
    float IntensityPercent;
    //
    // Assign noise penalty:
    for (PeakIndex = 0; PeakIndex < Spectrum->PeakCount; PeakIndex++)
    {
        IntensityPercent = Spectrum->Peaks[PeakIndex].Intensity / Spectrum->MaxIntensity;
        Spectrum->Peaks[PeakIndex].PercentIntensity = IntensityPercent;
        if (IntensityPercent < 0.05)
        {
            Spectrum->Peaks[PeakIndex].NoisePenalty = -921;//0.0001
        }
        else if (IntensityPercent < 0.3)
        {
            Spectrum->Peaks[PeakIndex].NoisePenalty = -1382; //0.000001
        }
        else if (IntensityPercent < 0.6)
        {
            Spectrum->Peaks[PeakIndex].NoisePenalty = -1842; //0.00000001
        }
        else 
        {
            Spectrum->Peaks[PeakIndex].NoisePenalty = -2303; //0.0000000001
        }
    }
    // First, look for isotope neighbors.  Scan downward from each peak:
    for (PeakIndex = 0; PeakIndex < Spectrum->PeakCount; PeakIndex++)
    {
        MaxMass = Spectrum->Peaks[PeakIndex].Mass - 79;
        MinMass = Spectrum->Peaks[PeakIndex].Mass - 121;
        IsotopeCount = 0;
        for (OldPeakIndex = max(0, PeakIndex - 1); OldPeakIndex; OldPeakIndex--)
        {
            if (Spectrum->Peaks[OldPeakIndex].Mass < MinMass)
            {
                break;
            }
            if (Spectrum->Peaks[OldPeakIndex].Mass > MaxMass)
            {
                continue;
            }
            Spectrum->Peaks[PeakIndex].IsotopeNeighbors[IsotopeCount++] = OldPeakIndex;
        }
    }
    // Now look for noise-neighbors (peaks which could be the same peak, but are split
    // due to limitations in recording).
    for (PeakIndex = 0; PeakIndex < Spectrum->PeakCount; PeakIndex++)
    {
        NoiseCount = 0;
        MaxMass = Spectrum->Peaks[PeakIndex].Mass + 21; // 0.2 Da radius
        MinMass = Spectrum->Peaks[PeakIndex].Mass - 21;
        for (OtherPeakIndex = PeakIndex + 1; OtherPeakIndex < min(Spectrum->PeakCount, PeakIndex + 5); OtherPeakIndex++)
        {
            if (Spectrum->Peaks[OtherPeakIndex].Mass > MaxMass)
            {
                break;
            }
            Spectrum->Peaks[PeakIndex].NoiseNeighbors[NoiseCount++] = OtherPeakIndex;
        }
        for (OtherPeakIndex = max(0, PeakIndex - 1); OtherPeakIndex > max(-1, PeakIndex - 5); OtherPeakIndex--)
        {
            if (Spectrum->Peaks[OtherPeakIndex].Mass < MinMass)
            {
                break;
            }
            Spectrum->Peaks[PeakIndex].NoiseNeighbors[NoiseCount++] = OtherPeakIndex;
        }
    }
}

void SpectrumSetCharge(MSSpectrum* Spectrum, int Charge)
{
    //MZ = ((Charge-1)*1.0078 + self->Spectrum->ParentMass) / self->Spectrum->Charge;
    Spectrum->Charge = Charge;
    Spectrum->PMCorrectedFlag = 0;
    Spectrum->ParentMass = (Spectrum->MZ * Charge) - (Charge - 1) * HYDROGEN_MASS;
}

// Compute the low/med/hi intensity cutoffs for the spectrum.
void ComputeSpectrumIntensityCutoffs(MSSpectrum* Spectrum)
{
    int PeakIndex;
    float GrassIntensity;
    float TotalIntensity;
    int CutoffRank;
    float SortedIntensity[200];
    int WeakPeakCount = 0;
    //
    TotalIntensity = 0;
    CutoffRank = (int)(Spectrum->ParentMass / (100 * DALTON));
    for (PeakIndex = 0; PeakIndex < Spectrum->PeakCount; PeakIndex++)
    {
        if (Spectrum->Peaks[PeakIndex].IntensityRank >= CutoffRank)
        {
            SortedIntensity[WeakPeakCount] = Spectrum->Peaks[PeakIndex].Intensity;
            WeakPeakCount++;
        }
        TotalIntensity += Spectrum->Peaks[PeakIndex].Intensity;
        if (WeakPeakCount == 200)
        {
            break;
        }
    }
    if (!WeakPeakCount)
    {
        //printf("** Error in ComputeSpectrumIntensityCutoffs: No weak peak ranks found?  Intensity ranking must be complete here.\n");
        if (!Spectrum->PeakCount)
        {
            return;
        }
        GrassIntensity = TotalIntensity / (2 * Spectrum->PeakCount);
    }
    else
    {
        qsort(SortedIntensity, WeakPeakCount, sizeof(float), (QSortCompare)CompareFloats);
        GrassIntensity = SortedIntensity[WeakPeakCount / 2];
    }
    Spectrum->IntensityCutoffLow = (float)0.25 * GrassIntensity;
    Spectrum->IntensityCutoffMedium = 3 * GrassIntensity;
    Spectrum->IntensityCutoffHigh = 10 * GrassIntensity;
}

//// Allocate and populate BinnedIntensities for the spectrum.  Assumes that ParentMass is set.
//void SpectrumComputeBinnedIntensities(SpectrumNode* Node) // OBSOLETE
//{
//    int MaxParentMass = 0;
//    MSSpectrum* Spectrum;
//    int PeakIndex;
//    int Bin;
//    int NearBin;
//    SpectralPeak* Peak;
//    float Intensity;
//    int BinScalingFactor = 100; // One bin per 0.1Da
//    
//    // A spectrum has at most this many "high" peaks (one per 100Da)
//    int SuperPeakCount;
//
//    static int* BestIntensityRank = NULL;
//    static int BestIntensityRankSize = 0;
//    //
//    Spectrum = Node->Spectrum;
//    if (!Spectrum)
//    {
//        return;
//    }
//    SuperPeakCount = Spectrum->ParentMass / (100 * DALTON);
//    MaxParentMass = Spectrum->MZ * 3;
//    Spectrum->IntensityBinCount = (MaxParentMass + DALTON) / BinScalingFactor; 
//    SafeFree(Spectrum->BinnedIntensities);
//    SafeFree(Spectrum->BinnedIntensitiesTight);
//    SafeFree(Spectrum->BinnedIntensityLevels);
//    SafeFree(Spectrum->BinPeakIndex);
//    Spectrum->BinnedIntensities = (float*)calloc(Spectrum->IntensityBinCount, sizeof(float));
//    Spectrum->BinnedIntensitiesTight = (float*)calloc(Spectrum->IntensityBinCount, sizeof(float));
//    Spectrum->BinnedIntensityLevels = (int*)calloc(Spectrum->IntensityBinCount, sizeof(int));
//    Spectrum->BinPeakIndex = (int*)calloc(Spectrum->IntensityBinCount, sizeof(int));
//
//    if (BestIntensityRankSize < Spectrum->IntensityBinCount)
//    {
//        SafeFree(BestIntensityRank);
//        BestIntensityRankSize = Spectrum->IntensityBinCount + 500;
//        BestIntensityRank = (int*)calloc(BestIntensityRankSize, sizeof(int));
//    }
//    for (Bin = 0; Bin < Spectrum->IntensityBinCount; Bin++)
//    {
//        Spectrum->BinPeakIndex[Bin] = -1;
//        BestIntensityRank[Bin] = 999;
//    }
//
//    // Iterate over spectral peaks, putting intensity into bins:
//    for (PeakIndex = 0; PeakIndex < Spectrum->PeakCount; PeakIndex++)
//    {
//        Peak = Spectrum->Peaks + PeakIndex;
//        Bin = (Peak->Mass + 50) / BinScalingFactor;
//        for (NearBin = Bin - 6; NearBin < Bin + 7; NearBin++)
//        {
//            if (NearBin < 0 || NearBin >= Spectrum->IntensityBinCount)
//            {
//                continue;
//            }
//            if (abs(Peak->Mass - (NearBin * BinScalingFactor)) > INTENSITY_BIN_RADIUS)
//            {
//                continue;
//            }
//            Spectrum->BinnedIntensities[NearBin] += Peak->Intensity;
//            BestIntensityRank[Bin] = min(BestIntensityRank[Bin], Peak->IntensityRank);
//            if (Spectrum->BinPeakIndex[NearBin] < 0)
//            {
//                Spectrum->BinPeakIndex[NearBin] = PeakIndex;
//            }
//            if (abs(Peak->Mass - (NearBin * BinScalingFactor)) <= INTENSITY_BIN_RADIUS_TIGHT)
//            {
//                Spectrum->BinnedIntensitiesTight[NearBin] += Peak->Intensity;
//            }
//        }
//    }
//    // Compute the intensity level (absent, lo, med, hi) for each bin:
//    ComputeSpectrumIntensityCutoffs(Spectrum);
//    for (Bin = 0; Bin < Spectrum->IntensityBinCount; Bin++)
//    {
//        Intensity = Spectrum->BinnedIntensities[Bin];
//        if (Intensity > Spectrum->IntensityCutoffHigh && BestIntensityRank[Bin] < SuperPeakCount)
//        {
//            Spectrum->BinnedIntensityLevels[Bin] = 3;
//        }
//        else if (Intensity > Spectrum->IntensityCutoffMedium)
//        {
//            Spectrum->BinnedIntensityLevels[Bin] = 2;
//        }
//        else if (Intensity > Spectrum->IntensityCutoffLow)
//        {
//            Spectrum->BinnedIntensityLevels[Bin] = 1;
//        }
//        else
//        {
//            Spectrum->BinnedIntensityLevels[Bin] = 0;
//        }
//    }
//}

void SpectrumComputeNoiseDistributions(SpectrumNode* Node)
{
    MSSpectrum* Spectrum;
    int BinCountA;
    int BinCountB;
    int BinCountC;
    int BinCountD;
    int BinCutoffA;
    int BinCutoffB;
    int Bin;
    int Index;
    int IntensityRank;
    SpectrumTweak* Tweak;
    int TweakIndex;
    // Compute the distributions of intensity-levels for the three sectors according
    // to each parent mass:
    Spectrum = Node->Spectrum;
    for (TweakIndex = 0; TweakIndex < TWEAK_COUNT; TweakIndex++)
    {
        Tweak = Node->Tweaks + TweakIndex;
        if (!Tweak->Charge)
        {
            continue;
        }
        BinCutoffA = (int)((Node->Tweaks[TweakIndex].ParentMass * 0.3333 + 5) / 100);
        BinCutoffB = (int)((Node->Tweaks[TweakIndex].ParentMass * 0.6667 + 5) / 100);
        BinCountA = 0;
        BinCountB = 0;
        BinCountC = 0;
        BinCountD = 0;
        // SECTOR_COUNT
        BinCutoffA = (int)((Node->Tweaks[TweakIndex].ParentMass * 0.5 + 5) / 100);
        for (Index = 0; Index < 8; Index++)
        {
            Node->Tweaks[TweakIndex].Intensities[Index] = 1; // padding-probability
        }
        for (Bin = 0; Bin < Spectrum->IntensityBinCount; Bin++)
        {
            if (Bin >= BinCutoffA)
            {
                BinCountB++;
                Tweak->Intensities[4 + Spectrum->BinnedIntensityLevels[Bin]] += 1.0;
            }
            else
            {
                BinCountA++;
                Tweak->Intensities[0 + Spectrum->BinnedIntensityLevels[Bin]] += 1.0;
            }
        }
        for (IntensityRank = 0; IntensityRank < 4; IntensityRank++)
        {
            Tweak->Intensities[0 + IntensityRank] = (float)log((Tweak->Intensities[0 + IntensityRank] + 2) / (BinCountA + 2));
            Tweak->Intensities[4 + IntensityRank] = (float)log((Tweak->Intensities[4 + IntensityRank] + 2) / (BinCountB + 2));
        }
        
    }
}

// Add a spectrum to the list of spectra to be searched. 
void AddSpectrumToList(InputFileNode* InputFile, int FilePos, int ScanNumber)
{
    SpectrumNode* NewNode;

    NewNode = (SpectrumNode*)calloc(1, sizeof(SpectrumNode));
    NewNode->InputFile = InputFile;
    if (GlobalOptions->LastSpectrum)
    {
        GlobalOptions->LastSpectrum->Next = NewNode;
    }
    else
    {
        GlobalOptions->FirstSpectrum = NewNode;
    }
    NewNode->FilePosition = FilePos;
    NewNode->ScanNumber = ScanNumber;
    GlobalOptions->LastSpectrum = NewNode;
    GlobalOptions->SpectrumCount++;
    InputFile->SpectrumCount++;
}


int GuessSpectrumFormatFromExtension(char* FileName)
{
    char* Extension;
    for (Extension = FileName + strlen(FileName); Extension > FileName; Extension--)
    {
        if (*Extension == '.')
        {
            break;
        }
    }
    if (!CompareStrings(Extension, ".out"))
    {
        // sequest gunk, ignore.
        return SPECTRUM_FORMAT_INVALID;
    }
    if (!CompareStrings(Extension, ".ms2"))
    {
        return SPECTRUM_FORMAT_MS2_COLONS; //SPECTRUM_FORMAT_MS2;
    }
    if (!CompareStrings(Extension, ".mzxml"))
    {
        return SPECTRUM_FORMAT_MZXML;
    }
    if (!CompareStrings(Extension, ".mzdata"))
    {
        return SPECTRUM_FORMAT_MZDATA;
    }
    if (!CompareStrings(Extension, ".mgf"))
    {
        return SPECTRUM_FORMAT_MGF;
    }
    if (!CompareStrings(Extension, ".dta"))
    {
        return SPECTRUM_FORMAT_DTA;
    }
    if (!CompareStrings(Extension, ".pkl"))
    {
        return SPECTRUM_FORMAT_PKL;
    }

    // Unexpected extension.  Let's ASSUME that it's a .dta file.
    REPORT_WARNING_S(30, FileName);
    return SPECTRUM_FORMAT_DTA;
}

void FreeSpectrumNode(SpectrumNode* Node)
{
    int TweakIndex;
    //
    if (!Node)
    {
        return;
    }
    for (TweakIndex = 0; TweakIndex < TWEAK_COUNT; TweakIndex++)
    {
        SafeFree(Node->Tweaks[TweakIndex].PRMScores);
        Node->Tweaks[TweakIndex].PRMScores = NULL;
    }
    if (Node->Spectrum)
    {
        FreeSpectrum(Node->Spectrum);
    }
    SafeFree(Node);
}
