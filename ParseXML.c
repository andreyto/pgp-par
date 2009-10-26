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
#include "CMemLeak.h"
#include "Tagger.h"
#include "base64.h"
#include "Errors.h"
#include "expat.h"

#define MZXML_BUFFER_SIZE 102400

typedef enum MZXMLScanState
{
    evMZXMLNone = 0,
    evMZXMLPrecursorMZ,
    evMZXMLPeaks
} MZXMLScanState;

typedef enum MZDataScanState
{
    evMZDataNone = 0,
    evMZDataIonSelection,
    evMZDataMZArray,
    evMZDataMZArrayBody,
    evMZDataIntensityArray,
    evMZDataIntensityArrayBody
} MZDataScanState;

#define MZXML_PARSE_LIST_SPECTRA 0
#define MZXML_PARSE_OBTAIN_PEAKS 1

// The MZXMLParseCursor is used for parsing MZXML and MZDATA formats.
// It holds the expat Parser object, and a limited amount of parse state
// (i.e. the current tag).
typedef struct MZXMLParseCursor
{
    int FirstScan;
    int LastScan;
    int ScanNumber;
    int ErrorFlag;
    XML_Parser Parser;
    int PeakCountAllocation;
    int PeakBufferSize;
    int PeakBufferPos;
    char* PeakBuffer;
    char* DecodedPeakBuffer;
    float* Peaks;
    char PrecursorMZBuffer[256];
    InputFileNode* InputFile;
    int State;
    MSSpectrum* Spectrum;
    int PeakCount;
    int ByteOrderLittle;
    int SpectrumPeaksCompleteFlag;
    void* XMLBuffer;
    int Mode;
    int MSLevel;
} MZXMLParseCursor;

typedef struct MZDataParseCursor
{
    int FirstScan;
    int LastScan;
    int ScanNumber;
    int ErrorFlag;
    XML_Parser Parser;
    int PeakCountAllocation;
    int PeakBufferSize;
    int PeakBufferPos;
    char* PeakBuffer;
    char* DecodedPeakBuffer;
    float* MZBuffer;
    float* IntensityBuffer;
    InputFileNode* InputFile;
    int State;
    MSSpectrum* Spectrum;
    int PeakCount;
    int ByteOrderLittle;
    int SpectrumPeaksCompleteFlag;
    void* XMLBuffer;
    int Mode;
    float* Peaks;
    int PrecursorMZ;
    int SpectrumStartFilePos;
    int MSLevel;
    int Precision;
} MZDataParseCursor;

// We build a single MZXMLCursor when needed, and free it when cleaning up.
MZXMLParseCursor* g_MZXMLParseCursor = NULL;
MZDataParseCursor* g_MZDataParseCursor = NULL;

void EndianByteSwap(char* Buffer, int EntrySize)
{
    char ByteSwap;
    int Pos;

    for (Pos = 0; Pos < (EntrySize / 2); Pos++)
    {
        ByteSwap = Buffer[Pos];
        Buffer[Pos] = Buffer[EntrySize - Pos - 1];
        Buffer[EntrySize - Pos - 1] = ByteSwap;
    }
}
// expat callback: Handle character data in the body of a tag. 
// The only mzxml body we care about is <peaks>
void MZXMLCharacterDataHandler(void* UserData, const XML_Char* String, int Length)
{
    MZXMLParseCursor* Cursor;
    int PeakCopySize;
    //
    Cursor = (MZXMLParseCursor*)UserData;
    if (Cursor->ErrorFlag)
    {
        return;
    }
    switch (Cursor->State)
    {
    case evMZXMLPrecursorMZ:
        strncat(Cursor->PrecursorMZBuffer, String, min(Length, 255));
        break;
    case evMZXMLPeaks:
        PeakCopySize = Length;
        if (Cursor->PeakBufferPos + PeakCopySize >= Cursor->PeakBufferSize)
        {
            REPORT_ERROR(0);
            PeakCopySize = max(0, Cursor->PeakBufferSize - Cursor->PeakBufferPos - 1);
        }
        memcpy(Cursor->PeakBuffer + Cursor->PeakBufferPos, String, PeakCopySize);
        Cursor->PeakBufferPos += PeakCopySize;
        Cursor->PeakBuffer[Cursor->PeakBufferPos] = '\0';
        break;
    // Default behavior is to ignore text:
    default:
        break;
    }
}

void MZXMLStartScan(MZXMLParseCursor* Cursor, const char** Attributes)
{
    int AttributeIndex;
    const char* Name;
    const char* Value;
    int MSLevel = 1;
    int ScanNumber = -1;
    int PeakCount = 0;
    int FilePos;
    //
    if (Cursor->ErrorFlag)
    {
        return;
    }
    for (AttributeIndex = 0; Attributes[AttributeIndex]; AttributeIndex += 2)
    {
        Name = Attributes[AttributeIndex];
        Value = Attributes[AttributeIndex + 1];
        if (!CompareStrings(Name, "msLevel"))
        {
            MSLevel = atoi(Value);
        }
        else if (!CompareStrings(Name, "peaksCount"))
        {
            PeakCount = atoi(Value);
        }
        else if (!CompareStrings(Name, "num"))
        {
            ScanNumber = atoi(Value);
        }
    }
    
    Cursor->ScanNumber = -1;
    Cursor->PeakBufferPos = 0;
    Cursor->PeakBuffer[0] = '\0';
    Cursor->PrecursorMZBuffer[0] = '\0';
    Cursor->PeakCount = PeakCount;
    Cursor->MSLevel = MSLevel;

    // If it's a level-2 scan with non-trivial peak count, then we should parse it:
    if (MSLevel >= 2 && PeakCount >= 10)
    {
        if (ScanNumber >= Cursor->FirstScan && (Cursor->LastScan < 0 || ScanNumber <= Cursor->LastScan))
        {
            FilePos = XML_GetCurrentByteIndex(Cursor->Parser);
            if (Cursor->Mode == MZXML_PARSE_LIST_SPECTRA)
            {
                AddSpectrumToList(Cursor->InputFile, FilePos, ScanNumber);
            }
            Cursor->ScanNumber = ScanNumber;
            // Allocate peak buffer, if necessary:
            if (PeakCount >= Cursor->PeakCountAllocation)
            {
                Cursor->PeakCountAllocation = PeakCount * 2;
                Cursor->PeakBufferSize = sizeof(double) * 4 * Cursor->PeakCountAllocation;
                SafeFree(Cursor->PeakBuffer);
                Cursor->PeakBuffer = (char*)malloc(Cursor->PeakBufferSize);
                SafeFree(Cursor->Peaks)
                Cursor->Peaks = (float*)malloc(sizeof(float) * 2 * Cursor->PeakCountAllocation);
                SafeFree(Cursor->DecodedPeakBuffer);
                Cursor->DecodedPeakBuffer = (char*)malloc(Cursor->PeakBufferSize);
            }
        }
    }
}

// Callback for reaching </peaks> in an mzXML parser - decode the peak array!
void MZXMLFinishPeaks(MZXMLParseCursor* Cursor, MSSpectrum* Spectrum)
{
    int Trail;
    int FloatIndex;
    int PeakCount;
    int PeakIndex;
    float Value;
    //

    PeakCount = Cursor->PeakCount;

	Trail = (PeakCount % 3);
    if (!(PeakCount % 3))
    {
        Cursor->PeakBuffer[PeakCount * 32/3] = '\0';
    }
    else
    {
        Cursor->PeakBuffer[(PeakCount * 32/3) + Trail + 1] = '\0';
    }
    b64_decode_mio(Cursor->PeakBuffer, Cursor->DecodedPeakBuffer);
    for (FloatIndex = 0; FloatIndex < (2 * PeakCount); FloatIndex++)
    {
#ifdef BYTEORDER_LITTLE_ENDIAN
        if (!Cursor->ByteOrderLittle)
        {
            EndianByteSwap(Cursor->DecodedPeakBuffer + (FloatIndex * 4), 4);
        }
#else
        if (Cursor->ByteOrderLittle)
        {
            EndianByteSwap(Cursor->DecodedPeakBuffer + (FloatIndex * 4), 4);
        }
#endif
        memcpy(Cursor->Peaks + FloatIndex, Cursor->DecodedPeakBuffer + FloatIndex * 4, 4);
    }

    Spectrum->PeakCount = PeakCount;
    Spectrum->PeakAllocation = PeakCount;
    Spectrum->Peaks = (SpectralPeak*)calloc(sizeof(SpectralPeak), PeakCount);
    if (!Spectrum->Peaks)
    {
        REPORT_ERROR_I(49, sizeof(SpectralPeak) * PeakCount);
    }

    for (PeakIndex = 0; PeakIndex < PeakCount; PeakIndex++)
    {
        Value = Cursor->Peaks[PeakIndex * 2];
        Spectrum->Peaks[PeakIndex].Mass = (int)(Value * MASS_SCALE + 0.5);
        Value = Cursor->Peaks[PeakIndex * 2 + 1];
        Spectrum->Peaks[PeakIndex].Intensity = Value;
    }
    if (Spectrum->Peaks[0].Mass < -1 || Spectrum->Peaks[0].Mass > (GlobalOptions->DynamicRangeMax + GlobalOptions->Epsilon))
    {
        REPORT_WARNING_SII(45, Spectrum->Node->InputFile->FileName, Spectrum->Node->ScanNumber,
            Spectrum->Peaks[0].Mass / MASS_SCALE);
    }
    if (Spectrum->Peaks[Spectrum->PeakCount - 1].Mass < -1 || Spectrum->Peaks[Spectrum->PeakCount - 1].Mass > (GlobalOptions->DynamicRangeMax + GlobalOptions->Epsilon))
    {
        REPORT_WARNING_SII(45, Spectrum->Node->InputFile->FileName, Spectrum->Node->ScanNumber,
            Spectrum->Peaks[Spectrum->PeakCount - 1].Mass / MASS_SCALE);
    }
    if (Spectrum->Peaks[0].Intensity < 0)
    {
        REPORT_WARNING_SIF(45, Spectrum->Node->InputFile->FileName, Spectrum->Node->ScanNumber,
            Spectrum->Peaks[0].Intensity);
    }
    if (Spectrum->Peaks[Spectrum->PeakCount - 1].Intensity < 0)
    {
        REPORT_WARNING_SIF(45, Spectrum->Node->InputFile->FileName, Spectrum->Node->ScanNumber,
            Spectrum->Peaks[0].Intensity);
    }

    Cursor->State = evMZXMLNone;
    Cursor->SpectrumPeaksCompleteFlag = 1;
    // After the end of the <peaks> flag, this scan ends.  
    // Nuke the handlers, so we can finish off the buffer in peace.
    XML_SetElementHandler(Cursor->Parser, NULL, NULL);
    XML_SetCharacterDataHandler(Cursor->Parser, NULL);
    XML_SetProcessingInstructionHandler(Cursor->Parser, NULL);
}

// expat callback: End a tag.
void MZXMLEndElement(void* UserData, const char* Tag)
{
    MZXMLParseCursor* Cursor;
    MSSpectrum* Spectrum;
    //
    Cursor = (MZXMLParseCursor*)UserData;
    Spectrum = Cursor->Spectrum;
    //printf("End '%s'\n", Tag); 
    // Set the precursor m/z, if appropriate:
    if (Cursor->State == evMZXMLPrecursorMZ)
    {
        if (Spectrum)
        {
            Spectrum->MZ = (int)(MASS_SCALE * atof(Cursor->PrecursorMZBuffer));
        }
    }

    // If we just finished <peaks>, and we have a spectrum, then set the peaks.
    if (Cursor->State != evMZXMLPeaks || !Spectrum)
    {
        Cursor->State = evMZXMLNone;
        return;
    }
    MZXMLFinishPeaks(Cursor, Spectrum);
}

void MZXMLStartPeaks(MZXMLParseCursor* Cursor, const char** Attributes)
{
    const char* Name;
    const char* Value;
    int ScanNumber = -1;
    int PeakCount = 0;
    int AttributeIndex;
    //
    if (Cursor->ErrorFlag)
    {
        return;
    }
    if (Cursor->MSLevel < 2)
    {
        return; // we don't care about peaks at level 1
    }
    
    Cursor->State = evMZXMLPeaks;
    Cursor->PeakBuffer[0] = '\0';
    Cursor->PeakBufferPos = 0;
    for (AttributeIndex = 0; Attributes[AttributeIndex]; AttributeIndex += 2)
    {
        Name = Attributes[AttributeIndex];
        Value = Attributes[AttributeIndex + 1];
        if (!CompareStrings(Name, "byteOrder"))
        {
            // Parse the byte ordering:
            if (!CompareStrings(Value, "network"))
            {
                Cursor->ByteOrderLittle = 0;
            }
            else if (!CompareStrings(Value, "little"))
            {
                Cursor->ByteOrderLittle = 0;
            }
            else if (!CompareStrings(Value, "big"))
            {
                Cursor->ByteOrderLittle = 0;
            }
        }
    }
}

// expat callback: Handle a tag and its attributes.
void MZXMLStartElement(void* UserData, const char* Tag, const char** Attributes)
{
    MZXMLParseCursor* Cursor;
    int ExpectedTag = 0;
    //
    Cursor = (MZXMLParseCursor*)UserData;
    if (Cursor->ErrorFlag)
    {
        return;
    }
    //printf("Start '%s'\n", Tag);
    // Switch on our current state, and handle the tags we expect to see in our current state.
    // Tags we don't expect are ignored (i.e. new tags can be added without breaking the parser)
    switch (Cursor->State)
    {
    default:
        // If we encounter <scan>, start the new scan:
        if (!strcmp(Tag, "scan"))
        {
            MZXMLStartScan(Cursor, Attributes);
            return;
        }
        if (!strcmp(Tag, "precursorMz"))
        {
            Cursor->State = evMZXMLPrecursorMZ;
            Cursor->PrecursorMZBuffer[0] = '\0';
            return;
        }
        if (!strcmp(Tag, "peaks"))
        {
            MZXMLStartPeaks(Cursor, Attributes);
            return;
        }
        break;
    }
}

void MZDataParseMSLevel(MZDataParseCursor* Cursor, const char** Attributes)
{
    const char* Name;
    const char* Value;
    int AttributeIndex;
    for (AttributeIndex = 0; Attributes[AttributeIndex]; AttributeIndex += 2)
    {
        Name = Attributes[AttributeIndex];
        Value = Attributes[AttributeIndex + 1];
        if (!strcmp(Name, "msLevel"))
        {
            Cursor->MSLevel = atoi(Value);
        }
    }
}

void MZDataGetPrecursorMZ(MZDataParseCursor* Cursor, const char** Attributes)
{
    const char* Name;
    const char* Value;
    int AttributeIndex;
    int MassChargeRatioFlag = 0;
    float FloatValue;
    for (AttributeIndex = 0; Attributes[AttributeIndex]; AttributeIndex += 2)
    {
        Name = Attributes[AttributeIndex];
        Value = Attributes[AttributeIndex + 1];
        if (!strcmp(Name, "name"))
        {
            if (!strcmp(Value, "MassToChargeRatio") || !strcmp(Value, "mz"))
            {
                MassChargeRatioFlag = 1;
            }
            continue;
        }
        if (!strcmp(Name, "value"))
        {
            FloatValue = (float)atof(Value);
        }
    }
    if (MassChargeRatioFlag)
    {
        Cursor->PrecursorMZ = (int)(MASS_SCALE * FloatValue);
    }
}

// Look up the number of peaks.  Re-allocate buffers if necessary.
void MZDataGetPeakCount(MZDataParseCursor* Cursor, const char** Attributes)
{
    const char* Name;
    const char* Value;
    int AttributeIndex;
    for (AttributeIndex = 0; Attributes[AttributeIndex]; AttributeIndex += 2)
    {
        Name = Attributes[AttributeIndex];
        Value = Attributes[AttributeIndex + 1];
        if (!strcmp(Name, "precision"))
        {
            Cursor->Precision = atoi(Value);
        }
        else if (!strcmp(Name, "endian"))
        {
            Cursor->ByteOrderLittle = 0; // default
            if (!strcmp(Value, "little"))
            {
                Cursor->ByteOrderLittle = 1;
            }
            else if (!strcmp(Value, "network"))
            {
                Cursor->ByteOrderLittle = 0;
            }
            else if (!strcmp(Value, "big"))
            {
                Cursor->ByteOrderLittle = 0;
            }
            continue;
        }
        if (!strcmp(Name, "length"))
        {
            Cursor->PeakCount = atoi(Value);
            // Is this more peaks than we can currently handle?
            if (Cursor->PeakCount >= Cursor->PeakCountAllocation)
            {
                Cursor->PeakCountAllocation = Cursor->PeakCount * 2;
                SafeFree(Cursor->PeakBuffer);
                Cursor->PeakBufferSize = sizeof(double) * 4 * Cursor->PeakCountAllocation;
                Cursor->PeakBuffer = (char*)malloc(Cursor->PeakBufferSize);
                SafeFree(Cursor->DecodedPeakBuffer);
                Cursor->DecodedPeakBuffer = (char*)malloc(sizeof(float) * Cursor->PeakCountAllocation);
                SafeFree(Cursor->IntensityBuffer);
                Cursor->IntensityBuffer = (float*)malloc(sizeof(float) * Cursor->PeakCountAllocation);
                SafeFree(Cursor->MZBuffer);
                Cursor->MZBuffer = (float*)malloc(sizeof(float) * Cursor->PeakCountAllocation);
            }
            continue;
        }
    }
}

// expat callback: Handle a tag and its attributes.
void MZDataStartElement(void* UserData, const char* Tag, const char** Attributes)
{
    MZDataParseCursor* Cursor;
    int ExpectedTag = 0;
    //
    Cursor = (MZDataParseCursor*)UserData;
    if (Cursor->ErrorFlag)
    {
        return;
    }

    // Switch on our current state, and handle the tags we expect to see in our current state.
    // Tags we don't expect are ignored (i.e. new tags can be added without breaking the parser)
    // If we encounter <spectrum>, start the new scan:
    if (!strcmp(Tag, "spectrum"))
    {
        //MZDataStartScan(Cursor, Attributes);
        Cursor->SpectrumStartFilePos = XML_GetCurrentByteIndex(Cursor->Parser);
        Cursor->ScanNumber = atoi(Attributes[1]);
        return;
    }
    // If we encounter <ionSelection>, update our state:
    if (!strcmp(Tag, "ionSelection"))
    {
        Cursor->State = evMZDataIonSelection;
        return;
    }
    // If we encounter <cvParam> within ionSelection, set precursor m/z if possible:
    if (!strcmp(Tag, "cvParam") && Cursor->State == evMZDataIonSelection)
    {
        MZDataGetPrecursorMZ(Cursor, Attributes);
        return;
    }
    
    if (!strcmp(Tag, "data"))
    {
        switch (Cursor->State)
        {
        case evMZDataMZArray:
            Cursor->State = evMZDataMZArrayBody;
            Cursor->PeakBufferPos = 0;
            MZDataGetPeakCount(Cursor, Attributes);
            break;
        case evMZDataIntensityArray:
            Cursor->State = evMZDataIntensityArrayBody;
            Cursor->PeakBufferPos = 0;
            MZDataGetPeakCount(Cursor, Attributes);
            break;
        default:
            REPORT_ERROR(0);
            break;
        }
        return;
    }
    if (!strcmp(Tag, "mzArrayBinary"))   
    {
        Cursor->State = evMZDataMZArray;
        return;
    }
    if (!strcmp(Tag, "intenArrayBinary"))   
    {
        Cursor->State = evMZDataIntensityArray;
        return;
    }
    if (!strcmp(Tag, "precursor"))
    {
        MZDataParseMSLevel(Cursor, Attributes);
        return;
    }
}

void MZDataCompleteSpectrum(MZDataParseCursor* Cursor, MSSpectrum* Spectrum)
{
    int PeakIndex;

    if (Cursor->Mode == MZXML_PARSE_LIST_SPECTRA)
    {
        if (Cursor->PeakCount >= 10 && Cursor->MSLevel > 1)
        {
            AddSpectrumToList(Cursor->InputFile, Cursor->SpectrumStartFilePos, Cursor->ScanNumber);
        }
        return;
    }
    
    Spectrum->PeakCount = Cursor->PeakCount;
    Spectrum->PeakAllocation = Cursor->PeakCount;
    Spectrum->Peaks = (SpectralPeak*)calloc(sizeof(SpectralPeak), Cursor->PeakCount);
    Spectrum->MZ = Cursor->PrecursorMZ;

    for (PeakIndex = 0; PeakIndex < Cursor->PeakCount; PeakIndex++)
    {
        Spectrum->Peaks[PeakIndex].Mass = (int)(Cursor->MZBuffer[PeakIndex] * MASS_SCALE + 0.5);
        Spectrum->Peaks[PeakIndex].Intensity = Cursor->IntensityBuffer[PeakIndex];
    }
    if (Spectrum->Peaks[0].Mass < -1 || Spectrum->Peaks[0].Mass > (GlobalOptions->DynamicRangeMax + GlobalOptions->Epsilon))
    {
        REPORT_WARNING_SII(45, Spectrum->Node->InputFile->FileName, Spectrum->Node->ScanNumber,
            Spectrum->Peaks[0].Mass);
    }
    if (Spectrum->Peaks[Spectrum->PeakCount - 1].Mass < -1 || Spectrum->Peaks[Spectrum->PeakCount - 1].Mass > (GlobalOptions->DynamicRangeMax + GlobalOptions->Epsilon))
    {
        REPORT_WARNING_SII(45, Spectrum->Node->InputFile->FileName, Spectrum->Node->ScanNumber,
            Spectrum->Peaks[0].Mass);
    }
    if (Spectrum->Peaks[0].Intensity < 0)
    {
        REPORT_WARNING_SIF(45, Spectrum->Node->InputFile->FileName, Spectrum->Node->ScanNumber,
            Spectrum->Peaks[0].Intensity);
    }
    if (Spectrum->Peaks[Spectrum->PeakCount - 1].Intensity < 0)
    {
        REPORT_WARNING_SIF(45, Spectrum->Node->InputFile->FileName, Spectrum->Node->ScanNumber,
            Spectrum->Peaks[0].Intensity);
    }

    Cursor->State = evMZDataNone;
    Cursor->SpectrumPeaksCompleteFlag = 1;
    // After the end of the <peaks> flag, this scan ends.  
    // Nuke the handlers, so we can finish off the buffer in peace.
    XML_SetElementHandler(Cursor->Parser, NULL, NULL);
    XML_SetCharacterDataHandler(Cursor->Parser, NULL);
    XML_SetProcessingInstructionHandler(Cursor->Parser, NULL);

}

// MZData callback for end </data> tag:
// - Decode the base64-encoded float array
// - Store the floats in the MZ or Intensity array
void MZDataProcessEncodedPeakData(MZDataParseCursor* Cursor, MSSpectrum* Spectrum)
{
    int PeakCount;
    int Trail;
    int FloatIndex;
    int EncodedRecordSize;
    //
    if (Cursor->State == evMZDataIntensityArrayBody)
    {
        Cursor->State = evMZDataIntensityArray;
    }
    else if (Cursor->State == evMZDataMZArrayBody)
    {
        Cursor->State = evMZDataMZArray;
    }
    else
    {
        REPORT_ERROR(0);
    }
    if (!Spectrum)
    {
        return;
    }
    PeakCount = Cursor->PeakCount;
    Trail = (PeakCount % 3);
    if (!(PeakCount % 3))
    {
        Cursor->PeakBuffer[PeakCount * 32/3] = '\0';
    }
    else
    {
        Cursor->PeakBuffer[(PeakCount * 32/3) + Trail + 1] = '\0';
    }
    b64_decode_mio(Cursor->PeakBuffer, Cursor->DecodedPeakBuffer);
    if (Cursor->Precision == 32)
    {
        EncodedRecordSize = 4;
    }
    else if (Cursor->Precision == 64)
    {
        EncodedRecordSize = 8;
    }
    else
    {
        // Default to 32bit:
        EncodedRecordSize = 4;
    }
    for (FloatIndex = 0; FloatIndex < PeakCount; FloatIndex++)
    {
#ifdef BYTEORDER_LITTLE_ENDIAN
        if (!Cursor->ByteOrderLittle)
        {
            EndianByteSwap(Cursor->DecodedPeakBuffer + (FloatIndex * EncodedRecordSize), EncodedRecordSize);
            //ByteSwap = Cursor->DecodedPeakBuffer[FloatIndex*4];
            //Cursor->DecodedPeakBuffer[FloatIndex*4] = Cursor->DecodedPeakBuffer[FloatIndex*4 + 3];
            //Cursor->DecodedPeakBuffer[FloatIndex*4 + 3] = ByteSwap;
            //ByteSwap = Cursor->DecodedPeakBuffer[FloatIndex*4 + 1];
            //Cursor->DecodedPeakBuffer[FloatIndex*4 + 1] = Cursor->DecodedPeakBuffer[FloatIndex*4 + 2];
            //Cursor->DecodedPeakBuffer[FloatIndex*4 + 2] = ByteSwap;
        }
#else
        if (Cursor->ByteOrderLittle)
        {
            EndianByteSwap(Cursor->DecodedPeakBuffer + (FloatIndex * EncodedRecordSize), EncodedRecordSize);
            //ByteSwap = Cursor->DecodedPeakBuffer[FloatIndex*4];
            //Cursor->DecodedPeakBuffer[FloatIndex*4] = Cursor->DecodedPeakBuffer[FloatIndex*4 + 3];
            //Cursor->DecodedPeakBuffer[FloatIndex*4 + 3] = ByteSwap;
            //ByteSwap = Cursor->DecodedPeakBuffer[FloatIndex*4 + 1];
            //Cursor->DecodedPeakBuffer[FloatIndex*4 + 1] = Cursor->DecodedPeakBuffer[FloatIndex*4 + 2];
            //Cursor->DecodedPeakBuffer[FloatIndex*4 + 2] = ByteSwap;
        }
#endif
        if (Cursor->State == evMZDataMZArrayBody || Cursor->State == evMZDataMZArray)
        {
            Cursor->MZBuffer[FloatIndex] = *((float*)(Cursor->DecodedPeakBuffer + FloatIndex * EncodedRecordSize));
            //memcpy(Cursor->MZBuffer + FloatIndex, Cursor->DecodedPeakBuffer + FloatIndex * 4, 4);
        }
        else
        {
            //memcpy(Cursor->IntensityBuffer + FloatIndex, Cursor->DecodedPeakBuffer + FloatIndex * 4, 4);
            Cursor->IntensityBuffer[FloatIndex] = *((float*)(Cursor->DecodedPeakBuffer + FloatIndex * EncodedRecordSize));
        }
    }
}

// expat callback: End a tag.
void MZDataEndElement(void* UserData, const char* Tag)
{
    MZDataParseCursor* Cursor;
    MSSpectrum* Spectrum;
    //
    Cursor = (MZDataParseCursor*)UserData;
    Spectrum = Cursor->Spectrum;

    if (!strcmp(Tag, "spectrum"))
    {
        MZDataCompleteSpectrum(Cursor, Spectrum);
        Cursor->SpectrumPeaksCompleteFlag = 1;
        return;
    }
    if (!strcmp(Tag, "data"))
    {
        MZDataProcessEncodedPeakData(Cursor, Spectrum);
        return;
    }
    if (!strcmp(Tag, "intenArrayBinary"))
    {
        Cursor->State = evMZDataNone;
        return;
    }
    if (!strcmp(Tag, "mzArrayBinary"))
    {
        Cursor->State = evMZDataNone;
        return;
    }
    if (!strcmp(Tag, "ionSelection"))
    {
        Cursor->State = evMZDataNone;
        return;
    }
}

// expat callback: Handle character data in the body of a tag. 
// The only mzdata body we care about is <data>
void MZDataCharacterDataHandler(void* UserData, const XML_Char* String, int Length)
{
    MZDataParseCursor* Cursor;
    int PeakCopySize;
    //
    Cursor = (MZDataParseCursor*)UserData;
    if (Cursor->ErrorFlag)
    {
        return;
    }
    switch (Cursor->State)
    {
    case evMZDataMZArrayBody:
    case evMZDataIntensityArrayBody: // deliberate fallthrough
        PeakCopySize = Length;
        if (Cursor->PeakBufferPos + PeakCopySize >= Cursor->PeakBufferSize)
        {
            REPORT_ERROR(0);
            PeakCopySize = max(0, Cursor->PeakBufferSize - Cursor->PeakBufferPos - 1);
        }
        memcpy(Cursor->PeakBuffer + Cursor->PeakBufferPos, String, PeakCopySize);
        Cursor->PeakBufferPos += PeakCopySize;
        Cursor->PeakBuffer[Cursor->PeakBufferPos] = '\0';
        break;
    // Default behavior is to ignore text:
    default:
        break;
    }
}

MZDataParseCursor* GetMZDataParseCursor()
{
    if (g_MZDataParseCursor)
    {
        return g_MZDataParseCursor;
    }
    g_MZDataParseCursor = (MZDataParseCursor*)calloc(1, sizeof(MZDataParseCursor));
    g_MZDataParseCursor->PeakCountAllocation = 1024;
    g_MZDataParseCursor->PeakBufferSize = sizeof(double) * 4 * g_MZDataParseCursor->PeakCountAllocation;
    g_MZDataParseCursor->PeakBuffer = (char*)malloc(g_MZDataParseCursor->PeakBufferSize);
    g_MZDataParseCursor->DecodedPeakBuffer = (char*)malloc(g_MZDataParseCursor->PeakCountAllocation * sizeof(float));
    g_MZDataParseCursor->MZBuffer = (float*)malloc(g_MZDataParseCursor->PeakCountAllocation * sizeof(float));
    g_MZDataParseCursor->IntensityBuffer = (float*)malloc(g_MZDataParseCursor->PeakCountAllocation * sizeof(float));
    g_MZDataParseCursor->Parser = XML_ParserCreate(NULL);
    //g_MZDataParseCursor->XMLBuffer = XML_GetBuffer(g_MZDataParseCursor->Parser, sizeof(char) * MZXML_BUFFER_SIZE);
    //if (!g_MZDataParseCursor->XMLBuffer)
    //{
    //    printf("* Error: Unable to get XML buffer of size %d\n", MZXML_BUFFER_SIZE);
    //}
    return g_MZDataParseCursor;
}

MZXMLParseCursor* GetMZXMLParseCursor()
{
    if (g_MZXMLParseCursor)
    {
        return g_MZXMLParseCursor;
    }
    g_MZXMLParseCursor = (MZXMLParseCursor*)calloc(1, sizeof(MZXMLParseCursor));
    g_MZXMLParseCursor->PeakCountAllocation = 1024;
    g_MZXMLParseCursor->PeakBufferSize = sizeof(double) * 4 * g_MZXMLParseCursor->PeakCountAllocation;
    g_MZXMLParseCursor->PeakBuffer = (char*)malloc(g_MZXMLParseCursor->PeakBufferSize);
    g_MZXMLParseCursor->DecodedPeakBuffer = (char*)malloc(g_MZXMLParseCursor->PeakBufferSize);
    g_MZXMLParseCursor->Peaks = (float*)malloc(sizeof(float) * 2 * g_MZXMLParseCursor->PeakCountAllocation);
    g_MZXMLParseCursor->Parser = XML_ParserCreate(NULL);
    //g_MZXMLParseCursor->XMLBuffer = XML_GetBuffer(g_MZXMLParseCursor->Parser, sizeof(char) * MZXML_BUFFER_SIZE);
    //if (!g_MZXMLParseCursor->XMLBuffer)
    //{
    //    printf("* Error: Unable to get XML buffer of size %d\n", MZXML_BUFFER_SIZE);
    //}

    return g_MZXMLParseCursor;
}

void FreeMZXMLParseCursor()
{
    if (!g_MZXMLParseCursor)
    {
        return;
    }
    SafeFree(g_MZXMLParseCursor->PeakBuffer);
    SafeFree(g_MZXMLParseCursor->Peaks);
    SafeFree(g_MZXMLParseCursor->DecodedPeakBuffer);
    if (g_MZXMLParseCursor->Parser)
    {
        XML_ParserFree(g_MZXMLParseCursor->Parser);
    }
    SafeFree(g_MZXMLParseCursor);
    g_MZXMLParseCursor = NULL;
}

void FreeMZDataParseCursor()
{
    if (!g_MZDataParseCursor)
    {
        return;
    }
    SafeFree(g_MZDataParseCursor->PeakBuffer);
    SafeFree(g_MZDataParseCursor);
    XML_ParserFree(g_MZDataParseCursor->Parser);
    g_MZDataParseCursor = NULL;
}

// Parse through an MZXML file to get a list of spectra and their byte offsets.
void ParseSpectraFromMZXML(char* FileName, InputFileNode* InputFile, int FirstScan, int LastScan)
{
    FILE* MZXMLFile;
    MZXMLParseCursor* Cursor;
    int FilePos = 0;
    int DoneFlag = 0;
    //void* XMLBuffer;
    int BytesRead;
    int XMLParseResult;
    int Error;
    //

    MZXMLFile = fopen(FileName, "rb");
    if (!MZXMLFile)
    {
        REPORT_ERROR_S(8, FileName);
        return;
    }
    printf("Parse spectra from '%s'...\n", FileName);
    Cursor = GetMZXMLParseCursor();
    Cursor->FirstScan = FirstScan;
    Cursor->LastScan = LastScan;
    Cursor->InputFile = InputFile;
    Cursor->ErrorFlag = 0;
    Cursor->Spectrum = NULL;
    Cursor->Mode = MZXML_PARSE_LIST_SPECTRA;
    XML_SetUserData(Cursor->Parser, Cursor);
    XML_SetElementHandler(Cursor->Parser, MZXMLStartElement, MZXMLEndElement);
    XML_SetCharacterDataHandler(Cursor->Parser, MZXMLCharacterDataHandler);
    while (!DoneFlag)
    {
        // Get a buffer (parser handles the memory):
        Cursor->XMLBuffer = XML_GetBuffer(Cursor->Parser, sizeof(char) * MZXML_BUFFER_SIZE);
        if (!Cursor->XMLBuffer)
        {
            printf("* ParseSpectraFromMZXML Error: Unable to get XML buffer of size %d\n", MZXML_BUFFER_SIZE);
            break;
        }

        // Read into the buffer:
        BytesRead = ReadBinary(Cursor->XMLBuffer, sizeof(char), MZXML_BUFFER_SIZE, MZXMLFile);
        if (!BytesRead)
        {
            // We'll call XML_Parse once more, this time with DoneFlag set to 1. 
            DoneFlag = 1;
        }

        // Parse this block o' text:
        XMLParseResult = XML_ParseBuffer(Cursor->Parser, BytesRead, DoneFlag);
        if (!XMLParseResult)
        {
            printf("XML Parse error - file position ~%d\n", XML_GetCurrentByteIndex(Cursor->Parser));
            Error = XML_GetErrorCode(Cursor->Parser);
            printf("Error code %d description '%s'\n", Error, XML_ErrorString(Error));
        }

        // If Cursor->ErrorFlag is set, then the file isn't valid!  Error out
        // now, since recovery could be difficult.
        if (Cursor->ErrorFlag)
        {
            break;
        }
        FilePos += BytesRead;
    }
    
    // Close file, free memory:
    fclose(MZXMLFile);
	FreeMZXMLParseCursor();
}

// Parse ONE spectrum from the file.  Return true on success.
int SpectrumLoadMZXML(MSSpectrum* Spectrum, FILE* MZXMLFile)
{
    MZXMLParseCursor* Cursor;
    int FilePos = 0;
    int DoneFlag = 0;
    //void* XMLBuffer;
    int BytesRead;
    int XMLParseResult;
    int ReturnResult = 1;
    //

    Cursor = GetMZXMLParseCursor();
    Cursor->Spectrum = Spectrum;
    Cursor->Mode = MZXML_PARSE_OBTAIN_PEAKS;
    Cursor->ErrorFlag = 0;
    Cursor->FirstScan = 0;
    Cursor->LastScan = -1;
    XML_ParserReset(Cursor->Parser, NULL);
    XML_SetUserData(Cursor->Parser, Cursor);
    XML_SetElementHandler(Cursor->Parser, MZXMLStartElement, MZXMLEndElement);
    XML_SetCharacterDataHandler(Cursor->Parser, MZXMLCharacterDataHandler);

    while (!DoneFlag)
    {
        // Get a buffer (parser handles the memory):
        Cursor->XMLBuffer = XML_GetBuffer(Cursor->Parser, sizeof(char) * MZXML_BUFFER_SIZE);
        if (!Cursor->XMLBuffer)
        {
            printf("* SpectrumLoadMZXML Error: Unable to get XML buffer of size %d\n", MZXML_BUFFER_SIZE);
            break;
        }

        // Read into the buffer:
        BytesRead = ReadBinary(Cursor->XMLBuffer, sizeof(char), MZXML_BUFFER_SIZE, MZXMLFile);
        if (!BytesRead)
        {
            // We'll call XML_Parse once more, this time with DoneFlag set to 1. 
            DoneFlag = 1;
        }

        // Parse this block o' text:
        XMLParseResult = XML_ParseBuffer(Cursor->Parser, BytesRead, DoneFlag);
        if (!XMLParseResult)
        {
            Cursor->ErrorFlag = 1;
            // If we have peaks...let's NOT report a warning...because we're parsing a sub-document, 
            // and we'll run off the edge and get well-formedness complaints.
            // Newer expat versions will have the ability to abort when we hit the </scan>
            // tag ending.
            if (!Cursor->Spectrum->PeakCount)
            {
                ReturnResult = 0;
            }
        }

        // If Cursor->ErrorFlag is set, then the file isn't valid!  Error out
        // now, since recovery could be difficult.
        if (Cursor->ErrorFlag)
        {
            break;
        }
        if (Cursor->SpectrumPeaksCompleteFlag)
        {
            break;
        }

        FilePos += BytesRead;
    }
    // Sanity check: We must have a precursor m/z!
    if (!Cursor->Spectrum->MZ)
    {
        ReturnResult = 0;
    }

    return ReturnResult;    
}

// Parse through an mzData file to get a list of spectra and their byte offsets.
void ParseSpectraFromMZData(char* FileName, InputFileNode* InputFile, int FirstScan, int LastScan)
{
    FILE* MZXMLFile;
    MZDataParseCursor* Cursor;
    int FilePos = 0;
    int DoneFlag = 0;
    //void* XMLBuffer;
    int BytesRead;
    int XMLParseResult;
    int Error;
    //

    MZXMLFile = fopen(FileName, "rb");
    if (!MZXMLFile)
    {
        REPORT_ERROR_S(8, FileName);
        return;
    }
    printf("Parse spectra from '%s'...\n", FileName);
    Cursor = GetMZDataParseCursor();
    Cursor->FirstScan = FirstScan;
    Cursor->LastScan = LastScan;
    Cursor->InputFile = InputFile;
    Cursor->ErrorFlag = 0;
    Cursor->Spectrum = NULL;
    Cursor->Mode = MZXML_PARSE_LIST_SPECTRA;
    XML_SetUserData(Cursor->Parser, Cursor);
    XML_SetElementHandler(Cursor->Parser, MZDataStartElement, MZDataEndElement);
    XML_SetCharacterDataHandler(Cursor->Parser, MZDataCharacterDataHandler);
    //XMLBuffer = Cursor->XMLBuffer;

    while (!DoneFlag)
    {
        // Get a buffer (parser handles the memory):
        Cursor->XMLBuffer = XML_GetBuffer(Cursor->Parser, sizeof(char) * MZXML_BUFFER_SIZE);
        if (!Cursor->XMLBuffer)
        {
            printf("* Error: Unable to get XML buffer of size %d\n", MZXML_BUFFER_SIZE);
            break;
        }

        // Read into the buffer:
        BytesRead = ReadBinary(Cursor->XMLBuffer, sizeof(char), MZXML_BUFFER_SIZE, MZXMLFile);
        if (!BytesRead)
        {
            // We'll call XML_Parse once more, this time with DoneFlag set to 1. 
            DoneFlag = 1;
        }

        // Parse this block o' text:
        XMLParseResult = XML_Parse(Cursor->Parser, Cursor->XMLBuffer, BytesRead, DoneFlag);
        if (!XMLParseResult)
        {
            printf("XML Parse error - file position ~%d\n", XML_GetCurrentByteIndex(Cursor->Parser));
            Error = XML_GetErrorCode(Cursor->Parser);
            printf("Error code %d description '%s'\n", Error, XML_ErrorString(Error));
        }

        // If Cursor->ErrorFlag is set, then the file isn't valid!  Error out
        // now, since recovery could be difficult.
        if (Cursor->ErrorFlag)
        {
            break;
        }
        FilePos += BytesRead;
    }
    
    // Close file, free memory:
    fclose(MZXMLFile);
	FreeMZDataParseCursor();
}


// Parse ONE spectrum from the file
void SpectrumLoadMZData(MSSpectrum* Spectrum, FILE* MZXMLFile)
{
    MZDataParseCursor* Cursor;
    int FilePos = 0;
    int DoneFlag = 0;
    //void* XMLBuffer;
    int BytesRead;
    int XMLParseResult;
    //

    Cursor = GetMZDataParseCursor();
    Cursor->Spectrum = Spectrum;
    Cursor->Mode = MZXML_PARSE_OBTAIN_PEAKS;
    Cursor->ErrorFlag = 0;
    XML_ParserReset(Cursor->Parser, NULL);
    XML_SetUserData(Cursor->Parser, Cursor);
    XML_SetElementHandler(Cursor->Parser, MZDataStartElement, MZDataEndElement);
    XML_SetCharacterDataHandler(Cursor->Parser, MZDataCharacterDataHandler);
    while (!DoneFlag)
    {
        // Get a buffer (parser handles the memory):
        Cursor->XMLBuffer = XML_GetBuffer(Cursor->Parser, sizeof(char) * MZXML_BUFFER_SIZE);
        if (!Cursor->XMLBuffer)
        {
            printf("* Error: Unable to get XML buffer of size %d\n", MZXML_BUFFER_SIZE);
            break;
        }

        // Read into the buffer:
        BytesRead = ReadBinary(Cursor->XMLBuffer, sizeof(char), MZXML_BUFFER_SIZE, MZXMLFile);
        if (!BytesRead)
        {
            // We'll call XML_Parse once more, this time with DoneFlag set to 1. 
            DoneFlag = 1;
        }

        // Parse this block o' text:
        XMLParseResult = XML_Parse(Cursor->Parser, Cursor->XMLBuffer, BytesRead, DoneFlag);
        if (!XMLParseResult)
        {
            // Let's NOT report a warning...because we're parsing a sub-document, 
            // and we'll run off the edge and get well-formedness complaints.
            // Newer expat versions will have the ability to abort when we hit the </scan>
            // tag ending.
        }

        // If Cursor->ErrorFlag is set, then the file isn't valid!  Error out
        // now, since recovery could be difficult.
        if (Cursor->ErrorFlag)
        {
            break;
        }
        if (Cursor->SpectrumPeaksCompleteFlag)
        {
            break;
        }

        FilePos += BytesRead;
    }
    
    // Close file, free memory:
    //fclose(MZXMLFile);
}

