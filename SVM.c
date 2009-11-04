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

// SVM support functions.
// We employ SVMs to distinguish between (1) true and false peptide classifications,
// and (2) true and false mutation assignments. 
// We also can use SVMs for charge state determination (+2 versus +3, currently; +1 is easy and for +4 and beyond we
// have no data) and parent mass correction.
#include "CMemLeak.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//#include <malloc.h>
#include "SVM.h"
#include "Utils.h"
#include "Inspect.h"
#include "Spectrum.h"
#include "Trie.h"
#include "Scorpion.h"
#include "BN.h"
#include "Score.h"
#include "Errors.h"
#include "IonScoring.h"

// Forward declarations:
float LDAClassify(float* Features);

// Global variables:
extern PRMBayesianModel* PRMModelCharge2;

SVMModel* PValueSVM = NULL;
float g_SVMToPValueMin;
int g_PValueBinCount;
float* g_SVMToPValue;

SVMModel* CCModel1SVM = NULL;
SVMModel* CCModel2SVM = NULL;

SVMModel* MQModel2SVM = NULL;
SVMModel* MQModel3SVM = NULL;


//SVMModel* PValueSVM = NULL;
extern float g_CutScores[];
extern float g_BAbsSkew[];
extern float g_YAbsSkew[];
extern float g_BSkew[];
extern float g_YSkew[];
extern float g_BIntensity[];
extern float g_YIntensity[];

float GetPValue(float MQScore) 
{
    int Bin;
    //
    Bin = rint((MQScore - g_SVMToPValueMin)*10);
    Bin = max(Bin, 0);
    Bin = min(Bin, g_PValueBinCount - 1);
    return g_SVMToPValue[Bin];
}

// Given a model and an array of feature-values, perform SVM classification.  
float SVMClassify(SVMModel* Model, float* Coords, int PreScaled)
{
    SupportVector* Vector;
    double Total = 0;
    double InnerProduct;
    int CoordIndex;
    double Diff;
    double ScaledCoords[64];
    //

    if (PreScaled)
    {
        for (CoordIndex = 0; CoordIndex < Model->Coords; CoordIndex++)
        {
            ScaledCoords[CoordIndex] = Coords[CoordIndex];
        }
    }
    else
    {
        // Scale coordinates to the range [-1, 1] based upon the extrema in the model:
        for (CoordIndex = 0; CoordIndex < Model->Coords; CoordIndex++)
        {
            ScaledCoords[CoordIndex] = (Coords[CoordIndex] - Model->ScaleMin[CoordIndex]) / Model->ScaleSize[CoordIndex] - 1.0;
            ScaledCoords[CoordIndex] = min(1, max(-1, ScaledCoords[CoordIndex]));
        }
    }

    // Compute the SVM value by taking weighted inner products with the support vectors ('border points')
    for (Vector = Model->FirstVector; Vector; Vector = Vector->Next)
    {
        InnerProduct = 0;
        for (CoordIndex = 0; CoordIndex < Model->Coords; CoordIndex++)
        {
            Diff = (ScaledCoords[CoordIndex] - Vector->Coords[CoordIndex]);
            InnerProduct += Diff * Diff;
        }
        InnerProduct = exp(-Model->Gamma * InnerProduct);
        Total += Vector->Weight * InnerProduct;
    }
    Total -= Model->Rho;
    return (float)Total;
}

// Free an SVMModel instance, including its list of vectors.
void FreeSVMModel(SVMModel* Model)
{
    SupportVector* Vector;
    SupportVector* Prev = NULL;
    //printf("Free SVM model.\n");
    if (Model)
    {
        for (Vector = Model->FirstVector; Vector; Vector = Vector->Next)
        {
            SafeFree(Prev);
            Prev = Vector;
        }
        SafeFree(Prev);
        SafeFree(Model);
    }
}

// Free all loaded SVM models.
void FreeSVMModels()
{
    FreeSVMModel(PValueSVM);
    FreeSVMModel(MQModel2SVM);
    FreeSVMModel(MQModel3SVM);
}

void InitPValueSVM()
{
    char FilePath[2048];

    // NEW models:
    if (!MQModel2SVM)
    {
        sprintf(FilePath, "%s%s.model", GlobalOptions->ResourceDir, "MQScoreSVM2");
        MQModel2SVM = ReadSVMModel(FilePath);
        sprintf(FilePath, "%s%s.range", GlobalOptions->ResourceDir, "MQScoreSVM2");
        ReadSVMScaling(MQModel2SVM, FilePath);
    }

    if (!MQModel3SVM)
    {
        sprintf(FilePath, "%s%s.model", GlobalOptions->ResourceDir, "MQScoreSVM3");
        MQModel3SVM = ReadSVMModel(FilePath);
        sprintf(FilePath, "%s%s.range", GlobalOptions->ResourceDir, "MQScoreSVM3");
        ReadSVMScaling(MQModel3SVM, FilePath);
    }
}

float SVMComputeMQScore(MSSpectrum* Spectrum, Peptide* Match, float* MQFeatures)
{
    SVMModel* Model;
    float Score;

    if (Spectrum->Charge < 3)
    {
        Model = MQModel2SVM;
    }
    else
    {
        Model = MQModel3SVM;
    }
    if (!Model)
    {
        return 0.0;
    }
    Score = SVMClassify(Model, MQFeatures, 0);
    Score = GetPenalizedScore(Spectrum, Match, Score);
    return Score;
}

int ReadSVMScalingCallback(int LineNumber, int FilePos, char* LineBuffer, void* UserData)
{
    SVMModel* Model;
    char* Value;
    int CoordIndex;
    float ScaleMin;
    float ScaleMax;
    //
    Model = (SVMModel*)UserData;
    //CoordIndex = LineNumber - 3;
    Value = strtok(LineBuffer, " \r\n\t\0");
    if (!Value)
    {
        return 1;
    }
    CoordIndex = atoi(Value) - 1;
    if (CoordIndex < 0)
    {
        return 1;
    }
    Value = strtok(NULL, " \r\n\t\0");
    if (!Value)
    {
        return 1;
    }
    ScaleMin = (float)atof(Value);
    Value = strtok(NULL, " \r\n\t\0");
    if (!Value)
    {
        return 1;
    }
    ScaleMax = (float)atof(Value);
    if (ScaleMax <= ScaleMin)
    {
        REPORT_ERROR(0);
    }
    Model->ScaleMin[CoordIndex] = ScaleMin;
    Model->ScaleMax[CoordIndex] = ScaleMax;
    Model->ScaleSize[CoordIndex] = (Model->ScaleMax[CoordIndex] - Model->ScaleMin[CoordIndex]) / 2.0;
    return 1;
}

// Read feature extrema (for scaling) for an SVM model.
void ReadSVMScaling(SVMModel* Model, char* ScaleFileName)
{
    FILE* File;
    //
    File = fopen(ScaleFileName, "r");
    if (!File)
    {
        REPORT_ERROR_S(8, ScaleFileName);
        return;
    }
    ParseFileByLines(File, ReadSVMScalingCallback, Model, 0);
    fclose(File);
}

typedef struct SVMParseInfo
{
    SVMModel* Model;
    int InVectors;
} SVMParseInfo;


int ReadSVMModelCallback(int LineNumber, int FilePos, char* LineBuffer, void* UserData)
{
    SVMModel* Model;
    SVMParseInfo* Info;
    SupportVector* Vector;
    char* Command;
    char* Value;
    int CoordIndex;
    char* CoordIndexStr;

    //
    Info = (SVMParseInfo*)UserData;
    Model = Info->Model;
    // Process either a support vector line, or a header command.
    if (Info->InVectors)
    {
        Vector = (SupportVector*)calloc(sizeof(SupportVector), 1);
        // Weight, then a list of values.
        Value = strtok(LineBuffer, " \r\n\t");
        if (!Value)
        {
            printf("* Critical error: strtok failed in ReadSVMModel\n");
            return 0;
        }
        Vector->Weight = atof(Value);
        CoordIndex = 0;
        while (1)
        {
            Value = strtok(NULL, " \r\n\t");
            if (!Value)
            {
                break;
            }
            CoordIndexStr = Value;
            while (*Value != ':' && *Value)
            {
                Value++;
            }
            *Value = '\0';
            CoordIndex = atoi(CoordIndexStr) - 1;
            Value++;
            if (CoordIndex >= SUPPORT_VECTOR_LENGTH)
            {
                printf("* Error: SVM vector too long!\n");
                break;
            }
            Vector->Coords[CoordIndex] = atof(Value);
            CoordIndex++;
        }
        Model->Coords = CoordIndex;
        if (Model->LastVector)
        {
            Model->LastVector->Next = Vector;
            Model->LastVector = Vector;
        }
        else
        {
            Model->FirstVector = Vector;
            Model->LastVector = Vector;
        }
    }
    else
    {
        // A header line.  We pay attention to parameters "gamma" and "rho".  The line "sv" marks the end of 
        // the header, and the start of the support vectors.
        Command = strtok(LineBuffer, " ");
        Value = strtok(NULL, " ");
        // First, handle commands that take no arguments:
        if (!CompareStrings(Command, "sv"))
        {
            Info->InVectors = 1;
        }
        else
        {
            // The remaining commands have a mandatory argument:
            if (!Value)
            {
                printf("* Invalid command line in ReadSVMModel\n");
                return 0;
            }
            if (!CompareStrings(Command, "gamma"))
            {
                Model->Gamma = atof(Value);
            }
            if (!CompareStrings(Command, "rho"))
            {
                Model->Rho = atof(Value);
            }
        }
    }
    return 1;
}

// Read an SVM model from a .model file.
SVMModel* ReadSVMModel(char* ModelFileName)
{
    SVMModel* Model;
    FILE* File;
    SVMParseInfo Info;
    //
    //printf("Reading SVM model.\n");
    Model = (SVMModel*)calloc(sizeof(SVMModel), 1);
    File = fopen(ModelFileName, "r");
    if (!File)
    {
        REPORT_ERROR_S(8, ModelFileName);
        return NULL;
    }
    Info.Model = Model;
    Info.InVectors = 0;
    ParseFileByLines(File, ReadSVMModelCallback, &Info, 0);   
    return Model;
}

void TestPValue(char* FeatureVectorFileName)
{
    FILE* FeatureVectorFile;
    int* HistogramFalse;
    int* HistogramTrue;
    float Coords[32];
    char* ValueString;
    int TrueFlag;
    int FeatureIndex;
    float Result;
    int HistogramBin;
    FILE* OutputFile;
    int FalseCount = 0;
    int TrueCount = 0;
    int TrueCumulative = 0;
    int FalseCumulative = 0;
    int BufferEnd = 0;
    int BufferPos = 0;
    int BytesRead;
    int BytesToRead;
    char* Buffer;
    char* LineBuffer;
    int LineNumber = 0;
    char* FieldString;
    char TextBuffer[BUFFER_SIZE * 2];
    //char* ValueString;
    //
    FeatureVectorFile = fopen(FeatureVectorFileName, "r");
    Buffer = (char*)malloc(sizeof(char) * 10240);
    LineBuffer = (char*)malloc(sizeof(char)*MAX_LINE_LENGTH);

    HistogramFalse = (int*)calloc(sizeof(int), 1000);
    HistogramTrue = (int*)calloc(sizeof(int), 1000);
    OutputFile = fopen("PValueTest.txt", "w");
    InitPValueSVM();
    while (1)
    {
        BytesToRead = BUFFER_SIZE - BufferEnd;
        BytesRead = ReadBinary(TextBuffer + BufferEnd, sizeof(char), BytesToRead, FeatureVectorFile);
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
        // Ok, it's a feature line.  Split into pieces...
        memset(Coords, 0, sizeof(float)*32);
        ValueString = strtok(LineBuffer, WHITESPACE);
        TrueFlag = atoi(ValueString);
        fprintf(OutputFile, "%d\t", TrueFlag);
        if (TrueFlag < 0)
        {
            TrueFlag = 0;
        }
        FeatureIndex = 0;
        while (1)
        {
            FieldString = strtok(NULL, WHITESPACE);
            ValueString = FieldString;
            if (!ValueString)
            {
                break;
            }
            while (*ValueString!=':')
            {
                ValueString++;
            }
            *ValueString = '\0';
            FeatureIndex = atoi(FieldString) - 1;
            ValueString++;
            Coords[FeatureIndex++] = (float)atof(ValueString);
            fprintf(OutputFile, "%s\t", ValueString);
        }
        Result = SVMClassify(PValueSVM, Coords, 1);
        fprintf(OutputFile, "%.4f\n", Result);
        HistogramBin = rint(Result*10) + 300;
        HistogramBin = max(0, min(999, HistogramBin));
        if (TrueFlag)
        {
            HistogramTrue[HistogramBin]++;
            TrueCount++;
        }
        else
        {
            HistogramFalse[HistogramBin]++;
            FalseCount++;
        }
    }
    FalseCount = max(FalseCount, 1); // avoid dividing by zero
    TrueCount = max(TrueCount, 1); // avoid dividing by zero
    
    for (HistogramBin = 0; HistogramBin < 1000; HistogramBin++)
    {
        TrueCumulative += HistogramTrue[HistogramBin];
        FalseCumulative += HistogramFalse[HistogramBin];
        fprintf(OutputFile, "%d\t%.2f\t%.2f\t%.2f\t\n",
            HistogramBin, (HistogramBin - 300) / 10.0,
            100*TrueCumulative/(float)TrueCount,
            100*FalseCumulative/(float)FalseCount);
    }
}

float LDAClassify(float* Features)
{
    float ScaledFeatures[6];
    double FeatureMin[] = {-1.88, 0, 0, 0};
    double FeatureMax[] = {3.81, 1, 1, 2};
    float HalfRange;
    int X;
    int Y;
    int FeatureCount;
    static double* CovInv[6];
    double* MeanVectorTrue;
    double* MeanVectorFalse;
    double SubProdTrue;
    double SubProdFalse;
    float ProdTemp[6];
    float ProdTrue;
    float ProdFalse;

    // Constants for TRYPTIC scoring:
    double TCovInvA[] = {6.037,-8.996,-7.351,-0.283};
    double TCovInvB[] = {-8.996,51.379,-3.536,2.428};
    double TCovInvC[] = {-7.351,-3.536,28.577,-0.271};
    double TCovInvD[] = {-0.283,2.428,-0.271,2.382};
    double TMeanVectorTrue[] = {2.048,0.022,0.187,0.622};
    double TMeanVectorFalse[] = {-0.352,-0.668,-0.629,0.102};
    double TSubProdTrue = (float)-10.052;
    double TSubProdFalse = (float)-12.146;

    // Constants for NON-TRYPTIC scoring:
    double NTCovInvA[] = {6.003,-8.708,-7.383};
    double NTCovInvB[] = {-8.708,48.904,-3.259};
    double NTCovInvC[] = {-7.383,-3.259,28.546};
    double NTMeanVectorTrue[] = {2.048,0.022,0.187};
    double NTMeanVectorFalse[] = {-0.352,-0.668,-0.629};
    double NTSubProdTrue = (float)-9.880;
    double NTSubProdFalse = (float)-11.888;

    // Choose the feature-set by digest type.
    if (GlobalOptions->DigestType == DIGEST_TYPE_TRYPSIN)
    {
        MeanVectorTrue = TMeanVectorTrue;
        MeanVectorFalse = TMeanVectorFalse;
        SubProdTrue = TSubProdTrue;
        SubProdFalse = TSubProdFalse;
        CovInv[0] = TCovInvA;
        CovInv[1] = TCovInvB;
        CovInv[2] = TCovInvC;
        CovInv[3] = TCovInvD;
        FeatureCount = 4;
    }
    else
    {
        MeanVectorTrue = NTMeanVectorTrue;
        MeanVectorFalse = NTMeanVectorFalse;
        SubProdTrue = NTSubProdTrue;
        SubProdFalse = NTSubProdFalse;
        CovInv[0] = NTCovInvA;
        CovInv[1] = NTCovInvB;
        CovInv[2] = NTCovInvC;
        FeatureCount = 3;
    }
    // Scale the features into [-1, 1]:
    for (X = 0; X < FeatureCount; X++)
    {
        HalfRange = (float)((FeatureMax[X] - FeatureMin[X]) / 2.0);
        ScaledFeatures[X] = (float)((Features[X] - FeatureMin[X]) / HalfRange - 1.0);
    }
    // Compute the product of the inverse covariance matrix with our feature vector:
    for (X = 0; X < FeatureCount; X++)
    {
        ProdTemp[X] = 0;
        for (Y = 0; Y < FeatureCount; Y++)
        {
            ProdTemp[X] += (float)(ScaledFeatures[Y] * CovInv[X][Y]);
        }
    }
    // Compute u0 * C-1 * X and u1 * C-1 * X
    ProdTrue = 0;
    ProdFalse = 0;
    for (X = 0; X < FeatureCount; X++)
    {
        ProdFalse += (float)(MeanVectorFalse[X] * ProdTemp[X]);
        ProdTrue += (float)(MeanVectorTrue[X] * ProdTemp[X]);
    }
    ProdTrue += (float)SubProdTrue;
    ProdFalse += (float)SubProdFalse;
    //printf("%.2f\t%.2f\t%.2f\t\n", (ProdTrue - ProdFalse), ProdTrue, ProdFalse);
    return (ProdTrue - ProdFalse);
}

void LoadCCModelSVM(int ForceRefresh)
{
    char FilePath[2048];
    if (CCModel1SVM)
    {
        if (ForceRefresh)
        {
            FreeSVMModel(CCModel1SVM);
            FreeSVMModel(CCModel2SVM);
        }
        else
        {
            return;
        }
    }

	if(GlobalOptions->PhosphorylationFlag)
	{ //separate model for charge 2 only.  not enuf training data for charge 1
		sprintf(FilePath, "%s%s", GlobalOptions->ResourceDir, "CCSVM1.model");
		CCModel1SVM = ReadSVMModel(FilePath);
		sprintf(FilePath, "%s%s", GlobalOptions->ResourceDir, "CCSVM1.range");
		ReadSVMScaling(CCModel1SVM, FilePath);
		sprintf(FilePath, "%s%s", GlobalOptions->ResourceDir, "CCSVM2Phos.model");
		CCModel2SVM = ReadSVMModel(FilePath);
		sprintf(FilePath, "%s%s", GlobalOptions->ResourceDir, "CCSVM2Phos.range");
		ReadSVMScaling(CCModel2SVM, FilePath);
	}
	else
	{
		sprintf(FilePath, "%s%s", GlobalOptions->ResourceDir, "CCSVM1.model");
		CCModel1SVM = ReadSVMModel(FilePath);
		sprintf(FilePath, "%s%s", GlobalOptions->ResourceDir, "CCSVM1.range");
		ReadSVMScaling(CCModel1SVM, FilePath);
		sprintf(FilePath, "%s%s", GlobalOptions->ResourceDir, "CCSVM2.model");
		CCModel2SVM = ReadSVMModel(FilePath);
		sprintf(FilePath, "%s%s", GlobalOptions->ResourceDir, "CCSVM2.range");
		ReadSVMScaling(CCModel2SVM, FilePath);
	}
}

void FreeCCModelSVM()
{
    FreeSVMModel(CCModel1SVM);
    CCModel1SVM = NULL;
    FreeSVMModel(CCModel2SVM);
    CCModel2SVM = NULL;

}
