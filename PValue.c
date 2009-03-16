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
#include <math.h>
#include "PValue.h"
#include "Score.h"
#include "Errors.h"

double GammaCof[] = {76.18009172947146, -86.50532032941677,
    24.01409824083091, -1.231739572450155, 
    0.1208650973866179e-2, -0.5395239384952e-5};

double Gamma(double Z)
{
    double X;
    double Y;
    double Temp;
    double Ser;
    int J;
    //////////
    X = Z;
    Y = Z;
    Temp = X + 5.5;
    Temp -= (X + 0.5) * log(Temp);
    Ser = 1.000000000190015;
    for (J = 0; J < 6; J++)
    {
        Y += 1;
        Ser += GammaCof[J] / Y;
    }
    Z = -Temp + log(2.5066282746310005 * Ser / X);
    return exp(Z);
}

#define F_BIN_MAX 511
#define F_BIN_OFFSET 40 //used because Fscores are negative, but array indexes cannot be

void DebugPrintPValueCurve(char* FileName, int* FScoreHistogram, double* OddsTrue)
{
    FILE* PValueFile;
    int FBin;
    float X;
    //
    PValueFile = fopen(FileName, "wb");
    if (!PValueFile)
    {
        printf("** Not debug-printing the p-value curve.\n");
        return;
    }
    for (FBin = 0; FBin <= F_BIN_MAX; FBin++)
    {
        X = (FBin - F_BIN_OFFSET) / (float)10.0;
        fprintf(PValueFile, "%d\t%d\t%.3f\t%.3f\t\n", FBin, FScoreHistogram[FBin], X, 1.0 - OddsTrue[FBin]);
    }
    fclose(PValueFile);
}

#define EM_CYCLE_COUNT 300
#define SQRT_2_PI (float)2.5066282746310002

#define MAX_BITS 30
#define PROCESS_COMPUTE_MEAN_DELTA 0
#define PROCESS_INITIALIZE_SCORE_HISTOGRAM 1
#define PROCESS_WRITE_PVALUES 2

#define BIT_INDEX_CHARGE 4
#define BIT_INDEX_MQSCORE 5
#define BIT_INDEX_PVALUE 13
#define BIT_INDEX_FSCORE 14
#define BIT_INDEX_DELTA_SCORE 16
#define DEFAULT_MQ_SCORE_WEIGHT (float)0.3
#define DEFAULT_DELTA_SCORE_WEIGHT (float)1.5
#define BLIND_MQ_SCORE_WEIGHT (float)0.3
#define BLIND_DELTA_SCORE_WEIGHT (float)2.25

float MQScoreWeight;
float DeltaScoreWeight;

typedef struct PValueParseInfo
{
    FILE* OutputFile;
    int TotalMatches;
    float MeanDeltaScore;
    int Action;
    char CurrentSpectrum[512];
} PValueParseInfo;

typedef struct PValueInfo
{
    float MeanDeltaScore2;
    float MeanDeltaScore3;
    int TotalMatches2;
    int TotalMatches3;
    char* FinalOutputPath;
    int FScoreHistogram2[F_BIN_MAX + 1];
    double OddsTrue2[F_BIN_MAX + 1];
    int FScoreHistogram3[F_BIN_MAX + 1];
    double OddsTrue3[F_BIN_MAX + 1];
    PValueParseInfo* ParseInfo;
} PValueInfo;

typedef struct PValueModel
{
    double MeanTrue;
    double VarianceTrue;
    double MeanFalse;
    double VarianceFalse;
    double PriorProbabilityTrue;
    double ThetaFalse;
    double KFalse;
    double GammaOffset;
    double StdDevTrue;
    double GammaDemonFalse;
    double CountTrue;
    double CountFalse;
} PValueModel;

// CustomTok is a variant of strtok: It returns once for every occurrence of a delimiter,
// rather than once for every contiguous block of delimiters.  Why the difference?  
// When processing tab-delimited files, we want to properly handle empty columns (corresponding to 
// occurrences of \t\t in the text).
static char* CustomTokEnd;
static char* CustomTokNext;
char* CustomTok(char* Buffer, char* Delimiters)
{
    char* CheckPos;
    char* StringStart;
    char* CheckDelimiter;
    //
    if (Buffer)
    {
        CustomTokEnd = Buffer + strlen(Buffer);
        CheckPos = Buffer;
        StringStart = Buffer;
    }
    else
    {
        CheckPos = CustomTokNext;
        StringStart = CustomTokNext;
    }

    // If we're out of bits, then say so:
    if (CheckPos >= CustomTokEnd)
    {
        return NULL;
    }
    // Scan forward until you see a delimiter, or until the end of the string:
    for (; CheckPos < CustomTokEnd; CheckPos++)
    {
        for (CheckDelimiter = Delimiters; *CheckDelimiter; CheckDelimiter++)
        {
            if (*CheckPos == *CheckDelimiter)
            {
                *CheckPos = '\0';
                CustomTokNext = CheckPos + 1;
                return StringStart;
            }
        }
    }
    // We didn't see the delimiter.  
    CustomTokNext = CustomTokEnd;
    return StringStart;
}

int PValueProcessResultsFileLine(int LineNumber, int FilePos, char* LineBuffer, void* UserData)
{
    PValueInfo* Info;
    PValueParseInfo* ParseInfo;
    char* Bits[MAX_BITS];
    int BitCount;
    char* Bit;
    char Spectrum[512];
    int TopMatchFlag;
    float MQScore;
    float DeltaScore;
    float FScore;
    int FBin;
    float PValue;
    char PValueBuffer[256];
    char FScoreBuffer[256];
    int BitIndex;
    int Charge;
    
    //
    Info = (PValueInfo*)UserData;
    ParseInfo = Info->ParseInfo;
    if (ParseInfo->Action == PROCESS_WRITE_PVALUES)
    {
        INSPECT_ASSERT(ParseInfo->OutputFile);
    }

    // Handle comments:
    if (LineBuffer[0] == '#')
    {
        // Comment lines are written out as-is:
        if (ParseInfo->Action == PROCESS_WRITE_PVALUES)
        {
            fprintf(ParseInfo->OutputFile, "%s\n", LineBuffer);
        }
        return 1;
    }
    // Split the line into tab-delimited bits:
    Bit = CustomTok(LineBuffer, "\t");
    Bits[0] = Bit;
    BitCount = 1;
    while (1)
    {
        Bit = CustomTok(NULL, "\t");
        if (!Bit)
        {
            break;
        }
        Bits[BitCount] = Bit;
        BitCount++;
        if (BitCount >= MAX_BITS)
        {
            break;
        }
    }

    // If we don't have enough tab-bits, then this isn't a valid line, and we skip it:
    if (BitCount < BIT_INDEX_PVALUE + 1)
    {
        return 1;
    }

    // Note whether this is the top-scoring match for the spectrum:
    sprintf(Spectrum, "%256s%50s", Bits[0], Bits[1]);
    if (strcmp(Spectrum, ParseInfo->CurrentSpectrum))
    {
        TopMatchFlag = 1;
        strncpy(ParseInfo->CurrentSpectrum, Spectrum, 512);
    }
    else
    {
        TopMatchFlag = 0;
    }
    Charge = atoi(Bits[BIT_INDEX_CHARGE]);
    
    // Now take various actions:
    switch (ParseInfo->Action)
    {
    case PROCESS_COMPUTE_MEAN_DELTA:
        if (TopMatchFlag)
        {
            if (Charge < 3)
            {
                Info->MeanDeltaScore2 += (float)atof(Bits[BIT_INDEX_DELTA_SCORE]);
                Info->TotalMatches2++;
            }
            else
            {
                Info->MeanDeltaScore3 += (float)atof(Bits[BIT_INDEX_DELTA_SCORE]);
                Info->TotalMatches3++;
            }
        }
        break;
    case PROCESS_INITIALIZE_SCORE_HISTOGRAM:
        if (TopMatchFlag)
        {
            MQScore = (float)atof(Bits[BIT_INDEX_MQSCORE]);
            DeltaScore = (float)atof(Bits[BIT_INDEX_DELTA_SCORE]) / (Charge < 3 ? Info->MeanDeltaScore2 : Info->MeanDeltaScore3);
            FBin = (int)floor((10 * (MQScoreWeight * MQScore + DeltaScoreWeight * DeltaScore)));
            FBin = min(max(FBin + F_BIN_OFFSET, 0), F_BIN_MAX);
            if (Charge < 3)
            {
                Info->FScoreHistogram2[FBin]++;
            }
            else
            {
                Info->FScoreHistogram3[FBin]++;
            }
        }
        break;
    case PROCESS_WRITE_PVALUES:
        MQScore = (float)atof(Bits[BIT_INDEX_MQSCORE]);
        DeltaScore = (float)atof(Bits[BIT_INDEX_DELTA_SCORE]) / (Charge < 3 ? Info->MeanDeltaScore2 : Info->MeanDeltaScore3);
        FScore = (MQScoreWeight * MQScore) + (DeltaScoreWeight * DeltaScore);
        sprintf(FScoreBuffer, "%.5f", FScore);
        Bits[BIT_INDEX_FSCORE] = FScoreBuffer;
        FBin = (int)(10 * FScore);
        FBin = min(max(FBin + F_BIN_OFFSET, 0), F_BIN_MAX);
        if (Charge < 3)
        {
            PValue = (float)(1.0 - Info->OddsTrue2[FBin]);
        }
        else
        {
            PValue = (float)(1.0 - Info->OddsTrue3[FBin]);
        }
        sprintf(PValueBuffer, "%.5f", PValue);
        Bits[BIT_INDEX_PVALUE] = PValueBuffer;
        for (BitIndex = 0; BitIndex < BitCount; BitIndex++)
        {
            fprintf(ParseInfo->OutputFile, "%s\t", Bits[BitIndex]);
        }
        fprintf(ParseInfo->OutputFile, "\n");
        break;
    default:
        printf("** Error: Unknown action '%d' in ProcessResultsFile!\n", ParseInfo->Action);
        return 0;
    }
    return 1;
}

void ProcessResultsFile(PValueInfo* Info, char* FilePath, int Action)
{
    FILE* File;
    PValueParseInfo ParseInfo;
    //
    Info->ParseInfo = &ParseInfo;
    File = fopen(FilePath, "rb");
    if (!File)
    {
        REPORT_ERROR_S(8, FilePath);
        return;
    }
    ParseInfo.Action = Action;
    if (Action == PROCESS_WRITE_PVALUES)
    {
        ParseInfo.OutputFile = fopen(Info->FinalOutputPath, "wb");
        if (!ParseInfo.OutputFile)
        {
            REPORT_ERROR_S(8, Info->FinalOutputPath);
            return;
        }
    }
    ParseFileByLines(File, PValueProcessResultsFileLine, Info, 1);
    fclose(File);
}

#define DEFAULT_MEAN_TRUE 4.48
#define DEFAULT_VARIANCE_TRUE 1.50
#define DEFAULT_MEAN_FALSE 0.19
#define DEFAULT_VARIANCE_FALSE 0.18
#define DEFAULT_PRIOR_PROBABILITY 0.25

PValueModel* InitPValueModel(int Charge3Flag)
{
    PValueModel* Model;
    //
    Model = (PValueModel*)calloc(1, sizeof(PValueModel));
    Model->MeanTrue = DEFAULT_MEAN_TRUE;
    Model->VarianceTrue = DEFAULT_VARIANCE_TRUE;
    Model->MeanFalse = DEFAULT_MEAN_FALSE;
    Model->VarianceFalse = DEFAULT_VARIANCE_FALSE;
    Model->PriorProbabilityTrue = DEFAULT_PRIOR_PROBABILITY;
    Model->GammaOffset = 0;
    Model->GammaOffset = max(Model->GammaOffset, -Model->MeanFalse + 0.1);
    Model->ThetaFalse = Model->VarianceFalse / (Model->MeanFalse + Model->GammaOffset);
    Model->KFalse = (Model->MeanFalse + Model->GammaOffset) / Model->ThetaFalse;
    Model->GammaDemonFalse = pow(Model->ThetaFalse, Model->KFalse) * Gamma(Model->KFalse);
    Model->StdDevTrue = sqrt(Model->VarianceTrue);
    return Model;
}

int FitPValueMixtureModel(PValueInfo* Info, PValueModel* Model, int Charge3Flag)
{
    int FBin;
    int* FScoreHistogram;
    double* OddsTrue;
    int LowestFScoreBin = -1;
    //double GammaOffset;
    int TotalMatches;
    int EMCycle;
    double X;
    int Count;
    int MaxBinPopulated = 0;
    float LowestFScore = 0;
    int CurveFitComplete;
    double Pow;
    double GX;
    double TrueNormal;
    double FalseGamma;
    //

    if (Charge3Flag)
    {
        FScoreHistogram = Info->FScoreHistogram3;
        OddsTrue = Info->OddsTrue3;
        TotalMatches = Info->TotalMatches3;
    }
    else
    {
        FScoreHistogram = Info->FScoreHistogram2;
        OddsTrue = Info->OddsTrue2;
        TotalMatches = Info->TotalMatches2;
    }

    // Note the highest score-bin that has any entries at all:
    for (FBin = 0; FBin <= F_BIN_MAX; FBin++)
    {
        if (FScoreHistogram[FBin] > 0 && LowestFScoreBin < 0)
        {
            LowestFScoreBin = FBin;
        }
        if (FScoreHistogram[FBin])
        {
            MaxBinPopulated = FBin;
        }
    }

    // Convert the lowest F-score bin# into the corresponding score:
    LowestFScore = (LowestFScoreBin - F_BIN_OFFSET) / (float)10.0;
    Model->GammaOffset = 0.0;
    if (LowestFScore <= 0)
    {
        Model->GammaOffset = max(Model->GammaOffset, -LowestFScore + 0.1);
    }
    if (Model->MeanFalse <= 0)
    {
        Model->GammaOffset = max(Model->GammaOffset, -Model->MeanFalse + 0.1);
    }

    ////////////////////////////////////////////////////////////////////
    // Fit the mixture model, using a gamma distribution for false match f-scores and a 
    // normal distribution for true match f-scores.
    Model->ThetaFalse = Model->VarianceFalse / (Model->MeanFalse + Model->GammaOffset);
    Model->KFalse = (Model->MeanFalse + Model->GammaOffset) / Model->ThetaFalse;
    Model->GammaDemonFalse = pow(Model->ThetaFalse, Model->KFalse) * Gamma(Model->KFalse);
    Model->StdDevTrue = sqrt(Model->VarianceTrue);
    if (TotalMatches < 200)
    {
        REPORT_WARNING_I(10, TotalMatches);
        CurveFitComplete = 0;
    }
    else
    {
        for (EMCycle = 0; EMCycle < EM_CYCLE_COUNT; EMCycle++)
        {
            // For each bin, compute the probability that it's true:
            for (FBin = 0; FBin <= F_BIN_MAX; FBin++)
            {
                // After the last histogram entry, just inherit the last true-probability.
                if (FBin > MaxBinPopulated)
                {
                    OddsTrue[FBin] = OddsTrue[FBin - 1];
                    continue;
                }
                X = (FBin - F_BIN_OFFSET) / 10.0;
                Pow = (X - Model->MeanTrue);
                Pow = - (Pow * Pow / (2 * Model->VarianceTrue));
                TrueNormal = exp(Pow) / (Model->StdDevTrue * SQRT_2_PI);
                GX = max(0.01, X + Model->GammaOffset);
                FalseGamma = pow(GX, Model->KFalse - 1) * exp(-GX / Model->ThetaFalse) / Model->GammaDemonFalse;
                // Avoid underflow:
                if (TrueNormal < 0.00001)
                {
                    if (X > 5)
                    {
                        OddsTrue[FBin] = (float)0.99999;
                    }
                    else
                    {
                        OddsTrue[FBin] = (float)0.0;
                    }
                }
                else
                {
                    OddsTrue[FBin] = (TrueNormal * Model->PriorProbabilityTrue) / (TrueNormal * Model->PriorProbabilityTrue + FalseGamma * (1 - Model->PriorProbabilityTrue));
                }
                
                //printf("%.8f\t%.8f\t%.8f\t%.8f\n", X, TrueNormal, FalseGamma, OddsTrue[FBin]);
                // Because the left tail of the normal distribution falls off slowly, the value of OddsTrue can be
                // high for negative values.  Cap it.
                if (FBin < F_BIN_OFFSET)
                {
                    OddsTrue[FBin] = min(OddsTrue[FBin], 1.0 / (F_BIN_OFFSET - FBin));
                }
            }
            /////////////////////////////////////////////////
            // Compute the mean of the true and the false distributions:
            Model->CountTrue = 0;
            Model->MeanTrue = 0;
            Model->CountFalse = 0;
            Model->MeanFalse = 0;
            for (FBin = 0; FBin <= F_BIN_MAX; FBin++)
            {
                X = (FBin - F_BIN_OFFSET) / 10.0;
                Count = FScoreHistogram[FBin];
                Model->MeanTrue += X * OddsTrue[FBin] * Count;
                Model->CountTrue += OddsTrue[FBin] * Count;
                Model->MeanFalse += X * (1.0 - OddsTrue[FBin]) * Count;
                Model->CountFalse += (1.0 - OddsTrue[FBin]) * Count;
            }
            Model->MeanTrue /= Model->CountTrue;
            Model->MeanFalse /= Model->CountFalse;
            Model->PriorProbabilityTrue = Model->CountTrue / (Model->CountTrue + Model->CountFalse);
            /////////////////////////////////////////////////
            // Compute the variance of the true and the false distributions:
            Model->VarianceTrue = 0;
            Model->VarianceFalse = 0;
            for (FBin = 0; FBin <= F_BIN_MAX; FBin++)
            {
                X = (FBin - F_BIN_OFFSET) / 10.0;
                Count = FScoreHistogram[FBin];
                Model->VarianceTrue += (X - Model->MeanTrue) * (X - Model->MeanTrue) * Count * OddsTrue[FBin];
                Model->VarianceFalse += (X - Model->MeanFalse) * (X - Model->MeanFalse) * Count * (1.0 - OddsTrue[FBin]);
            }
            Model->VarianceTrue /= Model->CountTrue;
            Model->StdDevTrue = sqrt(Model->VarianceTrue);
            Model->VarianceFalse /= Model->CountFalse;
            // Recompute other distribution parameters:
            Model->ThetaFalse = Model->VarianceFalse / (Model->MeanFalse + Model->GammaOffset);
            Model->KFalse = (Model->MeanFalse + Model->GammaOffset) / Model->ThetaFalse;
            Model->GammaDemonFalse = pow(Model->ThetaFalse, Model->KFalse) * Gamma(Model->KFalse);
            //printf("Cycle %d:\n", EMCycle);
            //printf("True: Count %.4f mean %.4f variance %.4f\n", CountTrue, MeanTrue, VarianceTrue);
            //printf("False: Count %.4f mean %.4f variance %.4f\n", CountFalse, MeanFalse, VarianceFalse);
            //printf("Prior probability true: %.4f\n", PriorProbabilityTrue);
        } // E-M cycle loop
        CurveFitComplete = 1;
    }

    ///////////////////////////////////////
    // Check to make sure the distribution is sensible.  If curve-fitting failed
    // due to underflow/overflow, then fall back to the default curve:
    if (Model->GammaDemonFalse <= 0 || Model->KFalse <= 0)
    {
        printf("** Error fitting p-value distribution; using default.  Consider running PValue.py\n");
        CurveFitComplete = 0;
    }

    if (!CurveFitComplete)
    {
        // COPY-PASTA: Fill in the OddsTrue array using all default parameters.
        Model = InitPValueModel(Charge3Flag);

        // For each bin, compute the probability that it's true:
        for (FBin = 0; FBin <= F_BIN_MAX; FBin++)
        {
            // After the last histogram entry, just inherit the last true-probability.
            if (FBin > MaxBinPopulated)
            {
                OddsTrue[FBin] = OddsTrue[FBin - 1];
                continue;
            }
            X = (FBin - F_BIN_OFFSET) / 10.0;
            Pow = (X - Model->MeanTrue);
            Pow = - (Pow * Pow / (2 * Model->VarianceTrue));
            TrueNormal = exp(Pow) / (Model->StdDevTrue * SQRT_2_PI);
            GX = max(0.01, X + Model->GammaOffset);
            FalseGamma = pow(GX, Model->KFalse - 1) * exp(-GX / Model->ThetaFalse) / Model->GammaDemonFalse;
            // Avoid underflow:
            if (TrueNormal < 0.00001)
            {
                if (X > 5)
                {
                    OddsTrue[FBin] = (float)0.99999;
                }
                else
                {
                    OddsTrue[FBin] = (float)0.0;
                }
            }
            else
            {
                OddsTrue[FBin] = (TrueNormal * Model->PriorProbabilityTrue) / (TrueNormal * Model->PriorProbabilityTrue + FalseGamma * (1 - Model->PriorProbabilityTrue));
            }
            
            //printf("%.8f\t%.8f\t%.8f\t%.8f\n", X, TrueNormal, FalseGamma, OddsTrue[FBin]);
            // Because the left tail of the normal distribution falls off slowly, the value of OddsTrue can be
            // high for negative values.  Cap it.
            if (FBin < F_BIN_OFFSET)
            {
                OddsTrue[FBin] = (float)min(OddsTrue[FBin], 1.0 / (F_BIN_OFFSET - FBin));
            }
        }
        // free the temp-model:
        SafeFree(Model);
    } // if curve fit didn't complete
    return 1;
}

// Iterate over all the matches, and get the distribution
// of F-scores; use those to derive p-values.  We compute
// one distribution for charge 1 and 2 spectra, a second 
// distribution for charge 3 spectra.
void CalculatePValues(char* ResultsFilePath, char* FinalOutputPath)
{
    PValueModel* Model2;
    PValueModel* Model3;
    PValueInfo* Info;
    //    
    Model2 = InitPValueModel(0);
    Model3 = InitPValueModel(1);
    Info = (PValueInfo*)calloc(1, sizeof(PValueInfo));

    if (GlobalOptions->RunMode & (RUN_MODE_MUTATION | RUN_MODE_BLIND))
    {
        MQScoreWeight = BLIND_MQ_SCORE_WEIGHT;
        DeltaScoreWeight = BLIND_DELTA_SCORE_WEIGHT;
    }
    else
    {
        MQScoreWeight = DEFAULT_MQ_SCORE_WEIGHT;
        DeltaScoreWeight = DEFAULT_DELTA_SCORE_WEIGHT;
    }

    //////////////////////////////////////////////////////////
    // Compute mean delta-score:
    ProcessResultsFile(Info, ResultsFilePath, PROCESS_COMPUTE_MEAN_DELTA);
    Info->MeanDeltaScore2 /= max(1, Info->TotalMatches2);
    Info->MeanDeltaScore2 = max(Info->MeanDeltaScore2, (float)0.01);
    Info->MeanDeltaScore3 /= max(1, Info->TotalMatches3);
    Info->MeanDeltaScore3 = max(Info->MeanDeltaScore3, (float)0.01);

    //////////////////////////////////////////////////////////
    // Initialze FScoreHistogram:
    memset(Info->FScoreHistogram2, 0, sizeof(int) * (F_BIN_MAX + 1));
    memset(Info->FScoreHistogram3, 0, sizeof(int) * (F_BIN_MAX + 1));
    ProcessResultsFile(Info, ResultsFilePath, PROCESS_INITIALIZE_SCORE_HISTOGRAM);

    //////////////////////////////////////////////////////////
    // Fit the mixture model, populating OddsTrue:
    FitPValueMixtureModel(Info, Model2, 0);
    FitPValueMixtureModel(Info, Model3, 1);

    // Verbose ouptut of the p-value curve:
    // (Disabled in production, especially for the web server!)
    //DebugPrintPValueCurve("PValueCurve2.txt", Info->FScoreHistogram2, Info->OddsTrue2);
    //DebugPrintPValueCurve("PValueCurve3.txt", Info->FScoreHistogram3, Info->OddsTrue3);
    
    // Write the p-values to the final output file:
    Info->FinalOutputPath = FinalOutputPath;
    ProcessResultsFile(Info, ResultsFilePath, PROCESS_WRITE_PVALUES);
    // Now we can erase the intermediate output file:
    remove(GlobalOptions->OutputFileName);
}
