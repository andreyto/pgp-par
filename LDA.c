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

// LDA support functions.

#include "CMemLeak.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "LDA.h"
#include "Utils.h"
#include "Inspect.h"
#include "Errors.h"
#include "Spectrum.h"
#include "Trie.h"
#include "Score.h"

// Global variables:
LDAModel* PMCCharge1LDA = NULL;
LDAModel* PMCCharge2LDA = NULL;
LDAModel* PMCCharge3LDA = NULL;

LDAModel* CCModel1LDA = NULL;
LDAModel* CCModel2LDA = NULL;

LDAModel* MQModel2LDA = NULL;
LDAModel* MQModel3LDA = NULL;


void LoadCCModelLDA(int ForceRefresh)
{
    char FilePath[2048];
    if (CCModel1LDA)
    {
        if (ForceRefresh)
        {
            FreeLDAModel(CCModel1LDA);
            FreeLDAModel(CCModel2LDA);
        }
        else
        {
            return;
        }
    }
    sprintf(FilePath, "%sCCLDA1.model", GlobalOptions->ResourceDir);
    CCModel1LDA = LoadLDAModel(FilePath);
    sprintf(FilePath, "%sCCLDA2.model", GlobalOptions->ResourceDir);
    CCModel2LDA = LoadLDAModel(FilePath);
}

void FreeLDAModels()
{
    FreeLDAModel(PMCCharge1LDA);
    PMCCharge1LDA = NULL;
    FreeLDAModel(PMCCharge2LDA);
    PMCCharge2LDA = NULL;
    FreeLDAModel(PMCCharge3LDA);
    PMCCharge3LDA = NULL;
    FreeLDAModel(MQModel2LDA);
    MQModel2LDA = NULL;
    FreeLDAModel(MQModel3LDA);
    MQModel3LDA = NULL;
}

// Load linear discriminant analysis (LDA) model for parent mass
// correction (PMC).  Special models for phosphorylation searches
void LoadPMCLDA(int ForceLoad)
{
    char FilePath[2048];
    if (PMCCharge1LDA)
    {
        if (ForceLoad)
        {
            FreeLDAModel(PMCCharge1LDA);
            FreeLDAModel(PMCCharge2LDA);
            FreeLDAModel(PMCCharge3LDA);
        }
        else
        {
            return;
        }
    }
    sprintf(FilePath, "%sPMCLDA1.model", GlobalOptions->ResourceDir);
    PMCCharge1LDA = LoadLDAModel(FilePath);
	if (GlobalOptions->PhosphorylationFlag)
	{//Load phosphorylation specific models, only different for charge 2 and 3
		sprintf(FilePath, "%sPMCLDA2Phos.model", GlobalOptions->ResourceDir);
		PMCCharge2LDA = LoadLDAModel(FilePath);
		sprintf(FilePath, "%sPMCLDA3Phos.model", GlobalOptions->ResourceDir);
		PMCCharge3LDA = LoadLDAModel(FilePath);
	}
	else
	{ 
		sprintf(FilePath, "%sPMCLDA2.model", GlobalOptions->ResourceDir);
		PMCCharge2LDA = LoadLDAModel(FilePath);
		sprintf(FilePath, "%sPMCLDA3.model", GlobalOptions->ResourceDir);
		PMCCharge3LDA = LoadLDAModel(FilePath);
	}
}

LDAModel* LoadLDAModel(char* LDAModelFileName)
{
    FILE* File;
    LDAModel* Model;
    float Value;
    int BytesRead;
    //
    File = fopen(LDAModelFileName, "rb");
    if (!File)
    {
        return NULL;
    }
    Model = (LDAModel*)calloc(1, sizeof(LDAModel));
    //ReadBinary(&Value, sizeof(float), 1, File);
    ReadBinary(&Model->FeatureCount, sizeof(int), 1, File);
    assert(Model->FeatureCount >= 1 && Model->FeatureCount < 100);
    Model->ScaledVector = (double*)calloc(Model->FeatureCount, sizeof(double));
    Model->TempProductVector = (double*)calloc(Model->FeatureCount, sizeof(double));
    
    // Read min and max values:
    Model->MinValues = (double*)calloc(Model->FeatureCount, sizeof(double));
    ReadBinary(Model->MinValues, sizeof(double), Model->FeatureCount, File);
    Model->MaxValues = (double*)calloc(Model->FeatureCount, sizeof(double));
    ReadBinary(Model->MaxValues, sizeof(double), Model->FeatureCount, File);
    // Read mean true vector and mean false vector:
    Model->MeanVectorTrue = (double*)calloc(Model->FeatureCount, sizeof(double));
    ReadBinary(Model->MeanVectorTrue, sizeof(double), Model->FeatureCount, File);
    Model->MeanVectorFalse = (double*)calloc(Model->FeatureCount, sizeof(double));
    ReadBinary(Model->MeanVectorFalse, sizeof(double), Model->FeatureCount, File);
    // Read constant ture and constant false:
    ReadBinary(&Model->ConstantTrue, sizeof(double), 1, File);
    ReadBinary(&Model->ConstantFalse, sizeof(double), 1, File);
    // Read inverted covariance matrix:
    Model->CovInv = (double*)calloc(Model->FeatureCount * Model->FeatureCount, sizeof(double));
    ReadBinary(Model->CovInv, sizeof(double), Model->FeatureCount * Model->FeatureCount, File);
    // Verify that we're at EOF:
    BytesRead = ReadBinary(&Value, sizeof(double), 1, File);
    assert(!BytesRead);
    //printf("\nLoading LDA from %s:\n", LDAModelFileName);
    //printf("%d features\n", Model->FeatureCount);
    //printf("MinValues: %.4f...%.4f\n", Model->MinValues[0], Model->MinValues[Model->FeatureCount - 1]);
    //printf("MaxValues: %.4f...%.4f\n", Model->MaxValues[0], Model->MaxValues[Model->FeatureCount - 1]);
    //printf("MeanVectorTrue: %.4f...%.4f\n", Model->MeanVectorTrue[0], Model->MeanVectorTrue[Model->FeatureCount - 1]);
    //printf("MeanVectorFalse: %.4f...%.4f\n", Model->MeanVectorFalse[0], Model->MeanVectorFalse[Model->FeatureCount - 1]);
    //printf("CovInv: %.4f, %.4f, ..., %.4f, %.4f\n",  Model->CovInv[0], Model->CovInv[1], 
    //    Model->CovInv[Model->FeatureCount * Model->FeatureCount - 2], 
    //    Model->CovInv[Model->FeatureCount * Model->FeatureCount - 1]);
    //printf("ConstantTrue %.4f, ConstantFalse %.4f\n", Model->ConstantTrue, Model->ConstantFalse);
    fclose(File);
    return Model;
}

void FreeLDAModel(LDAModel* Model)
{
    if (!Model)
    {
        return;
    }
    SafeFree(Model->MinValues);
    SafeFree(Model->MaxValues);
    SafeFree(Model->CovInv);
    SafeFree(Model->ScaledVector);
    SafeFree(Model->TempProductVector);
    SafeFree(Model->MeanVectorTrue);
    SafeFree(Model->MeanVectorFalse);
    SafeFree(Model);
}

float ApplyLDAModel(LDAModel* Model, float* Features)
{
    int FeatureIndex;
    double HalfRange;
    int ColumnIndex;
    double ProductTrue;
    double ProductFalse;
    //
    //printf("\nCFeatures %.4f...%.4f\n", Features[0], Features[Model->FeatureCount - 1]);
    // Scale the features into [-1, 1]:
    for (FeatureIndex = 0; FeatureIndex < Model->FeatureCount; FeatureIndex++)
    {
        HalfRange = (float)((Model->MaxValues[FeatureIndex] - Model->MinValues[FeatureIndex]) / 2.0);
        Model->ScaledVector[FeatureIndex] = (float)((Features[FeatureIndex] - Model->MinValues[FeatureIndex]) / HalfRange - 1.0);
    }
    //printf("Scaled vector %.4f...%.4f\n", Model->ScaledVector[0], Model->ScaledVector[Model->FeatureCount - 1]);
    // Compute the product of the inverse covariance matrix with our feature vector:
    for (FeatureIndex = 0; FeatureIndex < Model->FeatureCount; FeatureIndex++)
    {
        Model->TempProductVector[FeatureIndex] = 0;
        for (ColumnIndex = 0; ColumnIndex < Model->FeatureCount; ColumnIndex++)
        {
            Model->TempProductVector[FeatureIndex] += (float)(Model->ScaledVector[ColumnIndex] * Model->CovInv[FeatureIndex * Model->FeatureCount + ColumnIndex]);
        }
    }
    //printf("Temp product vector vector %.4f...%.4f\n", Model->TempProductVector[0], Model->TempProductVector[Model->FeatureCount - 1]);

    // Compute u0 * C-1 * X and u1 * C-1 * X
    ProductTrue = 0;
    ProductFalse = 0;
    for (FeatureIndex = 0; FeatureIndex < Model->FeatureCount; FeatureIndex++)
    {
        ProductTrue += (float)(Model->MeanVectorTrue[FeatureIndex] * Model->TempProductVector[FeatureIndex]);
        ProductFalse += (float)(Model->MeanVectorFalse[FeatureIndex] * Model->TempProductVector[FeatureIndex]);
    }
    ProductTrue += Model->ConstantTrue;
    ProductFalse += Model->ConstantFalse;
    //printf("ProdTrue %.4f ProdFalse %.4f result %.4f\n", ProductTrue, ProductFalse, ProductTrue - ProductFalse);
    //ProductTrue += (float)Model->Sub;
    //ProdFalse += (float)SubProdFalse;
    //printf("%.2f\t%.2f\t%.2f\t\n", (ProdTrue - ProdFalse), ProdTrue, ProdFalse);
    return (float)(ProductTrue - ProductFalse);
}

void InitPValueLDA()
{
    char FilePath[MAX_FILENAME_LEN];
    //
    if (!MQModel2LDA)
    {
        sprintf(FilePath, "%s%s.model", GlobalOptions->ResourceDir, "MQScoreLDA2");
        MQModel2LDA = LoadLDAModel(FilePath);
    }
    if (!MQModel3LDA)
    {
        sprintf(FilePath, "%s%s.model", GlobalOptions->ResourceDir, "MQScoreLDA3");
        MQModel3LDA = LoadLDAModel(FilePath);
    }
}

float LDAComputeMQScore(MSSpectrum* Spectrum, Peptide* Match, float* MQFeatures)
{
    LDAModel* Model;
    float Score;

    if (Spectrum->Charge < 3)
    {
        Model = MQModel2LDA;
    }
    else
    {
        Model = MQModel3LDA;
    }
    if (!Model)
    {
        return 0.0;
    }
    Score = ApplyLDAModel(Model, MQFeatures);
    Score = GetPenalizedScore(Spectrum, Match, Score);
    return Score;

}
