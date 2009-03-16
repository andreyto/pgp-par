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


// Bayesian network support functions.
// We employ a BN for scoring PRMs (prefix residue masses).  
#include "CMemLeak.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <malloc.h>
#include "Errors.h"
#include "BN.h"
#include "Utils.h"
#include "Inspect.h"
#include "Spectrum.h"
#include "Trie.h"
#include <math.h>

BayesianModel* BNCharge2ScoringBN = NULL;
BayesianModel* BNCharge3ScoringBN = NULL;
BayesianModel* BNCharge2TaggingBN = NULL;
BayesianModel* BNCharge3TaggingBN = NULL;

void FreeBayesianModel(BayesianModel* Model);

void OldInitBayesianModels()
{
    char FilePath[2048];
    if (GlobalOptions->InstrumentType == INSTRUMENT_TYPE_QTOF)
    {
        sprintf(FilePath, "%s%s", GlobalOptions->ResourceDir, "Ch2BNPEPQ.dat");
        BNCharge2ScoringBN = LoadBayesianModel(FilePath);
        sprintf(FilePath, "%s%s", GlobalOptions->ResourceDir, "PRMQ2.dat");
        BNCharge2TaggingBN = LoadBayesianModel(FilePath);
        sprintf(FilePath, "%s%s", GlobalOptions->ResourceDir, "Ch3BNPEPQ.dat");
        BNCharge3ScoringBN = LoadBayesianModel(FilePath);
        sprintf(FilePath, "%s%s", GlobalOptions->ResourceDir, "PRMQ3.dat");
        BNCharge3TaggingBN = LoadBayesianModel(FilePath);
    }
    else
    {
        sprintf(FilePath, "%s%s", GlobalOptions->ResourceDir, "Ch2BNPEP.dat");
        BNCharge2ScoringBN = LoadBayesianModel(FilePath);
        sprintf(FilePath, "%s%s", GlobalOptions->ResourceDir, "Ch3BNPEP.dat");
        BNCharge3ScoringBN = LoadBayesianModel(FilePath);
        if (GlobalOptions->DigestType == DIGEST_TYPE_TRYPSIN)
        {
            sprintf(FilePath, "%s%s", GlobalOptions->ResourceDir, "PRM2.dat");
            BNCharge2TaggingBN = LoadBayesianModel(FilePath);
            sprintf(FilePath, "%s%s", GlobalOptions->ResourceDir, "PRM3.dat");
            BNCharge3TaggingBN = LoadBayesianModel(FilePath);
        }
        else
        {
            sprintf(FilePath, "%s%s", GlobalOptions->ResourceDir, "PRM2.dat");
            BNCharge2TaggingBN = LoadBayesianModel(FilePath);
            sprintf(FilePath, "%s%s", GlobalOptions->ResourceDir, "PRM3.dat");
            BNCharge3TaggingBN = LoadBayesianModel(FilePath);
        }
    }
}

void OldFreeBayesianModels()
{
    if (BNCharge2ScoringBN)
    {
        FreeBayesianModel(BNCharge2ScoringBN);
        BNCharge2ScoringBN = NULL;
    }
    if (BNCharge2TaggingBN)
    {
        FreeBayesianModel(BNCharge2TaggingBN);
        BNCharge2TaggingBN = NULL;
    }
    if (BNCharge3ScoringBN)
    {
        FreeBayesianModel(BNCharge3ScoringBN);
        BNCharge3ScoringBN = NULL;
    }
    if (BNCharge3TaggingBN)
    {
        FreeBayesianModel(BNCharge3TaggingBN);
        BNCharge3TaggingBN = NULL;
    }

}

// Compute the probability of a basesian node:
// return ProbTable[ParentValue1*ParentBlock1 + ... + ParentValueN*ParentBlockN + FeatureValue]
float ComputeBNProbability(BayesianNode* BN, int NodeIndex, int* FeatureValues)
{
    int ProbTableIndex;
    int Parent;
    int ParentIndex;
    //
    ProbTableIndex = 0;
    for (ParentIndex = 0; ParentIndex < 4; ParentIndex++)
    {
        Parent = BN->Parents[ParentIndex];
        if (Parent >= 0)
        {
            //ProbTableIndex += BN->ParentBlocks[ParentIndex] * Model->Nodes[Parent].Value;
            ProbTableIndex += BN->ParentBlocks[ParentIndex] * FeatureValues[Parent];
        }
        else
        {
            break;
        }
    }
    ProbTableIndex += FeatureValues[NodeIndex];
    return BN->ProbTable[ProbTableIndex];
}

void FreeBayesianModel(BayesianModel* Model)
{
    int NodeIndex;
    BayesianNode* BN;
    if (Model)
    {
        for (NodeIndex = 0; NodeIndex < Model->NodeCount; NodeIndex++)
        {
            BN = Model->Nodes + NodeIndex;
            SafeFree(BN->ProbTable);
        }
        SafeFree(Model->Nodes);
        SafeFree(Model);
    }
}

BayesianModel* LoadBayesianModel(char* FileName)
{
    int FeatureCount;
    int FeatureIndex;
    FILE* File;
    BayesianNode* BN;
    BayesianModel* Model;
    //
    File = fopen(FileName, "rb");
    if (!File)
    {
        REPORT_ERROR_S(3, FileName);
        return NULL;
    }
    ReadBinary(&FeatureCount, sizeof(int), 1, File);
    if (FeatureCount < 1 || FeatureCount > 100)
    {
        REPORT_ERROR_I(6, FeatureCount);
        return NULL;
    }
    Model = (BayesianModel*)calloc(1, sizeof(BayesianModel));
    Model->NodeCount = FeatureCount;
    Model->Nodes = (BayesianNode*)calloc(FeatureCount, sizeof(BayesianNode));
    for (FeatureIndex = 0; FeatureIndex < FeatureCount; FeatureIndex++)
    {
        BN = Model->Nodes + FeatureIndex;
        ReadBinary(&BN->Flags, sizeof(int), 1, File);
        ReadBinary(&BN->ValueCount, sizeof(int), 1, File);
        ReadBinary(&BN->Name, sizeof(char), 64, File);
        if (BN->Flags & BNODE_HAS_PARENTS)
        {
            ReadBinary(BN->Parents, sizeof(int), 4, File);
            ReadBinary(BN->ParentBlocks, sizeof(int), 4, File);
            ReadBinary(&BN->ProbTableSize, sizeof(int), 1, File);
            if (BN->ProbTableSize <= 0 || BN->ProbTableSize > 1000)
            {
                REPORT_ERROR_II(7, BN->ProbTableSize, FeatureIndex);
            }
            BN->ProbTable = (float*)calloc(BN->ProbTableSize, sizeof(float));
            ReadBinary(BN->ProbTable, sizeof(float), BN->ProbTableSize, File);

        }
    }
    return Model;
}
