
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

#ifndef SCORPION_H
#define SCORPION_H

// Scorpion - Ion-based scoring of mass spectra.  Compute various features 
// for use by the match-scoring SVM.

#include "Inspect.h"
#include "Trie.h"
#include "BN.h"

#define CUT_FEATURE_COUNT 32

extern int g_CutFeatures[];

// Features used in scoring of cut points
typedef enum ScorpIons
{
    SIDynamicRange = 0,
    IonY,
    IonB,
    IonYI,
    IonBI,
    IonY2,
    IonB2,
    IonYH2O,
    IonA,
    IonBH2O,
    IonYNH3,
    IonBNH3,
    SICharge,
    SIFlankB, // Flanking amino acids that affect N-terminal fragments
    SIFlankY, // Flanking amino acids that affect C-terminal fragments
    SISector,
    //SIBasePrefix,
    //SIAcidPrefix,
    //SIBaseSuffix,
    //SIAcidSuffix,
    //SIFlankLeft,
    //SIFlankRight,
    SITestA,
    SITestB,
    SITestC,
    SITestD,
    SITestE,
    SITestF,
    SITestG,
    SITestH,
    SITestI,
    SITestJ,
    SIMax,
    IonParentLoss
} ScorpIons;

typedef enum CutFeature
{
    CFDynamic,
    CFCharge,
    CFFlank,
    CFSector,
    CFBasic,
    CFAcidic,
    CFY,
    CF
} CutFeature;

void GetCutFeatures(MSSpectrum* Spectrum, SpectrumTweak* Tweak, Peptide* Match, BayesianModel* Model);
void ScorpionSetPRMScores(MSSpectrum* Spectrum, SpectrumTweak* Tweak);
void TestPRMQuickScoring(char* OracleFile, char* OracleDir);
float GetExplainedIntensityPercent(MSSpectrum* Spectrum, int PeakCount, int BYOnly);
float GetExplainedPeakPercent(MSSpectrum* Spectrum, int PeakCount, int BYOnly);
int GetPeptideParentMass(Peptide* Match);
void TrainPepPRM(char* OracleFile, char* OracleDir);
void TestLDA(char* OracleFile, char* OracleDir);
void TestPepPRM(char* OracleFile, char* OracleDir);
char* GetScorpIonName(int IonType);
#endif // SCORPION_H
