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

#ifndef SCORE_H
#define SCORE_H
#include "Inspect.h"
#include "Trie.h"

// For marking forbidden paths in d.p. tables: 
#define FORBIDDEN_PATH -99999999

typedef enum evIonType
{
    evIonTypeNone = 0,
    evIonTypeB, 
    evIonTypeY,
    evIonTypeA,
    evIonTypeBH2O,
    evIonTypeAH2O,
    evIonTypeBNH3,
    evIonTypeANH3,
    evIonTypeYH2O,
    evIonTypeYNH3,
    evIonTypeB2,
    evIonTypeY2,
    evIonTypeNoise,
    evIonTypeBPhos, // B minus a phosphorylation
    evIonTypeYPhos, // Y minus a phosphorylation
    evIonTypeCount,
} evIonType;

typedef enum evLossType
{
    evLossNone = 0,
    evLossB,
    evLossY,
} evLossType;

// Workhorse function of Score.c: Compare the theoretical fragmentation pattern of a peptide to
// a spectrum, and assign a score.  
int ScoreMatch(MSSpectrum* Spectrum, Peptide* Match, int VerboseFlag);

// Apply a penalty if the peptide doesn't match GlobalOptions->DigestType.  We have a 
// special scoring model for tryptic peptides, but other less specific proteases - like 
// GluC - also should affect scoring based on endpoints.
int ApplyDigestBonus(Peptide* Match);

// Compute the p-value for a match, based upon explained intensity and explained peaks and b/y ladder
// and match score.
//void ComputeMatchConfidenceLevel(MSSpectrum* Spectrum, Peptide* Match);
void ScoreMatchTest(int VerboseFlag);
int InitPValue(char* FileName);
void SetMatchDeltaCN(SpectrumNode* Spectrum);
int GetPeptideMatchFeaturesFull(MSSpectrum* Spectrum, Peptide* Match, float* FeatureArray);
float GetPenalizedScore(MSSpectrum* Spectrum, Peptide* Match, float Score);
#endif // SCORE_H
