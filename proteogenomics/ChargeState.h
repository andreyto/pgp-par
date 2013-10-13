#ifndef CHARGE_STATE_H
#define CHARGE_STATE_H

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
#include "ParentMass.h"

// Code to support charge state determination.  Our plan is:
// - Organize a 'training+test' directory of spectra, half charge 2 and half charge 3.  The directory should 
// include some QTOF results.
// - Use an API in ChargeState.c to write out a set of features for these spectra.  +1 means charge 3, in this case.
// - Use libsvm to train a support vector machine on these features
// - Use the resulting model to guess charge states if the charge is unlisted, or if the MultiCharge option is set.  We 
// use easy heuristics to detect +1 spectra, then use the svm to separate +2 and +3.  If confidence is low, still search
// both charge states.

void TweakSpectrum(SpectrumNode* Node);

void GetChargeCorrectionFeatures1(PMCSpectrumInfo* SpectrumInfo1, PMCSpectrumInfo* SpectrumInfo2,
    PMCSpectrumInfo* SpectrumInfo3, float* Features);
void GetChargeCorrectionFeatures2(PMCSpectrumInfo* SpectrumInfo2, PMCSpectrumInfo* SpectrumInfo3,
    float* Features);
void GetChargeCorrectionFeatures2Phos(PMCSpectrumInfo* SpectrumInfo2, PMCSpectrumInfo* SpectrumInfo3,
    float* Features);
int ChargeCorrectSpectrum(SpectrumNode* Node, float* Model1Score, float* Model2Score);
#endif // CHARGE_STATE_H
