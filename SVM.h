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

#ifndef SVM_H
#define SVM_H
// Structs to support use of SVMs:
#include "Utils.h"
#include "Inspect.h"
#include "Spectrum.h"
#include "Trie.h"

// Support vectors are of this length (or shorter)
#define SUPPORT_VECTOR_LENGTH 32

typedef struct SupportVector
{
    //int Classification; // +1 or -1
    double Weight;
    double Coords[SUPPORT_VECTOR_LENGTH];
    struct SupportVector* Next;
} SupportVector;

typedef struct SVMModel
{
    SupportVector* FirstVector;
    SupportVector* LastVector;
    int Coords;
    double ScaleMin[SUPPORT_VECTOR_LENGTH];
    double ScaleMax[SUPPORT_VECTOR_LENGTH];
    double ScaleSize[SUPPORT_VECTOR_LENGTH];
    double Beta[SUPPORT_VECTOR_LENGTH]; // for computing classifier values
    double Beta0;
    double Rho;
    double Gamma; // for RBF kernel
} SVMModel;

extern SVMModel* PValueSVMModel;

float SVMComputeMQScore(MSSpectrum* Spectrum, Peptide* Match, float* MQFeatures);
float SVMClassify(SVMModel* Model, float* Coords, int PreScaled);
void FreeSVMModels();
SVMModel* ReadSVMModel(char* FileName);
void ReadSVMScaling(SVMModel* Model, char* ScaleFileName);
float GetPValue(float MQScore);
float LDAClassify(float* Features);
void TestPValue(char* FeatureVectorFileName);
void LoadCCModelSVM(int ForceRefresh);
void FreeCCModelSVM();
void InitPValueSVM();

#endif // SVM_H

