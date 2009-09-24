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

#ifndef BN_H
#define BN_H

// Structs to support use of Bayesian Networks.  We use bayesian 
// networks to score PRMs, for both tag generation and for final 
// scoring of matches.  Most of the nodes in the network correspond to
// fragment types.  The edges between nodes help capture the co-occurrence
// relations between peaks (e.g. b-h2o is more likely in presence of b),
// as well as other factors that predict peak strength (e.g. which spectrum
// sector the PRM lies in)
///
///
// The bayesian network file has the following format:
// There's one NodeRecord per bayesian network node.  The NodeRecord 
// has flags (int), ValueCount (int; the number of possible values for the node), 
// and a Name (char64).  If the node has parents, it then has:
// Parent indices (4 ints)
// Parent block-sizes (4 ints, used in computing positions in the probability table)
// ProbTableSize (equals the first parent block size * the ValueCount)
// Probability table (float array of size ProbTableSize)
// Note that the values in the probability table are log-probabilities, so that we can
// add them up.
#include "Utils.h"
#include "Inspect.h"
#include "Spectrum.h"
#include "Trie.h"

// Flags for BayesianNode.Flags:
// A node has the BNODE_USE_PROB flag set if it's a leaf node, whose
// probability is to be used.  A node has the BNODE_HAS_PARENTS flag
// set if it has one or more parent nodes.
#define BNODE_HAS_PARENTS 1
#define BNODE_USE_PROB 2

typedef struct BayesianNode
{
    int Flags;
    int Value;
    int ValueCount;
    char Name[64];
    int Parents[4];
    int ParentBlocks[4];
    int ProbTableSize; // redundant, but useful to keep around for sanity-checks
    float* ProbTable;
} BayesianNode;

typedef struct BayesianModel
{
    BayesianNode* Nodes;
    int NodeCount;
} BayesianModel;

extern BayesianModel* BNCharge2ScoringBN;
extern BayesianModel* BNCharge3ScoringBN;
extern BayesianModel* BNCharge2TaggingBN;
extern BayesianModel* BNCharge3TaggingBN;

BayesianModel* LoadBayesianModel(char* FileName);
float ComputeBNProbability(BayesianNode* BN, int NodeIndex, int* FeatureValues);
void OldFreeBayesianModels();
void OldInitBayesianModels();

#endif // BN_H

