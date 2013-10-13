#ifndef ION_SCORING_H
#define ION_SCORING_H

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

// New code to support scoring of PRMs, and of cuts.  
// Key ideas:
// - Feature meta-data is read from the bayesian network.  New network topologies,
// and even new ion types, can be used without the need to recompile code.
// - Ion types and network topologies will be learned in an objective, repeatable way.
#include "Spectrum.h"
#include "Tagger.h"
#include "Trie.h"

#define UNKNOWN_AMINO 'Z'

typedef enum FragmentTypes
{
    evFragmentTypeNone = 0,
    evFragmentY,
    evFragmentB,
    evFragmentYLoss,
    evFragmentBLoss,
    evParentLoss //used for parent-phosphate.  not really a fragment, but used to claim the peak.
} FragmentTypes;

typedef enum PRMBayesianNodeType
{
    evPRMBInvalid = 0,
    evPRMBPrefix,
    evPRMBSuffix,
    evPRMBPrefix2,
    evPRMBSuffix2,
    evSector,
    evFlank,
    evPrefixAA,
    evSuffixAA,
    evPrefixContain,
    evSuffixContain,
    evPrefixContainPhos,
    evSuffixContainPhos
} PRMBayesianNodeType;

// A node in a Bayesian Network for scoring PRMs or cuts:
typedef struct PRMBayesianNode
{
    // Index of the node in the bayesian network (0, 1, etc).  By convention, 
    // parents will always have lower indices than their children.
    int Index;
    PRMBayesianNodeType Type;
    // For ion type nodes:
    int MassOffset;
    int FragmentType; // from the FragmentTypes enum
    // The Flag on a BayesianNode of special type is
    // used to control how the node's values are computed.
    // Examples: For type evFlank, the flag tells us whether
    // we're looking for flanking aminos that affect b fragments,
    // or y fragments.
    int Flag;
    struct PRMBayesianNode* Next;
    // ValueCount is the number of distinct values this node can take on (2 or more).
    // ValueCount is determined by our Type and Flag, and possibly by the network's
    // intensity scheme.
    int ValueCount;
    int ParentCount;
    struct PRMBayesianNode** Parents;
    // ParentBlock[n] is the multiplier for parent n's value when indexing into 
    // the CountTable/ProbTable arrays.  For instance, if we have 4 possible values
    // and one parent, then ParentBlock[0] will be 4, and the index of a table entry
    // is ParentValue*4 + ChildValue.
    int* ParentBlocks;
    // Size of CountTable and ProbTable:
    int TableSize; 
    // Table counting the number of occurrences of a given value combination:
    int* CountTable;
    // Table giving natural logarithm of the probability of a given value combination:
    float* ProbTable;
    // Value is transiently set while scoring a PRM or cut point:
    int Value; 
    // An array of values for cut points:
    int Values[256]; 
    // Human-readable name of the node, mostly for debugging:
    char Name[256 + 1];
    // Flag to indicate whether this node, or an immediate parent, requires knowledge of 
    // flanking amino acids.  If this flag is set, then during tagging, we will delay
    // full scoring of this node until tag construction.
    int FlankFlag;
} PRMBayesianNode;

typedef struct PRMBayesianNodeHolder
{
    PRMBayesianNode* Node;
    struct PRMBayesianNodeHolder* Next;
} PRMBayesianNodeHolder;

typedef struct PRMBayesianModel
{
    PRMBayesianNode* Head;
    PRMBayesianNode* Tail;
    // Array of the nodes, for quickly looking them up by index:
    PRMBayesianNode** Nodes;
    int NodeCount;
    // Scheme for assigning intensity-levels to ion nodes:
    int IntensityScheme;
    // 0 is spectrum-specific, 1 is global
    int NoiseModel;
    // Radius (in daltons) of the window over which to sum intensities when
    // finding peaks:
    int IntensityRadius;
    int HalfIntensityRadius;
    // CutFlag is true if this model is used to score cut points.  A few operations
    // differ; in particular, we seize intensities for b and y peaks first, THEN consider
    // neutral losses.
    int CutFlag;
    // Intensity levels are sorted from highest (0) to lowest (MinIntensityLevel).
    int MinIntensityLevel;
    // RandomIntensityCounts and RandomIntensityScores track how often a *random* mass 
    // has a particular intensity level.  We'll try using a spectrum-specific noise 
    // model as well as this "global" noise model.
    int RandomIntensityCounts[10];
    float RandomIntensityScores[10];
    // Linked list of nodes which require flanking amino acid information (either directly,
    // or via parents):
    PRMBayesianNodeHolder* FirstFlank;
    PRMBayesianNodeHolder* LastFlank;
} PRMBayesianModel;

void AddPRMBayesianNode(PRMBayesianModel* Model, char* Name, int NodeType, int NodeFlag, float NodeMassOffset, int FragmentType);
void FreePRMBayesianModel(PRMBayesianModel* Model);
void FreePRMBayesianNode(PRMBayesianNode* Node);
void PrepareSpectrumForIonScoring(PRMBayesianModel* Model, MSSpectrum* Spectrum, int ForceRefresh);
int IonScoringGetNodeValue(PRMBayesianModel* Model, PRMBayesianNode* Node, MSSpectrum* Spectrum, int PRM,
    Peptide* Match, int AminoIndex);
void ComputePRMBayesianModelProbabilityTables(PRMBayesianModel* Model, int PaddingCount);
void SavePRMBayesianModel(PRMBayesianModel* Model, char* FileName);
PRMBayesianModel* LoadPRMBayesianModel(char* FileName);
void DebugPrintPRMBayesianModel(PRMBayesianModel* Model);
void TagGraphScorePRMNodes(PRMBayesianModel* Model, struct TagGraph* Graph, MSSpectrum* Spectrum, SpectrumTweak* Tweak);
float GetIonPRMFeatures(MSSpectrum* Spectrum, SpectrumTweak* Tweak, PRMBayesianModel* Model, int PRM, int VerboseFlag);
void BuildModelFlankList(PRMBayesianModel* Model);
void LoadFlankingAminoEffects();
int IonScoringGetFlank(PRMBayesianNode* Node, char Left, char Right);
float PRMBNGetCutScore(MSSpectrum* Spectrum, PRMBayesianModel* Model, int AminoIndex);
void InitBayesianModels();
int ReplacePRMScoringModel(int Charge, char* FileName);
int ReplaceTAGScoringModel(int Charge, char* FileName);
void SetSpectrumPRMScores(MSSpectrum* Spectrum, SpectrumTweak* Tweak);
void PopulateCutScores(PRMBayesianModel* Model, MSSpectrum* Spectrum, Peptide* Match, float* CutScores);
int CountTrypticTermini(Peptide* Match);
int ComputeMQScoreFeatures(MSSpectrum* Spectrum, Peptide* Match, float* MQFeatures, int VerboseFlag);
char* GetFragmentTypeName(int FragmentType);
void FreeBayesianModels();

extern PRMBayesianModel* PRMModelCharge2;
extern PRMBayesianModel* PRMModelCharge3;
extern PRMBayesianModel* TAGModelCharge2;
extern PRMBayesianModel* TAGModelCharge3;
extern PRMBayesianModel* PhosCutModelCharge2;
extern PRMBayesianModel* PhosCutModelCharge3;

#endif // ION_SCORING_H


