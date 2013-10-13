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

#ifndef TAGGER_H
#define TAGGER_H
// Tagger.h defines objects and functions related to the building of a tag graph,
// and the generation on short (usually tripeptide) tags from the graph.  Each 
// node in the graph represents a prefix residue mass (PRM), and so the nodes can 
// be thought of as points along an m/z axis.  A (directed) edge in the graph 
// represents a valid jump, where a "jump" is the mass of an amino acid (or modified
// amino acid).  The tag graph is used to construct tags (formally, paths of length
// three).  This is one approach to the local de novo interpretation problem.

#include "Inspect.h"
#include "Trie.h"
#include "Spectrum.h"
#include "IonScoring.h"

// Each graph node has some set of witness ions (e.g. b and y peaks).  We encode 
// the a set of ions as a bitfield.  
// For instance, B and Y is ION_FLAG_B | ION_FLAG_Y = 0x11
#define ION_FLAG_B 0x0001
#define ION_FLAG_BH2O 0x0002
#define ION_FLAG_BNH3 0x0004
#define ION_FLAG_A 0x0008
#define ION_FLAG_Y 0x0010
#define ION_FLAG_YH2O 0x0020
#define ION_FLAG_YNH3 0x0040
#define ION_FLAG_B2 0x0080
#define ION_FLAG_Y2 0x0100

// INTENSITY_RANK_COUNT is the number of entries in IntensityRankBOdds and IntensityRankYOdds
#define INTENSITY_RANK_COUNT 22

// Probability that a peak with this intensity-rank is a b peak:
extern float* IntensityRankBOdds;

// Probability that a peak with this intensity-rank is a y peak:
extern float* IntensityRankYOdds;

// Probability that a peak with reported m/z this far from the expected m/z is a true b or y peak:
extern double SkewHistoStep[100];

// The witness scores for a node are based upon the collection of ion
// types bearing witness to a particular break-point.  Similar to Danciek 
// scores.  Scores are empirically derived from a training dataset (currently ISB).  dataset.  We reckon these odds
// separately for low, medium and high mass peaks.  "Low" is "below 33% of precursor mass", 
// and "High" is "above 66% of precursor mass".

// A JumpNode captures the mass, amino acid, and PTM (if any) of a valid 
// edge length for the tag graph.
typedef struct JumpNode
{
    int Mass;
    struct JumpNode* Next;
    char Amino;
    // ASSUME: We only permit one modification per peptide in a tag. 
    MassDelta* Delta;
    float Score;
} JumpNode;

// The type of a graph node indicates whether it was created by interpreting a spectral peak as a b or y ion,
// or whether it is an endpoint (evGraphNodeLeft, evGraphNodeRight).  The special types evGraphNodeLeftMod,
// evGraphNodeRightMod are created when N- and C-terminal PTMs are allowed.  
typedef enum evGraphNodeType
{
    evGraphNodeB = 0,
    evGraphNodeY,
    evGraphNodeLeft,
    evGraphNodeLeftMod,
    evGraphNodeRight,
    evGraphNodeRightMod
} evGraphNodeType;

typedef struct TagGraphNode
{
#ifdef DEBUG_TAG_GENERATION
    char VerboseNodeInfo[2048];
#endif
    int OriginalPeakIndex;
    int BIndex;
    int YIndex;
    int IntensityRankB;
    int IntensityRankY;
    evGraphNodeType NodeType; 
    // A graph node is scored based upon its intensity score (intensity-rank of the b and y peak),
    // its isotope score (whether the b and y peaks are apparently secondary isotopic peaks, primary
    // peaks with children, or lone peaks) and its ion type score (the witness set).
    float IntensityScore;
    float IsotopeScore;
    float IonTypeScore;
    float Score;
    //float ScoreB;
    //float ScoreY;
    int IonTypeFlags;
    int Mass;
    // List of edges leading forward in the graph:
    struct TagGraphEdge* FirstEdge;
    struct TagGraphEdge* LastEdge;
    // Next and previous nodes (sorted by mass):
    struct TagGraphNode* Next;
    struct TagGraphNode* Prev;
    MassDelta* PTM; // Is non-null only if NodeType is LeftMod or RightMod
    // BackEdge, BackEdgeDouble, and BackEdgeTriple are set only when carrying out blind
    // mod search.  They speed up the big d.p. extension algorithm.
    struct TagGraphBackEdge** BackEdge; //[26]; // List of edges matching an unmodified aa
    struct TagGraphBackEdge** BackEdgeDouble; //[26*26]; // List of edges matching two unmodified aa's
    struct TagGraphBackEdge** BackEdgeTriple; //[26*26*26]; // List of edges matching three unmodified aa's
    int Index; // This is set AFTER all the graph nodes have been created and sorted.
} TagGraphNode;

// BackEdge points to a graph node whose mass is smaller by 1, 2, or 3 unmodified amino acid masses.
typedef struct TagGraphBackEdge
{
    TagGraphNode* FromNode;
    TagGraphNode* ToNode;
    int Score;
    int Skew;
    // If this edge is a double-amino-acid jump, then HalfMass is the mass after the first amino acid.
    int HalfMass; 
    int HalfMass2;  // For triples
    struct TagGraphBackEdge* Next;
} TagGraphBackEdge;

// Each NODE in the graph owns a list of EDGES.  Each edge joins to a higher-mass node
typedef struct TagGraphEdge
{
    TagGraphNode* FromNode;
    TagGraphNode* ToNode;
    JumpNode* Jump;
    float Score;
    struct TagGraphEdge* Next;
    int Skew;
} TagGraphEdge;

// A TagGraph has pointer to its first/last nodes, an index (for quickly finding nodes for a PRM),
// and a buffer of back edges (populated only in blind mode).
typedef struct TagGraph
{
    TagGraphNode* FirstNode;
    TagGraphNode* LastNode;
    // Index: Points to the first node that could match a given rounded-to-amu mass
    TagGraphNode** NodeIndex; 
    int NodeIndexSize;
    int NodeCount; // Number of nodes in the list FirstNode...LastNode.
    struct TagGraphBackEdge* BackEdgeBuffer;
} TagGraph;

void TagGraphBuildNodeIndex(TagGraph* Graph);
TrieNode* GenerateTagsFromSpectrum(MSSpectrum* Spectrum, TrieNode* Root, int MaximumTagCount, SpectrumTweak* Tweak);
void CorrectParentMass(MSSpectrum* Spectrum);
int LoadIntensityRankOdds(char* FileName);
int LoadWitnessScores(char* FileName);
void PopulateJumpingHash();
int FindIntensePeak(MSSpectrum* Spectrum, int Mass, float MaxIntensity, float* FoundIntensity);
void SpectrumFindIsotopicPeaks(MSSpectrum* Spectrum);
TagGraph* ConstructTagGraph(MSSpectrum* Spectrum);
void TagGraphAddNodes(TagGraph* Graph, MSSpectrum* Spectrum);
void TagGraphPopulateEdges(TagGraph* Graph);
void FreeTagGraph(TagGraph* Graph);
void FreeJumpingHash();
void FreeTagGraphNode(TagGraphNode* Node);
void TestTagging(char* OracleFile, char* OracleDir);
void TrainTagging(char* OracleFile, char* OracleDir);
int CompareTagScores(const TrieTag* TagA, const TrieTag* TagB);
TrieNode* BuildTrieFromTags(TrieTag* AllTags, int TagCount, TrieNode* Root, int MaximumTagCount);
void SetTagSkewScores();
void FreeTagSkewScores();
// declaration of TagGraphGenerateTags moved out, since it uses PRMBayesianModel
#endif // TAGGER_H
