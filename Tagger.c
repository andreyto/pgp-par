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

// Tag generation.  Given a spectrum, generate TrieTag objects.
// See TagTrainer.py for the generation of the tagging model based on
// empirical ion frequencies.

#include "CMemLeak.h"
#include "Inspect.h"
#include "Utils.h"
#include "Tagger.h"
#include "Spectrum.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Trie.h"
#include "Mods.h"
#include "Score.h"
#include "FreeMod.h"
#include "BN.h"
#include "Scorpion.h"
#include "Errors.h"

// For each mass (in daltons), we have an expected ratio of the +1 isotope to
// the main peak.  These ratios rise from near-0 to near 1 as mass goes from 0
// up to 1750.  (For larger masses, isotope sets are more complex!)
#define ISOTOPE_WEIGHT_COUNT 1750

// How far can the isotope-peak-to-main-peak ratio differ from what we expect?
// (This controls whether the members Is and HasPlausibleIsotopicPeak are set
// for spectral peaks.  In practice, isotope ratios can vary quite a bit, so 
// we're fairly permissive)
#define TIGHT_ISOTOPE_RATIO_DIFFERENCE 0.5
#define MAX_ISOTOPE_RATIO_DIFFERENCE 0.8

// How many intensity-ranks do we track?  (We're granular for top-10 peaks, then less so 
// for crappier peaks)
#define INTENSITY_RANK_COUNT 22

// Number of charges and sectors for which we have a tagging model.  The size
// of the array g_TruthByNodeType will equal CHARGE_COUNT * CHUNK_COUNT, similarly for other
// model params.
#define CHARGE_COUNT 3
#define CHUNK_COUNT 3

// JumpingHash[n] has a list of all amino acids (and modifications) that have a mass
// rounding off to n daltons.  To find, e.g., a jump matching 80.3, we'd check 
// JumpingHash[79] and JumpingHash[80] and JumpingHash[81].
#define MAX_JUMPING_HASH 1024
JumpNode** JumpingHash = NULL;

// 
StringNode* FirstTagCheckNode;
StringNode* LastTagCheckNode;

// JumpsByAA is a 2D array containing pointers to all jumps for each amino acid.
// We iterate over JumpsByAA when setting jump scores.
JumpNode** JumpsByAA; //[AMINO_ACIDS * GlobalOptions->DeltasPerAA];

// Skew histo is used in scoring edges - SkewHistoStep[n] is the score penalty we
// apply to an edge that deviates from the ideal jump size by floor(n/100).
// Hard-coded, and should probably be part of the tagging model!
double SkewHistoStep[] = {0.0000, 0.0000, -0.0037, -0.0501, -0.1061, -0.2106, -0.2418, 
    -0.3466, -0.4198, -0.4800, -0.4861, -0.4863, -0.4926, -0.5421, -0.5893, -0.6851, 
    -0.7173, -0.8108, -0.8811, -0.9275, -0.9302, -0.9302, -0.9353, -0.9802, -1.0259, 
    -1.1027, -1.1395, -1.2253, -1.2640, -1.2921, -1.2931, -1.2931, -1.2960, -1.3207, 
    -1.3587, -1.4249, -1.4567, -1.5158, -1.5444, -1.5696, -1.5710, -1.5710, -1.5710, 
    -1.5989, -1.6153, -1.6641, -1.6884, -1.7345, -1.7544, -1.7708, -1.7714, -1.7714, 
    -1.7714, -1.7903, -1.8150, -1.8556, -1.8819, -1.9274, -1.9462, -1.9716, -1.9723, 
    -1.9723, -1.9737, -1.9976, -2.0244, -2.0640, -2.0957, -2.1574, -2.1966, -2.2151, 
    -2.2151, -2.2151, -2.2178, -2.2512, -2.2848, -2.3753, -2.4070, -2.4806, -2.5463, 
    -2.5879, -2.5891, -2.5891, -2.5917, -2.6571, -2.6951, -2.8598, -2.9352, -3.0991, 
    -3.1702, -3.2468, -3.2493, -3.2493, -3.2669, -3.4412, -3.5462, -3.9350, -4.1434, 
    -4.8738, -5.7776, -6.5501};

// IsotopeWeights[n] is the expected ratio for the +1 isotope to the +0 isotope, given
// a peptide that whose weight is n daltons.
float IsotopeWeights[ISOTOPE_WEIGHT_COUNT];

// Forward declarations:
void TagGraphAddEndpointNodes(TagGraph* Graph, MSSpectrum* Spectrum);
void DebugPrintTagList(MSSpectrum* Spectrum, TrieTag* Tags, int TagCount);
void PrintTagToLine(FILE* OutputFile, TrieTag* Tag);
void DebugPrintTagGraph(MSSpectrum* Spectrum, TagGraph* Graph);
void DebugPrintTagsForPeptide(MSSpectrum* Spectrum, TagGraph* Graph, TrieTag* Tags, int TagCount);
TrieTag* TagGraphGenerateTags(TagGraph* Graph, MSSpectrum* Spectrum, int* TagCount, 
    int MaximumTagCount, SpectrumTweak* Tweak, float TagEdgeScoreMultiplier,
    PRMBayesianModel* Model);


// New (as of 3/2005) tagging model.  We load one model in for each NodeType.  (A NodeType is a combination of charge
// and sector - e.g. charge 3+ and middle sector).  This model is used in scoring tag graph nodes, and in scoring
// PRMs.

#define BY_RANK_TINY 20
#define BY_RANK_MISSING 21
#define BY_RANK_COUNT 22
#define BY_SKEW_COUNT 5
typedef struct TaggingModel
{
    int BRank[BY_RANK_COUNT];
    int SisterBRank[BY_RANK_COUNT];
    int SisterBSkew[BY_SKEW_COUNT];
    int BSkew[BY_SKEW_COUNT];
    int SkewableBRank[BY_RANK_COUNT];
    int YRank[BY_RANK_COUNT];
    int SisterYRank[BY_RANK_COUNT];
    int SisterYSkew[BY_SKEW_COUNT];
    int YSkew[BY_SKEW_COUNT];
    int SkewableYRank[BY_RANK_COUNT];    
    int Witness[512];
    int BIsotope[4];
    int YIsotope[4];
} TaggingModel;

typedef struct TagMaster
{
    TaggingModel Models[CHARGE_COUNT * CHUNK_COUNT];
    float PTMPenalty;
} TagMaster;

TagMaster MasterTaggingModel;

// Constructor for TagGraph
TagGraph* ConstructTagGraph(MSSpectrum* Spectrum)
{
    TagGraph* Graph;
    Graph = (TagGraph*)calloc(sizeof(TagGraph), 1);
    return Graph;
}

// Destructor for a node from a TagGraph, as well as the nodes edges.
void FreeTagGraphNode(TagGraphNode* Node)
{
    TagGraphEdge* Edge;
    TagGraphEdge* PrevEdge = NULL;
    if (!Node)
    {
        return;
    }
    for (Edge = Node->FirstEdge; Edge; Edge = Edge->Next)
    {
        SafeFree(PrevEdge);
        PrevEdge = Edge;
    }
    SafeFree(PrevEdge);
    // Back edges:
    SafeFree(Node->BackEdge);
    SafeFree(Node->BackEdgeDouble);
    SafeFree(Node->BackEdgeTriple);
    SafeFree(Node);
}

// Destructor for a TagGraph.
void FreeTagGraph(TagGraph* Graph)
{
    TagGraphNode* TagNode;
    TagGraphNode* PrevTagNode = NULL;

    //
    if (!Graph)
    {
        return;
    }
    SafeFree(Graph->BackEdgeBuffer);
    for (TagNode = Graph->FirstNode; TagNode; TagNode = TagNode->Next)
    {
        if (PrevTagNode)
        {
            FreeTagGraphNode(PrevTagNode);
        }
        PrevTagNode = TagNode;
    }
    if (PrevTagNode)
    {
        FreeTagGraphNode(PrevTagNode);
    }
    SafeFree(Graph->NodeIndex);
    Graph->FirstNode = NULL;
    Graph->LastNode = NULL;
    SafeFree(Graph);
}

// Somewhat ugly macro for inserting a tag graph node into the list.
// (Note that since we're not inserting phosphate loss peaks, this is
// overkill - we will always be inserting at the end of the list in practice!)
#define INSERT_TAGNODE_ASC(First, Last, Node)\
{\
InsertAfter = (Last);\
while ((InsertAfter) && (InsertAfter)->Mass > (Node)->Mass)\
{\
    (InsertAfter) = (InsertAfter)->Prev;\
}\
if (InsertAfter)\
{\
    if ((InsertAfter)->Next)\
    {\
        InsertAfter->Next->Prev = Node;\
    }\
    Node->Next = InsertAfter->Next;\
    InsertAfter->Next = Node;\
    Node->Prev = InsertAfter;\
}\
else\
{\
    Node->Next = First;\
    if (First)\
    {\
        First->Prev = Node;\
    }\
    First = Node;\
}\
if (InsertAfter == Last)\
{\
    Last = Node;\
}\
}

#define INSERT_TAGNODE_DESC(First, Last, Node)\
{\
InsertAfter = Last;\
while (InsertAfter && InsertAfter->Mass < Node->Mass)\
{\
    InsertAfter = InsertAfter->Prev;\
}\
if (InsertAfter)\
{\
    if (InsertAfter->Next)\
    {\
        InsertAfter->Next->Prev = Node;\
    }\
    Node->Next = InsertAfter->Next;\
    InsertAfter->Next = Node;\
    Node->Prev = InsertAfter;\
}\
else\
{\
    Node->Next = First;\
    if (First)\
    {\
        First->Prev = Node;\
    }\
    First = Node;\
}\
if (InsertAfter == Last)\
{\
    Last = Node;\
}\
}

// Take a new (empty) tag graph, and add nodes to it.  For each spectral peak, add 2 nodes (one b, one y).
// Also add endpoint nodes. 
void TagGraphAddNodes(TagGraph* Graph, MSSpectrum* Spectrum)
{
    int PeakIndex;
    TagGraphNode* FirstBNode = NULL;
    TagGraphNode* LastBNode = NULL;
    TagGraphNode* FirstYNode = NULL;
    TagGraphNode* LastYNode = NULL;
    TagGraphNode* InsertAfter = NULL;
    TagGraphNode* Node;
    TagGraphNode* MergingBNode;
    TagGraphNode* MergingYNode;
    int NodeMass;
    int NodeIndex;
    int MinPRMMass = 50 * DALTON;
    // Iterate over peaks.  For each peak, add a b and a y node.  We'll build two lists,
    // b node list (FirstBNode / LastBNode) and y node list (FirstYNode / LastYNode), then
    // merge the lists.
    for (PeakIndex = 0; PeakIndex < Spectrum->PeakCount; PeakIndex++)
    {
        NodeMass = Spectrum->Peaks[PeakIndex].Mass - HYDROGEN_MASS;
        // Filter any nodes whose mass is negative, or <50 but not zero, or larger than the precursor mass.  
        // (A node at, say, PRM 19 couldn't possibly be part of a true peptide, since the smallest mass-jump 
        // is 57)
        if (NodeMass > -GlobalOptions->Epsilon && (NodeMass < GlobalOptions->Epsilon || NodeMass > MinPRMMass) && (NodeMass < Spectrum->ParentMass + GlobalOptions->ParentMassEpsilon))
        {
            Node = (TagGraphNode*)calloc(sizeof(TagGraphNode), 1);
            Node->NodeType = evGraphNodeB;
            Node->OriginalPeakIndex = PeakIndex;
            Node->IntensityRankB = Spectrum->Peaks[PeakIndex].IntensityRank;
            Node->BIndex = PeakIndex;
            Node->YIndex = -1;
            Node->IntensityRankY = -1;
            Node->Mass = NodeMass;
            Node->IonTypeFlags = ION_FLAG_B;
            INSERT_TAGNODE_ASC(FirstBNode, LastBNode, Node);
			if (0)//(Spectrum->Charge > 2 && Spectrum->Peaks[PeakIndex].IntensityRank < 16)
			{ //charge 3 spectra have high intensity doubly-charged peaks. I need to put those into the graph
				Node = (TagGraphNode*)calloc(sizeof(TagGraphNode), 1);
				Node->NodeType = evGraphNodeB;
				Node->OriginalPeakIndex = PeakIndex;
				Node->IntensityRankB = Spectrum->Peaks[PeakIndex].IntensityRank;
				Node->BIndex = PeakIndex;
				Node->YIndex = -1;
				Node->IntensityRankY = -1;
				Node->Mass = (NodeMass * 2) - HYDROGEN_MASS; //the single charge mass, if the peak was doubly charged
				Node->IonTypeFlags = ION_FLAG_B;
				INSERT_TAGNODE_ASC(FirstBNode, LastBNode, Node);
			}
			//printf("Peak %d intensity %f rank %d\n",Spectrum->Peaks[PeakIndex].Mass, Spectrum->Peaks[PeakIndex].Intensity, Spectrum->Peaks[PeakIndex].IntensityRank);
        }
        NodeMass = Spectrum->ParentMass - Spectrum->Peaks[PeakIndex].Mass;
        if (NodeMass > -GlobalOptions->Epsilon && (NodeMass < GlobalOptions->Epsilon || NodeMass > MinPRMMass) && (NodeMass < Spectrum->ParentMass + GlobalOptions->ParentMassEpsilon))
        {
            Node = (TagGraphNode*)calloc(sizeof(TagGraphNode), 1);
            Node->NodeType = evGraphNodeY;
            Node->OriginalPeakIndex = PeakIndex;
            Node->IntensityRankY = Spectrum->Peaks[PeakIndex].IntensityRank;
            Node->IntensityRankB = -1;
            Node->YIndex = PeakIndex;
            Node->BIndex = -1;
            Node->Mass = NodeMass;
            Node->IonTypeFlags = ION_FLAG_Y;
            INSERT_TAGNODE_DESC(FirstYNode, LastYNode, Node);
			if (0)//(Spectrum->Charge > 2 && Spectrum->Peaks[PeakIndex].IntensityRank < 16)
			{ //charge 3 spectra have high intensity doubly-charged peaks. I need to put those into the graph
				Node = (TagGraphNode*)calloc(sizeof(TagGraphNode), 1);
				Node->NodeType = evGraphNodeY;
				Node->OriginalPeakIndex = PeakIndex;
				Node->IntensityRankY = Spectrum->Peaks[PeakIndex].IntensityRank;
				Node->IntensityRankB = -1;
				Node->YIndex = PeakIndex;
				Node->BIndex = -1;
				Node->Mass = (NodeMass * 2) - HYDROGEN_MASS; //the single charge mass, if the peak was doubly charged
				Node->IonTypeFlags = ION_FLAG_Y;
				INSERT_TAGNODE_DESC(FirstYNode, LastYNode, Node);
			}
        }
        // We could insert phosphate-loss peaks for b and y nodes at this point.
        // There are cases (particularly for breaks next to the phosphorylation site,
        // and for phosphoserines) when the phosphate-loss peak is MORE LIKELY than
        // the original peak.  However, if we insert phosphate-loss peaks, we end up
        // with 4 nodes per peak rather than 2...that slows the speed of tag generation 
        // down *considerably*, and probably lowers selectivity quite a bit.
    }
    ///////////////////////////////////////////
    // Merge b and y node lists into the list Graph->FirstNode...Graph->LastNode
    MergingBNode = FirstBNode;
    MergingYNode = LastYNode;
    while (1)
    {
        if (!MergingBNode && !MergingYNode)
        {
            break;
        }
        if (!MergingBNode || (MergingYNode && (MergingBNode->Mass > MergingYNode->Mass)))
        {
            // Insert the y node into the master list:
            Node = MergingYNode->Prev; // temp
            if (!Graph->FirstNode)
            {
                Graph->FirstNode = MergingYNode;
                Graph->LastNode = MergingYNode;
                MergingYNode->Next = NULL;
                MergingYNode->Prev = NULL;
            }
            else
            {
                MergingYNode->Prev = Graph->LastNode;
                Graph->LastNode->Next = MergingYNode;
                Graph->LastNode = MergingYNode;
            }
            MergingYNode->Next = NULL;
            MergingYNode = Node;            
        }
        else
        {
            // Insert the b node into the master list:
            Node = MergingBNode->Next; // temp
            if (!Graph->FirstNode)
            {
                Graph->FirstNode = MergingBNode;
                Graph->LastNode = MergingBNode;
                MergingBNode->Next = NULL;
                MergingBNode->Prev = NULL;
            }
            else
            {
                MergingBNode->Prev = Graph->LastNode;
                Graph->LastNode->Next = MergingBNode;
                Graph->LastNode = MergingBNode;
            }
            MergingBNode->Next = NULL;
            MergingBNode = Node;
        }
        Graph->NodeCount++;
    }

    TagGraphAddEndpointNodes(Graph, Spectrum);
    for (Node = Graph->FirstNode, NodeIndex = 0; Node; Node = Node->Next, NodeIndex++)
    {
        Node->Index = NodeIndex;
    }
}

// Insert another node into the tag-graph.  (Used only for a few nodes, as this isn't super fast)
void InsertTagGraphNode(TagGraph* Graph, TagGraphNode* Node)
{
    TagGraphNode* TempNode;
    // Iterate backwards, until either TempNode points to a smaller PRM or we fall of the edge of the list.
    for (TempNode = Graph->LastNode; TempNode && TempNode->Mass > Node->Mass; TempNode = TempNode->Prev)
    {
        ;;
    }
    if (!TempNode)
    { 
        // This new node is smaller than any we've seen.
        if (Graph->FirstNode)
        {
            Graph->FirstNode->Prev = Node;
        }
        else
        {
            Graph->LastNode = Node;
        }
        Node->Next = Graph->FirstNode;
        Graph->FirstNode = Node;
    }
    else if (TempNode->Next)
    {
        TempNode->Next->Prev = Node;
        Node->Next = TempNode->Next;
        Node->Prev = TempNode;
        TempNode->Next = Node;
    }
    else
    {
        Node->Prev = Graph->LastNode;
        Graph->LastNode->Next = Node;
        if (Graph->LastNode == Graph->FirstNode)
        {
            Graph->FirstNode = Node;
        }
        Graph->LastNode = Node;
    }
    Graph->NodeCount++;
}

// Add the "goalpost nodes" to our tag graph, mass 0 and parent mass:
void TagGraphAddEndpointNodes(TagGraph* Graph, MSSpectrum* Spectrum)
{
    TagGraphNode* Node;
    int ModType;

    // LEFT edge:
    Node = (TagGraphNode*)calloc(sizeof(TagGraphNode), 1);
    Node->Mass = 0;
    Node->NodeType = evGraphNodeLeft;
    Node->IonTypeFlags = ION_FLAG_B;
    Node->OriginalPeakIndex = -1;
    
    InsertTagGraphNode(Graph, Node);

    // LEFT EDGE plus N-terminal mod:
    for (ModType = 0; ModType < AllPTModCount; ModType++)
    {
        if (AllKnownPTMods[ModType].Flags & DELTA_FLAG_N_TERMINAL)
        {
            Node = (TagGraphNode*)calloc(sizeof(TagGraphNode), 1);
            Node->Mass = AllKnownPTMods[ModType].Mass;
            Node->NodeType = evGraphNodeLeftMod;
            Node->IonTypeFlags = ION_FLAG_B;
            Node->OriginalPeakIndex = -1;
            // The node stores a pointer to the MassDelta, so that the tag
            // will also include the MassDelta:
            Node->PTM = MassDeltaByIndex[MAX_PT_MODTYPE * MDBI_ALL_MODS + ModType];
            //Node->PTM = ModType;
            InsertTagGraphNode(Graph, Node);
        }
    }

    Node = (TagGraphNode*)calloc(sizeof(TagGraphNode), 1);
    Node->Mass = Spectrum->ParentMass - PARENT_MASS_BOOST;
    Node->NodeType = evGraphNodeRight;
    Node->IonTypeFlags = ION_FLAG_Y;
    Node->OriginalPeakIndex = -1;
    InsertTagGraphNode(Graph, Node);

    // RIGHT EDGE minus C-terminal PTM:
    for (ModType = 0; ModType < AllPTModCount; ModType++)
    {
        if (AllKnownPTMods[ModType].Flags & DELTA_FLAG_C_TERMINAL)
        {
            Node = (TagGraphNode*)calloc(sizeof(TagGraphNode), 1);
            Node->Mass = Spectrum->ParentMass - PARENT_MASS_BOOST - AllKnownPTMods[ModType].Mass;
            Node->NodeType = evGraphNodeRightMod;
            Node->IonTypeFlags = ION_FLAG_Y;
            Node->OriginalPeakIndex = -1;
            Node->PTM = MassDeltaByIndex[MAX_PT_MODTYPE * MDBI_ALL_MODS + ModType];
            //Node->PTM = ModType;
            InsertTagGraphNode(Graph, Node);
        }
    }
}

// Print all the PRM nodes from a tag graph.  (Handy for debugging tagging)
void DebugPrintTagGraph(MSSpectrum* Spectrum, TagGraph* Graph)
{
    TagGraphNode* Node;
    //
    printf(">->Printing tag graph...\n");
    for (Node = Graph->FirstNode; Node; Node = Node->Next)
    {
#ifdef DEBUG_TAG_GENERATION
        printf("%s\n", Node->VerboseNodeInfo);
#else
        printf("%.2f %.2f\n", Node->Mass / (float)MASS_SCALE, Node->Score);
        //printf("At %.2f node %d ion types %d score %.3f:\n", Node->Mass/100.0, Node->NodeType, Node->IonTypeFlags, Node->Score);
        //printf("  b%d y%d in%.2f io%.2f is%.2f\n", Node->IntensityRankB, Node->IntensityRankY, Node->IntensityScore, Node->IonTypeScore, Node->IsotopeScore);

        if (Node->BIndex > -1)
        {
        }
        if (Node->YIndex > -1)
        {
        }
#endif
        
    }
    printf("<-<End of tag graph.\n");
}

// The JumpingHash stores, for each mass (rounded to nearest integer), a list
// of amino acids (or modified amino acids) matching the mass.  When constructing
// tags, we allow a move from node A to node B if we find a jump whose size matches
// the mass difference between nodes A and B.  (We check three hash buckets, to
// compensate for roundoff screwery)
JumpNode* JumpingHashAddJump(int Mass, char Amino, MassDelta* Delta)
{
    int HashBucket;
    JumpNode* Node;
    JumpNode* NewNode;
    JumpNode* Prev;
    // HashBucket = Mass, rounded to nearest int
    FAST_ROUND(Mass / (float)MASS_SCALE, HashBucket); 
    if (HashBucket < 0 || HashBucket >= MAX_JUMPING_HASH)
    {
        printf("** ERROR: Bad mass in JumpingHashAddJump\n");
        printf("Mass %d amino %c delta %s\n", Mass, Amino, Delta->Name);
        return NULL;
    }
    // HashBucket = (int)Mass;
    //if (Mass > HashBucket + 0.5)
    //{
    //    HashBucket += 1;
    //}
    NewNode = (JumpNode*)calloc(sizeof(JumpNode), 1);
    NewNode->Amino = Amino;
    NewNode->Mass = Mass;
    NewNode->Delta = Delta;  // The PTM for this jump, or -1 if there's no mod.
    if (NewNode->Delta)
    {
        NewNode->Score = Delta->Score;
    }
    Node = JumpingHash[HashBucket];

    if (!Node)
    {
        // Add a brand new entry to the hash:
        JumpingHash[HashBucket] = NewNode;
    }
    else
    {
        // Add this jump to the end of the list.  
        // (Lists are short, so we don't bother keeping a tail pointer or two-way links)
        for (; Node; Node = Node->Next)
        {
            Prev = Node;
        }
        Prev->Next = NewNode;
    }
    return NewNode;
}

// Populate the jumping hash with each amino acid, and each modified amino acid.
void PopulateJumpingHash()
{
    int Amino;
    int PTModIndex;
    MassDelta* Delta;
    int ModForAAIndex;
    JumpNode* JNode;

    FreeJumpingHash(); // free any old stuff

    // Allocate memory:
    JumpingHash = (JumpNode**)calloc(sizeof(JumpNode*), MAX_JUMPING_HASH);
    SafeFree(JumpsByAA);
    JumpsByAA = (JumpNode**)calloc(sizeof(JumpNode*), AMINO_ACIDS * GlobalOptions->DeltasPerAA);

    memset(JumpsByAA, 0, sizeof(JumpNode*) * AMINO_ACIDS * GlobalOptions->DeltasPerAA);
    for (Amino = 'A'; Amino<='Y'; Amino++)
    {
        if (PeptideMass[Amino]<0.01)
        {
            continue; // Not an amino acid ("O" or "U" or "J" or somesuch)
        }
        ModForAAIndex = 0;
        // Don't build a jump node for unmodified Q or unmodified I, because they are accounted
        // for by the jumps for unmodified K and L.
        if (Amino != 'Q' && Amino != 'I')
        {
            JNode = JumpingHashAddJump(PeptideMass[Amino], (char)Amino, NULL);
            JumpsByAA[(Amino-'A') * GlobalOptions->DeltasPerAA] = JNode;
            ModForAAIndex = 1;
        }
        for (PTModIndex = 0; PTModIndex < GlobalOptions->DeltasPerAA; PTModIndex++)
        {
            Delta = &MassDeltas[Amino - 'A'][PTModIndex];
            
            if (Delta->Flags)
            {
                if (!(Delta->Flags & (DELTA_FLAG_C_TERMINAL | DELTA_FLAG_N_TERMINAL)))
                {
                    JNode = JumpingHashAddJump(Delta->RealDelta + PeptideMass[Amino], (char)Amino, Delta);
                    JumpsByAA[(Amino-'A') * GlobalOptions->DeltasPerAA + ModForAAIndex] = JNode;
                    ModForAAIndex++;
                }
            }
            else
            {
                // There are no more PTMs in MassDeltas[Amino], so stop iterating:
                break;
            }
        }
    }
}

// Destructor for the JumpingHash contents.
void FreeJumpingHash()
{
    int HashBucket;
    JumpNode* Node;
    JumpNode* Prev = NULL;
    //
    if (JumpingHash)
    {
        for (HashBucket = 0; HashBucket < MAX_JUMPING_HASH; HashBucket++)
        {
            Prev = NULL;
            for (Node = JumpingHash[HashBucket]; Node; Node = Node->Next)
            {
                SafeFree(Prev);
                Prev = Node;
            }
            SafeFree(Prev);
            JumpingHash[HashBucket] = NULL;
        }
        SafeFree(JumpingHash);
        JumpingHash = NULL;
    }
    SafeFree(JumpsByAA);
    JumpsByAA = NULL;
}

void DebugPrintTagGraphEdges(TagGraph* Graph)
{
    TagGraphNode* Node;
    TagGraphEdge* Edge;
    for (Node = Graph->FirstNode; Node; Node = Node->Next)
    {
        printf("Node at %.2f: (score %.2f)\n", Node->Mass / (float)MASS_SCALE, Node->Score);
        for (Edge = Node->FirstEdge; Edge; Edge = Edge->Next)
        {
            printf("-->Add '%c' (%.2f, skew %.2f) to reach %.2f\n", Edge->Jump->Amino, Edge->Jump->Mass / (float)MASS_SCALE,
                ((Edge->ToNode->Mass - Node->Mass) - (Edge->Jump->Mass))/(float)MASS_SCALE,
                Edge->ToNode->Mass / (float)MASS_SCALE);
        }
    }
}

// Called after populating the tag graph with nodes.
// Now we add edges between any two nodes that can be linked by a JUMP (an amino acid, or 
// an amino acid plus a decoration)
void TagGraphPopulateEdges(TagGraph* Graph)
{
    TagGraphNode* Node;
    TagGraphNode* OtherNode;
    int JumpSize;
    // For efficiency, we *never* consider a jump smaller or larger than these boundaries.
    // Note that glycine has mass 57.02, and tryptophan has size 186
    // (If there are PTMs, we do consider jumps of MaxJumpSize + MaxPTMMass; that's probably overkill)
    int MinJumpSize = GLYCINE_MASS - (DALTON * 2);
    int MaxJumpSize;
    int MaxAA;
    int ModIndex;
    int MaxSkew;
    TagGraphEdge* Edge;
    int IntSkew;
    int EdgeCount = 0;
    int ModJumpCount = 0;
    int AA;
    JumpNode* JNode;
    int Bucket;
    int HashBucket;
    MaxAA = 0;
    for (AA = 'A'; AA < 'X'; AA++)
    {
        MaxAA = max(MaxAA, PeptideMass[AA]);
    }
    MaxJumpSize = MaxAA;
    for (ModIndex = 0; ModIndex < AllPTModCount; ModIndex++)
    {
        MaxJumpSize = (int)max(MassDeltaByIndex[MAX_PT_MODTYPE * MDBI_ALL_MODS + ModIndex]->RealDelta + MaxAA, MaxJumpSize);
        for (AA = 0; AA < 26; AA++)
        {
            ModJumpCount += AllKnownPTMods[ModIndex].Allowed[AA];
        }
    }
    MaxJumpSize += GlobalOptions->ParentMassEpsilon;

    MaxSkew = sizeof(SkewHistoStep) / sizeof(double) - 1;

    // We do a double-loop over the graph to find all legal edges.
    for (Node = Graph->FirstNode; Node; Node = Node->Next)
    {
        if (Node->NodeType == evGraphNodeRight || Node->NodeType == evGraphNodeRightMod)
        {
            // This is a right-edge, so no edges emit from it:
            continue;
        }

        for (OtherNode = Node->Next; OtherNode; OtherNode = OtherNode->Next)
        {
            if (OtherNode->NodeType == evGraphNodeLeft || OtherNode->NodeType == evGraphNodeLeftMod)
            {
                // This is a left-edge, so no edges enter it:
                continue;
            }
            JumpSize = OtherNode->Mass - Node->Mass;
            if (JumpSize < MinJumpSize)
            {
                continue;
            }
            if (JumpSize > MaxJumpSize)
            {
                break;
            }
            FAST_ROUND(JumpSize / (float)MASS_SCALE, HashBucket);
            for (Bucket = HashBucket - 1; Bucket < HashBucket + 2; Bucket++)
            {
                if (Bucket < 0 || Bucket >= MAX_JUMPING_HASH)
                {
                    continue;
                }
                for (JNode = JumpingHash[Bucket]; JNode; JNode = JNode->Next)
                {
                    IntSkew = JumpSize - JNode->Mass;
                    if (abs(IntSkew) > GlobalOptions->Epsilon)
                    {
                        continue;
                    }
                    if (JNode->Delta)
                    {
                        if (GlobalOptions->TagPTMMode == 1 || GlobalOptions->MaxPTMods == 0 || JNode->Delta->Score < -5)
                        {
                            continue;
                        }
                    }
                    // Allocate a TagGraphEdge, initialize it, and add it to this node's list of edges:
                    Edge = (TagGraphEdge*)calloc(sizeof(TagGraphEdge), 1);
                    Edge->Jump = JNode;
                    Edge->FromNode = Node;
                    Edge->ToNode = OtherNode;
                    Edge->Skew = IntSkew;
                    // For now, no skew scoring:
                    //if (IntSkew > MaxSkew)
                    //{
                    //    Edge->Score = (float)SkewHistoStep[MaxSkew];
                    //}
                    //else
                    //{
                    //    Edge->Score = (float)SkewHistoStep[IntSkew];
                    //}
                    Edge->Score = JNode->Score;
                    if (!Node->FirstEdge)
                    {
                        Node->FirstEdge = Edge;
                        Node->LastEdge = Edge;
                    }
                    else
                    {
                        Node->LastEdge->Next = Edge;
                        Node->LastEdge = Edge;
                    }
                    if (Edge->Jump->Delta)
                    {
                        Edge->Score += MasterTaggingModel.PTMPenalty;
                    }
                    //GlobalStats->TagGraphEdges++;
                    EdgeCount++;
                } // Jnode loop
            } // bucket loop
        }
    }
}

// For quick-sort of tags - list from highest to lowest score.
int CompareTagScores(const TrieTag* TagA, const TrieTag* TagB)
{
    if (TagA->Score > TagB->Score)
    {
        return -1;
    }
    if (TagA->Score < TagB->Score)
    {
        return 1;
    }
    if (TagA->PrefixMass < TagB->PrefixMass)
    {
        return -1;
    }
    if (TagA->PrefixMass > TagB->PrefixMass)
    {
        return 1;
    }
    if (TagA < TagB)
    {
        return -1;
    }
    if (TagA > TagB)
    {
        return -1;
    }
    return 0;
}

int TagSkewBinCount;
float* TagSkewScore = NULL;
float* TagTotalAbsSkewScore = NULL;

void FreeTagSkewScores()
{
    SafeFree(TagSkewScore);
    TagSkewScore = NULL;
    SafeFree(TagTotalAbsSkewScore);
    TagTotalAbsSkewScore = NULL;
}

void SetTagSkewScores()
{
    char FilePath[2048];
    FILE* TagSkewFile;
    //
    if (TagSkewScore)
    {
        return;
    }
    sprintf(FilePath, "%s%s", GlobalOptions->ResourceDir, "TagSkewScores.dat");
    TagSkewFile = fopen(FilePath, "rb");
    if (!TagSkewFile)
    {
        REPORT_ERROR_S(3, FilePath);
        // To avoid crashing later, set up a length-1 array:
        TagSkewBinCount = 1;
        TagSkewScore = (float*)calloc(1, sizeof(float));
        TagSkewScore[0] = 0;
        TagTotalAbsSkewScore = (float*)calloc(1, sizeof(float));
        TagTotalAbsSkewScore[0] = 0;
        return;
    }
    // Read the number of entries:
    ReadBinary(&TagSkewBinCount, sizeof(int), 1, TagSkewFile);
    // Allocate arrays:
    TagSkewScore = (float*)calloc(TagSkewBinCount, sizeof(float));
    TagTotalAbsSkewScore = (float*)calloc(TagSkewBinCount, sizeof(float));
    // Populate arrays:
    ReadBinary(TagSkewScore, sizeof(float), TagSkewBinCount, TagSkewFile);
    ReadBinary(TagTotalAbsSkewScore, sizeof(float), TagSkewBinCount, TagSkewFile);
    fclose(TagSkewFile);
}

static TrieTag* AllTags = NULL;
//// New tag generation function: Generates tags of a (more-or-less) arbitrary length!
//TrieTag* TagGraphGenerateTagsOld(TagGraph* Graph, MSSpectrum* Spectrum, int* TagCount, 
//    int MaximumTagCount, SpectrumTweak* Tweak, float TagEdgeScoreMultiplier)
//{
//    TagGraphNode* TagNodes[12];
//    TagGraphEdge* TagEdges[12];
//    int NodeIndex;
//    int EdgeIndex;
//    TagGraphNode* Node;
//    TagGraphEdge* Edge;
//    TagGraphNode* LeftNode;
//    TagGraphNode* RightNode;
//    int CurrentDepth;
//    int InternalNodes;
//    float NodeScore;
//    TrieTag* Tag;
//    int TagAllocation;
//    int BacktrackFlag;
//    int AllTagCount = 0;
//    int Bin;
//    float ScoreToBeat = -9999;
//    //
//    *TagCount = 0;
//    TagAllocation = 1024;
//    if (!AllTags)
//    {
//        AllTags = (TrieTag*)calloc(sizeof(TrieTag), TagAllocation);
//    }
//    NodeIndex = 0;
//    EdgeIndex = -1;
//    BacktrackFlag = 0;
//    CurrentDepth = 0;
//    TagNodes[0] = Graph->FirstNode;
//    while (1)
//    {
//        // If we're BACKTRACKING, then move to a sibling or parent:
//        if (BacktrackFlag)
//        {
//            // Move the root of the subtree, if necessary:
//            if (CurrentDepth == 0)
//            {
//                // Move to the next 'first' node:
//                TagNodes[0] = TagNodes[0]->Next;
//                if (!TagNodes[0])
//                {
//                    break;
//                }
//                BacktrackFlag = 0;
//                continue;
//            }
//            // Move to a sibling, if we can:
//            TagEdges[CurrentDepth - 1] = TagEdges[CurrentDepth - 1]->Next;
//            if (TagEdges[CurrentDepth - 1])
//            {
//                TagNodes[CurrentDepth] = TagEdges[CurrentDepth - 1]->ToNode;
//                BacktrackFlag = 0;
//                continue;
//            }
//            // No more siblings - move up one level.
//            CurrentDepth--;
//            continue;
//        }
//
//        // Special case for level 1: Skip tag nodes with silly masses like 20Da.
//        if (CurrentDepth == 0)
//        {
//            Node = TagNodes[0];
//            if (Node->Mass > GlobalOptions->ParentMassEpsilon && Node->Mass < GLYCINE_MASS - GlobalOptions->Epsilon)
//            {
//                BacktrackFlag = 1;
//                continue;
//            }
//        }
//
//        // If we're deep enough, report a tag and start backtracking:
//        if (CurrentDepth >= GlobalOptions->GenerateTagLength)
//        {
//            BacktrackFlag = 1;
//            LeftNode = TagNodes[0];
//            RightNode = TagNodes[CurrentDepth];
//            Tag = AllTags + (*TagCount);
//            InternalNodes = 0;
//            NodeScore = 0;
//            for (NodeIndex = 0; NodeIndex <= CurrentDepth; NodeIndex++)
//            {
//                Node = TagNodes[NodeIndex];
//                if (Node->OriginalPeakIndex > 0)
//                {
//                    NodeScore += Node->Score;
//                    InternalNodes++;
//                }
//                Tag->Nodes[NodeIndex] = TagNodes[NodeIndex];
//            }
//            NodeScore *= (GlobalOptions->GenerateTagLength + 1) / (float)max(1, InternalNodes);
//            Tag->Score = NodeScore;
//            Tag->ModsUsed = 0;
//            memset(Tag->ModType, 0, sizeof(MassDelta*) * MAX_PT_MODS);
//            memset(Tag->AminoIndex, -1, sizeof(int) * MAX_PT_MODS);
//            if (LeftNode->NodeType == evGraphNodeLeftMod)
//            {
//                Tag->AminoIndex[Tag->ModsUsed] = 0;
//                Tag->ModType[Tag->ModsUsed] = LeftNode->PTM;
//                Tag->ModsUsed++;
//            }
//            for (EdgeIndex = 0; EdgeIndex < CurrentDepth; EdgeIndex++)
//            {
//                Edge = TagEdges[EdgeIndex];
//                Tag->Score += TagEdges[EdgeIndex]->Score;
//                Tag->Tag[EdgeIndex] = Edge->Jump->Amino;
//                if (Edge->Jump->Delta)
//                {
//                    Tag->AminoIndex[Tag->ModsUsed] = EdgeIndex;
//                    Tag->ModType[Tag->ModsUsed] = Edge->Jump->Delta;
//                    Tag->ModsUsed++;
//                }
//            }
//            // Set skew info:
//            Tag->TotalSkew = 0;
//            Tag->TotalAbsSkew = 0;
//            for (EdgeIndex = 0; EdgeIndex < CurrentDepth; EdgeIndex++)
//            {
//                Edge = TagEdges[EdgeIndex];
//                Tag->TotalSkew += Edge->Skew;
//                Tag->TotalAbsSkew += abs(Edge->Skew);
//            }
//            ////////////////////////////////////////////////////
//            // If the total skew is large, penalize the tag's score:
//            Bin = (int)fabs((Tag->TotalSkew / 50.0) + 0.5);
//            if (Bin >= TagSkewBinCount)
//            {
//                Bin = TagSkewBinCount - 1;
//            }
//            Tag->Score += TagSkewScore[Bin] * TagEdgeScoreMultiplier;
//            Bin = (int)fabs((Tag->TotalAbsSkew / 50.0) + 0.5);
//            if (Bin >= TagSkewBinCount)
//            {
//                Bin = TagSkewBinCount - 1;
//            }
//            Tag->Score += TagTotalAbsSkewScore[Bin] * TagEdgeScoreMultiplier;
//            ////////////////////////////////////////////////////
//            Tag->Tag[EdgeIndex] = '\0';
//            if (Tag->Score < ScoreToBeat)
//            {
//                // Abort the tag - it's not good enough!
//                continue; 
//            }
//            if (RightNode->NodeType == evGraphNodeRightMod)
//            {
//                Tag->AminoIndex[Tag->ModsUsed] = CurrentDepth;
//                Tag->ModType[Tag->ModsUsed] = RightNode->PTM;
//                Tag->ModsUsed++;
//            }
//            Tag->PSpectrum = Spectrum;
//            Tag->Tweak = Tweak;
//            Tag->TagLength = CurrentDepth;
//            Tag->ParentMass = Spectrum->ParentMass;
//            Tag->Charge = Spectrum->Charge;
//            Tag->PrefixMass = TagNodes[0]->Mass;
//            Tag->SuffixMass = Spectrum->ParentMass - PARENT_MASS_BOOST - TagNodes[CurrentDepth]->Mass;
//            (*TagCount)++;
//            AllTagCount++;
//            // If we've got as many tags as we can handle, drop all but the best.  (Don't
//            // just reallocate; we could end up with a *lot*!)
//            if ((*TagCount) + 5 >= TagAllocation)
//            {
//                qsort(AllTags, *TagCount, sizeof(TrieTag), (QSortCompare)CompareTagScores);
//                *TagCount = TagAllocation / 2;
//                if (MaximumTagCount >= 0)
//                {
//                    ScoreToBeat = AllTags[min(TagAllocation - 5, MaximumTagCount)].Score;
//                }
//                else
//                {
//                    ScoreToBeat = AllTags[*TagCount].Score;
//                }
//            }
//            continue;
//        } // If we're at tag depth
//
//        // We're not at tag depth yet. 
//        // Move to our first available child:
//        TagEdges[CurrentDepth] = TagNodes[CurrentDepth]->FirstEdge;
//        if (!TagEdges[CurrentDepth])
//        {
//            BacktrackFlag = 1;
//            continue;
//        }
//        else
//        {
//            CurrentDepth++;
//            TagNodes[CurrentDepth] = TagEdges[CurrentDepth - 1]->ToNode;
//        }
//    }
//    // Sort the tags, by score:
//    qsort(AllTags, *TagCount, sizeof(TrieTag), (QSortCompare)CompareTagScores);
//    return AllTags;
//
//}

// Build a trie from a list of tags.  Returns the trie root. 
// AllTags is the tag array, TagCount its size.
// Since we construct one big trie for many spectra, we take Root as an
// argument; Root is NULL on the first call.
// If GlobalOptions->GenerateCount is >= 0 and < TagCount, then we stop after adding
// that many tags.
TrieNode* BuildTrieFromTags(TrieTag* AllTags, int TagCount, TrieNode* Root, int MaximumTagCount)
{
    int DuplicateFlag;
    int TagsInTrie = 0;
    int TagIndex;
    TrieTag* Tag;

    // Construct a root, if we don't have one already.  
    if (!Root)
    {
        Root = NewTrieNode();
        Root->FailureNode = Root;
    }
    for (TagIndex = 0; TagIndex < TagCount; TagIndex++)
    {
        AddTagToTrie(Root, AllTags + TagIndex, &DuplicateFlag);
        if (!DuplicateFlag)
        {
            TagsInTrie++;
            Tag = AllTags + TagIndex;
            if (MaximumTagCount >= 0 && TagsInTrie >= MaximumTagCount)
            {
                break;
            }
        }
    }
    //DebugPrintTrieTags(Root);
    return Root;
}

void DebugPrintTagList(MSSpectrum* Spectrum, TrieTag* Tags, int TagCount)
{
    int TagIndex;
    TrieTag* Tag;
    int Index;
    for (TagIndex = 0; TagIndex < TagCount; TagIndex++)
    {
        Tag = Tags + TagIndex;
#ifdef DEBUG_TAG_GENERATION
        printf("%s\n", Tag->TagScoreDetails);
#endif
        printf("%d: %.2f: %s %.2f %.2f\n", TagIndex, Tag->Score, Tag->Tag, Tag->PrefixMass / (float)MASS_SCALE, (Spectrum->ParentMass - PARENT_MASS_BOOST - Tag->SuffixMass) / (float)MASS_SCALE);
        for (Index = 0; Index < Tag->TagLength; Index++)
        {
            printf("%c", Tag->Tag[Index]);
            fflush(stdout);
            if (Tag->AminoIndex[Index]>-1)
            {
                printf("%s", Tag->ModType[Index]->Name);
                fflush(stdout);
            }
        }
        printf("\n");
    }
}

// Called when searching in tagless mode.  (Tagless mode performs *no* database filtering; it's 
// appropriate for searching a small database, typically a database formed by an initial search run) 
// The trie, in this case, will have a child for each amino acid (prefix 0, suffix = parent mass - amino mass)
TrieNode* GenerateDummyTags(MSSpectrum* Spectrum, TrieNode* Root)
{
    TrieTag* Tag;
    char* Aminos = "ACDEFGHKLMNPRSTVWY"; // skip I and Q, because they're synonymous with L and K
    char* Amino;
    int DuplicateFlag;
    int ModIndex;
    int TweakIndex;
    SpectrumTweak* Tweak;
    // Set up the root, if it doesn't exist already:
    if (!Root)
    {
        Root = NewTrieNode();
        Root->FailureNode = Root;
    }
    for (TweakIndex = 0; TweakIndex < TWEAK_COUNT; TweakIndex++)
    {
        Tweak = Spectrum->Node->Tweaks + TweakIndex;
        if (!Tweak->Charge)
        {
            continue;
        }
        // Loop over alphabet soup, add one tag per amino:
        for (Amino = Aminos; *Amino; Amino++)
        {
            Tag = NewTrieTag();
            Tag->Tag[0] = *Amino;
            Tag->Tag[1] = '\0';
            memset(Tag->ModType, 0, sizeof(MassDelta*) * MAX_PT_MODS);
            memset(Tag->AminoIndex, -1, sizeof(int) * MAX_PT_MODS);
            Tag->PSpectrum = Spectrum;
            Tag->Charge = Tweak->Charge;
            Tag->ParentMass = Tweak->ParentMass;
            Tag->Tweak = Tweak;
            Tag->PrefixMass = 0;
            Tag->SuffixMass = Tweak->ParentMass - PeptideMass[*Amino] - PARENT_MASS_BOOST;
            Tag->TagLength = 1;
            //GlobalStats->TagsGenerated++;
            AddTagToTrie(Root, Tag, &DuplicateFlag);
            // ...ok, also allow mods on this first amino
            for (ModIndex = 0; ModIndex < AllPTModCount; ModIndex++)
            {
                if (AllKnownPTMods[ModIndex].Allowed[*Amino - 'A'])
                {
                    Tag = NewTrieTag();
                    Tag->Tag[0] = *Amino;
                    Tag->Tag[1] = '\0';
                    memset(Tag->ModType, 0, sizeof(MassDelta*) * MAX_PT_MODS);
                    memset(Tag->AminoIndex, -1, sizeof(int) * MAX_PT_MODS);
                    Tag->ModType[0] = MassDeltaByIndex[(*Amino-'A') * MAX_PT_MODTYPE + ModIndex];
                    Tag->AminoIndex[0] = 0;
                    Tag->ModsUsed = 1;
                    Tag->PSpectrum = Spectrum;
                    Tag->Charge = Tweak->Charge;
                    Tag->Tweak = Tweak;
                    Tag->ParentMass = Tweak->ParentMass;
                    Tag->PrefixMass = 0;
                    Tag->SuffixMass = Tweak->ParentMass - PeptideMass[*Amino] - PARENT_MASS_BOOST - AllKnownPTMods[ModIndex].Mass;
                    Tag->TagLength = 1;
                    //GlobalStats->TagsGenerated++;
                    AddTagToTrie(Root, Tag, &DuplicateFlag);
                }
            }
        }
    }
    return Root;
}

TrieNode* GenerateTagsFromSpectrum(MSSpectrum* Spectrum, TrieNode* Root, int MaximumTagCount, SpectrumTweak* Tweak)
{
    TrieTag* Tags;
    int TagCount;
    
    // Note: Spectrum load and preprocessing methods need to be called before calling this function.
    // Call these: 
    //SpectrumFindIsotopicPeaks(Spectrum);
    //IntensityRankPeaks(Spectrum);
    //SpectrumCorrectParentMass(Spectrum);

    if (Spectrum->Graph)
    {
        FreeTagGraph(Spectrum->Graph);
        Spectrum->Graph = NULL;
    }
    if (GlobalOptions->TaglessSearchFlag)
    {
        return GenerateDummyTags(Spectrum, Root);
    }
    Spectrum->Graph = ConstructTagGraph(Spectrum);
    TagGraphAddNodes(Spectrum->Graph, Spectrum);
    TagGraphScorePRMNodes(NULL, Spectrum->Graph, Spectrum, Tweak);
    //DebugPrintTagGraph(Spectrum, Spectrum->Graph); 
    TagGraphPopulateEdges(Spectrum->Graph);
    
#ifdef DEBUG_TAG_GENERATION
    DebugPrintTagGraph(Spectrum, Spectrum->Graph);
    DebugPrintTagGraphEdges(Spectrum->Graph); ////
#endif
    Tags = TagGraphGenerateTags(Spectrum->Graph, Spectrum, &TagCount, MaximumTagCount, Tweak, TAG_EDGE_SCORE_MULTIPLIER, NULL);

#ifdef DEBUG_TAG_GENERATION
    DebugPrintTagList(Spectrum, Tags, 300);
#endif
    DebugPrintTagsForPeptide(Spectrum, Spectrum->Graph, Tags, TagCount);
    Root = BuildTrieFromTags(Tags, TagCount, Root, MaximumTagCount);
    if (0)
    {
        DebugPrintTrieTags(Root);
    }

    // The caller should usually invoke InitializeTrieFailureNodes next.  When doing a batch of 
    // spectra, however, we do InitializeTrieFailureNodes once at the end.

    return Root;
}

// Build a hash (Graph->NodeIndex) for quick lookup of nodes based on mass.  This is used
// in GetBYScore, when choosing a PTM attachment point.  
void TagGraphBuildNodeIndex(TagGraph* Graph)
{
    TagGraphNode* Node;
    int Bucket;
    int BucketMax;
    SafeFree(Graph->NodeIndex);
    Graph->NodeIndexSize = ((int)(Graph->LastNode->Mass / DALTON)) + 1;
    Graph->NodeIndex = (TagGraphNode**)calloc(sizeof(TagGraphNode*), Graph->NodeIndexSize);
    for (Node = Graph->FirstNode; Node; Node = Node->Next)
    {
        BucketMax = (int)min(Graph->NodeIndexSize - 1, Node->Mass / DALTON + 1);
        for (Bucket = max(0, ((int)Node->Mass / DALTON) - 1); Bucket <= BucketMax; Bucket++)
        {
            if (!Graph->NodeIndex[Bucket])
            {
                Graph->NodeIndex[Bucket] = Node;
            }
        }
    }
}

int NiceCheckAA(char AA1, char AA2)
{
    if (AA1 == 'I')
    {
        AA1= 'L';
    }
    if (AA2 == 'I')
    {
        AA2= 'L';
    }
    if (AA1 == 'Q')
    {
        AA1= 'K';
    }
    if (AA2 == 'Q')
    {
        AA2= 'K';
    }
    return (AA1 == AA2);
}

void DebugCheckTagMatch(int TagIndex, TrieTag* Tag, int* Masses, int MassCount, char* Peptide)
{
    int Pos;
    int Diff;
    int ParentMass;
    int TagAAIndex;
    //
    ParentMass = Masses[MassCount-1];
    for (Pos = 0; Pos < MassCount - 3; Pos++)
    {
        Diff = abs(Tag->PrefixMass - Masses[Pos]);
        if (Diff > 2 * DALTON)
        {
            continue;
        }
        Diff = abs((ParentMass - Masses[Pos+3]) - Tag->SuffixMass);
        if (Diff > 2 * DALTON)
        {
            continue;
        }
        for (TagAAIndex = 0; TagAAIndex < Tag->TagLength; TagAAIndex++)
        {
            if (!NiceCheckAA(Peptide[Pos + TagAAIndex + 1], Tag->Tag[TagAAIndex]))
            {
                continue;
            }
        }
        printf("Matched by tag #%d: '%s', prefix %.2f, suffix %.2f\n", TagIndex, Tag->Tag, Tag->PrefixMass / (float)MASS_SCALE, Tag->SuffixMass / (float)MASS_SCALE);
    }
}

// Sometimes we don't generate tags for a peptide, and it's not obvious why.
// In such casses, include a line of the form "tagcheck,PEPTIDE" in the input file.
// Then, this function will compare the tags for this peptide against the actual tags
// (and actual tag graph).
void DebugPrintTagsForPeptide(MSSpectrum* Spectrum, TagGraph* Graph, TrieTag* Tags, int TagCount)
{
    StringNode* Node;
    int MassCount;
    int Masses[64];
    char* Amino;
    int AminoMass;
    int AccumMass;
    char Peptide[64];
    int PeptideLength;
    int TagIndex;
    TrieTag* Tag;
    int MassIndex;
    TagGraphNode* GraphNode;
    //
    for (Node = FirstTagCheckNode; Node; Node = Node->Next)
    {
        printf("--- Check tagging results for %s Charge %d\n", Node->String, Spectrum->Charge);
        MassCount = 0;
        AccumMass = 0;
        PeptideLength = 0;
        // Parse the peptide string.  For now, DROP all mods.
        for (Amino = Node->String; *Amino; Amino++)
        {
            AminoMass = PeptideMass[*Amino];
            if (AminoMass)
            {
                AccumMass += AminoMass;
                Masses[MassCount++] = AccumMass;
                Peptide[PeptideLength++] = *Amino;
            }
        }
        Peptide[PeptideLength] = '\0';
        ///////////////////////////////////////////////////////////
        // Loop over tags, and see whether any tag matches the peptide:
        for (TagIndex = 0; TagIndex < TagCount; TagIndex++)
        {   
            Tag = Tags + TagIndex;
            DebugCheckTagMatch(TagIndex, Tag, Masses, MassCount, Peptide);
        }
        ///////////////////////////////////////////////////////////
        // Loop over PRMs in the peptide, and see how taggable they are:
        for (MassIndex = 0; MassIndex < MassCount; MassIndex++)
        {
            AccumMass = Masses[MassIndex];
            printf("Mass %d (%.2f):\n", MassIndex, AccumMass / (float)MASS_SCALE);
            for (GraphNode = Graph->FirstNode; GraphNode; GraphNode = GraphNode->Next)
            {
                if (GraphNode->Mass > AccumMass + DALTON)
                {
                    break;
                }
                if (GraphNode->Mass < AccumMass - DALTON)
                {
                    continue;
                }
                printf("  Node at %.2f (%.2f) score %.2f\n", GraphNode->Mass / (float)MASS_SCALE, 
                    (GraphNode->Mass - AccumMass) / (float)MASS_SCALE, GraphNode->Score);
            }
        }
    }
}
TagGraphNode* TagTestGetBestNode(TagGraph* Graph, int PRM)
{
    int MinMass;
    int MaxMass;
    TagGraphNode* TGNode;
    TagGraphNode* BestNode = NULL;
    //
    MinMass = PRM - 50;
    MaxMass = PRM + 50;
    for (TGNode = Graph->FirstNode; TGNode; TGNode = TGNode->Next)
    {
        if (TGNode->Mass > MaxMass)
        {
            break;
        }
        if (TGNode->Mass < MinMass)
        {
            continue;
        }
        if (!BestNode || BestNode->Score < TGNode->Score)
        {
            BestNode = TGNode;
        }
    }
    return BestNode;
}

void TestTaggingCallback(SpectrumNode* Node, int Charge, int ParentMass, Peptide* Annotation)
{
    static int* TrueTagRankHistogram = NULL;
    static int SpectrumCount = 0;
    int Rank;
    int Cumulative;
    FILE* ResultsFile;
    int TagIndex;
    int TagCount;
    TrieTag* Tags;
    int FoundFlag;
    TrieTag* TestTag;
    int PRM[64];
    int Mass;
    int AminoIndex;
    int ModIndex;
    int TrieTagCount;
    int MatchLength;
    BayesianModel* Model;
    TagGraphNode* Node0;
    TagGraphNode* Node1;
    TagGraphNode* Node2;
    TagGraphNode* Node3;
    float TagScore;
    TrieNode* Root;
    int DuplicateFlag;
    int VerboseFlag = 0;
    //
    Root = NULL;
    if (!Node)
    {
        if (!Charge)
        {
            // Initialization call:
            TrueTagRankHistogram = (int*)calloc(sizeof(int), 512);
        }
        else
        {
            // Completion call:
            ResultsFile = fopen("TagTestingResults.txt", "w");
            Cumulative = 0;
            fprintf(ResultsFile, "Tagging results on %d spectra\n", SpectrumCount);
            for (Rank = 0; Rank < 512; Rank++)
            {
                Cumulative += TrueTagRankHistogram[Rank];
                fprintf(ResultsFile, "%d\t%d\t%.2f\t%.2f\t\n",
                    Rank, TrueTagRankHistogram[Rank], TrueTagRankHistogram[Rank] / (float)SpectrumCount,
                    Cumulative / (float)SpectrumCount);
            }
            //SafeFree(TrueTagRankHistogram);
        }
        return;
    }
    // Standard call: Given a spectrum, generate some tags.  Remember the rank of the first true tag.

    Root = NewTrieNode();
    Root->FailureNode = Root;

    Node->Tweaks[0].Charge = Charge;
    Node->Tweaks[0].ParentMass = Annotation->ParentMass;
    Node->Spectrum->Charge = Charge;
    Node->Spectrum->ParentMass = Annotation->ParentMass;
    WindowFilterPeaks(Node->Spectrum, 0, 0);
    PrepareSpectrumForIonScoring(PRMModelCharge2, Node->Spectrum, 0);
    //SpectrumComputeBinnedIntensities(Node);
    //SpectrumComputeNoiseDistributions(Node);
    //SpectrumAssignIsotopeNeighbors(Node->Spectrum);
    //SpectrumFindIsotopicPeaks(Node->Spectrum);
    Node->Spectrum->Graph = ConstructTagGraph(Node->Spectrum);
    TagGraphAddNodes(Node->Spectrum->Graph, Node->Spectrum);
    TagGraphPopulateEdges(Node->Spectrum->Graph);
    TagGraphScorePRMNodes(NULL, Node->Spectrum->Graph, Node->Spectrum, Node->Tweaks);
    Tags = TagGraphGenerateTags(Node->Spectrum->Graph, Node->Spectrum, &TagCount, 1024, Node->Tweaks, TAG_EDGE_SCORE_MULTIPLIER, NULL);
    if (Charge > 2)
    {
        Model = BNCharge3TaggingBN;
    }
    else
    {
        Model = BNCharge2TaggingBN;
    }
    // Set our PRM array, so we can check tag prefix masses:
    //printf("\nTags for: %s\n", Annotation->Bases); 
    Mass = 0;
    MatchLength = strlen(Annotation->Bases);
    for (AminoIndex = 0; AminoIndex < MatchLength-1; AminoIndex++)
    {
        switch (Annotation->Bases[AminoIndex])
        {
        case 'I':
            Annotation->Bases[AminoIndex] = 'L';
            break;
        case 'Q':
            Annotation->Bases[AminoIndex] = 'K';
            break;
        default:
            break;
        }
        PRM[AminoIndex] = Mass;
        Mass += PeptideMass[Annotation->Bases[AminoIndex]];
        for (ModIndex = 0; ModIndex < MAX_PT_MODS; ModIndex++)
        {
            if (Annotation->AminoIndex[ModIndex] == AminoIndex && Annotation->ModType[ModIndex])
            {
                Mass +=Annotation->ModType[ModIndex]->RealDelta;
            }
        }
    }
    ///////////////////////////////////////////////////////////////
    // Optional verbose debugging:
    // For each theoretical tag, look for the best tag that can be generated.
    if (VerboseFlag)
    {
        //GetPRMFeatures(Node->Spectrum, Node->Tweaks, Model, 97870);
        //DebugPrintScorpPRMScores(Node->Spectrum, Node->Tweaks);
        for (AminoIndex = 0; AminoIndex < MatchLength-2; AminoIndex++)
        {
            Node0 = TagTestGetBestNode(Node->Spectrum->Graph, PRM[AminoIndex]);
            Node1 = TagTestGetBestNode(Node->Spectrum->Graph, PRM[AminoIndex + 1]);
            Node2 = TagTestGetBestNode(Node->Spectrum->Graph, PRM[AminoIndex + 2]);
            Node3 = TagTestGetBestNode(Node->Spectrum->Graph, PRM[AminoIndex + 3]);
            printf("Theoretical tag %.2f %s:\n", PRM[AminoIndex] / (float)MASS_SCALE, Annotation->Bases + AminoIndex);
            if (Node0)
            {
                TagScore = Node0->Score;
                printf(" Node0 %.2f score %.2f\n", Node0->Mass / (float)MASS_SCALE, Node0->Score);
            }
            else
            {
                printf(" <Node0 missing>\n");
                TagScore = -9999;
            }
            if (Node1)
            {
                TagScore += Node1->Score;
                printf(" Node1 %.2f score %.2f\n", Node1->Mass / (float)MASS_SCALE, Node1->Score);
            }
            else
            {
                printf(" <Node1 missing>\n");
                TagScore = -9999;
            }
            if (Node2)
            {
                TagScore += Node2->Score;
                printf(" Node2 %.2f score %.2f\n", Node2->Mass / (float)MASS_SCALE, Node2->Score);
            }
            else
            {
                printf(" <Node2 missing>\n");
                TagScore = -9999;
            }
            if (Node3)
            {
                TagScore += Node3->Score;
                printf(" Node3 %.2f score %.2f\n", Node3->Mass / (float)MASS_SCALE, Node3->Score);
            }
            else
            {
                printf(" <Node3 missing>\n");
                TagScore = -9999;
            }
            if (Node0 && Node0->OriginalPeakIndex < 0)
            {
                TagScore *= (float)1.3333;
            }
            if (Node3 && Node3->OriginalPeakIndex < 0)
            {
                TagScore *= (float)1.3333;
            }

            printf("overall: %.2f\n", TagScore);
        }
    }
    ///////////////////////////////////////////////////////////////
    // Check each tag to see whether it's correct:
    TagCount = min(TagCount, 512);
    FoundFlag = 0;
    TrieTagCount = 0;
    for (TagIndex = 0; TagIndex < TagCount; TagIndex++)
    {
        TestTag = Tags + TagIndex;
        DuplicateFlag = 0;
        AddTagToTrie(Root, TestTag, &DuplicateFlag);
        if (!DuplicateFlag)
        {
            TrieTagCount++;
            //if (TrieTagCount <= 10)
            {
                if (VerboseFlag)
                {
                    printf("%.2f\t%s\t%.2f\n", TestTag->PrefixMass / (float)MASS_SCALE, TestTag->Tag, TestTag->Score);
                }
                
            }
            for (AminoIndex = 0; AminoIndex < MatchLength-2; AminoIndex++)
            {
                if (abs(TestTag->PrefixMass - PRM[AminoIndex]) < GlobalOptions->ParentMassEpsilon)
                {
                    if (TestTag->Tag[0] == Annotation->Bases[AminoIndex] &&
                        TestTag->Tag[1] == Annotation->Bases[AminoIndex + 1] &&
                        TestTag->Tag[2] == Annotation->Bases[AminoIndex + 2])
                    {
                        TrueTagRankHistogram[TrieTagCount]++;
                        FoundFlag = 1;
                    }
                }
            }
        }
        if (FoundFlag)
        {
            break;
        }
    }
    if (!FoundFlag)
    {
        // we missed, too bad.  don't need to poke the histogram.
        //TrueTagRankHistogram[511]++;
    }
    FreeTrieNode(Root);
    Root = NULL;

    SpectrumCount++;

}

void TestTagging(char* OracleFile, char* OracleDir)
{
    InitBayesianModels(); // to use new PRM scoring
    InitStats();
    TestTaggingCallback(NULL, 0, 0, NULL); // initialization
    TrainOnOracleFile(OracleFile, OracleDir, TestTaggingCallback);
    TestTaggingCallback(NULL, 1, 0, NULL); // completion
}

void TrainTaggingCallback(SpectrumNode* Node, int Charge, int ParentMass, Peptide* Annotation)
{
    static int SpectrumCount = 0;
    int TagIndex;
    int TagCount;
    TrieTag* Tags;
    int FoundFlag;
    TrieTag* TestTag;
    int PRM[64];
    int Mass;
    int AminoIndex;
    int ModIndex;
    int TrieTagCount;
    int MatchLength;
    int FeatureIndex;
    BayesianModel* Model;
    TrieNode* Root;
    int DuplicateFlag;
    int TrueTagFlag;
    static FILE* TagTrainingFile = NULL;
    //
    Root = NULL;

    if (!TagTrainingFile)
    {
        TagTrainingFile = fopen("TagTraining.txt", "w");
    }
    // Standard call: Given a spectrum, generate some tags.  Test the first n tags, and write
    // out a feature-vector for each.

    // ** skip modded peptides:
    if (Annotation->ModType[0])
    {
        return;
    }
    Root = NewTrieNode();
    Root->FailureNode = Root;

    Node->Tweaks[0].Charge = Charge;
    Node->Tweaks[0].ParentMass = Annotation->ParentMass;
    Node->Spectrum->Charge = Charge;
    Node->Spectrum->ParentMass = Annotation->ParentMass;
    WindowFilterPeaks(Node->Spectrum, 0, 0);
    PrepareSpectrumForIonScoring(PRMModelCharge2, Node->Spectrum, 0);
    //SpectrumComputeBinnedIntensities(Node);
    //SpectrumComputeNoiseDistributions(Node);
    //SpectrumAssignIsotopeNeighbors(Node->Spectrum);
    //SpectrumFindIsotopicPeaks(Node->Spectrum);
    Node->Spectrum->Graph = ConstructTagGraph(Node->Spectrum);
    TagGraphAddNodes(Node->Spectrum->Graph, Node->Spectrum);
    TagGraphPopulateEdges(Node->Spectrum->Graph);
    TagGraphScorePRMNodes(NULL, Node->Spectrum->Graph, Node->Spectrum, Node->Tweaks);
    Tags = TagGraphGenerateTags(Node->Spectrum->Graph, Node->Spectrum, &TagCount, 1024, Node->Tweaks, TAG_EDGE_SCORE_MULTIPLIER, NULL);
    if (Charge > 2)
    {
        Model = BNCharge3TaggingBN;
    }
    else
    {
        Model = BNCharge2TaggingBN;
    }
    // Set our PRM array, so we can check tag prefix masses:
    Mass = 0;
    MatchLength = strlen(Annotation->Bases);
    for (AminoIndex = 0; AminoIndex < MatchLength-1; AminoIndex++)
    {
        switch (Annotation->Bases[AminoIndex])
        {
        case 'I':
            Annotation->Bases[AminoIndex] = 'L';
            break;
        case 'Q':
            Annotation->Bases[AminoIndex] = 'K';
            break;
        default:
            break;
        }
        PRM[AminoIndex] = Mass;
        Mass += PeptideMass[Annotation->Bases[AminoIndex]];
        for (ModIndex = 0; ModIndex < MAX_PT_MODS; ModIndex++)
        {
            if (Annotation->AminoIndex[ModIndex] == AminoIndex && Annotation->ModType[ModIndex])
            {
                Mass +=Annotation->ModType[ModIndex]->RealDelta;
            }
        }
    }
    ///////////////////////////////////////////////////////////////
    // Check each tag to see whether it's correct:
    TagCount = min(TagCount, 512);
    FoundFlag = 0;
    TrieTagCount = 0;
    for (TagIndex = 0; TagIndex < min(10, TagCount); TagIndex++)
    {
        TrueTagFlag = 0;
        TestTag = Tags + TagIndex;
        DuplicateFlag = 0;
        Root = AddTagToTrie(Root, TestTag, &DuplicateFlag);
        if (!DuplicateFlag)
        {
            TrieTagCount++;
            for (AminoIndex = 0; AminoIndex < MatchLength-2; AminoIndex++)
            {
                if (abs(TestTag->PrefixMass - PRM[AminoIndex]) < GlobalOptions->ParentMassEpsilon)
                {
                    if (TestTag->Tag[0] == Annotation->Bases[AminoIndex] &&
                        TestTag->Tag[1] == Annotation->Bases[AminoIndex + 1] &&
                        TestTag->Tag[2] == Annotation->Bases[AminoIndex + 2])
                    {
                        TrueTagFlag = 1;
                    }
                }
            }
        }
        if (TrueTagFlag)
        {
            fprintf(TagTrainingFile, "+1 ");
        }
        else
        {
            fprintf(TagTrainingFile, "-1 ");
        }
        FeatureIndex = 1;
        fprintf(TagTrainingFile, "%d:%.3f ", FeatureIndex++, TestTag->Score);
        fprintf(TagTrainingFile, "\n");
    }
    if (!FoundFlag)
    {
        // we missed, too bad.  don't need to poke the histogram.
        //TrueTagRankHistogram[511]++;
    }
    FreeTrieNode(Root);
    Root = NULL;

    SpectrumCount++;

}

void TrainTagging(char* OracleFile, char* OracleDir)
{
    InitBayesianModels(); // to use new PRM scoring
    InitStats();
    TrainOnOracleFile(OracleFile, OracleDir, TrainTaggingCallback);
}

// Using flanking amino acid info, score the remaining nodes in Model.
float SetTaggingFlankScore(PRMBayesianModel* Model, TagGraphNode** TagNodes, TagGraphEdge** TagEdges, int Depth, int RightEndpointFlag)
{
    PRMBayesianNode* Node;
    PRMBayesianNodeHolder* Holder;
    char PrefixAA = UNKNOWN_AMINO;
    char SuffixAA = UNKNOWN_AMINO;
    float Score = 0;
    int TableIndex;
    int ParentIndex;
    //
    if (Depth)
    {
        PrefixAA = TagEdges[Depth - 1]->Jump->Amino;
    }
    if (!RightEndpointFlag)
    {
        SuffixAA = TagEdges[Depth]->Jump->Amino;
    }
    for (Holder = Model->FirstFlank; Holder; Holder = Holder->Next)
    {
        Node = Holder->Node;
        switch (Node->Type)
        {
        case evFlank:
            Node->Value = IonScoringGetFlank(Node, PrefixAA, SuffixAA);
            break;
        case evPrefixAA:
            if ((PrefixAA - 'A') == Node->Flag)
            {
                Node->Value = 1;
            }
            else
            {
                Node->Value = 0;
            }
            break;
        case evSuffixAA:
            if ((SuffixAA - 'A') == Node->Flag)
            {
                Node->Value = 1;
            }
            else
            {
                Node->Value = 0;
            }
            break;
        default:
            // We already knew this node's value (based on intensity).  Now we know its parents' values (based 
            // in part on flanking amino acids).  ASSUME that all parents are in the FlankList.
            TableIndex = Node->Value;
            for (ParentIndex = 0; ParentIndex < Node->ParentCount; ParentIndex++)
            {
                TableIndex += Node->Parents[ParentIndex]->Value * Node->ParentBlocks[ParentIndex];
            }
            Score = Node->ProbTable[TableIndex];
            // The score from the NOISE MODEL has already been integrated.  So, we're done.
            break;
        }
    }
    return Score;
}


// New tag generation function: Generates tags of a (more-or-less) arbitrary length!
// Incorporates a more sophisticated intensity scoring function that considers
// amino acid effects.  
TrieTag* TagGraphGenerateTags(TagGraph* Graph, MSSpectrum* Spectrum, int* TagCount, 
    int MaximumTagCount, SpectrumTweak* Tweak, float TagEdgeScoreMultiplier,
    struct PRMBayesianModel* Model)
{
    TagGraphNode* TagNodes[12];
    TagGraphEdge* TagEdges[12];
    int NodeIndex;
    int EdgeIndex;
    TagGraphNode* Node;
    TagGraphEdge* Edge;
    TagGraphNode* LeftNode;
    TagGraphNode* RightNode;
    int CurrentDepth;
    int InternalNodes;
    float NodeScore;
    TrieTag* Tag;
    int TagAllocation;
    int BacktrackFlag;
    int AllTagCount = 0;
    int Bin;
    float FlankScore[12];
    float ScoreToBeat = -9999;
    //
    if (!Model)
    {
        if (Tweak->Charge < 3)
        {
            Model = TAGModelCharge2;
        }
        else
        {
            Model = TAGModelCharge3;
        }
    }
    *TagCount = 0;
    TagAllocation = 1024;
    if (!AllTags)
    {
        AllTags = (TrieTag*)calloc(sizeof(TrieTag), TagAllocation);
    }
    NodeIndex = 0;
    EdgeIndex = -1;
    BacktrackFlag = 0;
    CurrentDepth = 0;
    TagNodes[0] = Graph->FirstNode;
    // Main iteration: Depth-first traversal through the DAG, up to a maximum depth of 
    // GlobalOptions->GenerateTagLength, and with each possible root (TagNodes[0]).
    while (1)
    {
        // If we're BACKTRACKING, then move to a sibling or parent:
        if (BacktrackFlag)
        {
            // Move the root of the subtree, if necessary:
            if (CurrentDepth == 0)
            {
                // Move to the next 'first' node:
                TagNodes[0] = TagNodes[0]->Next;
                if (!TagNodes[0])
                {
                    break;
                }
                BacktrackFlag = 0;
                continue;
            }
            // Move to a sibling, if we can:
            TagEdges[CurrentDepth - 1] = TagEdges[CurrentDepth - 1]->Next;
            if (TagEdges[CurrentDepth - 1])
            {
                TagNodes[CurrentDepth] = TagEdges[CurrentDepth - 1]->ToNode;
                BacktrackFlag = 0;
                FlankScore[CurrentDepth - 1] = SetTaggingFlankScore(Model, TagNodes, TagEdges, CurrentDepth - 1, 0);
                continue;
            }
            // No more siblings - move up one level.
            CurrentDepth--;
            continue;
        }

        // Special case for level 1: Skip tag nodes with silly masses like 20Da.
        if (CurrentDepth == 0)
        {
            Node = TagNodes[0];
            if (Node->Mass > GlobalOptions->ParentMassEpsilon && Node->Mass < GLYCINE_MASS - GlobalOptions->Epsilon)
            {
                BacktrackFlag = 1;
                continue;
            }
        }

        // If we're deep enough, report a tag and start backtracking:
        if (CurrentDepth >= GlobalOptions->GenerateTagLength)
        {
            FlankScore[CurrentDepth] = SetTaggingFlankScore(Model, TagNodes, TagEdges, CurrentDepth, 1);
            BacktrackFlag = 1;
            LeftNode = TagNodes[0];
            RightNode = TagNodes[CurrentDepth];
            Tag = AllTags + (*TagCount);
            InternalNodes = 0;
            NodeScore = 0;
            for (NodeIndex = 0; NodeIndex <= CurrentDepth; NodeIndex++)
            {
                Node = TagNodes[NodeIndex];
                if (Node->OriginalPeakIndex > 0)
                {
                    NodeScore += Node->Score;
                    NodeScore += FlankScore[NodeIndex];
                    InternalNodes++;
                }
                Tag->Nodes[NodeIndex] = TagNodes[NodeIndex];
            }
            NodeScore *= (GlobalOptions->GenerateTagLength + 1) / (float)max(1, InternalNodes);
            Tag->Score = NodeScore;
            Tag->ModsUsed = 0;
            memset(Tag->ModType, 0, sizeof(MassDelta*) * MAX_PT_MODS);
            memset(Tag->AminoIndex, -1, sizeof(int) * MAX_PT_MODS);
            if (LeftNode->NodeType == evGraphNodeLeftMod)
            {
                // Sanity check: The first AA must be one where this mod can
                // be attached!
                if (!AllKnownPTMods[LeftNode->PTM->Index].Allowed[TagEdges[0]->Jump->Amino - 'A'])
                {
                    continue;
                }
                Tag->AminoIndex[Tag->ModsUsed] = 0;
                Tag->ModType[Tag->ModsUsed] = LeftNode->PTM;
                Tag->ModsUsed++;
            }
            for (EdgeIndex = 0; EdgeIndex < CurrentDepth; EdgeIndex++)
            {
                Edge = TagEdges[EdgeIndex];
                Tag->Score += TagEdges[EdgeIndex]->Score;
                Tag->Tag[EdgeIndex] = Edge->Jump->Amino;
                if (Edge->Jump->Delta)
                {
                    Tag->AminoIndex[Tag->ModsUsed] = EdgeIndex;
                    Tag->ModType[Tag->ModsUsed] = Edge->Jump->Delta;
                    Tag->ModsUsed++;
                }
            }
            // Set skew info:
            Tag->TotalSkew = 0;
            Tag->TotalAbsSkew = 0;
            for (EdgeIndex = 0; EdgeIndex < CurrentDepth; EdgeIndex++)
            {
                Edge = TagEdges[EdgeIndex];
                Tag->TotalSkew += Edge->Skew;
                Tag->TotalAbsSkew += abs(Edge->Skew);
            }
            ////////////////////////////////////////////////////
            // If the total skew is large, penalize the tag's score:
            Bin = (int)fabs((Tag->TotalSkew / 50.0) + 0.5);
            if (Bin >= TagSkewBinCount)
            {
                Bin = TagSkewBinCount - 1;
            }
            Tag->Score += TagSkewScore[Bin] * TagEdgeScoreMultiplier;
            Bin = (int)fabs((Tag->TotalAbsSkew / 50.0) + 0.5);
            if (Bin >= TagSkewBinCount)
            {
                Bin = TagSkewBinCount - 1;
            }
            Tag->Score += TagTotalAbsSkewScore[Bin] * TagEdgeScoreMultiplier;
            ////////////////////////////////////////////////////
            Tag->Tag[EdgeIndex] = '\0';
            if (Tag->Score < ScoreToBeat)
            {
                // Abort the tag - it's not good enough!
                continue; 
            }
            if (RightNode->NodeType == evGraphNodeRightMod)
            {
                // Sanity check: The first AA must be one where this mod can
                // be attached!
                if (!AllKnownPTMods[RightNode->PTM->Index].Allowed[TagEdges[CurrentDepth - 1]->Jump->Amino - 'A'])
                {
                    continue;
                }
                Tag->AminoIndex[Tag->ModsUsed] = CurrentDepth;
                Tag->ModType[Tag->ModsUsed] = RightNode->PTM;
                Tag->ModsUsed++;
            }
            Tag->PSpectrum = Spectrum;
            Tag->Tweak = Tweak;
            Tag->TagLength = CurrentDepth;
            Tag->ParentMass = Spectrum->ParentMass;
            Tag->Charge = Spectrum->Charge;
            Tag->PrefixMass = TagNodes[0]->Mass;
            Tag->SuffixMass = Spectrum->ParentMass - PARENT_MASS_BOOST - TagNodes[CurrentDepth]->Mass;
            (*TagCount)++;
            AllTagCount++;
            // If we've got as many tags as we can handle, drop all but the best.  (Don't
            // just reallocate; we could end up with a *lot*!)
            if ((*TagCount) + 5 >= TagAllocation)
            {
                qsort(AllTags, *TagCount, sizeof(TrieTag), (QSortCompare)CompareTagScores);
                *TagCount = TagAllocation / 2;
                if (MaximumTagCount >= 0)
                {
                    ScoreToBeat = AllTags[min(TagAllocation - 5, MaximumTagCount)].Score;
                }
                else
                {
                    ScoreToBeat = AllTags[*TagCount].Score;
                }
            }
            continue;
        } // If we're at tag depth

        // We're not at tag depth yet. 
        // Move to our first available child:
        TagEdges[CurrentDepth] = TagNodes[CurrentDepth]->FirstEdge;
        if (!TagEdges[CurrentDepth])
        {
            BacktrackFlag = 1;
            continue;
        }
        else
        {
            CurrentDepth++;
            TagNodes[CurrentDepth] = TagEdges[CurrentDepth - 1]->ToNode;
            FlankScore[CurrentDepth - 1] = SetTaggingFlankScore(Model, TagNodes, TagEdges, CurrentDepth - 1, 0);
        }
    }
    // Sort the tags, by score:
    qsort(AllTags, *TagCount, sizeof(TrieTag), (QSortCompare)CompareTagScores);
    return AllTags;

}
