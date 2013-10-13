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

#ifndef TRIE_H
#define TRIE_H

// Implementation of the Aho-Corasic algorithm for string search using 
// trie automaton.  Also, implementation of the d.p. tag extension algorithm
// in presence of PTMs.

#include "Inspect.h"
#include <stdio.h>
#include "Utils.h"
#include "Spectrum.h"
#include "Mods.h"

// Tags, produced from mass-spec analysis.  A tag consists of sequence of bases
// (e.g. "QVL"), a prefix mass, and a Suffix mass.  Tags are stored in the nodes
// of a trie.  
// In the simple case, leaf nodes in the trie have just one tag.  But
// there may be multiple tags with the same bases but different prefix/Suffix
// masses, these two tags will be children of the same trie node.  (We use 
// the TrieTagHanger struct to hold lists of tags.  The second class adds some
// overhead; one advantage is that tags can be in more than one list)
typedef struct TrieTag
{
#ifdef DEBUG_TAG_GENERATION
    char TagScoreDetails[2048];
#endif
    int PrefixMass;
    int SuffixMass;
    int Charge;
    int ParentMass;
    int TagLength;
    int CandidatesScored;
    // How far this tag's mass is from the ACTUAL mass of the peptide
    float Error; 
    // Rank of this tag (0 being best)
    int Rank;
    // Score of this tag (higher is better):
    float Score;
    // PTMods used up in generating this tag.  These count against our allowable total.
    // AminoIndex is set to -1 for unused entries.
    int AminoIndex[MAX_PT_MODS]; 
    MassDelta* ModType[MAX_PT_MODS]; 
    int ModsUsed;
    int MandatoryModUsed; // if GlobalOptions->MandatoryModIndex is set.
    // The tag itself:
    char Tag[MAX_TAG_LENGTH + 1];
    struct MSSpectrum* PSpectrum;
    SpectrumTweak* Tweak;
    int PrefixExtends;
    int SuffixExtends;
    int DBTagMatches;
    // Some members for training edge skew measures:
    int TotalSkew;
    int TotalAbsSkew;
    struct TagGraphNode* Nodes[4];
} TrieTag;

// A trie (from 'reTRIEval') is a tree where each node corresponds to a word.  The root 
// corresponds to an empty string, and a node's children correspond to that node's word 
// extended by one letter.  The trie data structure allows fast searches for any of the
// words in the trie.  In this case, the 'words' are short peptides, and the database
// is swiss-prot or some other collection of peptides.
//
// During the search, the ANCHOR is the start of our current match (if any).  
typedef struct TrieNode
{
    struct TrieNode* Children[TRIE_CHILD_COUNT]; 
    // Depth == length of our string.
    // Root has depth 0, its children have depth 1...
    int Depth; 

    // The failure node is an optimization which makes trie searching fast.
    // Suppose we just finished matching the tag PANTS.  Naively, we would move on
    // to tags starting with A.  But if we have a node ANT, maybe we can jump there 
    // directly.  If we have no nodes starting with A, we can jump to the N.
    // The FailureNode 'precomputes our jumps' - we move the anchor FailureLetterSkip 
    // letters forward from the old anchor, and switch to the given failure node.
    // If FailureNode is equal to the trie root, then we CLEAR THE ANCHOR.  
    int FailureLetterSkip;
    struct TrieNode* FailureNode; 

    // Our (most recently added) letter:
    char Letter;
    // Our list of tags:
    struct TrieTagHanger* FirstTag;
    struct TrieTagHanger* LastTag;
} TrieNode;

typedef struct PeptideSpliceNode
{
    int ChromosomeNumber;
    int DonorPos;
    int AcceptorPos;
    struct PeptideSpliceNode* Next;
} PeptideSpliceNode;

// A Peptide struct represents the peptide we use to annotate a tandem mass spectrum - possibly 
// with PTMs, and prefix and suffix residues.
typedef struct Peptide
{
    int ParentMassError;
    char Bases[256];

    char PrefixAmino; // The base BEFORE the peptide starts.  (Useful for checking trypsin fragmentation)

    char SuffixAmino; // The base AFTER the peptide starts.

    // FilePos is the byte-offset in the database where the peptide starts.  If the peptide is found
    // multiple times, FilePos the position within the file of the FIRST occurrence of the peptide.
    int FilePos; 

    // RecordNumber is the protein record # where the peptide is found
    int RecordNumber; 

    int InitialScore;

    float MatchQualityScore;

    // For the best match, DeltaCN is the difference in score between this match and its runner-up.  For
    // other matches, DeltaCN is the difference in score between them and the best match (i.e. DeltaCN is 
    // negative for them).  We compute DeltaCN because a larger DeltaCN value generally indicates 
    // a better match.  
    float DeltaCN;

    // DeltaCNOther is the difference in score between this peptide and the best runner-up that's NOT 
    // the same peptide.  "Same" means "file-pos at most 2 steps away, or sequence has at most two diffs".
    float DeltaCNOther;
    struct Peptide* FamilyLeader;
    float FScore;
    float PValue;

    // We may "own" our own mass delta, in which case we must free it when we dealloc:
    MassDelta* PetDelta; 

    // Track the nth post-translational modification by setting AminoIndex[n] to the index of the
    // modified amino acid, and ModType[n] to the modification type.  Set AminoIndex to -1 for
    // the extra records.
    int AminoIndex[MAX_PT_MODS]; 
    MassDelta* ModType[MAX_PT_MODS];
    struct PeptideMatch* First;
    struct PeptideMatch* Last;
    struct Peptide* Next;
    struct Peptide* Prev;
    int PrefixMass; // Used only if this is a tag
    int SuffixMass; // Used only if this is a tag 

    // DB is a pointer to the database which this match comes from.
    DatabaseFile* DB;
    PeptideSpliceNode* SpliceHead; // if splice-tolerant
    int GenomicLocationStart; // if splice-tolerant
    int GenomicLocationEnd; // if splice-tolerant
    int ChromosomeNumber; // if splice-tolerant
    int ChromosomeForwardFlag; // if splice-tolerant
    char* SplicedBases; // if splice-tolerant
    int ParentMass;
    SpectrumTweak* Tweak;
    float ScoreFeatures[16];
    int SpecialFragmentation;
    int SpecialModPosition;
} Peptide;

typedef struct PeptideMatch
{
    int FilePos;
    int RecordNumber;
    TrieTag* Tag; 
    struct PeptideMatch* Next;
} PeptideMatch;

typedef float (*ScoringFunction)(MSSpectrum* Spectrum, Peptide* Match, int VerboseFlag);

typedef struct SearchInfo
{
    DatabaseFile* DB;
    int RecordNumber;
    //ScoringFunction Scorer;
    MSSpectrum* Spectrum;
    TrieNode* Root;
    int VerboseFlag;
} SearchInfo;

//container for information about the blind tag match.
//These matches extend on only one side, which I called Anchored
typedef struct BlindTagMatch
{
    TrieTag* Tag;
    struct BlindTagMatch* Next;
    struct BlindTagMatch* Prev;
    //denotes the direction of the matched (modless) extension
    int ExtendLR; // -1 for Left, 1 for Right
    int ExtendDBLoc; // the location in the database where the extension matches
    int TagDBLoc; // the location in the DB where the first letter of the tag matches.
    int ExtendLength; //length of the anchored extension XXXTAG---- means an extend len of 3
} BlindTagMatch;

// PTModCount lists how many post-translational mods exist for each peptide:
extern int PTModCount[TRIE_CHILD_COUNT];

// PTMods lists the mass of each post-translational mod for each peptide:
extern float PTMods[TRIE_CHILD_COUNT][MAX_PT_MODTYPE];

// Table of prefix and suffix peptide masses
extern int PeptideMass[256];
extern int StandardPeptideMass[256];

#define IS_ROOT(node) ((node)->Depth == 0)


// For constructing lists of TrieTags.  (A single TrieTag can be part of more than
// one list, by using more than one TrieTagHanger)
typedef struct TrieTagHanger
{
    struct TrieTagHanger* Prev;
    struct TrieTagHanger* Next;
    TrieTag* Tag;
} TrieTagHanger;

// Constructor: TrieNode
TrieNode* NewTrieNode();

// Destructor: TrieNode
void FreeTrieNode(TrieNode* This);

// Constructor: TrieTag
TrieTag* NewTrieTag();

// Destructor: TrieTag
void FreeTrieTag(TrieTag* This);

// Add a new tag to the trie.  New trie nodes will be added, if necessary, in order
// to hold the tag.  (For instance, adding "CAT" to a root node with no children would
// add three nodes: C, CA, and CAT).
TrieNode* AddTagToTrie(TrieNode* Root, TrieTag* Tag, int* DuplicateFlag);

// Constructor: TrieTagHanger
TrieTagHanger* NewTrieTagHanger();

// Destructor: TrieTagHanger
void FreeTrieTagHanger(TrieTagHanger* This);

// Debug: Print a trie to stdout
void DebugPrintTrie(TrieNode* Root);

// Print all matches 
void PrintMatches(MSSpectrum* Spectrum, char* IndexFileName);

// Load the masses for amino acids: n-terminal (left) and c-terminal (right) masses
int LoadPeptideMasses(char* FileName);

// Initialize GlobalOptions, a global variable storing configurable options.
void InitOptions();

void InitStats();

// Important main method: Use a trie to search a data-file.
int ScanFileWithTrie(SearchInfo* Info);

int GetMaxTagRank(TrieNode* Root);
//int ComparePeptideScores(const Peptide* A, const Peptide* B);
void PrintMatch(Peptide* Match, FILE* IndexFile);
void DebugPrintTrieTags(TrieNode* Node);
void InitializeTrieFailureNodes(TrieNode* Root, TrieNode* Node, char* Tag);
void FreePeptideNode(Peptide* Pep);
void GetProteinID(int RecordNumber, DatabaseFile* DB, char* Name);
void FlagMandatoryModUsage(TrieNode* Node);
void WriteMatchToString(Peptide* Match, char* Buffer, int IncludeMods);
Peptide* StoreSpectralMatch(MSSpectrum* Spectrum, Peptide* Match, int PeptideLength, int MQScoreFlag);
Peptide* NewPeptideNode();
int CheckForPTAttachmentPoints(int DecorationMassIndex, char* Buffer, int Start, int End, int BufferDir);
Peptide* AddNewMatch(SearchInfo* Info, int FilePos, TrieTag* Tag, char* Peptide, 
                 int PeptideLength, int TagPosition, int PrefixDecoration, int SuffixDecoration, 
                 int GenomicStart, int GenomicEnd);
Peptide* GetPeptideFromAnnotation(char* Annotation);
int GetPeptideParentMass(Peptide* Match);
void WriteTagToString(TrieTag* Tag, char* Buffer, int IncludeMods);

#endif //TRIE_H
