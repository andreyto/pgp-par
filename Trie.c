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

#include "CMemLeak.h"
#include "Trie.h"
#include "Utils.h"
#include <memory.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h> // for fabs
#include <ctype.h>
#include "Spectrum.h"
#include "Mods.h"
#include "Score.h"
#include "Tagger.h"
#include "BN.h"
#include "Scorpion.h"
#include "Errors.h"
#include "SVM.h"
#include "LDA.h"

// If two tags have the same peptides, and are within this amount on prefix/Suffix masses, then
// consider them identical and only use the top scorer. 
// 1.5 daltons:
#define IDENTICAL_TAG_EPSILON 1500

// Number of chars allowed in a post-translational modification name.  ('phosphorylation' is a typical name)
#define MAX_PTMOD_NAME 256

#define TRIE_INDEX_RECORD_SIZE (LONG_LONG_SIZE + sizeof(int) + 80*sizeof(char))
#define SPLICEDB_INDEX_RECORD_SIZE (LONG_LONG_SIZE + sizeof(int) + 80*sizeof(char))

// Global variable storing configurable options:
Options* GlobalOptions;
InspectStats* GlobalStats;

////////////////////////////////////////////////////////////////////////////////////////
// Forward declarations:
void FreeTrieTagHangerList(TrieTagHanger* Head, int FreeTags);
void FlagMandatoryModUsage(TrieNode* Node);
int ExtendTagMatchBlind(SearchInfo* Info, TrieNode* Node, char* Buffer, int BufferPos, int BufferEnd, int FilePos);
int ProcessGeneHitsBlindTag();
int InsertBlindTagMatch(BlindTagMatch* Match);
void FreeBlindTagMatch(BlindTagMatch* This);
void FreeAllBlindTagMatches(BlindTagMatch* This);
int IsIdenticalBlindTagMatches(BlindTagMatch* NodeA, BlindTagMatch* NodeB);

Peptide* FindMatchedPeptide(char* Bases);
void AddPTMassesToTagTable(int TagTableSize, char* CurrentTag, float Mass, char FirstAllowedPep, int CharsLeft, int ModsLeft,
                           int Peptide, int MinMod);


// Indexed by characters.  (Use upper-case amino acid codes)
int StandardPeptideMass[256];
// PeptidMass may be different from StandardPeptideMass if a fixed modification
// (e.g. +57 to all cysteine residues) has been applied.  
int PeptideMass[256];

// A decoration is a collection of post-translational modification.  This includes the
// 'empty decoration', with no modifications, and mass 0.  Each decoration has an index;
// they are ordered from smallest mass to largest.

// Size of the decoration array:
int DecorationMassCount;

// Mass of each decoration:
float* DecorationMasses;

// Largest mass over all our decorations
float DecorationMaxMass;

// DecorationMassMods[DecorationIndex][n] is the index of the nth mod used in a particular decoration.
// For decorations that use fewer than the maximum allowed number of mods, we store an index of -1.
int** DecorationMassMods;

// DecorationModCount[DecorationIndex] is the number of mods used in a decoration
int* DecorationModCount;

// PTModCount lists how many post-translational mods exist for each peptide.  (Faster than iterating
// over a full 2D table of flags).  Indexed by peptide-char (entry #0 is alanine)
int PTModCount[TRIE_CHILD_COUNT];

// How many modifications are there, in all?
int TotalPTMods;

// SubDecorations tells how to get to a sub-decoration (a decoration containing fewer post-translational
// modificaitons) from a parent decoration.  SubDecorations[DecorIndex][Modification] is the index 
// of the decoration Decor with one such modification removed.  SubDecorations entries are -1 if the specified
// mod isn't part of the specified decoration.
int** SubDecorations;

// PTMods lists the mass of each post-translational mod for each peptide.  (Redundant storage,
// for fast lookups)
float PTModMass[TRIE_CHILD_COUNT][MAX_PT_MODTYPE];

// PTMods lists the index of each post-translational mod for each peptide.  (So, PTMods[0][0] is the 
// first modification available to alanine
int PTModIndex[TRIE_CHILD_COUNT][MAX_PT_MODTYPE];

// Names of all known PTMods.  
char PTModName[MAX_PT_MODTYPE][MAX_PTMOD_NAME];
float ModMasses[MAX_PT_MODTYPE];

//BlindTagMatch Pointers for the list of matches to a single Gene
BlindTagMatch* FirstBlindTag = NULL;
BlindTagMatch* LastBlindTag = NULL;

void InitStats()
{
    if (GlobalStats)
    {
        memset(GlobalStats, 0, sizeof(InspectStats));
    }
    else
    {
        GlobalStats = (InspectStats*)calloc(sizeof(InspectStats), 1);
    }
}

// Set global options to reasonable default values:
void InitOptions()
{
    GlobalOptions = (Options*)calloc(sizeof(Options), 1);
    GlobalOptions->MaxPTMods = 0; 
    GlobalOptions->Epsilon = DEFAULT_EPSILON;
    GlobalOptions->FlankingMassEpsilon = DEFAULT_FLANKING_MASS_EPSILON; 
    GlobalOptions->OutputFile = stdout;
    sprintf(GlobalOptions->ErrorFileName, "inspect.err");
    GlobalOptions->ErrorCount = 0;
    GlobalOptions->WarningCount = 0;
    GlobalOptions->ReportAllMatches = 1;
    GlobalOptions->ParentMassEpsilon = DEFAULT_PARENT_MASS_EPSILON;
    GlobalOptions->ParentMassPPM = DEFAULT_PARENT_MASS_PPM;
    GlobalOptions->ReportMatchCount = 10; // Don't report more than 10 to the page!
    GlobalOptions->StoreMatchCount = 100;  // 
    GlobalOptions->MandatoryModIndex = -1; // By default, there is no mandatory modification
    GlobalOptions->GenerateTagCount = 100; 
    GlobalOptions->GenerateTagLength = DEFAULT_TAG_LENGTH;
    GlobalOptions->DynamicRangeMin = 105 * DALTON; 
    GlobalOptions->DynamicRangeMax = 2000 * DALTON; 
    GlobalOptions->TrieBlockSize = 100;
    GlobalOptions->TagPTMMode = 2;
    strcpy(GlobalOptions->AminoFileName, FILENAME_AMINO_ACID_MASSES);
    sprintf(GlobalOptions->InputFileName, "Input.txt");
    GlobalOptions->MinPTMDelta = -200;
    // Default maxptmdelta is 250, this allows us to find GlcNac (203) and biotin (226)
    GlobalOptions->MaxPTMDelta = 250;
    GlobalOptions->DeltaBinCount = (GlobalOptions->MaxPTMDelta - GlobalOptions->MinPTMDelta) * 10 + 1;
    GlobalOptions->DeltasPerAA = max(512, GlobalOptions->DeltaBinCount * 2);
    GlobalOptions->NewScoring = 0;
    
}

//constructor for a new BlindTagMatch
BlindTagMatch* NewBlindTagMatch()
{
    BlindTagMatch* This;
    This = (BlindTagMatch*)calloc(sizeof(BlindTagMatch),1);
    This->Next = NULL; //set a few pointers up for tidyness
    This->Prev = NULL;
    return This;
}
//destructor for BlindTagMatch.  This frees ALL connected nodes
//following the next pointer
void FreeAllBlindTagMatches(BlindTagMatch* This)
{
    This->Tag = NULL;  //free pointer to tag, but KEEP TAG
    This->Prev = NULL;
    if (This->Next){
        FreeAllBlindTagMatches(This->Next);
    }
    This->Next = NULL;
    SafeFree(This);
}
//This destructor assumes that the links have been
//previously nullified, and the linked list fixed
//and ready for this node to be wipedout
void FreeBlindTagMatch(BlindTagMatch* This)
{
    This->Tag = NULL;//free pointer to tag, but KEEP TAG, it's part of the trie
    This->Next = NULL;
    This->Prev = NULL;
    SafeFree(This);
}
// Constructor for a new TrieNode
TrieNode* NewTrieNode()
{
    TrieNode* This;
    This = (TrieNode*)calloc(sizeof(TrieNode), 1);
    return This;
}


// Free a trie node.  Also frees its tag-nodes (if any), and recursively frees its children.
void FreeTrieNode(TrieNode* This)
{
    int Letter;
    if (!This)
    {
        return;
    }
    // Free our tag nodes:
    FreeTrieTagHangerList(This->FirstTag, 1);

    // Free our children too!  Free them only AFTER
    // we iterate over them, since Node->Next must be
    // valid at the end of each loop-cycle.
    for (Letter = 0; Letter < TRIE_CHILD_COUNT; Letter++)
    {
        // Nodes I and K always point to the same child as L and Q, respectively.
        // So...don't free them twice!
        if (Letter == ('I'-'A') || Letter == ('K'-'A'))
        {
            continue;
        }
        if (This->Children[Letter])
        {
            FreeTrieNode(This->Children[Letter]);
        }
    }

    // Ok, now free ourselves:
    SafeFree(This);
}

// Constructor for a TrieTag
TrieTag* NewTrieTag()
{
    TrieTag* This;
    int ModIndex;
    //
    This = (TrieTag*)calloc(sizeof(TrieTag), 1);
    for (ModIndex = 0; ModIndex < MAX_PT_MODS; ModIndex++)
    {
        This->AminoIndex[ModIndex] = -1;
    }
    This->DBTagMatches = 0;
    This->PrefixExtends = 0;
    This->SuffixExtends = 0;
    return This;
}

// Destructor for a TrieTag
void FreeTrieTag(TrieTag* This)
{
    SafeFree(This);
}

// Trie construction helper function.
// We've got a node on the Trie which completely matches a tag.  So,
// add this tag to the list of TrieTagNodes on this TrieNode.
TrieTagHanger* TrieNodeAddTagNode(TrieNode* Node, TrieTag* Tag, int* DuplicateFlag)
{
    TrieTagHanger* Hanger;
    TrieTag* LocalTag;
    // Look at our current list of tags for the node.
    // DON'T add this tag if we already have the same pre-and-post masses.
    // (We normally add tags in order from best to worst, so in the event that
    // two tags are quit similar, we keep the one with the higher score)
    //printf("Adding tag '%s' %.2f...\n", Tag->Tag, Tag->PrefixMass); 
    for (Hanger = Node->FirstTag; Hanger; Hanger = Hanger->Next)
    {
        if (Hanger->Tag->PSpectrum == Tag->PSpectrum)
        {
            if ((abs(Hanger->Tag->SuffixMass - Tag->SuffixMass) < IDENTICAL_TAG_EPSILON) &&
                (abs(Hanger->Tag->PrefixMass - Tag->PrefixMass) < IDENTICAL_TAG_EPSILON))
            {
                // Prefer the new prefix/suffix, if this new tag scores higher:
                if (Hanger->Tag->Score < Tag->Score)
                {
                    Hanger->Tag->PrefixMass = Tag->PrefixMass;
                    Hanger->Tag->SuffixMass = Tag->SuffixMass;
                }
                *DuplicateFlag = 1;
                return Hanger;
            }
        }
    }

    Hanger = NewTrieTagHanger();
    LocalTag = NewTrieTag();
    memcpy(LocalTag, Tag, sizeof(TrieTag));
    Hanger->Tag = LocalTag;
    //Hanger->Tag = Tag;
    if (Node->LastTag)
    {
        Node->LastTag->Next = Hanger;
        Hanger->Prev = Node->LastTag;
    }
    else
    {
        Node->FirstTag = Hanger;
    }
    Node->LastTag = Hanger;
    *DuplicateFlag = 0;
    return Hanger;
}

// We've got a tag - add it to the trie.
TrieNode* AddTagToTrie(TrieNode* Root, TrieTag* Tag, int* DuplicateFlag)
{
    TrieNode* Node;
    TrieNode* NextNode;
    int Index;
    int TagLength;
    char TagChar;
    // 
    // First, travel down the trie, matching the specified tag as far as possible (perhaps completely):
    Index = 0;
    TagLength = Tag->TagLength; 
    Node = Root;
    while (1)
    {
        TagChar = Tag->Tag[Index];
        if (TagChar == 'I')
        {
            TagChar = 'L';
        }
        else if (TagChar == 'K')
        {
            TagChar = 'Q';
        }
        // Look up our child for this letter:
        NextNode = Node->Children[TagChar - 'A'];
        if (!NextNode)
        {
            // Ok, we matched as far as possbile - next we'll add children to match the remainder of the tag.
            break;
        }
        Node = NextNode;
        Index++;
        // Did we match the tag completely?
        if (Index == TagLength)
        {
            // Aha - this tag is in the trie!  Add the tag to the list:
            TrieNodeAddTagNode(Node, Tag, DuplicateFlag);
            return Node;  
        }
    }
    // Ok, we didn't match the entire tag...so, start adding child nodes now!
    while (Index < TagLength)
    {
        NextNode = NewTrieNode();
        NextNode->Letter = Tag->Tag[Index];
        NextNode->Depth = Index + 1;
        Node->Children[Tag->Tag[Index] - 'A'] = NextNode;
        // Extra child I and L to same place, K and Q to same place.
        switch (Tag->Tag[Index])
        {
        // Special case for aminos with same mass (I and L equal, K and Q are off by <.1):
        // Child pointers for I and L, and for K and Q both point to the same place.
        case 'I':
            Node->Children['L'-'A'] = NextNode;
            break;
        case 'L':
            Node->Children['I'-'A'] = NextNode;
            break;
        case 'K':
            Node->Children['Q'-'A'] = NextNode;
            break;
        case 'Q':
            Node->Children['K'-'A'] = NextNode;
            break;
        default:
            break;
        }
        Node = NextNode;
        Index++;
    }
    TrieNodeAddTagNode(Node, Tag, DuplicateFlag);

    return Node;
}

// Constructor for a TagHanger
TrieTagHanger* NewTrieTagHanger()
{
    TrieTagHanger* This;
    //
    This = (TrieTagHanger*)calloc(sizeof(TrieTagHanger), 1);
    return This;
}

// Destructor for a TagHanger
void FreeTrieTagHanger(TrieTagHanger* This)
{
    SafeFree(This);
}

// Destructor for a TagHanger list
void FreeTrieTagHangerList(TrieTagHanger* Head, int FreeTags)
{
    TrieTagHanger* Prev = NULL;
    //
    for (; Head; Head = Head->Next)
    {
        if (Prev)
        {
            if (FreeTags)
            {
                FreeTrieTag(Prev->Tag);
            }
            FreeTrieTagHanger(Prev);
        }
        Prev = Head;
    }
    if (Prev)
    {
        if (FreeTags)
        {
            FreeTrieTag(Prev->Tag);
        }
        FreeTrieTagHanger(Prev);
    }
}

// Prints a Trie node, using indentation to denote depth.
// The entry point is DebugPrintTrie, which calls this function.
int DebugPrintTrieHelper(TrieNode* Root, char* TagSoFar)
{
    char Buffer[1024];
    char TagBuffer[1024];
    int Index;
    int BufferPos;
    int TagLength;
    TrieTagHanger* Node;
    int Letter;
    int TagCount = 0;
    //
    TagLength = strlen(TagSoFar);
    BufferPos = 0;
    for (Index = 0; Index < TagLength; Index++)
    {
        Buffer[BufferPos++]=' ';
        Buffer[BufferPos++]=' ';
        Buffer[BufferPos++]=' ';
    }
    Buffer[BufferPos] = '\0';
    //
    strcpy(TagBuffer, TagSoFar);
    if (Root->Letter)
    {
        BufferPos = strlen(TagBuffer);
        TagBuffer[BufferPos++] = Root->Letter;
        TagBuffer[BufferPos++] = '\0';
    }
    printf("%s%s\n", Buffer, TagBuffer);

    // Print attached tags:
    for (Node = Root->FirstTag; Node; Node = Node->Next)
    {
        printf("%s%s: prefix %d,suffix %d mods %d\n", Buffer, Node->Tag->Tag, Node->Tag->PrefixMass,Node->Tag->SuffixMass, Node->Tag->ModsUsed);
        TagCount++;
    }
    if (Root->FailureNode)
    {
        printf("%s Failure: Skip %d, depth %d\n", Buffer, Root->FailureLetterSkip, Root->FailureNode->Depth);
    }
    else
    {
        printf("%s (no failure node set)\n", Buffer);
    }

    // Print children:
    for (Letter = 0; Letter < TRIE_CHILD_COUNT; Letter++)
    {
        if (Root->Children[Letter])
        {
            
            TagCount += DebugPrintTrieHelper(Root->Children[Letter], TagBuffer);
        }
    }
    return TagCount;
}

// Print out a trie and all its nodes to stdout.
void DebugPrintTrie(TrieNode* Root)
{
    int TagCount;

    printf("-->Trie:\n");
    TagCount = DebugPrintTrieHelper(Root, "");
    printf("Total tags: %d\n", TagCount);
    printf("(end of trie nodes)\n");
}

// Set up all the failure nodes for our trie
// The failure node for a trie node is the node you jump to when you're currently
// matching that trie node, but then break a match because none of your children
// match the *next* character.  If another trie node matches a substring of this
// node (not a prefix, but any other substring), we must try matching that node as 
// well.
// Example: Suppose we have nodes ABCDE and BCD, and we're scanning text ABCDF.
// when we reach F, we jump to ABCDE's failure node BCD, and move the anchor to B.
void InitializeTrieFailureNodes(TrieNode* Root, TrieNode* Node, char* Tag)
{
    TrieNode* FailureNode;
    int Letter;
    int TagLength;
    int StartIndex = 0;
    int EndIndex = 0;
    //
    TagLength = strlen(Tag);
    if (!Root)
    {
        return;
    }
    if (Node == Root)
    {
        // Failure on the root means the letter can't start a tag; no speedup, just step forward:
        Root->FailureNode = Root;
        Root->FailureLetterSkip = 1;
    }
    else
    {
        // There's a real tag.  Navigate to the SHORTEST node-with-tags which matches a suffix of our tag.
        // Try knocking off one letter, then two, and so on:
        for (StartIndex=1; StartIndex<TagLength; StartIndex++)
        {
            FailureNode = Root;
            for (EndIndex = StartIndex; EndIndex<TagLength; EndIndex++)
            {
                Letter = Tag[EndIndex];
                if (Letter == 'I')
                {
                    Letter = 'L';
                }
                if (Letter == 'Q')
                {
                    Letter = 'K';
                }
                if (!FailureNode->Children[Letter - 'A'])
                {
                    // We can't go deeper in the trie...and we saw NO TAGS!  So, we needn't go here.
                    // Suppose you have tags GGGROOVY and GROOVY.  After matching GGGROOVY, you can
                    // jump to the 3rd G (you needn't handle the 2nd, even though it matches
                    // partway down the trie).  Of course, if you had GGGROOVY, GROOVY and GGROO, you would only
                    // jump to the 2nd G - because we'd have found a tag in here..
                    FailureNode = Root;
                    break; 
                }
                FailureNode = FailureNode->Children[Letter - 'A']; // move down the tree.
                // If there are tags, STOP NOW.  (Don't jump from PANTS->ANTS if ANT has a tag)
                if (FailureNode->FirstTag)
                {
                    break;
                }
            }
            // If we're not pointing at root, then we're pointing at a good failure node:
            if (FailureNode != Root)
            {
                Node->FailureNode = FailureNode;
                Node->FailureLetterSkip = StartIndex;
                break;
            }
        }
        if (!Node->FailureNode)
        {
            // Hmm...no good failure nodes found?  That means we can jump forward over our full tag!
            Node->FailureNode = Root;
            Node->FailureLetterSkip = strlen(Tag);
        }
    }
    // Now, handle all our children:
    for (Letter = 0; Letter < TRIE_CHILD_COUNT; Letter++)
    {
        if (Node->Children[Letter])
        {
            Tag[TagLength] = 'A'+Letter;
            Tag[TagLength+1] = '\0';
            InitializeTrieFailureNodes(Root, Node->Children[Letter], Tag);
        }
    }
}

// Constructor for a Peptide
Peptide* NewPeptideNode()
{
    Peptide* This;
    This = (Peptide*)calloc(sizeof(Peptide), 1);
    if (!This)
    {
        printf("** Fatal error: Unable to allocate a new peptide!\n");
        return NULL;
    }
    memset(This->AminoIndex, -1, sizeof(int) * MAX_PT_MODS);
    return This;
}

// Destructor for a Peptide
void FreePeptideNode(Peptide* Pep)
{
    PeptideMatch* Node;
    PeptideMatch* Prev = NULL;
    PeptideSpliceNode* PSNode;
    PeptideSpliceNode* PSPrev;
    if (!Pep)
    {
        return;
    }
    // Free the list of PeptideSpliceNodes, starting with SpliceHead:
    PSPrev = NULL;
    for (PSNode = Pep->SpliceHead; PSNode; PSNode = PSNode->Next)
    {
        SafeFree(PSPrev);
        PSPrev = PSNode;
    }
    SafeFree(PSPrev);
    // Free the list of PeptideMatch instances, starting with First:
    for (Node = Pep->First; Node; Node = Node->Next)
    {
        SafeFree(Prev);
        Prev = Node;
    }
    SafeFree(Prev);
    SafeFree(Pep->SplicedBases);
    SafeFree(Pep->PetDelta);
    SafeFree(Pep);
}

MassDelta* GetPeptideModFromAnnotation(Peptide* Match, char* ModBuffer, int ModCount, int AminoIndex)
{
    int MaxModsFromParsedPeptide = 10;
    MassDelta* Delta;
    //
    if (!Match->PetDelta)
    {
        Match->PetDelta = (MassDelta*)calloc(sizeof(MassDelta), MaxModsFromParsedPeptide);
    }
    if (ModCount >= MaxModsFromParsedPeptide)
    {
        return NULL;
    }
    Delta = Match->PetDelta + ModCount;
    Delta->Flags = DELTA_FLAG_VALID;
    if(!CompareStrings(ModBuffer,"phos"))
    { //necessary switch, cannot do atoi("phos") and expect real numbers
        Delta->RealDelta = 80 * DALTON;
        Delta->Flags |= DELTA_FLAG_PHOSPHORYLATION;
        Match->SpecialFragmentation = FRAGMENTATION_PHOSPHO; // special flags we need
        Match->SpecialModPosition = Match->AminoIndex[ModCount];
    }
    else
    {
        Delta->RealDelta = atoi(ModBuffer) * DALTON;
    }
    Delta->Delta = Delta->RealDelta / 100; // tenth-of-a-dalton
    Match->AminoIndex[ModCount] = AminoIndex - 1;
    Match->ModType[ModCount] = Delta;
    return Delta;
}

// Produce a peptide from an annotation string.  The annotation string
// consists of amino acids, plus - possibly - some modification masses.
// Valid examples: 
// GPLLVQDVVFTDEMAHFDR
// VLVLDTDY+16KK
// SVTDC-2TSNFCLFQSNSK
Peptide* GetPeptideFromAnnotation(char* Annotation)
{
    char ModBuffer[32];
    int AminoIndex = 0;
    int ModCount = 0;
    int ModBufferPos;
    Peptide* Match;
    MassDelta* Delta;
    int PRM = 0;
    char* BaseAnnotation;
    int MaxModsFromParsedPeptide = 10;
    //
    if (!Annotation)
    {
        return NULL;
    }
    ModBufferPos = 0;
    Match = NewPeptideNode();
    BaseAnnotation = Annotation;
    if (BaseAnnotation[1] == '.')
    {
        Match->PrefixAmino = BaseAnnotation[0];
        Annotation += 2;
    }
    while (*Annotation)
    {
        if ((*Annotation >= 'A' && *Annotation <= 'Z') || *Annotation == '.')
        {
            // It's an amino acid, or period.  
            // Finish any pending mod:
            if (ModBufferPos)
            {
                ModBuffer[ModBufferPos] = '\0';
                Delta = GetPeptideModFromAnnotation(Match, ModBuffer, ModCount, AminoIndex);
                if (!Delta)
                {
                    printf("*** Warning: Invalid modifications in '%s', not parsing\n", Annotation);
                    FreePeptideNode(Match);
                    return NULL;
                }
                PRM += Delta->RealDelta;
                ModBufferPos = 0;
                ModCount += 1;
                // Bail out if we have too many PTMs to cope with:
                if (ModCount == MAX_PT_MODS)
                {
                    return NULL;
                }
            }
            // It's a dot - set the prefix and break:
            if (*Annotation == '.')
            {
                Match->SuffixAmino = *(Annotation + 1);
                break;
            }
            // It's an amino acid - add the AA mass:
            Match->Bases[AminoIndex++] = *Annotation;
            PRM += PeptideMass[*Annotation];
        }
        else
        {
            ModBuffer[ModBufferPos++] = *Annotation;
        }
        Annotation++;
    }
    Match->Bases[AminoIndex] = '\0';
    // Finish any pending mod:
    if (ModBufferPos)
    {
        ModBuffer[ModBufferPos] = '\0';
        Delta = GetPeptideModFromAnnotation(Match, ModBuffer, ModCount, AminoIndex);
        if (!Delta)
        {
            printf("*** Warning: Invalid modifications in '%s', not parsing\n", Annotation);
            FreePeptideNode(Match);
            return NULL;
        }
        PRM += Delta->RealDelta;
        ModBufferPos = 0;
        ModCount += 1;
        // Bail out if we have too many PTMs to cope with:
        if (ModCount == MAX_PT_MODS)
        {
            return NULL;
        }
    }
    Match->ParentMass = PRM + PARENT_MASS_BOOST;
    return Match;
}


int GetPeptideParentMass(Peptide* Match)
{
    int Mass = PARENT_MASS_BOOST;
    char* Amino;
    int ModIndex;
    for (Amino = Match->Bases; *Amino; Amino++)
    {
        Mass += PeptideMass[*Amino];
    }
    for (ModIndex = 0; ModIndex < MAX_PT_MODS; ModIndex++)
    {
        if (Match->AminoIndex[ModIndex] >= 0)
        {
            Mass += Match->ModType[ModIndex]->RealDelta;
        }
    }
    Match->ParentMass = Mass;
    return Mass;
}

int LoadPeptideMassesCallback(int LineNumber, int FilePos, char* LineBuffer, void* UserData)
{
    float Mass;
    char* Str;
    char Letter;
    // Name:
    Str = strtok(LineBuffer, " ");
    // 3-letter abbrev:
    Str = strtok(NULL, " ");
    if (!Str)
    {
        REPORT_ERROR(0);
        return 0;
    }
    // 1-letter abbrev:
    Str = strtok(NULL, " ");
    if (!Str)
    {
        REPORT_ERROR(0);
        return 0;
    }
    Letter = Str[0];
    // mass:
    Str = strtok(NULL, " ");
    if (!Str)
    {
        REPORT_ERROR(0);
        return 0;
    }
    Mass = (float)atof(Str);
    ROUND_MASS(Mass, StandardPeptideMass[Letter]);
    return 1;
}

// Read peptide masses from a file.
int LoadPeptideMasses(char* FileName)
{
    FILE* AAFile;
    //
    if (!FileName)
    {
        // Use a sensible default:
        FileName = FILENAME_AMINO_ACID_MASSES;
    }
    AAFile = fopen(FileName, "r");
    if (!AAFile)
    {
        REPORT_ERROR_S(8, FileName);
        return 0;
    }
    ParseFileByLines(AAFile, LoadPeptideMassesCallback, NULL, 0);

    // This absurdly high mass for the delimiter, *, ensures that it won't be part of a match:
    StandardPeptideMass[42] = 9999999;
    memcpy(PeptideMass, StandardPeptideMass, sizeof(int) * 256);
    return 1;
}


// We read in large chunks of the file at once.  When we get past SCAN_BUFFER_A, it's time to shunt what we've got
// to the front of the buffer.  And if our buffer ends before SCAN_BUFFER_B, we try to read more data (until
// we reach eof)
#define SCAN_BUFFER_SIZE 5242880
#define SCAN_BUFFER_A 5232680
#define SCAN_BUFFER_B 5242680
#define RECORD_END '*'
 
// We have matched a tag in the peptide database, and the flanking series (plus some PTMods) matches our
// flanking mass.  Check to be sure these PTMods can be attached to the flanking series.  (For example:
// if the flanking sequence AAG plus a phosphate mass matches our prefix mass, that's NOT a match, because
// neither glycine nor alanine is phosphorylatable)
// Simplification: Assume that multiple PTMods can be attached to one base.  (This assumption isn't always valid, but
// it's nontrivial to know when it is; the user can toss out any unreasonable constructs later)
// From Start to End, INCLUSIVE.
int CheckForPTAttachmentPoints(int DecorationMassIndex, char* Buffer, int Start, int End, int BufferDir)
{
    int ModIndex;
    int ModsLeft[MAX_PT_MODTYPE];
    int BufferPos;
    int Done;
    int PeptideIndex;
    int Legal;

    memcpy(ModsLeft, AllDecorations[DecorationMassIndex].Mods, sizeof(int)*MAX_PT_MODTYPE);
    for (BufferPos = Start; BufferPos <= End; BufferPos++)
    {
        Done = 1; //by default
        PeptideIndex = Buffer[BufferPos] - 'A';
        for (ModIndex = 0; ModIndex < AllPTModCount; ModIndex++)
        {
            if (ModsLeft[ModIndex])
            {
                Legal = 1;
                // Avoid attaching a C-terminal PTM, if we're not at the C terminus:
                if (AllKnownPTMods[ModIndex].Flags & DELTA_FLAG_C_TERMINAL)
                {
                    if (BufferDir < 0 || BufferPos != End)
                    {
                        Legal = 0;
                    }
                }
                // Avoid attaching an N-terminal PTM, if we're not at the N terminus:
                if (AllKnownPTMods[ModIndex].Flags & DELTA_FLAG_N_TERMINAL)
                {
                    if (BufferDir > 0 || BufferPos != Start)
                    {
                        Legal = 0;
                    }
                }
                if (Legal)
                {
                    ModsLeft[ModIndex] = max(0, ModsLeft[ModIndex] - AllKnownPTMods[ModIndex].Allowed[PeptideIndex]);
                }
                if (ModsLeft[ModIndex])
                {
                    Done = 0;
                }
            }
        }
        if (Done)
        {
            return 1;
        }
    }
    return 0;
}

// MAX_SIDE_MODS is how many flanking matches we're allowed for an initial tag match.
// (For instance: The preceding aminos may match with no PTMs, or we may be able to match
// with one fewer amino and a PTM)
#define MAX_SIDE_MODS 10
int LeftMatchPos[MAX_SIDE_MODS];
int LeftMatchDecoration[MAX_SIDE_MODS];
int RightMatchPos[MAX_SIDE_MODS];
int RightMatchDecoration[MAX_SIDE_MODS];

// MatchFlankingMass is called when we matched a trie tag, and we are checking whether the
// flanking amino acids match our prefix or suffix mass.
// WARNING: If there are two or more decorations with the same mass, this method will FAIL, because we'll only
// consider ONE such decoration.
int MatchFlankingMass(MSSpectrum* Spectrum, TrieTag* Tag, char* Buffer, int StartPos, int BufferDir, int BufferEnd, int MatchMass, int ModsRemaining)
{
    int MatchCount = 0;
    int Pos;
    int Mass;
    int Diff;
    int AbsDiff;
    int FlankingMass;
    int MandatoryDecorationChange = 0;
    int DecorationMassIndex;
    int Verbose = 0;
    int* MatchPos;
    int* MatchDecoration;
    int MinMatchMass = MatchMass - GlobalOptions->FlankingMassEpsilon;
    //
    if (BufferDir<0)
    {
        MatchPos = LeftMatchPos;
        MatchDecoration = LeftMatchDecoration;
    }
    else
    {
        MatchPos = RightMatchPos;
        MatchDecoration = RightMatchDecoration;
    }

    /////////////////////////////////////////////////////////
    // If prefix mass is zero, that qualifies as a match always.
    if (MatchMass < GlobalOptions->FlankingMassEpsilon) 
    {
        MatchPos[MatchCount] = StartPos - BufferDir;
        MatchDecoration[MatchCount] = PlainOldDecorationIndex; 
        return 1;
    }
    DecorationMassIndex = AllDecorationCount - 1;
    // Skip over any decorations that use up too many pt-mods:
    while (1)
    {
        if (AllDecorations[DecorationMassIndex].TotalMods > ModsRemaining)
        {
            DecorationMassIndex--;
            continue;
        }
        break;        
    }
    FlankingMass = 0;
    for (Pos = StartPos; Pos >= 0; Pos += BufferDir)
    {
        if (Pos >= BufferEnd)
        {
            break;
        }
        if (Buffer[Pos] == '>' || Buffer[Pos] == '*')
        {
            break;
        }
        Mass = PeptideMass[Buffer[Pos]];
        if (Mass == 0)
        {
            // Invalid peptide!
            break;
        }
        FlankingMass += Mass;
        Diff = MatchMass  - (FlankingMass + AllDecorations[DecorationMassIndex].Mass);
        AbsDiff = abs(Diff);
        if (AbsDiff < GlobalOptions->FlankingMassEpsilon)
        {
            // Aha!  This is *probably* a match.  Check to be sure we have the bases we need:
            if (CheckForPTAttachmentPoints(DecorationMassIndex, Buffer, min(Pos, StartPos), max(Pos, StartPos), BufferDir))
            {
                if (Verbose)
                {
                    printf("Side is match!  Dec-index %d, flank %.2f.\n", DecorationMassIndex, FlankingMass / (float)MASS_SCALE);
                }
                MatchPos[MatchCount] = Pos;
                MatchDecoration[MatchCount] = DecorationMassIndex;
                MatchCount++;
                if (MatchCount == MAX_SIDE_MODS)
                {
                    return MatchCount;
                }

            }
        }
        // Move the DecorationMassIndex, if needed.
        while (MandatoryDecorationChange || FlankingMass + AllDecorations[DecorationMassIndex].Mass > MinMatchMass)
        {
            // The flanking sequence's mass is significantly bigger than our (decorated) target mass.
            // Move to a smaller decoration:
            MandatoryDecorationChange = 0;
            DecorationMassIndex--;
            if (DecorationMassIndex<0)
            {
                break;
            }
            // Skip any decorations that include phosphorylation, if we're not on phospho mode:
            if (!GlobalOptions->PhosphorylationFlag && g_PhosphorylationMod>-1 && AllDecorations[DecorationMassIndex].Mods[g_PhosphorylationMod])
            {
                MandatoryDecorationChange = 1;
                continue;
            }
            if (AllDecorations[DecorationMassIndex].TotalMods > ModsRemaining)
            {
                continue;
            }
            // And, check for a match:
            Diff = MatchMass  - (FlankingMass + AllDecorations[DecorationMassIndex].Mass);
            AbsDiff = abs(Diff);
            if (AbsDiff < GlobalOptions->FlankingMassEpsilon) 
            {
                // Aha!  This is *probably* a match.  Check to be sure we have the bases we need:
                if (CheckForPTAttachmentPoints(DecorationMassIndex, Buffer, min(Pos, StartPos), max(Pos, StartPos), BufferDir))
                {
                    if (Verbose)
                    {
                        printf("Left is match!  Dec-index %d, flank %.2f.\n", DecorationMassIndex, FlankingMass / (float)MASS_SCALE);
                    }
                    MatchPos[MatchCount] = Pos;
                    MatchDecoration[MatchCount] = DecorationMassIndex;
                    MatchCount++;
                    if (MatchCount == MAX_SIDE_MODS)
                    {
                        return MatchCount;
                    }
                    MandatoryDecorationChange = 1;
                }
            }
        }
        if (DecorationMassIndex<0)
        {
            break;
        }
    }
    return MatchCount;
}
// We extend LEFT and RIGHT from the match region (running from BufferPos to BufferEnd, INCLUSIVE),
// attempting to match our tag's prefix mass and Suffix mass.  Extension works like this:
// - DecoratedMassIndex starts out pointing at our largest decoration, FlankingMass starts at 0
// - At each iteration step:
// --  move one base further along, and add its mass to FlankingMass
// --  If FlankingMass plus the mass of our decoration matches our tag, we have a match.
// --  If FlankingMass plus the mass of our decoration is too LARGE, decrement DecoratedMassIndex
//        until we have a match, run out of decorations, or the mass again becomes too SMALL.
// --  At some point, we'll run out of decorations (FlankingMass becomes larger than the tag mass), and stop.
void GetMatches(SearchInfo* Info, TrieNode* Node, char* Buffer, int BufferPos, int BufferEnd, int FilePos)
{
    TrieTagHanger* TagNode;
    int ModsRemaining;

    int LeftMatchCount;
    int RightMatchCount;
    int LeftMatchIndex;
    int RightMatchIndex;
    int ModIndex;
    int UsedTooMany;
    static int PTMLimit[MAX_PT_MODTYPE];
    // To avoid repeated scoring:
    int ExtensionIndex = 0;
    int ExtensionCount = 0;
    static int StartingPoints[512];
    static int EndingPoints[512];
    static int ExtensionLeftDecorations[512]; 
    static int ExtensionRightDecorations[512];
    static MSSpectrum* ExtensionSpectra[512];
    int ExtensionFound;
    MSSpectrum* Spectrum;
    //
    if (!Node->FirstTag)
    {
        return;
    }
    //GlobalStats->TagMatches++;
    //Log("Extend matches of '%s' at position %d\n", Node->FirstTag->Tag->Tag, FilePos);
    // Try each tag corresponding to this TrieNode.
    for (TagNode = Node->FirstTag; TagNode; TagNode = TagNode->Next)
    {
        Spectrum = TagNode->Tag->PSpectrum;
        Info->Spectrum = Spectrum;
        memcpy(PTMLimit, g_PTMLimit, sizeof(int) * AllPTModCount);
        for (ModIndex = 0; ModIndex < MAX_PT_MODS; ModIndex++)
        {
            if (TagNode->Tag->AminoIndex[ModIndex] < 0)
            {
                break;
            }
            PTMLimit[TagNode->Tag->ModType[ModIndex]->Index] -= 1;
        }
        ModsRemaining = GlobalOptions->MaxPTMods - TagNode->Tag->ModsUsed;
        if (ModsRemaining < 0)
        {
            continue;
        }
        LeftMatchCount = MatchFlankingMass(Spectrum, TagNode->Tag, Buffer, BufferPos - TagNode->Tag->TagLength, -1, BufferEnd, TagNode->Tag->PrefixMass, ModsRemaining);
        if (LeftMatchCount == 0)
        {
            continue;
        }
        RightMatchCount = MatchFlankingMass(Spectrum, TagNode->Tag, Buffer, BufferPos + 1, 1, BufferEnd, TagNode->Tag->SuffixMass, ModsRemaining);
        if (RightMatchCount == 0)
        {
            continue;
        }
        // Consider each combination of left-extension and right-extension:
        for (LeftMatchIndex = 0; LeftMatchIndex<LeftMatchCount; LeftMatchIndex++)
        {
            for (RightMatchIndex = 0; RightMatchIndex<RightMatchCount; RightMatchIndex++)
            {
                UsedTooMany = 0;
                for (ModIndex = 0; ModIndex < AllPTModCount; ModIndex++)
                {
                    if (AllDecorations[LeftMatchDecoration[LeftMatchIndex]].Mods[ModIndex] + 
                        AllDecorations[RightMatchDecoration[RightMatchIndex]].Mods[ModIndex] > PTMLimit[ModIndex])
                    {
                        UsedTooMany = 1;
                        break;
                    }
                }
                if (UsedTooMany)
                {
                    continue;
                }
                if (AllDecorations[LeftMatchDecoration[LeftMatchIndex]].TotalMods + 
                    AllDecorations[RightMatchDecoration[RightMatchIndex]].TotalMods > ModsRemaining)
                {
                    continue;
                }
                if (GlobalOptions->MandatoryModIndex > -1 && 
                    !TagNode->Tag->MandatoryModUsed &&
                    AllDecorations[LeftMatchDecoration[LeftMatchIndex]].Mods[GlobalOptions->MandatoryModIndex] == 0 &&
                    AllDecorations[RightMatchDecoration[RightMatchIndex]].Mods[GlobalOptions->MandatoryModIndex] == 0)
                {
                    continue; // We don't have our mandatory PTM (biotin, or whatever)
                }
                ExtensionFound = 0;
                for (ExtensionIndex = 0; ExtensionIndex < ExtensionCount; ExtensionIndex++)
                {
                    if (StartingPoints[ExtensionIndex] == LeftMatchPos[LeftMatchIndex] && EndingPoints[ExtensionIndex] == RightMatchPos[RightMatchIndex]
                    && ExtensionLeftDecorations[ExtensionIndex] == LeftMatchDecoration[LeftMatchIndex]
                    && ExtensionRightDecorations[ExtensionIndex] == RightMatchDecoration[RightMatchIndex]
                    && ExtensionSpectra[ExtensionIndex] == TagNode->Tag->PSpectrum)
                    {
                        ExtensionFound = 1;
                        break;
                    }
                }
                if (ExtensionFound)
                {
                    continue;
                }
                StartingPoints[ExtensionCount] = LeftMatchPos[LeftMatchIndex];
                EndingPoints[ExtensionCount] = RightMatchPos[RightMatchIndex];
                ExtensionLeftDecorations[ExtensionCount] = LeftMatchDecoration[LeftMatchIndex];
                ExtensionRightDecorations[ExtensionCount] = RightMatchDecoration[RightMatchIndex];
                ExtensionSpectra[ExtensionCount] = TagNode->Tag->PSpectrum;
                Info->Spectrum = TagNode->Tag->PSpectrum;
                AddNewMatch(Info, FilePos, TagNode->Tag, 
                    Buffer + LeftMatchPos[LeftMatchIndex], 
                    RightMatchPos[RightMatchIndex] - LeftMatchPos[LeftMatchIndex] + 1,
                    (BufferPos - Node->Depth + 1) - LeftMatchPos[LeftMatchIndex],
                    LeftMatchDecoration[LeftMatchIndex], RightMatchDecoration[RightMatchIndex],
                    0, 0);
                ExtensionCount = min(511, ExtensionCount);
            }
        }
    }
    return;
}

//Extending Tags for a blind search requires a separate function.  
//We keep Tags where only one side (suffix or prefix) is extendable.  It is a simple 
//extension, because no PTMs are allowed.  If both sides are extendable
//then it is a nomod match, and sent to the regular scorer. 
int ExtendTagMatchBlind(SearchInfo* Info, TrieNode* Node, char* Buffer, int BufferPos, int BufferEnd, int FilePos)
{
    TrieTagHanger* Hanger;
    MSSpectrum* Spectrum;
    int LeftMatchCount;
    int RightMatchCount;
    int ModsRemaining = 0; //always zero for this simple extension
    int Extensions = 0;

    for (Hanger = Node->FirstTag; Hanger; Hanger = Hanger->Next)
    {
        Spectrum = Hanger->Tag->PSpectrum;
        Info->Spectrum = Spectrum;
        //by virtue of getting here, we know that this TAG (tripeptide) has matched the database
        Hanger->Tag->DBTagMatches++;
        LeftMatchCount = MatchFlankingMass(Spectrum, Hanger->Tag, Buffer, BufferPos - Hanger->Tag->TagLength, -1, BufferEnd, Hanger->Tag->PrefixMass, ModsRemaining);
        RightMatchCount = MatchFlankingMass(Spectrum, Hanger->Tag, Buffer, BufferPos + 1, 1, BufferEnd, Hanger->Tag->SuffixMass, ModsRemaining);
        if (LeftMatchCount + RightMatchCount == 1)
        {
            //set up the BlindTagMatchObject, representing this match.
            //Match = NewBlindTagMatch();
            //Match->Tag = Hanger->Tag;
            //Match->TagDBLoc = BufferPos - Hanger->Tag->TagLength; //pos of the first char
            if (LeftMatchCount)
            {
                Hanger->Tag->PrefixExtends ++;
            //    Match->ExtendLR = -1;
            //    Match->ExtendDBLoc = LeftMatchPos[0]; //only one match position possible, bc no mods
            //    Match->ExtendLength = Match->TagDBLoc - Match->ExtendDBLoc;
            }
            else
            {
                Hanger->Tag->SuffixExtends ++;
            //    Match->ExtendLR = 1; //right extension
            //    Match->ExtendDBLoc = RightMatchPos[0];
            //    Match->ExtendLength = Match->ExtendDBLoc - Match->TagDBLoc;
            }
            //InsertBlindTagMatch(Match);
            Hanger->Tag->PrefixExtends += LeftMatchCount;
            Hanger->Tag->SuffixExtends += RightMatchCount;
            Extensions++;
        }
        else if (LeftMatchCount + RightMatchCount == 2)
        {
            //send to regular scorer, it's a two sided hit
        }
        //printf ("Extend matches of '%s' at position %d\n", Node->FirstTag->Tag->Tag, FilePos);
        //printf ("Returned RightMatch %d, returned LeftMatch %d\n",RightMatchCount,LeftMatchCount);
    }
    return Extensions;
}

void GetProteinID(int RecordNumber, DatabaseFile* DB, char* Name)
{
    int Dummy[16];
    int RecordSize;
    if (!DB || !DB->IndexFile)
    {
        Name[0] = '?';
        Name[1] = '\0';
        return;
    }
    if (DB->Type == evDBTypeSpliceDB)
    {
        RecordSize = SPLICEDB_INDEX_RECORD_SIZE;
    }
    else
    {
        RecordSize = TRIE_INDEX_RECORD_SIZE;
    }
    fseek(DB->IndexFile, TRIE_INDEX_RECORD_SIZE * RecordNumber, SEEK_SET);
    ReadBinary(&Dummy, LONG_LONG_SIZE, 1, DB->IndexFile);
    ReadBinary(&Dummy, sizeof(int), 1, DB->IndexFile);
    ReadBinary(Name, sizeof(char), 80, DB->IndexFile);
    Name[80] = '\0';
    //Log("Record %d has ID %s\n", Pep->RecordNumber, Pep->Name);
}


void SortModifications(int* AminoIndex, MassDelta** ModType)
{
    int AminoIndexSorted[MAX_PT_MODS];
    MassDelta* ModTypeSorted[MAX_PT_MODS];
    int MinAminoIndex = 0;
    int NextSortedPosition = 0;
    int Index;
    int MinAminoPos = 0;
    //
    memset(AminoIndexSorted, -1, sizeof(int)*MAX_PT_MODS);
    memset(ModTypeSorted, 0, sizeof(MassDelta*)*MAX_PT_MODS);
    while (1)
    {
        // Find the smallest amino acid index in AminoIndex, skipping
        // over entries of -1 (which are empty)
        MinAminoIndex = -1;
        for (Index = 0; Index < MAX_PT_MODS; Index++)
        {
            if (AminoIndex[Index]>-1 && (MinAminoIndex<0 || AminoIndex[Index]<MinAminoIndex))
            {
                MinAminoIndex = AminoIndex[Index];
                MinAminoPos = Index;
            }
        }
        if (MinAminoIndex==-1)
        {
            // Everything's been moved to the sorted list.  Jane, stop this crazy thing!
            break;
        }
        // MOVE these entries into the sorted list:
        AminoIndexSorted[NextSortedPosition] = AminoIndex[MinAminoPos];
        AminoIndex[MinAminoPos] = -1;
        ModTypeSorted[NextSortedPosition] = ModType[MinAminoPos];
        ModType[MinAminoPos] = NULL;
        NextSortedPosition++;
    }
    // Move the sorted shadows back into the real arrays:
    memcpy(AminoIndex, AminoIndexSorted, sizeof(int)*MAX_PT_MODS);
    memcpy(ModType, ModTypeSorted, sizeof(MassDelta*)*MAX_PT_MODS);
    // Hooray!
}

#define SCORE_PTM_ATTACH_IMPOSSIBLE (float)-999999999.0

// Diagram of the dynamic programming table for optimal mod positioning:
// Suppose we have three decorations (zero, one or two attachments of the same PTM),
// and the PTMs should be attached at B and C in prefix ABCDE.  Then the grid
// looks like this:
//      A  B  C  D  E
//  0   x--x
//         |
//  1      x--x 
//            |
//  2         x--x--x
//
// (Columns for amino acids, rows for decorations, vertical moves mean an attachment)
//
// Find the optimal way to place modifications (from FullDecoration) on a polypeptide
// (Peptide) with length PeptideLength; store the mod-placements in AminoIndex / ModType
void FindOptimalPTModPositions(MSSpectrum* Spectrum, char* Peptide, 
    int PeptideLength, int FullDecoration, int BaseMass, int* AminoIndex, 
    MassDelta** ModType, int VerboseFlag, SpectrumTweak* Tweak)
{
    float* ScoreMatrix = NULL;
    int* SubDecorationMatrix = NULL;
    int PeptidePos;
    int DecorationIndex;
    float BestScore;
    char Amino;
    int ModIndex;
    int Mass;
    float Score;
    int ModCount;
    int CanBridge;
    int ModsNeeded;
    float BYScore;
    int ModAdder;
    int AminoAcidIndex;
    int BestSubDecoration;
    int SubDecoration;
    PRMBayesianModel* Model;
    /// 
    VerboseFlag  = 0;
    memset(AminoIndex, -1, sizeof(int) * MAX_PT_MODS);
    memset(ModType, 0, sizeof(MassDelta*) * MAX_PT_MODS);
    if (FullDecoration == PlainOldDecorationIndex)
    {
        return; // No mods to place!
    }
    if (Spectrum->Charge > 2)
    {
        Model = PRMModelCharge3;
    }
    else
    {
        Model = PRMModelCharge2;
    }
    // D.P. tables.  ScoreMatrix holds the score at each cell; SubDecorationMatrix tells
    // the previous subdecoration (and hence, how to backtrack through the table)
    ScoreMatrix = (float*)calloc(sizeof(float), PeptideLength * AllDecorationCount);
    SubDecorationMatrix = (int*)calloc(sizeof(int), PeptideLength * AllDecorationCount);

    // Fill the dynamic programming table.  Outer loop over amino acids,
    // inner loop over decorations.
    Mass = BaseMass;
    for (PeptidePos = 0; PeptidePos < PeptideLength; PeptidePos++)
    {
        Amino = Peptide[PeptidePos];
        AminoAcidIndex = Amino - 'A';
        Mass += PeptideMass[Amino];
        BestScore = 0;
        for (DecorationIndex = 0; DecorationIndex < AllDecorationCount; DecorationIndex++)
        {
            if (!IsSubDecoration[DecorationIndex][FullDecoration])
            {
                continue;
            }
            BestScore = 0; 
            BestSubDecoration = DecorationIndex;
            BYScore = GetIonPRMFeatures(Spectrum, Tweak, Model, Mass + AllDecorations[DecorationIndex].Mass, 0);
            //BYScore = (int)(100 * GetPRMFeatures(Spectrum, Tweak, Model, Mass + AllDecorations[DecorationIndex].Mass, 0));
            if (PeptidePos)
            {
                // Consider attaching nothing at this peptide:
                BestScore += ScoreMatrix[(PeptidePos - 1) * AllDecorationCount + DecorationIndex];
            }
            else
            {
                if (DecorationIndex != PlainOldDecorationIndex)
                {
                    BestScore += SCORE_PTM_ATTACH_IMPOSSIBLE;
                }
            }
            BestScore += BYScore;
            BestSubDecoration = DecorationIndex;
			//Log printf("    No mod here: Score %.2f\n", BestScore);

            // Consider attaching a modification at this peptide:
            for (SubDecoration = 0; SubDecoration < AllDecorationCount; SubDecoration++)
            {
                if (SubDecoration == DecorationIndex)
                {
                    continue;
                }
                if (!IsSubDecoration[SubDecoration][DecorationIndex])
                {
                    continue;
                }
                CanBridge = 1;
                for (ModIndex = 0; ModIndex < AllPTModCount; ModIndex++)
                {
                    // This decoration must contain all the mods from the subdecoration:
                    ModsNeeded = AllDecorations[DecorationIndex].Mods[ModIndex] - AllDecorations[SubDecoration].Mods[ModIndex];
                    if (ModsNeeded < 0)
                    {
                        CanBridge = 0;
                        break;
                    }
                    // This amino acid must be able to support the modification(s):
                    if (ModsNeeded > AllKnownPTMods[ModIndex].Allowed[AminoAcidIndex])
                    {
                        CanBridge = 0;
                        break;
                    }
                    // If the decoration is terminal, then this attachment position must be terminal:
                    if (ModsNeeded)
                    {
                        if ((AllKnownPTMods[ModIndex].Flags & DELTA_FLAG_C_TERMINAL) && PeptidePos < (PeptideLength - 1))
                        {
                            CanBridge = 0;
                            break;
                        }
                        if ((AllKnownPTMods[ModIndex].Flags & DELTA_FLAG_N_TERMINAL) && PeptidePos)
                        {
                            CanBridge = 0;
                            break;
                        }
                    }
                }
                if (CanBridge)
                {
                    if (PeptidePos)
                    {
                        Score = ScoreMatrix[(PeptidePos - 1) * AllDecorationCount + SubDecoration];
                    }
                    else
                    {
                        if (SubDecoration != PlainOldDecorationIndex)
                        {
                            Score = SCORE_PTM_ATTACH_IMPOSSIBLE; // Impossible!
                        }
                        else
                        {
                            Score = 0;
                        }
                    }
                    //Log printf("    To Sub-decoration %d: Score %d\n", SubDecoration, Score);
                    Score += BYScore;
                    if (Score >= BestScore)
                    {
                        BestScore = Score;
                        BestSubDecoration = SubDecoration;
                    }
                }
            }
            if (VerboseFlag)
            {
                //Log printf("    PeptidePos %d decoration %d: \n      Mass %d BYscore %.2f, best score %.2f, sub decoration %d\n", PeptidePos, DecorationIndex, 
                    //(Mass + AllDecorations[DecorationIndex].Mass), BYScore, BestScore, BestSubDecoration);
            }
            ScoreMatrix[PeptidePos * AllDecorationCount + DecorationIndex] = BestScore;
            SubDecorationMatrix[PeptidePos * AllDecorationCount + DecorationIndex] = BestSubDecoration;
        }
    }
    // Fill in AminoIndex, ModType.  Start at the bottom right of the DP table (last amino acid,
    // and full decoration), work back to the top row (first amino acid, no more decorations)
    ModCount = 0;
    DecorationIndex = FullDecoration;
    PeptidePos = PeptideLength - 1;
    while (PeptidePos >= 0)
    {
        SubDecoration = SubDecorationMatrix[PeptidePos * AllDecorationCount + DecorationIndex];
        if (SubDecoration != DecorationIndex)
        {
            for (ModIndex = 0; ModIndex < MAX_PT_MODTYPE; ModIndex++)
            {
                ModsNeeded = AllDecorations[DecorationIndex].Mods[ModIndex] - AllDecorations[SubDecoration].Mods[ModIndex];
                for (ModAdder = 0; ModAdder<ModsNeeded; ModAdder++)
                {
                    AminoIndex[ModCount] = PeptidePos;
                    AminoAcidIndex = Peptide[PeptidePos] - 'A';
                    ModType[ModCount] = MassDeltaByIndex[AminoAcidIndex * MAX_PT_MODTYPE + ModIndex];
                    ModCount++;
                }
            }
        }
        PeptidePos--;
        DecorationIndex = SubDecoration;
    }
    // Free temp storage:
    SafeFree(ScoreMatrix);
    SafeFree(SubDecorationMatrix);
}

// Return TRUE if two matches are the same.
// If we're performing an exon-graph search, then we only consider matches
// to be the same if they have the same sequence AND genomic coordinates.
int IsMatchDuplicate(Peptide* Match, Peptide* OldMatch, int PeptideLength)
{
    int CompareGenomicLocation = 1;

    if (Match->DB && Match->DB->Type == evDBTypeTrie && OldMatch->DB && OldMatch->DB->Type == evDBTypeTrie)
    {
        CompareGenomicLocation = 0;
    }
    if (!CompareGenomicLocation)
    {
        if (!strncmp(Match->Bases, OldMatch->Bases, PeptideLength) && 
            !memcmp(Match->AminoIndex, OldMatch->AminoIndex, sizeof(int)*MAX_PT_MODS) && 
            !memcmp(Match->ModType, OldMatch->ModType, sizeof(int)*MAX_PT_MODS))
        {
            return 1;
        }
    }
    else
    {
        // For exon graph search, we consider a match to be different if it has a different
        // genomic location.  We may see the same peptide inside two different exons, and 
        // we may have different options for splicing.
        if (!strncmp(Match->Bases, OldMatch->Bases, PeptideLength) && 
            !memcmp(Match->AminoIndex, OldMatch->AminoIndex, sizeof(int)*MAX_PT_MODS) && 
            !memcmp(Match->ModType, OldMatch->ModType, sizeof(int)*MAX_PT_MODS) &&
            Match->GenomicLocationStart == OldMatch->GenomicLocationStart &&
            Match->GenomicLocationEnd == OldMatch->GenomicLocationEnd)
        {
            return 1;
        }
    }
    return 0;
}

// Store a match in this Spectrum's Node's match list.  Don't store duplicate matches.
// Don't store more than GlobalOptions->StoreMatchCount matches.  Keep matches sorted
// by InitialScore (or, if MQScoreFlag is set, by MatchQualityScore)
Peptide* StoreSpectralMatch(MSSpectrum* Spectrum, Peptide* Match, int PeptideLength, int MQScoreFlag)
{
    Peptide* OldMatch;
    Peptide* CrummyScoreOldMatch;
    int VerboseFlag = 0;
    int SameFlag = 0;
    SpectrumNode* Node = Spectrum->Node;
    int NTT;
    //
    if (GlobalOptions->RequireTermini)
    {
        NTT = CountTrypticTermini(Match);
        if (NTT < GlobalOptions->RequireTermini)
        {
            FreePeptideNode(Match);
            return NULL;
        }
    }
    //printf("Store match %d '%s'\n", Match->Score, Match->Bases); 
    if (!Node->FirstMatch)
    {
        Node->FirstMatch = Match;
        Node->LastMatch = Match;
        Node->MatchCount++;
    }
    else
    {
        OldMatch = Node->FirstMatch;
        while (1)
        {
            SameFlag = IsMatchDuplicate(Match, OldMatch, PeptideLength);
            // Check to see whether it's the SAME as an existing match:
            if (SameFlag)
            {
                // Old match is the same as our new peptide.  Free the new guy, and break:
                OldMatch->MatchQualityScore = max(OldMatch->MatchQualityScore, Match->MatchQualityScore); 
                OldMatch->InitialScore = max(OldMatch->InitialScore, Match->InitialScore); 
                SafeFree(Match);
                Match = OldMatch;
                //OldMatch->SeenCount++;
                break;
            }
            if ((MQScoreFlag && Match->MatchQualityScore > OldMatch->MatchQualityScore) || (!MQScoreFlag && Match->InitialScore > OldMatch->InitialScore))
            {
                if (Node->FirstMatch == OldMatch)
                {
                    Node->FirstMatch = Match;
                }
                Match->Next = OldMatch;
                Match->Prev = OldMatch->Prev;
                if (OldMatch->Prev)
                {
                    OldMatch->Prev->Next = Match;
                }
                OldMatch->Prev = Match;
                Node->MatchCount++;
                // It's possible that we've already seen this peptide, but with a lower score.  (Why a lower score?
                // probably because we searched with the WRONG parent mass before, and the RIGHT parent mass now!) So, iterate over
                // the rest of the old matches, and if any is the same as this match, free it.
                for (CrummyScoreOldMatch = Match->Next; CrummyScoreOldMatch; CrummyScoreOldMatch = CrummyScoreOldMatch->Next)
                {
                    SameFlag = IsMatchDuplicate(Match, CrummyScoreOldMatch, PeptideLength);
                    if (SameFlag)
                    {
                        if (Node->LastMatch == CrummyScoreOldMatch)
                        {
                            Node->LastMatch = Node->LastMatch->Prev;
                        }
                        if (CrummyScoreOldMatch->Next)
                        {
                            CrummyScoreOldMatch->Next->Prev = CrummyScoreOldMatch->Prev;
                        }
                        if (CrummyScoreOldMatch->Prev)
                        {
                            CrummyScoreOldMatch->Prev->Next = CrummyScoreOldMatch->Next;
                        }
                        FreePeptideNode(CrummyScoreOldMatch);
                        break;
                    }
                }
                break;
            }
            OldMatch = OldMatch->Next;
            if (!OldMatch)
            {
                // Save our new match at the end of the list.
                Node->LastMatch->Next = Match;
                Match->Prev = Node->LastMatch;
                Node->LastMatch = Match;
                Node->MatchCount++;
                break;
            }
        }
    }
    if (Node->MatchCount > GlobalOptions->StoreMatchCount)
    {
        if (Match == Node->LastMatch)
        {
            Match = NULL;
        }
        OldMatch = Node->LastMatch->Prev;
        FreePeptideNode(Node->LastMatch);
        Node->LastMatch = OldMatch;
        if (OldMatch)
        {
            OldMatch->Next = NULL;
        }
        Node->MatchCount--;
    }
    return Match;
}


// Record a new match in the global match list.  If it's a duplicate peptide, then
// don't add it again.
Peptide* AddNewMatch(SearchInfo* Info, int FilePos, TrieTag* Tag, char* MatchedBases, 
                 int PeptideLength, int TagPosition, int PrefixDecoration, int SuffixDecoration, 
                 int GenomicStart, int GenomicEnd)
{
    Peptide* Match;
    PeptideMatch* PepInfo;
    int AminoIndex[MAX_PT_MODS];
    MassDelta* ModType[MAX_PT_MODS];
    int PrefixAminoIndex[MAX_PT_MODS];
    MassDelta* PrefixModType[MAX_PT_MODS];
    int SuffixAminoIndex[MAX_PT_MODS];
    MassDelta* SuffixModType[MAX_PT_MODS];
    int Mass;
    int SuffixStart;
    int AminoPos;
    int ModIndex;
    int TotalMods = 0;
    float ScoreToBeat;
    //int Score;
    int VerboseFlag;
    char* Amino;
    //int PrecursorMass;
    int ParentMassError;
    MSSpectrum* Spectrum = Info->Spectrum;
    //
    memset(PrefixAminoIndex, -1, sizeof(int) * MAX_PT_MODS);
    memset(SuffixAminoIndex, -1, sizeof(int) * MAX_PT_MODS);
    memset(PrefixModType, 0, sizeof(MassDelta*) * MAX_PT_MODS);
    memset(SuffixModType, 0, sizeof(MassDelta*) * MAX_PT_MODS);
    memset(ModType, 0, sizeof(MassDelta*) * MAX_PT_MODS);
    memset(AminoIndex, -1, sizeof(int) * MAX_PT_MODS);

    SuffixStart = TagPosition + strlen(Tag->Tag);
    //Log("***Tag extension complete!  Prefix '%s' tag '%s' suffix '%s'\n", StrPrefix, Tag->Tag, StrSuffix); 
    //Log("Prefix mods %d, suffix mods %d\n", PrefixDecoration, SuffixDecoration); 
    TotalMods = Tag->ModsUsed + AllDecorations[PrefixDecoration].TotalMods + AllDecorations[SuffixDecoration].TotalMods;
    //////////////////////////////////////////////////
    // Optimally place the prefix and suffix PTMs:
    VerboseFlag = 0;

    ////////////////////////////////////////////////////////////////////////////////////////
    // Temporarily adjust the charge and parent mass to reflect this candidate:
    Spectrum->Charge = Tag->Charge;
    Spectrum->ParentMass = PARENT_MASS_BOOST;
    for (AminoPos = 0, Amino = MatchedBases; AminoPos<PeptideLength; AminoPos++,Amino++)
    {
        Spectrum->ParentMass += PeptideMass[*Amino];
    }
    Spectrum->ParentMass += AllDecorations[PrefixDecoration].Mass;
    Spectrum->ParentMass += AllDecorations[SuffixDecoration].Mass;
    for (ModIndex = 0; ModIndex < MAX_PT_MODS; ModIndex++)
    {
        if (!Tag->ModType[ModIndex])
        {
            break;
        }
        Spectrum->ParentMass += Tag->ModType[ModIndex]->RealDelta; 
    }
    ////////////////////////////////////////////////////////////////////////
    // Reject this parent mass, if it's too far from the theoretical mass.
	//Use the corrected parent mass from the tweak, not the file mass.
    //PrecursorMass = Spectrum->MZ * Spectrum->Charge - (HYDROGEN_MASS * (Spectrum->Charge - 1));
    //ParentMassError = PrecursorMass - Spectrum->ParentMass;
    ParentMassError = Tag->Tweak->ParentMass - Spectrum->ParentMass;
    if (abs(ParentMassError) > GlobalOptions->ParentMassEpsilon)
    {
        // *** Reject this match, it doesn't match the parent mass!
        return NULL;
    } 

    ////////////////////////////////////////////////////////////////////////////////////////
    FindOptimalPTModPositions(Spectrum, MatchedBases, TagPosition, PrefixDecoration, 0, 
        PrefixAminoIndex, PrefixModType, VerboseFlag, Tag->Tweak);
    // Get the starting mass for our suffix match:
    Mass = 0;
    for (AminoPos = 0; AminoPos < TagPosition + Tag->TagLength; AminoPos++)
    {
        Mass += PeptideMass[MatchedBases[AminoPos]];
    }
    Mass += AllDecorations[PrefixDecoration].Mass;
    for (ModIndex = 0; ModIndex < MAX_PT_MODS; ModIndex++)
    {
        if (!Tag->ModType[ModIndex])
        {
            break;
        }
        Mass += Tag->ModType[ModIndex]->RealDelta;
    }
    FindOptimalPTModPositions(Spectrum, 
        MatchedBases + TagPosition + Tag->TagLength, 
        PeptideLength - TagPosition - Tag->TagLength, 
        SuffixDecoration, 
        Mass, 
        SuffixAminoIndex, 
        SuffixModType, 0, Tag->Tweak);
    /////////////////////////////////////////////////////////
    // Merge all the mods into one array, then sort it:
    TotalMods = 0;
    for (ModIndex = 0; ModIndex < MAX_PT_MODS; ModIndex++)
    {
        if (!PrefixModType[ModIndex])
        {
            break;
        }
        ModType[TotalMods] = PrefixModType[ModIndex];
        AminoIndex[TotalMods] = PrefixAminoIndex[ModIndex];
        TotalMods++;
    }
    for (ModIndex = 0; ModIndex < MAX_PT_MODS; ModIndex++)
    {
        if (!Tag->ModType[ModIndex])
        {
            break;
        }
        ModType[TotalMods] = Tag->ModType[ModIndex];
        AminoIndex[TotalMods] = Tag->AminoIndex[ModIndex] + TagPosition;
        TotalMods++;
    }
    for (ModIndex = 0; ModIndex < MAX_PT_MODS; ModIndex++)
    {
        if (!SuffixModType[ModIndex])
        {
            break;
        }
        ModType[TotalMods] = SuffixModType[ModIndex];
        AminoIndex[TotalMods] = SuffixAminoIndex[ModIndex] + TagPosition + Tag->TagLength;
        TotalMods++;
    }
    SortModifications(AminoIndex, ModType);

    /////////////////////////////////////////////////////////
    // Score the match.  If the score's not good enough, then toss it:
    ScoreToBeat = -999999;
    if (Spectrum->Node->MatchCount >= GlobalOptions->StoreMatchCount)
    {
        ScoreToBeat = Spectrum->Node->LastMatch->MatchQualityScore;
    }
    Match = NewPeptideNode();
    Match->ParentMassError = ParentMassError;
    Match->Tweak = Tag->Tweak;
    Match->DB = Info->DB;
    memcpy(Match->AminoIndex, AminoIndex, sizeof(int)*MAX_PT_MODS);
    memcpy(Match->ModType, ModType, sizeof(int)*MAX_PT_MODS);
    if (FilePos)
    {
        Match->PrefixAmino = *(MatchedBases - 1);
    }
    Match->SuffixAmino = *(MatchedBases + PeptideLength);
    strncpy(Match->Bases, MatchedBases, PeptideLength);
    Match->FilePos = FilePos;
    Match->RecordNumber = Info->RecordNumber;
    VerboseFlag = 0;
    Spectrum->CandidatesScored++;
    Tag->CandidatesScored++;
    GetPeptideParentMass(Match);
    // Invoke the scoring function now:
    ComputeMQScoreFeatures(Spectrum, Match, Match->ScoreFeatures, 0);
#ifdef MQSCORE_USE_SVM
    Match->MatchQualityScore = SVMComputeMQScore(Spectrum, Match, Match->ScoreFeatures);
#else
    Match->MatchQualityScore = LDAComputeMQScore(Spectrum, Match, Match->ScoreFeatures);
#endif
    Match->InitialScore = (int)(1000 * Match->MatchQualityScore);
    Match->GenomicLocationEnd = GenomicEnd;
    Match->GenomicLocationStart = GenomicStart;
    if (Match->MatchQualityScore < ScoreToBeat)
    {
        // Not good enough - forget it!
        SafeFree(Match);
        return NULL;
    }
    // It's good enough to add to the list:
    Match = StoreSpectralMatch(Spectrum, Match, PeptideLength, 1);
    if (!Match)
    {
        return NULL;
    }

    // Store the match details, if requested:
    if (GlobalOptions->ReportAllMatches)
    {
        PepInfo = (PeptideMatch*)calloc(sizeof(PeptideMatch), 1);
        PepInfo->FilePos = FilePos;
        PepInfo->RecordNumber = Info->RecordNumber;
        PepInfo->Tag = Tag;
        if (Match->Last)
        {
            Match->Last->Next = PepInfo;
        }
        else
        {
            Match->First = PepInfo;
        }
        Match->Last = PepInfo;
    }
    return Match;
}

// Print a list (one per line) of all the decorations we generated for
// the available post-translational modifications.
void DebugPrintDecoratedMassList()
{
    int Index;
    int ModIndex;
    //
    printf("Decorated masses:  (%d in all)\n", DecorationMassCount);
    for (Index = 0; Index < DecorationMassCount; Index++)
    {
        printf("  %.2f: ",DecorationMasses[Index]);
        for (ModIndex = 0; ModIndex <= GlobalOptions->MaxPTMods; ModIndex++)
        {
            if (DecorationMassMods[Index][ModIndex]<0)
            {
                // That's all the modifications in this one.
                printf("(end)\n");
                break;
            }
            printf("%d: %s (%.2f), ", DecorationMassMods[Index][ModIndex], PTModName[DecorationMassMods[Index][ModIndex]], ModMasses[DecorationMassMods[Index][ModIndex]]);
        }
    }
    printf("End of decorated mass list.\n");
}

// Helper macro for quick-sort
#define DECO_SWAP(a,b) \
{ \
fSwap = Masses[(a)]; \
memcpy(TempSpace, Mods[(a)], sizeof(int) * GlobalOptions->MaxPTMods); \
Masses[(a)] = Masses[b]; \
memcpy(Mods[(a)], Mods[(b)], sizeof(int) * GlobalOptions->MaxPTMods); \
Masses[(b)] = fSwap; \
memcpy(Mods[(b)], TempSpace, sizeof(int) * GlobalOptions->MaxPTMods); \
}

// Sort decorations using QuickSort.  We're sorting the array Masses, but we'll 
// also make the corresponding
// changes to the 2D array Mods, to keep the arrays in synch.
// Reminder: Quick-sort is done recursively.  Take the first element of the array as a pivot, then 
// 'pseudo-sort' the remaining elements so that all the EARLY elements (those less than the pivot) 
// come before all the LATE elements (those larger than the pivot).  The 'pseudo-sort' is done by 
// moving a left-index and right-index in from the edges of the array until they meet).  
// Then - here's the recursion-part - use quick-sort to sort the early and late elements.
void QuickSortDecoratedMasses(float* Masses, int** Mods, int Count)
{
    float fSwap;
    int TempSpace[1024];
    float Pivot;
    int LeftIndex;
    int RightIndex;
    // Sorting a list of one element is easy:
    if (Count<2)
    {
        return;
    }
    // Sorting a list of two elements is easy:
    if (Count == 2)
    {
        if (Masses[0] > Masses[1])
        {
            DECO_SWAP(0,1);
        }
        return;
    }
    // Now the REAL case begins:
    Pivot = Masses[0];
    LeftIndex = 1;
    RightIndex = Count-1;
    while (LeftIndex < RightIndex)
    {
        while (Masses[LeftIndex] <= Pivot) 
        {
            LeftIndex++;
            if (LeftIndex == Count)
            {
                // Pivot element is the biggest of all!
                DECO_SWAP(0, Count-1);
                QuickSortDecoratedMasses(Masses, Mods, Count-1);
                return;
            }
        }
        while (Masses[RightIndex] > Pivot)
        {
            RightIndex--;
        }
        if (RightIndex == 0)
        {
            // Pivot element is the smallest of all!
            QuickSortDecoratedMasses(Masses+1, Mods+1, Count-1);
            return;
        }
        if (RightIndex > LeftIndex)
        {
            DECO_SWAP(LeftIndex, RightIndex);
        }
    }
    DECO_SWAP(0, RightIndex);
    QuickSortDecoratedMasses(Masses, Mods, RightIndex);
    QuickSortDecoratedMasses(Masses+RightIndex+1, Mods+RightIndex+1, Count-RightIndex-1);

}

int PopulateDecoratedMassList(float* TotalMass, int** Mods, 
                              float MassSoFar, int* UsedMods, int UsedModCount)
{
    int Index;
    int MinModIndex;
    int RecordsBuilt = 0;
    //
    // If our prefix is mod #1, don't do 1,0; just to 1,1 and onward.  (Decorations
    // are listed from lowest PTM-index to largest)
    if (UsedModCount)
    {
        MinModIndex = UsedMods[UsedModCount-1];
    }
    else
    {
        MinModIndex = 0;
    }
    // Consider adding no mods at all:
    for (Index = 0; Index < UsedModCount; Index++)
    {
        Mods[0][Index] = UsedMods[Index];
    }
    TotalMass[0] = MassSoFar;
    RecordsBuilt++;
    if (UsedModCount == GlobalOptions->MaxPTMods)
    {
        return 1;
    }
    // Ok: Extend with each legal (lexigraphically subsequent) modificaiton!
    for (Index = MinModIndex; Index < TotalPTMods; Index++)
    {
        UsedMods[UsedModCount] = Index;
        RecordsBuilt += PopulateDecoratedMassList(TotalMass + RecordsBuilt, Mods + RecordsBuilt, 
            MassSoFar + ModMasses[Index], UsedMods, UsedModCount+1);
    }
    return RecordsBuilt;

}

int GetDecoratedMassCount(int AvailableMods, int PermissibleModCount)
{
    int ModIndex;
    int Total;
    if (PermissibleModCount == 0)
    {
        return 1;
    }
    Total = 1; // If we add no more
    for (ModIndex = 0; ModIndex < AvailableMods; ModIndex++)
    {
        Total += GetDecoratedMassCount(AvailableMods - ModIndex, PermissibleModCount - 1);
    }
    return Total;
}

//Trie.c::ProcessGeneHitsBlindTag
//This function processes all the onesided hits to a single gene, from the
//blind tagging option, send to a function to find the PTM site, and the scores.
//1. Tags (or the container) is sent to the function SeekMatch1PTM
int ProcessGeneHitsBlindTag()
{
    BlindTagMatch* Match;
    int Counter = 0;
    
    for (Match = FirstBlindTag; Match; Match = Match->Next)
    {
        Counter ++;
    }
    printf ("Processing a gene with %d matches\n",Counter);
    return 1;
}
//Trie.c::IsIdenticalBlindTagMatches
//returns true (1) if the two tag matches are identical
//else false
//Conditions for Identity
//1. Tags from the same spectra and Tweak
//2. Tags have identical DBAnchorPoint
//3. Tags extend in the same direction.

int IsIdenticalBlindTagMatches(BlindTagMatch* NodeA, BlindTagMatch* NodeB)
{
    if (NodeA->Tag->PSpectrum != NodeB->Tag->PSpectrum)
    {
        return 0;
    }
    if (NodeA->Tag->Tweak != NodeB->Tag->Tweak)
    {
        return 0;
    }
    if (NodeA->ExtendDBLoc != NodeB->ExtendDBLoc)
    {
        return 0;
    }
    if (NodeA->ExtendLR != NodeB->ExtendLR)
    { //this one may be unnecessary but it is in there for completeness
        return 0;
    }
    return 1;
}
//Trie.c :: InsertBlindTagMatch
//Inserts an object into a linked list, first testing
//if the object is identical to an already existing entry.
//If an object is not inserted, then it is freed (because the calling
//function expects us to deal with this type of thing).  Similarly
//with an object in the list which must be replaced.
int InsertBlindTagMatch(BlindTagMatch* Match)
{
    BlindTagMatch* NodeA;
    BlindTagMatch* Prev = NULL; //in case we do some swapping in the list
    BlindTagMatch* Next = NULL;

    if (FirstBlindTag == NULL) //just started
    {
        FirstBlindTag = Match;
        LastBlindTag = Match;
        return 1;
    }

    //cycle through the list, and see if there are any identical tags.
    // if identical tags exist, then we keep only the one with the 
    //longer extension.  In the absence of any twin, we put it at the end
    for (NodeA = FirstBlindTag; NodeA; NodeA = NodeA->Next)
    {
        if (IsIdenticalBlindTagMatches(Match, NodeA))
        {  //decide who to keep
            if (NodeA->ExtendLength >= Match->ExtendLength)
            { //winner is already in the list
                FreeBlindTagMatch(Match);
                return 1;
            }
            //have to remove item in the list. swap in Match
            if (NodeA->Prev == NULL) //first Item
            {
                FirstBlindTag = Match;
                Next = NodeA->Next; //temp
                Match->Next = Next;
                Next->Prev = Match;
                FreeBlindTagMatch(NodeA);
                return 1;
            }
            if (NodeA->Next == NULL) //last item
            {
                LastBlindTag = Match;
                Prev = NodeA->Prev;
                Match->Prev = Prev;
                Prev->Next = Match;
                FreeBlindTagMatch(NodeA);
                return 1;
            }
            //default else, nodeA in the middle
            Prev = NodeA->Prev;
            Next = NodeA->Next;
            Match->Prev = Prev;
            Match->Next = Next;
            Prev->Next = Match;
            Next->Prev = Match;
            FreeBlindTagMatch(NodeA);
            return 1;
        }
    }
    LastBlindTag->Next = Match; //add onto the end
    Match->Prev = LastBlindTag; //point back
    LastBlindTag = Match; //move end
    return 1;
}
// Main method: Use a trie to search a data-file.  Return the number of proteins searched.
int ScanFileWithTrie(SearchInfo* Info)
{
    FILE* File;
    int FilePos = 0;
    char* Buffer;
    int BufferPos = 0;
    int BufferEnd = 0;
    int AnchorPos = -1; // -1 means that no anchor is set
    TrieNode* Node;
    TrieNode* NextNode;
    int IsEOF = 0;
    int BytesRead;
    int OldPos;
    int PaddingDistance = 50;
    int Verbose = 0;
    //
    Info->RecordNumber = 0;
    File = Info->DB->DBFile;
    if (!File)
    {
        return 0;
    }
    fseek(File, 0, 0);
    if (!Info->Root)
    {
        return 0;
    }
    Buffer = (char*)calloc(SCAN_BUFFER_SIZE, sizeof(char));
    Node = Info->Root;
    // We'll scan in chunks of the file, and scan across them.  We try to always keep a buffer of 50 characters
    // before and after the current position, so that we can look forward and back to get masses.  (When we match
    // a tag-string, we look at surrounding masses).
    while (1)
    {
        //printf("Anc%d Buf%d BufEnd%d F%d Char%c\n", AnchorPos, BufferPos, BufferEnd, FilePos, Buffer[BufferPos]);
        // Periodically shunt data toward the front of the buffer:
        if (BufferPos > SCAN_BUFFER_A && AnchorPos==-1)
        {
            // ......ppppBbbbbbbbbbE... <- diagram (p = pad, B = buffer start, E = buffer end)
            // ppppBbbbbbbbbbE....      <- after move
            memmove(Buffer, Buffer + BufferPos - PaddingDistance, BufferEnd - (BufferPos - PaddingDistance));
            BufferEnd -= (BufferPos - PaddingDistance);
            BufferPos = PaddingDistance;   
        }

        // Read more data, if we have room and we can:
        if (BufferEnd < SCAN_BUFFER_B && !IsEOF)
        {
            BytesRead = ReadBinary(Buffer + BufferEnd, sizeof(char), SCAN_BUFFER_SIZE - BufferEnd, File);
            if (!BytesRead)
            {
                IsEOF = 1;
            }
            BufferEnd += BytesRead;
        }

        if (AnchorPos!=-1)
        {
            // If we're anchored: Attempt to extend the current match.
            if (Buffer[BufferPos] >= 'A' && Buffer[BufferPos] <= 'Z')
            {
                NextNode = Node->Children[Buffer[BufferPos] - 'A'];
            }
            else
            {
                NextNode = NULL;
            }
            // If we can extend the current match...
            if (NextNode)
            {
                // Note any new matches:
                if (NextNode->FirstTag)
                {
                    if(GlobalOptions->RunMode & RUN_MODE_BLINDTAG)
                    {
                        ExtendTagMatchBlind(Info, NextNode, Buffer, BufferPos, BufferEnd, FilePos);
                    }
                    else
                    {
                        GetMatches(Info, NextNode, Buffer, BufferPos, BufferEnd, FilePos);
                    }
                }
                // Travel down the trie:
                Node = NextNode;
                BufferPos++;
                FilePos++;
            }
            else
            {
                // We could NOT extend the match.
                // We're done with this anchor.  Clear the anchor, and use our FailureNode to jump
                // forward.  (AnchorPos moves forward by FailureLetterSkip chars, and the BufferPos
                // moves to the correct distance ahead of the anchor)
                if (IS_ROOT(Node->FailureNode))
                {
                    AnchorPos = -1;
                }
                else
                {
                    AnchorPos = AnchorPos + Node->FailureLetterSkip;
                    OldPos = BufferPos;
                    BufferPos = AnchorPos + Node->FailureNode->Depth - 1;
                    FilePos += (BufferPos - OldPos);
                    // Process matches immediately:
                    if (Node->FailureNode->FirstTag)
                    {
                        if (GlobalOptions->RunMode & RUN_MODE_BLINDTAG)
                        {
                            ExtendTagMatchBlind(Info, NextNode, Buffer, BufferPos, BufferEnd, FilePos);
                        }
                        else
                        {
                            GetMatches(Info, Node->FailureNode, Buffer, BufferPos, BufferEnd, FilePos);
                        }
                    }
                    BufferPos++;
                    FilePos++;
                }
                Node = Node->FailureNode;
            }
        }
        else
        {
            // We're not currently anchored.  Process end-of-record tags, or attempt to start a 
            // brand new match.
            if (BufferPos>=BufferEnd || Buffer[BufferPos] == RECORD_END || !Buffer[BufferPos])
            {
                // END of a protein.
                if (GlobalOptions->RunMode & RUN_MODE_BLINDTAG)
                {
                    //ProcessGeneHitsBlindTag(); // Process the blind tags a gene at a time.
                    //FreeAllBlindTagMatches(FirstBlindTag); //free up the hits
                    //FirstBlindTag = NULL; //reset the pointers
                    //LastBlindTag = NULL;
                }
                Info->RecordNumber++;   
                AnchorPos = -1;
            }
            else 
            {
                // Now: Start a new match, if possible:
                if (Buffer[BufferPos] >= 'A' && Buffer[BufferPos] <= 'Z')
                {
                    NextNode = Node->Children[Buffer[BufferPos] - 'A'];
                }
                else
                {
                    NextNode = NULL;
                }
                if (NextNode)
                {
                    // Note any new matches.  (Not likely, because 
                    // at this point in the code, we're only at depth 1 in the
                    // tree; tags of length 1 aren't very good)
                    if (NextNode->FirstTag)
                    {
                        if (GlobalOptions->RunMode & RUN_MODE_BLINDTAG)
                        {
                            ExtendTagMatchBlind(Info, NextNode, Buffer, BufferPos, BufferEnd, FilePos);
                        }
                        else
                        {
                            GetMatches(Info, NextNode, Buffer, BufferPos, BufferEnd, FilePos);
                        }
                    }
                    Node = NextNode;
                    AnchorPos = BufferPos;
                }
            }
            BufferPos++;
            FilePos++;
            if (BufferPos >= BufferEnd)
            {
                break;
            }
        } // if not anchored
    } // Master while-loop

    SafeFree(Buffer);
    
    return Info->RecordNumber + 1;
}


// Print just the tags from our trie:
void DebugPrintTrieTags(TrieNode* Node)
{
    TrieTagHanger* Hanger;
    int ChildIndex;
    TrieNode* Failure;
    char TagBuffer[256];
    int Len;
    int ModIndex;
    if (!Node)
    {
      return;
    }
    for (Hanger = Node->FirstTag; Hanger; Hanger = Hanger->Next)
    {
        ModIndex = 0;
        TagBuffer[0] = Hanger->Tag->Tag[0];
        TagBuffer[1] = '\0';
        if (Hanger->Tag->ModType[ModIndex] && Hanger->Tag->AminoIndex[ModIndex] == 0)
        {
            strcat(TagBuffer, Hanger->Tag->ModType[ModIndex]->Name);
            ModIndex++;
        }
        Len = strlen(TagBuffer);
        TagBuffer[Len] = Hanger->Tag->Tag[1];
        TagBuffer[Len+1] = '\0';
        if (Hanger->Tag->ModType[ModIndex] && Hanger->Tag->AminoIndex[ModIndex] == 0)
        {
            strcat(TagBuffer, Hanger->Tag->ModType[ModIndex]->Name);
            ModIndex++;
        }
        Len = strlen(TagBuffer);
        TagBuffer[Len] = Hanger->Tag->Tag[2];
        TagBuffer[Len+1] = '\0';

        if (Hanger->Tag->ModType[ModIndex] && Hanger->Tag->AminoIndex[ModIndex] == 0)
        {
            strcat(TagBuffer, Hanger->Tag->ModType[ModIndex]->Name);
            ModIndex++;
        }

        printf("Tag '%s' (prefix %.2f, Suffix %.2f) %.2f hits %d\n", TagBuffer, 
            Hanger->Tag->PrefixMass / (float)MASS_SCALE, 
            Hanger->Tag->SuffixMass / (float)MASS_SCALE,
            Hanger->Tag->Score, Hanger->Tag->CandidatesScored);
#ifdef DEBUG_TAG_GENERATION
        printf("%s\n", Hanger->Tag->TagScoreDetails);
#endif
    }
    Failure = Node->FailureNode;
    if (Node->FirstTag && Failure)
    {
        printf("  Node %s has failure node depth %d letter %c.\n", Node->FirstTag->Tag->Tag, Failure->Depth, Failure->Letter);
    }

    for (ChildIndex = 0; ChildIndex < 26; ChildIndex++)
    {
        if (ChildIndex == 'I'-'A' || ChildIndex == 'Q'-'A')
        {
            continue; 
        }
        if (Node->Children[ChildIndex])
        {
            DebugPrintTrieTags(Node->Children[ChildIndex]);
        }
    }
}

void FlagMandatoryModUsage(TrieNode* Node)
{
    TrieTagHanger* Hanger;
    int CharIndex;
    int ModIndex;
    //
    if (!Node)
    {
        return;
    }
    for (Hanger = Node->FirstTag; Hanger; Hanger = Hanger->Next)
    {
        for (ModIndex = 0; ModIndex < MAX_PT_MODS; ModIndex++)
        {
            if (Hanger->Tag->ModType[ModIndex]->Index == GlobalOptions->MandatoryModIndex)
            {
                Hanger->Tag->MandatoryModUsed = 1;
            }
        }
    }

    for (CharIndex = 0; CharIndex < TRIE_CHILD_COUNT; CharIndex++)
    {
        FlagMandatoryModUsage(Node->Children[CharIndex]);
    }
}

// COPYPASTA from WriteMatchToString.
void WriteTagToString(TrieTag* Tag, char* Buffer, int IncludeMods)
{
    char* Stuff;
    int AminoPos;
    char NameChar;
    int ModIndex;
    int NameIndex;
    //
    Stuff = Buffer;

    for (AminoPos = 0; AminoPos < (int)strlen(Tag->Tag); AminoPos++)
    {
        *Stuff++ = Tag->Tag[AminoPos];
        if (IncludeMods)
        {
            for (ModIndex = 0; ModIndex < MAX_PT_MODS; ModIndex++)
            {
                if (Tag->AminoIndex[ModIndex] == AminoPos && Tag->ModType[ModIndex])
                {
                    // Write the modification:
                    for (NameIndex = 0; NameIndex < 4; NameIndex++)
                    {
                        NameChar = Tag->ModType[ModIndex]->Name[NameIndex];
                        if (!NameChar)
                        {
                            break;
                        }
                        *Stuff++ = ConvertToLower(NameChar);
                    }
                }
            }
        }
    }
    *Stuff = '\0';
}

// Write (to a char buffer) the string version of a peptide, including modifications.
// For example: "EAM+16APK".  Similar to the method PeptideClass.GetModdedName
void WriteMatchToString(Peptide* Match, char* Buffer, int IncludeMods)
{
    char* Stuff;
    int AminoPos;
    char NameChar;
    int ModIndex;
    int NameIndex;
    //
    Stuff = Buffer;
    
    if (Match->PrefixAmino)
    {
        *Stuff++ = Match->PrefixAmino;
    }
    else
    {
        *Stuff++ = '-';
    }
    *Stuff++ = '.';
    for (AminoPos = 0; AminoPos < (int)strlen(Match->Bases); AminoPos++)
    {
        *Stuff++ = Match->Bases[AminoPos];
        if (IncludeMods)
        {
            for (ModIndex = 0; ModIndex < MAX_PT_MODS; ModIndex++)
            {
                if (Match->AminoIndex[ModIndex]==AminoPos)
                {
                    // Write the modification:
                    for (NameIndex = 0; NameIndex < 4; NameIndex++)
                    {
                        NameChar = Match->ModType[ModIndex]->Name[NameIndex];
                        if (!NameChar)
                        {
                            break;
                        }
                        *Stuff++ = ConvertToLower(NameChar);
                    }
                }
            }
        }
    }
    *Stuff++ = '.';
    if (Match->SuffixAmino)
    {
        *Stuff++ = Match->SuffixAmino;
    }
    else
    {
        *Stuff++ = '-';
    }
    *Stuff = '\0';
}
