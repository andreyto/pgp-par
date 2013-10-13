//
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
#include "Utils.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "Errors.h"
#include "Trie.h"
#include "Inspect.h"
#include "Spliced.h"
#include "SpliceDB.h"

// SpliceScan.c is not used in production.  Its job is to take a standard
// database (.trie format), and check whether the proteins in that database
// are present in a splice-tolerant database.  Because EST coverage is 
// incomplete, we may be missing some exons from a gene, or missing some 
// genes entirely.  And due to polymorphisms, errors
// in gene and protein sequencing, there will be some minor differences.
// We want to quantify how many genes are missing, and how extensive the
// differences are.

// We do the following for each protein:
// - Take all 8-mers from the protein, put them in a trie.  
//  ASSUMPTION: If the protein is present at all, there will be 8 consecutive residues 
//   present without error
//  ASSUMPTION: any proteins with length <8aa can be ignored (yes, IPI has records of length <8...)
// - Use the trie to search each gene in the splice-tolerant database
// - Count the number of leaves that were matched.  If the rate is above (near?) our best so far,
//   flag all the characters that were matched.  Remember the percentage of characters matched,
//   span from the first to the last, and which gene record was so good.
//  ASSUMPTION: 8mers aren't repeated within a protein
//  ASSUMPTION: Chance matches are very rare.  Our numbers will be distorted only slightly
//   by considering 8mers that appear at the wrong position.  
// - When you run out of genes, or when you get 95% or better match, report the best match.  
//  Free the trie and start the next record.

//#define SS_TRIE_DEBUG 1

#define PROTEIN_NAME_BLOCK 81
#define SS_BLOCK_SIZE 1000
#define MAX_TRIE_NODE_COUNT 2000000
int g_TrieNodeMatches[SS_BLOCK_SIZE];
int TrieNodeHitFlags[MAX_TRIE_NODE_COUNT];
int g_NextTrieLeafIndex;

typedef struct SSTrieNode
{
    void* Children[26]; // SSTrieNodes or SSTrieLeafs
    struct SSTrieNode* FailureNode;
    int FailureDepth;
#ifdef SS_TRIE_DEBUG
    int Depth; // for debugging only!
    char Buffer[16]; // for debugging only!
#endif // SS_TRIE_DEBUG
} SSTrieNode;

typedef struct SSTrieLeafNode
{
    int ProteinNumber;
    int ProteinPos;
    struct SSTrieLeafNode* Next;
} SSTrieLeafNode;

typedef struct SSTrieLeaf
{
    struct SSTrieNode* FailureNode;
    int FailureDepth;
    // Index is a pointer into the array MAX_TRIE_NODE_COUNT.  We use an index
    // instead of storing the hit flag in the leaf.  Why?  Because we have 
    // to reset all the flags to zero after every gene record.  Doing so with 
    // memset is much faster than traversing the trie!
    int Index; 
    SSTrieLeafNode* Head;
#ifdef SS_TRIE_DEBUG
    char Buffer[16]; // for debugging only!
#endif // SS_TRIE_DEBUG
} SSTrieLeaf;

void SSTrieFailureNodeHelper(SSTrieNode* Root, char* Buffer, SSTrieNode* FailedNode, int Depth)
{
    int SuffixStart;
    int BufferPos;
    int FailureDepth;
    SSTrieNode* FailureNode = NULL;
    SSTrieNode* Node;
    SSTrieLeaf* Leaf;
    int AA;
    //
    ////////////////////////////////////////////////////////////////////////////
    // Set this node's failure-node, by finding the longest proper suffix of Buffer
    // that reaches a node:
    if (Depth > 1)
    {
        for (SuffixStart = 1; SuffixStart < Depth; SuffixStart++)
        {
            Node = Root;
            for (BufferPos = SuffixStart; BufferPos < Depth; BufferPos++)
            {
                Node = Node->Children[Buffer[BufferPos]];
                if (!Node)
                {
                    break;
                }
            }
            if (Node)
            {
                // The suffix matched!
                FailureDepth = Depth - SuffixStart;
                FailureNode = Node;
                break;
            }
        }
        if (!FailureNode)
        {
            FailureNode = Root;
            FailureDepth = 0;
        }
        if (Depth == 8)
        {
            Leaf = (SSTrieLeaf*)FailedNode;
            Leaf->FailureDepth = FailureDepth;
            Leaf->FailureNode = FailureNode;
        }
        else
        {
            FailedNode->FailureDepth = FailureDepth;
            FailedNode->FailureNode = FailureNode;
        }
    }
    else
    {
        // A depth-1 node.  Always gets the root as its failure node:
        FailedNode->FailureDepth = 0;
        FailedNode->FailureNode = Root;
    }
    ////////////////////////////////////////////////////////////////////////////
    // Set our children's failure-nodes:
    if (Depth < 8)
    {
        for (AA = 0; AA < 26; AA++)
        {
            Node = FailedNode->Children[AA];
            if (Node)
            {
                Buffer[Depth] = AA;
                SSTrieFailureNodeHelper(Root, Buffer, Node, Depth + 1);
            }
        }
    }
}

// Initialize all the failure nodes for the trie.
void SetSSTrieFailureNodes(SSTrieNode* Root)
{
    char Buffer[16];
    int AA;
    SSTrieNode* Child;

    // Root never gets a failure node:
    Root->FailureDepth = 0;
    Root->FailureNode = NULL;

    // For other nodes: Populate a Buffer with your string.  Then find the 
    // longest proper suffix of Buffer that reaches a node.
    for (AA = 0; AA < 26; AA++)
    {
        Child = Root->Children[AA];
        if (Child)
        {
            Buffer[0] = AA;
            SSTrieFailureNodeHelper(Root, Buffer, Child, 1);
        }
    }
}

void FreeSSTrieNode(SSTrieNode* Root, int Depth)
{
    SSTrieLeaf* Leaf;
    SSTrieLeafNode* Node;
    SSTrieLeafNode* Prev = NULL;
    int AA;
    //
    if (!Root)
    {
        return;
    }
    if (Depth == 8)
    {
        Leaf = (SSTrieLeaf*)Root;
        for (Node = Leaf->Head; Node; Node = Node->Next)
        {
            SafeFree(Prev);
            Prev = Node;
        }
        SafeFree(Prev);
        SafeFree(Leaf);
        return;
    }
    for (AA = 0; AA < 26; AA++)
    {
        if (Root->Children[AA])
        {
            FreeSSTrieNode(Root->Children[AA], Depth + 1);
        }
    }
    SafeFree(Root);
}

SSTrieNode* ConstructSSTrie(char** Sequences, int BlockSize)
{
    SSTrieNode* Root;
    SSTrieNode* CurrentNode;
    SSTrieNode* NextNode;
    SSTrieLeaf* Leaf;
    SSTrieLeafNode* Node;
    int Len;
    int StartPos;
    int PeptidePos;
    int SequencePos;
    int AA;
    int ProteinNumber;
    char* Sequence;
    //
    Root = (SSTrieNode*)calloc(1, sizeof(SSTrieNode));
    g_NextTrieLeafIndex = 0;
    for (ProteinNumber = 0; ProteinNumber < BlockSize; ProteinNumber++)
    {
        Sequence = Sequences[ProteinNumber];
        Len = strlen(Sequence);
        for (StartPos = 0; StartPos <= Len-8; StartPos++)
        {
            SequencePos = StartPos;
            PeptidePos = 0;
            CurrentNode = Root;
            // Add nodes for the first n-1 positions:
            for (PeptidePos = 0; PeptidePos < 7; PeptidePos++)
            {
                AA = Sequence[StartPos + PeptidePos] - 'A';
                if (AA < 0 || AA > 25)
                {
                    break; // invalid character in protein sequence!
                }
                NextNode = CurrentNode->Children[AA];
                if (!NextNode)
                {
                    NextNode = (SSTrieNode*)calloc(1, sizeof(SSTrieNode));
                    CurrentNode->Children[AA] = NextNode;
#ifdef SS_TRIE_DEBUG
                    NextNode->Depth = PeptidePos + 1;
                    strncpy(NextNode->Buffer, Sequence + StartPos, PeptidePos + 1);
                    NextNode->Buffer[PeptidePos + 1] = '\0';
#endif //SS_TRIE_DEBUG
                }
                CurrentNode = NextNode;
            }
            // Add a leaf node for the nth position:
            AA = Sequence[StartPos + PeptidePos] - 'A';
            if (AA < 0 || AA > 25)
            {
                continue; // invalid character in protein sequence!
            }
            Leaf = CurrentNode->Children[AA];
            if (!Leaf)
            {
                Leaf = (SSTrieLeaf*)calloc(1, sizeof(SSTrieLeaf));
                //Leaf->ProteinPos = StartPos + PeptidePos;
                Leaf->Index = g_NextTrieLeafIndex++;
                CurrentNode->Children[AA] = Leaf;
                Leaf->Head = (SSTrieLeafNode*)calloc(1, sizeof(SSTrieLeafNode));
                Leaf->Head->ProteinNumber = ProteinNumber;
                Leaf->Head->ProteinPos = StartPos;
#ifdef SS_TRIE_DEBUG
                strncpy(Leaf->Buffer, Sequence + StartPos, 8);
                Leaf->Buffer[8] = '\0';
#endif
            }
            else
            {
                for (Node = Leaf->Head; Node->Next; Node = Node->Next)
                {
                    ;;
                }
                Node->Next = (SSTrieLeafNode*)calloc(1, sizeof(SSTrieLeafNode));
                Node->Next->ProteinNumber = ProteinNumber;
                Node->Next->ProteinPos = StartPos;
            }
        } // Loop on start positions
    } // Loop on proteins
    return Root;
}

int SSTrieCoverSequence(SSTrieNode* Root, char* MatchFlags, int Depth, int ProteinNumber)
{
    SSTrieLeaf* Leaf;
    int AA;
    int Sum = 0;
    SSTrieNode* Child;
    SSTrieLeafNode* Node;
    int X;
    //
    if (Depth == 8)
    {
        Leaf = (SSTrieLeaf*)Root;
        if (TrieNodeHitFlags[Leaf->Index])
        {
            for (Node = Leaf->Head; Node; Node = Node->Next)
            {
                if (Node->ProteinNumber == ProteinNumber)
                {
                    for (X = 0; X < 8; X++)
                    {
                        if (!MatchFlags[Node->ProteinPos + X])
                        {
                            Sum += 1;
                            MatchFlags[Node->ProteinPos + X] = 1;
                        }
                    }
                }
            }
        }
    }
    else
    {
        for (AA = 0; AA < 26; AA++)
        {
            Child = Root->Children[AA];
            if (Child)
            {
                Sum += SSTrieCoverSequence(Child, MatchFlags, Depth + 1, ProteinNumber);
            }
        }
    }
    return Sum;
}


// Recursive main function for scanning through splice-tolerant database with a trie.  
void SSDatabaseScanHelper(ExonStruct* Exon, int Len, int Pos, SSTrieNode* Node, int Depth)
{
    SSTrieLeaf* Leaf;
    int NextLen = 0;
    int EdgeIndex;
    SSTrieNode* NextNode;
    SSTrieLeafNode* LeafNode;
    ExonEdge* Edge;
    //
    if (!Node)
    {
        return;
    }
    if (Depth == 8)
    {
        Leaf = (SSTrieLeaf*)Node;
        if (!TrieNodeHitFlags[Leaf->Index])
        {
            //Leaf->HitFlag = 1;
            TrieNodeHitFlags[Leaf->Index] = 1;
            for (LeafNode = Leaf->Head; LeafNode; LeafNode = LeafNode->Next)
            {        
                g_TrieNodeMatches[LeafNode->ProteinNumber]++;
            }
        }
        return;
    }
    Len = Exon->Length;
    if (Pos >= Len)
    {
        for (EdgeIndex = 0; EdgeIndex < Exon->ForwardEdgeCount; EdgeIndex++)
        {
            Edge = Exon->ForwardEdges + EdgeIndex;
            if (Edge->AA)
            {
                NextNode = Node->Children[Edge->AA - 'A'];
                if (NextNode)
                {
                    NextLen = Edge->Exon->Length;
                    SSDatabaseScanHelper(Edge->Exon, NextLen, 0, NextNode, Depth + 1);
                }
            }
            else
            {
                NextLen = Edge->Exon->Length;
                SSDatabaseScanHelper(Edge->Exon, NextLen, 0, Node, Depth);
            }
        }
    }
    else
    {
        NextNode = Node->Children[Exon->Sequence[Pos] - 'A'];
        if (NextNode)
        {
            SSDatabaseScanHelper(Exon, Len, Pos + 1, NextNode, Depth + 1);
        }
    }
}

void DebugPrintSSTrie(SSTrieNode* Node, int Depth, char* Buffer)
{
    int AA;
    SSTrieNode* Child;
#ifdef SS_TRIE_DEBUG
    SSTrieNode* FailureNode;
#endif
    SSTrieLeafNode* LeafNode;
    SSTrieLeaf* Leaf;
    //
#ifdef SS_TRIE_DEBUG
    for (AA = 0; AA < Depth; AA++)
    {
        printf(" ");
    }
    if (Depth == 8)
    {
        Leaf = (SSTrieLeaf*)Node;
        FailureNode = Leaf->FailureNode;
        if (FailureNode)
        {
            printf("Leaf '%s' failure '%s' (depth %d)\n", Leaf->Buffer, FailureNode->Buffer, FailureNode->Depth);
        }
        else
        {
            printf("Leaf '%s' (NO FAILURE NODE)\n", Leaf->Buffer);
        }
    }
    else
    {
        FailureNode = Node->FailureNode;
        if (FailureNode)
        {
            printf("Node '%s' d%d failure node '%s' (depth %d)\n", Node->Buffer, Node->Depth, 
                FailureNode->Buffer, FailureNode->Depth);
        }
        else
        {
            printf("Node '%s' d%d NO FAILURE NODE\n", Node->Buffer, Node->Depth);
        }
        for (AA = 0; AA < 26; AA++)
        {
            Child = Node->Children[AA];
            if (Child)
            {
                Buffer[Depth] = 'A' + AA;
                DebugPrintSSTrie(Child, Depth + 1, Buffer);
            }
        }
    }
    return;
#endif
    if (Depth < 8)
    {
        for (AA = 0; AA < 26; AA++)
        {
            Child = Node->Children[AA];
            if (Child)
            {
                Buffer[Depth] = 'A' + AA;
                DebugPrintSSTrie(Child, Depth + 1, Buffer);
            }
        }
    }
    Buffer[8] = '\0';
    Leaf = (SSTrieLeaf*)Node;
    for (LeafNode = Leaf->Head; LeafNode; LeafNode = LeafNode->Next)
    {
        printf("%s at pos %d in record #%d\n", Buffer, LeafNode->ProteinPos, LeafNode->ProteinNumber);
    }
}

void SSDatabaseScanExon(ExonStruct* Exon, SSTrieNode* Root, int StartDepth)
{
    SSTrieNode* CurrentNode;
    SSTrieNode* Child;
    int StartPos;
    int Pos;
    int Depth;
    int AA;
    SSTrieLeaf* Leaf;
    SSTrieLeafNode* LeafNode;
    int LinkIndex;
    ExonEdge* Edge;
    //
    CurrentNode = Root;
    Pos = 0;
    Depth = StartDepth;
    //printf("\n--->Start exon scan: Exon len %d, start depth is %d\n", Exon->Length, StartDepth);

    // It's possible that we started at a leaf, if an edge-AA finished us off.  If so,
    // flag the match and return:
    if (StartDepth == 8)
    {
        // We're at a leaf!  Flag this match:
        Leaf = (SSTrieLeaf*)CurrentNode;
        if (!TrieNodeHitFlags[Leaf->Index])
        {
            TrieNodeHitFlags[Leaf->Index] = 1;
            for (LeafNode = Leaf->Head; LeafNode; LeafNode = LeafNode->Next)
            {        
                g_TrieNodeMatches[LeafNode->ProteinNumber]++;
                //printf("At start of exon %d hit word %d\n", Exon->Index, LeafNode->ProteinPos); 
            }
            
        }
        return;
    }
    StartPos = 0;
    while (1)
    {
        if (Pos >= Exon->Length)
        {
            // We've reached the end of the exon.  Follow all outgoing edges:
            for (LinkIndex = 0; LinkIndex < Exon->ForwardEdgeCount; LinkIndex++)
            {
                Edge = Exon->ForwardEdges + LinkIndex;
                AA = Edge->AA;
                if (AA)
                {
                    Child = CurrentNode->Children[AA - 'A'];
                    if (Child)
                    {
                        SSDatabaseScanExon(Edge->Exon, Child, Depth + 1);
                    }
                }
                else
                {
                    SSDatabaseScanExon(Edge->Exon, CurrentNode, Depth);
                }
            }
            // If we were already partway down the trie when we started this exon, return:
            if (StartDepth)
            {
                return;
            }
            else
            {
                // Advance to the next starting point, and jump back to the root:
                StartPos++;
                if (StartPos >= Exon->Length)
                {
                    return; // done with all peptides that begin in this exon!
                }
                CurrentNode = Root;
                Pos = StartPos;
                Depth = 0;
            }
        }
        //printf("%c Pos %d, current node '%s' depth %d=%d\n", Exon->Sequence[Pos], Pos, CurrentNode->Buffer, CurrentNode->Depth, Depth);
        AA = Exon->Sequence[Pos] - 'A';
        if (CurrentNode->Children[AA])
        {
            CurrentNode = CurrentNode->Children[AA];
            Depth++;
            if (Depth == 8)
            {
                // We're at a leaf!  Flag this match:
                Leaf = (SSTrieLeaf*)CurrentNode;
                if (!TrieNodeHitFlags[Leaf->Index])
                {
                    TrieNodeHitFlags[Leaf->Index] = 1;
                    for (LeafNode = Leaf->Head; LeafNode; LeafNode = LeafNode->Next)
                    {        
                        g_TrieNodeMatches[LeafNode->ProteinNumber]++;
                        //printf("At position %d in exon %d hit word %d\n", Pos, Exon->Index, LeafNode->ProteinPos); 
                    }
                    
                }
                // If our starting depth is >0, then don't go to a failure node,
                // just return.  (We only use failure nodes when doing "normal"
                // linear parsing along an exon)
                if (StartDepth)
                {
                    return;
                }
                //Depth = Leaf->FailureDepth;
                //CurrentNode = Leaf->FailureNode;
                StartPos++;
                Depth = 0;
                CurrentNode = Root;
                Pos = StartPos;
                continue;
            }
            Pos++;
        }
        else
        {
            // Our match has ended.  Stop now, or use a failure node:
            if (StartDepth)
            {
                return;
            }
            // If we are in the root now, then we should just advance by
            // one character:
            if (!Depth)
            {
                Pos++;
                StartPos++;
            }
            else
            {
                //// We're not in the root, so we can use a failure node.
                //// Pos stays where it is.
                //Depth = CurrentNode->FailureDepth;
                //CurrentNode = CurrentNode->FailureNode;
                StartPos++;
                CurrentNode = Root;
                Pos = StartPos;
                Depth = 0;
            }
        }
    }
}

typedef struct SSMatchInfo
{
    int Coverage;
    int RecordNumber;
    int ChromosomeNumber;
    int Strand;
    int ApproximatePosition;
    int CoverageStart;
    int CoverageEnd;
} SSMatchInfo;

// From high to low coverage
int CompareSSMatchInfo(const SSMatchInfo* a, const SSMatchInfo* b)
{
    if (a->Coverage > b->Coverage)
    {
        return -1;
    }
    if (a->Coverage < b->Coverage)
    {
        return 1;
    }
    return 0;
}

// Keep the top n matches for each protein, where n is this number:
#define MATCHES_PER_PROTEIN 5
#define LAST_MATCH_FOR_PROTEIN 4

// Main function: Given an array of protein Sequences (with names in NameBuffer) and the splicedb file name,
// scan through the genes in the splice-db to find words from the proteins
void SSDatabaseScanProteins(int FirstRecordNumber, char** Sequences, char* NameBuffer, char* SpliceDBFileName, 
    int BlockSize, FILE* OutputFile)
{
    int ProteinNumber;
    int MaxSequenceLength = 0;
    SSTrieNode* Root;
    SSMatchInfo* AllMatchInfo;
    SSMatchInfo* MatchInfo;
    int SequenceLengths[SS_BLOCK_SIZE];
    int RecordNumber = 0;
    char* MatchFlags;
    FILE* SpliceDBFile;
    GeneStruct* CurrentGene;
    ExonStruct* Exon;
    int EdgeIndex;
    int AA;
    int ExonIndex;
    int Len;
    int Coverage;
    SSTrieNode* Node;
    int GeneNumber;
    ExonEdge* Edge;
    int MatchIndex;
    //
    AllMatchInfo = (SSMatchInfo*)calloc(SS_BLOCK_SIZE * MATCHES_PER_PROTEIN, sizeof(SSMatchInfo));
    for (ProteinNumber = 0; ProteinNumber < BlockSize; ProteinNumber++)
    {
        SequenceLengths[ProteinNumber] = strlen(Sequences[ProteinNumber]);
        MaxSequenceLength = max(MaxSequenceLength, SequenceLengths[ProteinNumber]);
    }
    MatchFlags = (char*)calloc(MaxSequenceLength, sizeof(char));
    Root = ConstructSSTrie(Sequences, BlockSize);
    SetSSTrieFailureNodes(Root);
    //DebugPrintSSTrie(Root, 0, DebugBuffer);
    SpliceDBFile = fopen(SpliceDBFileName, "rb");
    if (!SpliceDBFile)
    {
        REPORT_ERROR_S(8, SpliceDBFileName);
        return;
    }
    memset(TrieNodeHitFlags, 0, sizeof(int) * MAX_TRIE_NODE_COUNT);
    GeneNumber = 0;
    while (1)
    {
        GeneNumber++;
        if (GeneNumber%100 == 0)
        {
            printf("%d ", GeneNumber);
        }
        //ResetSSTrieFlags(Root, 0);
        memset(TrieNodeHitFlags, 0, sizeof(int) * g_NextTrieLeafIndex);
        memset(g_TrieNodeMatches, 0, sizeof(int) * SS_BLOCK_SIZE);
        CurrentGene = LoadGene(SpliceDBFile);
        if (!CurrentGene)
        {
            break;
        }
        // Iterate over exons:
        for (ExonIndex = 0; ExonIndex < CurrentGene->ExonCount; ExonIndex++)
        {
            Exon = CurrentGene->Exons + ExonIndex;
            Len = Exon->Length;

            // Try starting a match with the incoming edge:
            for (EdgeIndex = 0; EdgeIndex < Exon->BackEdgeCount; EdgeIndex++)
            {
                Edge = Exon->BackwardEdges + EdgeIndex;
                AA = Edge->AA - 'A';
                if (AA >= 0 && AA < 26)
                {
                    Node = Root->Children[AA];
                    //SSDatabaseScanHelper(Exon, Len, 0, Node, 1);
                    if (Node)
                    {
                        SSDatabaseScanExon(Exon, Node, 1);
                    }
                }
            }
            SSDatabaseScanExon(Exon, Root, 0);
            //for (Pos = 0; Pos < Len; Pos++)
            //{
            //    SSDatabaseScanHelper(Exon, Len, Pos, Root, 0);
            //}
        }
        // Rate the quality of the match, saving it if it's good:
        for (ProteinNumber = 0; ProteinNumber < BlockSize; ProteinNumber++)
        {
            Coverage = g_TrieNodeMatches[ProteinNumber];
            if (Coverage > SequenceLengths[ProteinNumber] * 0.1)
            {
                //printf("\nProtein #%d: %s\n", ProteinNumber, NameBuffer + (ProteinNumber * PROTEIN_NAME_BLOCK));
                //DebugPrintGene(CurrentGene);
                memset(MatchFlags, 0, sizeof(char) * MaxSequenceLength);
                Coverage = SSTrieCoverSequence(Root, MatchFlags, 0, ProteinNumber);
                //for (AA = 0; AA < SequenceLengths[ProteinNumber]; AA++)
                //{
                //    printf("%d\t%c\t%d\t\n", AA, Sequences[ProteinNumber][AA], MatchFlags[AA]);
                //}
                // If this coverage is better than the lowest-saved-coverage for the protein,
                // then replace the lowest-saved-coverage and sort the list:
                if (Coverage > AllMatchInfo[ProteinNumber * MATCHES_PER_PROTEIN + LAST_MATCH_FOR_PROTEIN].Coverage)
                { 
                    MatchInfo = AllMatchInfo + ProteinNumber * MATCHES_PER_PROTEIN + LAST_MATCH_FOR_PROTEIN;
                    MatchInfo->Coverage = Coverage;
                    MatchInfo->RecordNumber = RecordNumber;
                    MatchInfo->ChromosomeNumber = CurrentGene->ChromosomeNumber;
                    MatchInfo->Strand = CurrentGene->ForwardFlag;
                    MatchInfo->ApproximatePosition = CurrentGene->Exons[0].Start;

                    for (AA = 0; AA < SequenceLengths[ProteinNumber]; AA++)
                    {
                        if (MatchFlags[AA])
                        {
                            MatchInfo->CoverageStart = AA;
                            //BestStartPos[ProteinNumber] = AA;
                            break;
                        }
                    }
                    for (AA = SequenceLengths[ProteinNumber] - 1; AA >= 0; AA--)
                    {
                        if (MatchFlags[AA])
                        {
                            MatchInfo->CoverageEnd = AA;
                            //BestEndPos[ProteinNumber] = AA;
                            break;
                        }
                    }
                    qsort(AllMatchInfo + ProteinNumber * MATCHES_PER_PROTEIN, MATCHES_PER_PROTEIN, sizeof(SSMatchInfo), (QSortCompare)CompareSSMatchInfo);
                }
            }
        }
        RecordNumber += 1;
        FreeGene(CurrentGene);
        // If we've got 95% of the protein, then stop now - we probably won't 
        // get any more!
        //if (BestCoverage > 0.95*SequenceLength)
        //{
        //    break;
        //}
    }
    // Print the match:
    for (ProteinNumber = 0; ProteinNumber < BlockSize; ProteinNumber++)
    {
        fprintf(OutputFile, "%d\t", FirstRecordNumber + ProteinNumber);
        fprintf(OutputFile, "%s\t", NameBuffer + (ProteinNumber * PROTEIN_NAME_BLOCK));
        fprintf(OutputFile, "%d\t", SequenceLengths[ProteinNumber]);
        for (MatchIndex = 0; MatchIndex < MATCHES_PER_PROTEIN; MatchIndex++)
        {
            MatchInfo = AllMatchInfo + ProteinNumber * MATCHES_PER_PROTEIN + MatchIndex;
            if (MatchInfo->Coverage)
            {
                fprintf(OutputFile, "%d\t", MatchInfo->ChromosomeNumber);
                fprintf(OutputFile, "%d\t", MatchInfo->Strand);
                fprintf(OutputFile, "%d\t", MatchInfo->ApproximatePosition);
                fprintf(OutputFile, "%d\t", MatchInfo->Coverage);
            }
        }
        fprintf(OutputFile, "\n");
    }
    FreeSSTrieNode(Root, 0);

    // Cleanup:
    SafeFree(MatchFlags);
    SafeFree(AllMatchInfo);
}

typedef struct SSHashNode
{
    char TrueSequence[8];
    int ProteinIndex;
    int ProteinPos;
    int MatchFlag;
    struct SSHashNode* Next;
} SSHashNode;

#define SS_HASH_MAX 5000000
// Big hash:
SSHashNode* SSHash[SS_HASH_MAX];

void ClearSSHash()
{
    int HashIndex;
    SSHashNode* Node;
    SSHashNode* Prev;
    // 
    for (HashIndex = 0; HashIndex <  SS_HASH_MAX; HashIndex++)
    {
        Prev = NULL;
        Node = SSHash[HashIndex];
        while (Node)
        {
            SafeFree(Prev);
            Prev = Node;
            Node = Node->Next;
        }
        SafeFree(Prev);
        SSHash[HashIndex] = NULL;
    }
}

#define HASH_SEQUENCE(Buffer)\
HashValue = 0;\
for (X = 0; X < 8; X++)\
{\
    HashValue += Buffer[X] * X * X;\
    HashValue %= SS_HASH_MAX;\
}\
    
// Hashing *may* be faster than trie; it hasn't been implemented yet.
void PopulateSSHash(char** SequenceBuffer, int BlockSize)
{
    int ProteinIndex;
    int Pos;
    int Len;
    //
    for (ProteinIndex = 0; ProteinIndex < BlockSize; ProteinIndex++)
    {
        Len = strlen(SequenceBuffer[ProteinIndex]);
        for (Pos = 0; Pos < Len - 7; Pos++)
        {
        }
    }
    
}

// For more rapid scanning of proteins...let's use a hash instead of a trie.  
void SSQDatabaseScanProteins(char** SequenceBuffer, char* NameBuffer, char* SpliceDBFileName, int BlockSize)
{
    ClearSSHash();
    PopulateSSHash(SequenceBuffer, BlockSize);
}

// Main method:
void SSDatabaseScan(char* TrieFileName, char* IndexFileName, char* SpliceDBFileName,
    int FirstRecord, int LastRecord)
{
    //GeneStruct* CurrentGene;
    //GeneStruct* LoadGene(FILE* File)
    int DummyInt;
    int LastFilePos = -1;
    int FilePos;
    FILE* TrieFile;
    FILE* IndexFile;
    char* SequenceBuffer[SS_BLOCK_SIZE];
    char NameBuffer[PROTEIN_NAME_BLOCK * SS_BLOCK_SIZE];
    int BytesRead;
    int RecordLength;
    int BlockIndex = 0;
    FILE* OutputFile;
    int RecordNumber;
    int BlockFirstRecordNumber = 0;
    //
    TrieFile = fopen(TrieFileName, "rb");
    if (!TrieFile)
    {
        REPORT_ERROR_S(8, TrieFileName);
        return;
    }
    IndexFile = fopen(IndexFileName, "rb");
    if (!IndexFile)
    {
        REPORT_ERROR_S(8, IndexFileName);
        return;
    }
    OutputFile = fopen("SSDatabaseScan.txt", "wb");
    if (!OutputFile)
    {
        printf("** Error: Failed to open SSDatabaseScan.txt\n");
        return;
    }
    // Header:
    fprintf(OutputFile, "RecordNumber\tProtein\tLength\tChromosome\tForwardFlag\tApproxPos\tCoverage\t\n");
    // Read protein records from the trie database.  Once you accumulate a block
    // of them in the trie, launch a scan through the exon graph with SSDatabaseScanProteins.
    RecordNumber = 0;
    while (1)
    {
        BytesRead = ReadBinary(&DummyInt, sizeof(int), 1, IndexFile);
        if (!BytesRead)
        {
            // End of file.  Scan our last block, if we have anything in the block:
            if (LastFilePos >= 0 && BlockIndex)
            {
                fseek(TrieFile, FilePos, SEEK_SET);
                SequenceBuffer[BlockIndex - 1] = (char*)calloc(30000, sizeof(char));
                ReadBinary(SequenceBuffer[BlockIndex - 1], sizeof(char), 30000, TrieFile);
                SSDatabaseScanProteins(BlockFirstRecordNumber, SequenceBuffer, NameBuffer, SpliceDBFileName, BlockIndex, OutputFile);
            }
            break;
        }
        BytesRead = ReadBinary(&DummyInt, sizeof(int), 1, IndexFile);
        BytesRead = ReadBinary(&FilePos, sizeof(int), 1, IndexFile);
        // 
        if (LastFilePos >= 0 && RecordNumber >= FirstRecord)
        {
            RecordLength = FilePos - LastFilePos - 1;
            SequenceBuffer[BlockIndex - 1] = (char*)calloc(RecordLength + 1, sizeof(char));
            fseek(TrieFile, LastFilePos, SEEK_SET);
            ReadBinary(SequenceBuffer[BlockIndex - 1], sizeof(char), RecordLength, TrieFile);
            if (BlockIndex == SS_BLOCK_SIZE || (LastRecord >= 0 && RecordNumber >= LastRecord))
            {
                SSDatabaseScanProteins(BlockFirstRecordNumber, SequenceBuffer, NameBuffer, SpliceDBFileName, BlockIndex, OutputFile);
                for (BlockIndex = 0; BlockIndex < BlockIndex; BlockIndex++)
                {
                    SafeFree(SequenceBuffer[BlockIndex]);
                }
                BlockIndex = 0;
                // If we hit the last record, then stop now.
                if (LastRecord >= 0 && RecordNumber >= LastRecord)
                {
                    break;
                }
            }
        }
        LastFilePos = FilePos;
        ReadBinary(NameBuffer + BlockIndex*PROTEIN_NAME_BLOCK, sizeof(char), 80, IndexFile);
        NameBuffer[BlockIndex*PROTEIN_NAME_BLOCK + 80] = '\0';
        if (RecordNumber >= FirstRecord)
        {
            if (BlockIndex == 0)
            {
                BlockFirstRecordNumber = RecordNumber;
            }
            BlockIndex++;
        }
        RecordNumber++;
    }
    fclose(IndexFile);
    fclose(TrieFile);
}
