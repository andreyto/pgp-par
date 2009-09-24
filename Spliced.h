
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

#ifndef SPLICED_H
#define SPLICED_H

#include <stdio.h>
#include "Utils.h"
#include "Trie.h"

#define GENE_NAME_LENGTH 256
// Maximum length (in AA) of an extension:
#define MAX_EXTENSION_LENGTH 64
// max length (in chars) of extension with splice chars included
#define MAX_SEXTENSION_LENGTH 192
// No gene can have this many exons or more:
#define MAX_GENE_EXONS 50000 

#define MAX_EXTENSION_EXONS 128

// Splice-aware database search code.
// Here is our basic approach:
// - Construct an exon-only nucleotide database, stored as a graph.  Each node is an exon (in some reading frame).  Edges
// may contain one additional amino acid (to ''glue'' the extra nucleotides at the edges of the exons).  This construction
// is performed offline by the script SplicePrepare.py
// - Using a trie of tags (built in Tagger.c, just as they are for ordinary search), search the graph.  
// Tags (and their extensions) may follow edges between nodes.

typedef struct ExonEdge
{
    char AA; // can be NULL
    int Power; // if zero, this is an adjacent-edge and not a splice junction
    struct ExonStruct* Source; // the source exon (DONOR, if this is a forward edge)
    struct ExonStruct* Exon; // the target exon (ACCEPTOR, if this is a forward edge)

    // We construct a linked list of exon edges when parsing from an XML file.  Then, when
    // the gene is complete, we convert the linked list to an array.  The linked list
    // uses the Next member; the finished array does not.
    struct ExonEdge* Next;
} ExonEdge;

typedef struct ExonStruct
{
    int Start;
    int End;
    int Index;
    char* Sequence;
    int Length; // length of our sequence, in amino acids. 
    char Prefix[3];
    char Suffix[3];
    int BackEdgeCount;
    int ForwardEdgeCount;
    ExonEdge* ForwardEdges;
    ExonEdge* BackwardEdges;
    ExonEdge* BackEdgeHead; // used during XML parse only
    ExonEdge* BackEdgeTail; // used during XML parse only
    int Occurrences;
    struct GeneStruct* Gene;
} ExonStruct;

// GeneStructs can be stored in a doubly-linked list.
typedef struct GeneStruct
{
    char Name[GENE_NAME_LENGTH + 1];
    char SprotName[GENE_NAME_LENGTH + 1];
    int ChromosomeNumber;
    int ForwardFlag;
    int ExonCount; // Size of the Exons arrays
    struct ExonStruct* Exons;
    struct GeneStruct* Next; 
    struct GeneStruct* Prev;
} GeneStruct;

void TestSplicing(); // internal testing junk!

// Main method for Spliced.c: Given a collection of tags (Root) and a binary splicedb (FileName), search
// for matches to the current list of spectra (list head GlobalOptions->FirstSpectrum, but we get to them
// via back-links from tags).  Score matches with Scorer.  If GeneNames is not null, then it's an array
// of gene names to be searched, and we skip any gene whose name isn't on the list.
void SearchSplicableGenes(SearchInfo* Info);
//void SearchSplicableGenes(TrieNode* Root, char* FileName, char** GeneNames, ScoringFunction Scorer,
//    int DBNumber);
GeneStruct* LoadGene(FILE* File);
void FreeGene(GeneStruct* Gene);
void DebugPrintGene(GeneStruct* Gene);
//void SearchSplicableGene(TrieNode* Root, GeneStruct* Gene, ScoringFunction Scorer, int DBNumber);
void SearchSplicableGene(SearchInfo* Info, GeneStruct* Gene);
void AllocSpliceStructures();
void SetExonForwardEdges(GeneStruct* Gene);

#endif // SPLICED_H
