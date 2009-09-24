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

#ifndef SPLICEDB_H
#define SPLICEDB_H

#include <stdio.h>
#include "Utils.h"
#include "Trie.h"

typedef struct GenomeDAGLink
{
    struct GenomeDAGNode* Node;
    struct GenomeDAGLink* Next;
    int Count;
} GenomeDAGLink;

#define MAX_DAG_NODE_LINKS 3

typedef struct GenomeDAGNode
{
    int Start;
    int End;
    char* Sequence;
    GenomeDAGLink* FirstForward;
    GenomeDAGLink* FirstBack;
    //GenomeDAGNode** Next[MAX_DAG_NODE_LINKS];
    //GenomeDAGNode** Prev[MAX_DAG_NODE_LINKS];
    struct ExonNode** Exons;
    //GenomeDAGLink* FirstForward;
    //GenomeDAGLink* FirstBack;
} GenomeDAGNode;

typedef struct IntervalNode
{
    int Start;
    int End; // exclusive
    int Occurrences;
    int Satisfied;
    unsigned int OriginalFilePos;
    struct EdgeNode* FirstForward;
    struct EdgeNode* LastForward;
    struct EdgeNode* FirstBack;
    struct EdgeNode* LastBack;
    struct IntervalNode* Prev;
    struct IntervalNode* Next;
    struct ExonNode* FirstExon;
    struct ExonNode* LastExon;
    struct GeneNode* GNode; // non-null while this interval is in a pending gene.
    int DAGNodeCount;
    GenomeDAGNode* DAGNodes;
    int Flags; // for keeping track of which reading frames we permit!
} IntervalNode;

typedef struct GeneNode
{
    IntervalNode* Interval;
    // RX is the minimum covered length of any path originating at this interval and extending forward.
    // We set RX during the 'satisfaction' process, so that we can note the (partial) satisfaction
    // of intervals other than the seed.
    // The 'covered' section of a path is the portion that's already part of the gene.
    // RX is initialized to 0.  
    // During satisfaction procedure: 
    // If RX is big enough already, return without recursing.  Otherwise:
    // RX is set to 9999 if we have no forward edges.
    // Otherwise, RX is set to the minimum value of the outgoing edge's interval's length (plus return value of the recursive 
    // satisfaction call, if any).
    int RX;
    int LX;
    struct GeneNode* Prev;
    struct GeneNode* Next;
} GeneNode;

// ExonNode is used while constructing the database.
// In production, use ExonStruct from Spliced.h instead.
typedef struct ExonNode
{
    IntervalNode* Interval;
    struct ExonLink* FirstForward;
    struct ExonLink* LastForward;
    struct ExonLink* FirstBack;
    struct ExonLink* LastBack;
    struct ExonNode* Next;
    char Prefix[3]; // one or two characters, and null-terminator
    char Suffix[3];
    int Index;
    int Start; // start (in genomic coordinates, usually same as the parent interval's start)
    int End; //end (in genomic coordinates)
    int Length; // length in amino acids (not in genomic coordinates)
    char* Sequence;
    GenomeDAGNode* DAGNode;
    int MaxForwardOverall;
    int MaxBackOverall;

} ExonNode;

typedef struct ExonLink
{
    char AA;
    ExonNode* Exon;
    int Power;
    // maximum peptide length achievable with this amino acid (if any) and the next exon,
    // until stop codon or edge of graph.
    int MaxLength; 
    struct ExonLink* Next;
} ExonLink;

typedef struct EdgeNode
{
    int Count;
    float Score;
    IntervalNode* Interval;
    struct EdgeNode* Prev;
    struct EdgeNode* Next;
} EdgeNode;

void PrepareSpliceDB(int ChromosomeNumber, int ReverseFlag, int MinORFLength);
void PrepareOneGeneSpliceDB(int ChromosomeNumber, int ReverseFlag, int IntervalStart, int IntervalEnd,
    char* CustomFileName, char* GeneName, int MinORFLength);
void TestSpliceDB(int argc, char** argv);
#endif // SPLICEDB_H
