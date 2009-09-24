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


// BuildMS2DB.c is responsible for building a splice-tolerant database
// from an input file in .gff format.  Our overall plan is as follows:

//- Parse all lines of the .gff file.  Discard those lines which don't come from
//  the correct chromosome and strand.  Build up a linked list of exons, 
//  which is indexed by a hash based on (start, end, reading-frame).  
//  NOTE: Reading frame is defined as "the first base pair in a codon, modulo 3".  
//  NOTE: Records of type EST give rise to three exons (one for each reading frame),
//        which will be pruned later.
//  NOTE: We may parse SEVERAL gff files!
//- Merge and split exons, to produce a minimal disjoint list.  
//- Prune exons with short ORFs.
//- ITERATE, until all exons are covered:
//  - Take the first uncovered exon.  Grab all the exons which it links to.  Build 
//    a corresponding gene record.  Flag these exons as covered.
//- Write cross-reference records for the GFFGenes.
//  
//  IMPORTANT NOTE: The data structures built in this file aren't used during a search.
//  When it comes to searching, look in MS2DB.c, which uses different (simpler) data
//  structures to track this stuff.

#include "Utils.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "Trie.h"
#include "Inspect.h"
#include "Spliced.h"
#include "SpliceDB.h"
#include "SNP.h"
#include "Errors.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Macros:
#define MAX_NAME 256
#define EXON_HASH_SIZE 10000

#define MS2_EXON_NONE 0
#define MS2_EXON_CUSTOMAA 1

#define CODON_LENGTH 3

//#define GFF_QUICKPARSE

// Use this flag for a single-base exon which covers the middle
// base of an exon, and which has only one link forward to an
// adjacent exon.  
// The flag reflects this scenario (where B is the one-base exon):
//  AAAA---(X)---BCCCC
// Rather than these scenarios:
//  AAAAB---(X)---CCC
//  AAAA---B---(X)---CCC
#define MS2_EXON_CUSTOMAA_HEAD 2

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Structs:
typedef struct MS2Exon
{
    int Start;
    int End;
    int ReadingFrame;
    struct MS2Exon* Next;
    struct MS2Exon* Prev;
    char Prefix[CODON_LENGTH];
    char Suffix[CODON_LENGTH];
    struct MS2Edge* FirstForward;
    struct MS2Edge* LastForward;
    struct MS2Edge* FirstBackward;
    struct MS2Edge* LastBackward;
    struct MS2Gene* Gene;
    int Index; // Assigned to a completed gene
    char* Sequence;
    // Flags for special bookkeeping.  "bookkeeping" and its derivatives are perhaps the 
    // only English words with three consecutive double-letters.
    int Flags; 
    // CustomAA is set only for length-1 exons that are part of special bridges.
    char CustomAA;
} MS2Exon;

typedef struct MS2Edge
{
    struct MS2Exon* LinkTo;
    struct MS2Exon* LinkFrom;
    struct MS2Edge* Next;
    // The amino acid for this edge can be set specially, up-front:
    char SpecialAA;
} MS2Edge;

// A wrapper for an exon.  This lets us stick an exon into multiple linked lists.
typedef struct MS2ExonNode
{
    struct MS2ExonNode* Next;
    struct MS2ExonNode* Prev;
    MS2Exon* Exon;
} MS2ExonNode;

typedef struct IntNode
{
    struct IntNode* Next;
    int Value;
} IntNode;

typedef struct MS2CrossReference
{
    struct MS2Gene* Gene;
    struct GFFGeneClass* GFFGene;
    IntNode* FirstExonID;
    IntNode* LastExonID;
    struct MS2CrossReference* Next;
} MS2CrossReference;

typedef struct MS2Gene
{
    MS2ExonNode* FirstExon;
    MS2ExonNode* LastExon;
    struct MS2Gene* Next;
    int Index; // for debugging, mostly!
    MS2CrossReference* FirstCrossReference;
    MS2CrossReference* LastCrossReference;
} MS2Gene;

// Singleton class tracking high-level data like the exon hashes.
typedef struct MS2Builder
{
    FILE* GenomeFile;
    int ForwardFlag;
    MS2ExonNode** ExonHash;
    MS2Exon* FirstExon;
    MS2Exon* LastExon;
    struct GFFGeneClass* FirstGFFGene;
    struct GFFGeneClass* LastGFFGene;
    MS2Gene* FirstGene;
    MS2Gene* LastGene;
    char ChromosomeName[MAX_NAME];
    int ExonCount;
    int GeneCount;
    int VerboseFlag;
} MS2Builder;

typedef struct GFFExonClass
{
    int Start;
    int End;
    int ReadingFrame;
    struct GFFExonClass* Next;
} GFFExonClass;

typedef struct GFFGeneClass
{
    char Name[MAX_NAME + 1];
    char DatabaseName[MAX_NAME + 1];
    GFFExonClass* FirstExon;
    GFFExonClass* LastExon;
    struct GFFGeneClass* Next;
    struct MS2CrossReference* CrossReference;
} GFFGeneClass;

typedef struct GFFParser
{
    // Link to our builder, where the REAL data (not transient parse-state) lives:
    MS2Builder* Builder;
    // Keep a link to the current gene, so we can add exons to it:
    struct GFFGeneClass* CurrentGene;
    // Remember our filename (mostly for error reporting)
    char* CurrentFileName;
    // Keep a link to the last exon of the current gene, so that we 
    // can add edges between exons as needed:
    MS2Exon* PrevExon;
} GFFParser;

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Forward declarations:
void DebugPrintMS2Builder(MS2Builder* Builder, char* Notes);
void ExonInheritOneForwardEdge(MS2Exon* Exon, MS2Edge* OldEdge);
void ExonInheritOneBackwardEdge(MS2Exon* Exon, MS2Edge* OldEdge);
void FreeMS2CrossReference(MS2CrossReference* CR);

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor functions:
MS2Exon* NewExon(int Start, int End, int ReadingFrame)
{
    MS2Exon* Exon;
    Exon = (MS2Exon*)calloc(1, sizeof(MS2Exon));
    Exon->Start = Start;
    Exon->End = End;
    Exon->ReadingFrame = ReadingFrame;
    return Exon;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Destructor functions:

void FreeGFFGene(GFFGeneClass* Gene)
{
    GFFExonClass* Exon;
    GFFExonClass* Prev;
    if (!Gene)
    {
        return;
    }
    FreeMS2CrossReference(Gene->CrossReference);
    Prev = NULL;
    for (Exon = Gene->FirstExon; Exon; Exon = Exon->Next)
    {
        SafeFree(Prev);
        Prev = Exon;
    }
    SafeFree(Prev);
    SafeFree(Gene);
}

// Free an MS2Exon, and its associated edges.
void FreeMS2Exon(MS2Exon* Exon)
{
    MS2Edge* Edge;
    MS2Edge* Prev;
    //
    if (!Exon)
    {
        return;
    }
    // Free the SEQUENE:
    SafeFree(Exon->Sequence);
    Exon->Sequence = NULL;
    // Free the list of FORWARD edges:
    Prev = NULL;
    for (Edge = Exon->FirstForward; Edge; Edge = Edge->Next)
    {
        SafeFree(Prev);
        Prev = Edge;
    }
    SafeFree(Prev);
    // Free the list of BACKWARD edges:
    Prev = NULL;
    for (Edge = Exon->FirstBackward; Edge; Edge = Edge->Next)
    {
        SafeFree(Prev);
        Prev = Edge;
    }
    SafeFree(Prev);
    SafeFree(Exon);
}

void FreeExonHash(MS2Builder* Builder)
{
    int HashIndex;
    MS2ExonNode* Node;
    MS2ExonNode* Prev;
    //
    if (!Builder->ExonHash)
    {
        return;
    }
    for (HashIndex = 0; HashIndex < EXON_HASH_SIZE; HashIndex++)
    {
        Prev = NULL;
        for (Node = Builder->ExonHash[HashIndex]; Node; Node = Node->Next)
        {
            SafeFree(Prev);
            Prev = Node;
        }
        SafeFree(Prev);
    }
    SafeFree(Builder->ExonHash);
    Builder->ExonHash = NULL;
}

// Free an MS2CrossReference, and its associated list of integers.
void FreeMS2CrossReference(MS2CrossReference* CR)
{
    IntNode* Node;
    IntNode* Prev;
    if (!CR)
    {
        return;
    }
    Prev = NULL;
    for (Node = CR->FirstExonID; Node; Node = Node->Next)
    {
        SafeFree(Prev);
        Prev = Node;
    }
    SafeFree(Prev);
    SafeFree(CR);
}

// Free an MS2Gene, and its associated list of MS2ExonNodes
void FreeMS2Gene(MS2Gene* Gene)
{
    MS2ExonNode* Node;
    MS2ExonNode* Prev = NULL; 
    //
    if (!Gene)
    {
        return;
    }
    for (Node = Gene->FirstExon; Node; Node = Node->Next)
    {
        SafeFree(Prev);
        Prev = Node;
    }
    SafeFree(Prev);
}

void FreeGFFGenes(MS2Builder* Builder)
{
    GFFGeneClass* Gene;
    GFFGeneClass* Prev;
    //
    Prev = NULL;
    for (Gene = Builder->FirstGFFGene; Gene; Gene = Gene->Next)
    {
        FreeGFFGene(Prev);
        Prev = Gene;
    }
    FreeGFFGene(Prev);
    Builder->FirstGFFGene = NULL;
    Builder->LastGFFGene = NULL;
}

void FreeMS2Genes(MS2Builder* Builder)
{
    MS2Gene* Gene;
    MS2Gene* Prev = NULL;
    for (Gene = Builder->FirstGene; Gene; Gene = Gene->Next)
    {
        FreeMS2Gene(Prev);
        Prev = Gene;
    }
    FreeMS2Gene(Prev);
    Builder->FirstGene = NULL;
    Builder->LastGene = NULL;
}

void FreeMS2Exons(MS2Builder* Builder)
{
    MS2Exon* Exon;
    MS2Exon* Prev = NULL;
    //
    for (Exon = Builder->FirstExon; Exon; Exon = Exon->Next)
    {
        FreeMS2Exon(Prev);
        Prev = Exon;
    }
    FreeMS2Exon(Prev);
    Builder->FirstExon = NULL;
    Builder->LastExon = NULL;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Other functions:

void MS2CrossReferenceAddID(MS2CrossReference* CR, int ID)
{
    IntNode* Node;
    //
    Node = (IntNode*)calloc(1, sizeof(IntNode));
    Node->Value = ID;
    if (CR->LastExonID)
    {
        CR->LastExonID->Next = Node;
    }
    else
    {
        CR->FirstExonID = Node;
    }
    CR->LastExonID = Node;
}

// Return true if we successfully cover everything from start...end
int BuildGFFCrossReference(MS2Builder* Builder, MS2CrossReference* CR, 
                           MS2Gene* Gene, int Start, int End, int ReadingFrame)
{
    MS2Exon* Exon;
    MS2ExonNode* Node;
    int Result = 0;
    //
    for (Node = Gene->FirstExon; Node; Node = Node->Next)
    {
        Exon = Node->Exon;
        // Ignore exons with wrong reading frame:
        if (Exon->ReadingFrame != ReadingFrame)
        {
            continue;
        }
        // Handle perfect-overlap:
        if (Exon->Start == Start && Exon->End == End)
        {
            MS2CrossReferenceAddID(CR, Exon->Index);
            return 1;
        }
        // Handle partial-overlap:
        if (Builder->ForwardFlag && Exon->Start == Start && Exon->End <= End)
        {
            MS2CrossReferenceAddID(CR, Exon->Index);
            return BuildGFFCrossReference(Builder, CR, Gene, Exon->End, End, ReadingFrame);
        }
        if (!Builder->ForwardFlag && Exon->End == End && Exon->Start >= Start)
        {
            MS2CrossReferenceAddID(CR, Exon->Index);
            return BuildGFFCrossReference(Builder, CR, Gene, Start, Exon->Start, ReadingFrame);
        }
    }
    return Result;
}

void BuildGFFCrossReferences(MS2Builder* Builder)
{
    // Iterate over GFFGenes.  Inner loop over MS2Genes.  When you find an MS2Gene which covers
    // any of the GFFGene, build a cross-reference.
    GFFGeneClass* GFFGene;
    GFFExonClass* GFFExon;
    MS2Gene* Gene;
    MS2Exon* Exon;
    MS2ExonNode* Node;
    int OverlapFlag;
    int CoveredFlag;
    int CoverageComplete;
    int Result;
    //
    for (GFFGene = Builder->FirstGFFGene; GFFGene; GFFGene = GFFGene->Next)
    {
        CoveredFlag = 0;
        for (Gene = Builder->FirstGene; Gene; Gene = Gene->Next)
        {
            // First, check whether the gene overlaps this GFF gene's first exon:
            OverlapFlag = 0;
            GFFExon = GFFGene->FirstExon;
            for (Node = Gene->FirstExon; Node; Node = Node->Next)
            {
                Exon = Node->Exon;
                if (Exon->Start < GFFExon->End && Exon->End > GFFExon->Start && Exon->ReadingFrame == GFFExon->ReadingFrame)
                {
                    OverlapFlag = 1;
                }
            }
            if (!OverlapFlag)
            {
                continue;
            }
            // They overlap - create a cross-reference and start adding exon IDs to the list!
            GFFGene->CrossReference = (MS2CrossReference*)calloc(1, sizeof(MS2CrossReference));
            GFFGene->CrossReference->GFFGene = GFFGene;
            GFFGene->CrossReference->Gene = Gene;
            // the MS2Gene keeps a list of its cross references:
            if (Gene->LastCrossReference)
            {
                Gene->LastCrossReference->Next = GFFGene->CrossReference;
            }
            else
            {
                Gene->FirstCrossReference = GFFGene->CrossReference;
            }
            Gene->LastCrossReference = GFFGene->CrossReference;
            CoverageComplete = 1;
            for (GFFExon = GFFGene->FirstExon; GFFExon; GFFExon = GFFExon->Next)
            {
                Result = BuildGFFCrossReference(Builder, GFFGene->CrossReference, Gene, GFFExon->Start, GFFExon->End, GFFExon->ReadingFrame);
                if (!Result)
                {
                    CoverageComplete = 0;
                }
            }
            if (!CoverageComplete)
            {
                REPORT_ERROR_S(23, GFFGene->Name);
            }
            CoveredFlag = 1;
            break;
        }
        // Sanity check: The GFF gene MUST be covered (since we did, after all, create exons for all these 
        // GFF exons!)
        if (!CoveredFlag)
        {
            REPORT_ERROR_S(22, GFFGene->Name);
        }
    }
}

void AddMS2Exon(MS2Builder* Builder, MS2Exon* Exon)
{
    if (!Builder->FirstExon)
    {
        Builder->FirstExon = Exon;
    }
    else
    {
        Builder->LastExon->Next = Exon;
        Exon->Prev = Builder->LastExon;
    }
    Builder->LastExon = Exon;
    Builder->ExonCount++;
}

void LinkExonForward(MS2Exon* FromExon, MS2Exon* ToExon)
{
    // Add a link forward from FromExon to ToExon, as well as a reciprocal link back.
    MS2Edge* Edge;

    // Sanity checking: ToExon starts after FromExon
    INSPECT_ASSERT(ToExon->Start >= FromExon->End);
    
    // If the two exons are ALREADY linked, then return immediately:
    for (Edge = FromExon->FirstForward; Edge; Edge = Edge->Next)
    {
        if (Edge->LinkTo == ToExon)
        {
            return;
        }
    }
    // Sanity check: There's no forward link, therefore there must be no reciprocal link:
    for (Edge = ToExon->FirstBackward; Edge; Edge = Edge->Next)
    {
        INSPECT_ASSERT(Edge->LinkTo != FromExon);
    }
    // Add an edge linking FORWARD:
    Edge = (MS2Edge*)calloc(1, sizeof(MS2Edge));
    Edge->LinkFrom = FromExon;
    Edge->LinkTo = ToExon;
    if (!FromExon->FirstForward)
    {
        FromExon->FirstForward = Edge;
    }
    else
    {
        FromExon->LastForward->Next = Edge;
    }
    FromExon->LastForward = Edge;

    // Add a reciprocal edge linking BACKWARD:
    Edge = (MS2Edge*)calloc(1, sizeof(MS2Edge));
    Edge->LinkFrom = ToExon;
    Edge->LinkTo = FromExon;
    if (!ToExon->FirstBackward)
    {
        ToExon->FirstBackward = Edge;
    }
    else
    {
        ToExon->LastBackward->Next = Edge;
    }
    ToExon->LastBackward = Edge;
}
void RemoveBackwardEdge(MS2Exon* Exon, MS2Exon* LinkedExon)
{
    MS2Edge* Prev = NULL;
    MS2Edge* Edge;
    for (Edge = Exon->FirstBackward; Edge; Edge = Edge->Next)
    {
        if (Edge->LinkTo == LinkedExon)
        {
            // Remove this edge!
            if (Prev)
            {
                Prev->Next = Edge->Next;
            }
            else
            {
                Exon->FirstBackward = Edge->Next;
            }
            if (Exon->LastBackward == Edge)
            {
                Exon->LastBackward = Prev;
            }
            SafeFree(Edge);
            break;
        }
        Prev = Edge;
    }
}

void RemoveForwardEdge(MS2Exon* Exon, MS2Exon* LinkedExon)
{
    MS2Edge* Prev = NULL;
    MS2Edge* Edge;
    for (Edge = Exon->FirstForward; Edge; Edge = Edge->Next)
    {
        if (Edge->LinkTo == LinkedExon)
        {
            // Remove this edge!
            if (Prev)
            {
                Prev->Next = Edge->Next;
            }
            else
            {
                Exon->FirstForward = Edge->Next;
            }
            if (Exon->LastForward == Edge)
            {
                Exon->LastForward = Prev;
            }
            SafeFree(Edge);
            break;
        }
        Prev = Edge;
    }
}

void ExonInheritForwardEdges(MS2Exon* Exon, MS2Exon* DeadExon)
{
    MS2Edge* Edge;
    MS2Edge* Prev = NULL;
    //
    // Sanity checking: Exon and DeadExon share their right endpoint.
    INSPECT_ASSERT(Exon->End == DeadExon->End);
    for (Edge = DeadExon->FirstForward; Edge; Edge = Edge->Next)
    {
        // Add a link forward from Exon->LinkToExon:
        LinkExonForward(Exon, Edge->LinkTo);
        // Remove reciprocal exon link from LinkTo back to DeadExon:
        RemoveBackwardEdge(Edge->LinkTo, DeadExon);
        SafeFree(Prev);
        Prev = Edge;
    }
    SafeFree(Prev);
    DeadExon->FirstForward = NULL;
    DeadExon->LastForward = NULL;
}

void ExonInheritBackwardEdges(MS2Exon* Exon, MS2Exon* DeadExon)
{
    MS2Edge* Edge;
    MS2Edge* Prev = NULL;
    //
    // Sanity checking: Exon and DeadExon share their right endpoint.
    INSPECT_ASSERT(Exon->Start == DeadExon->Start);
    for (Edge = DeadExon->FirstBackward; Edge; Edge = Edge->Next)
    {
        // Add a link forward from LinkToExon->Exon:
        LinkExonForward(Edge->LinkTo, Exon);
        // Remove reciprocal exon link from LinkTo to DeadExon:
        RemoveForwardEdge(Edge->LinkTo, DeadExon);
        SafeFree(Prev);
        Prev = Edge;
    }
    SafeFree(Prev);
    DeadExon->FirstBackward = NULL;
    DeadExon->LastBackward = NULL;
}

// Given this start, end, and reading frame, look up the corresponding 
// exon in Builder->ExonHash.  If the exon doesn't exist yet, then create it.
// Return the exon.
MS2Exon* HashExon(MS2Builder* Builder, int Start, int End, int ReadingFrame)
{
    int HashValue;
    MS2ExonNode* Node;
    MS2ExonNode* Prev = NULL;
    MS2Exon* Exon;
    //
    HashValue = (Start + CODON_LENGTH * End + ReadingFrame) % EXON_HASH_SIZE;
    for (Node = Builder->ExonHash[HashValue]; Node; Node = Node->Next)
    {
        if (Node->Exon->Start == Start && Node->Exon->End == End && Node->Exon->ReadingFrame == ReadingFrame)
        {
            return (Node->Exon);
        }
        Prev = Node;
    }
    // There's no node for this exon yet.  Add one:
    Exon = (MS2Exon*)calloc(1, sizeof(MS2Exon));
    Exon->Start = Start;
    Exon->End = End;
    Exon->ReadingFrame = ReadingFrame;
    Node = (MS2ExonNode*)calloc(1, sizeof(MS2ExonNode));
    Node->Exon = Exon;
    if (Prev)
    {
        Prev->Next = Node;
        Node->Prev = Prev;
    }
    else
    {
        Builder->ExonHash[HashValue] = Node;
    }
    AddMS2Exon(Builder, Exon);
    return Exon;
}

// Parse one line of a .gff file.  Callback for ParseFileByLines
int HandleGFFFileLine(int LineNumber, int FilePos, char* LineBuffer, void* ParseData)
{
    GFFParser* Parser;
    MS2Builder* Builder;
    GFFGeneClass* GFFGene;
    GFFExonClass* GFFExon;
    char* GeneName;
    char* DatabaseName;
    char* DummyStr;
    char* IntervalType;
    int SyntaxErrorFlag = 0;
    int Start;
    int End;
    int ReadingFrame;
    MS2Exon* Exon;
    Parser = (GFFParser*)ParseData;
    Builder = Parser->Builder;
    // Break the line up by tabs.
    // Bit 0: Gene name
    GeneName = strtok(LineBuffer, "\t");
    // Bit 1: Source (used to populate the Database field)
    DatabaseName = strtok(NULL, "\t");
    if (!DatabaseName)
    {
        SyntaxErrorFlag = 1;
        goto cleanup;
    }
    // Debugging option:
#ifdef GFF_QUICKPARSE
    if (LineNumber > 1000)
    {
        return 0;
    }
#endif
    // Bit 2: interval type (est or exon)
    IntervalType = strtok(NULL, "\t");
    if (!IntervalType)
    {
        SyntaxErrorFlag = 1;
        goto cleanup;
    }
    // Bit 3: start
    DummyStr = strtok(NULL, "\t");
    if (!DummyStr)
    {
        SyntaxErrorFlag = 1;
        goto cleanup;
    }
    Start = atoi(DummyStr) - 1; // fix one-based numbering!

    // Bit 4: end
    DummyStr = strtok(NULL, "\t");
    if (!DummyStr)
    {
        SyntaxErrorFlag = 1;
        goto cleanup;
    }
    End = atoi(DummyStr); 

    // Error checking:
    if (Start < 0 || End < 0 || Start >= End)
    {
        REPORT_ERROR_IIIS(20, Start, End, LineNumber, Parser->CurrentFileName);
    }

    // Bit 5: score (ignored)
    DummyStr = strtok(NULL, "\t");
    if (!DummyStr)
    {
        SyntaxErrorFlag = 1;
        goto cleanup;
    }
    // Bit 6: forward flag
    DummyStr = strtok(NULL, "\t");
    if (!DummyStr || !*DummyStr)
    {
        SyntaxErrorFlag = 1;
        goto cleanup;
    }
    // Skip over this exon, if it comes from the wrong strand:
    if (*DummyStr == '+')
    {
        if (!Builder->ForwardFlag)
        {
            goto cleanup;
        }
    }
    else if (*DummyStr == '-')
    {
        if (Builder->ForwardFlag)
        {
            goto cleanup;
        }
    }
    else
    {
        SyntaxErrorFlag = 1;
        goto cleanup;
    }
    // Bit 6: reading frame
    DummyStr = strtok(NULL, "\t");
    if (!DummyStr)
    {
        SyntaxErrorFlag = 1;
        goto cleanup;
    }
    ReadingFrame = atoi(DummyStr);
    if (ReadingFrame < 0 || ReadingFrame > 2)
    {
        SyntaxErrorFlag = 1;
        goto cleanup;
    }
    /////////////////////////////////////////////////////////////////////////////////
    // We've parsed a valid gff file line.  Create a new GFFGene (if necessary), and 
    // add a new GFFExon to our current GFFGene.
    // Fix up the reading frame.  As always, reading frame is the modulus of the first base 
    // pair of a codon.  In GFF format, reading frame is the number of bases to be skipped over
    // before the first base pair of a codon.
    if (Builder->ForwardFlag)
    {
        ReadingFrame = (Start + ReadingFrame) % CODON_LENGTH;
    }
    else
    {
        ReadingFrame = (End - 1 - ReadingFrame) % CODON_LENGTH;
    }
    // Create a new GFFGene, if necessary:
    if (!Parser->CurrentGene || CompareStrings(GeneName, Parser->CurrentGene->Name))
    {
        GFFGene = (GFFGeneClass*)calloc(1, sizeof(GFFGeneClass));
        strncpy(GFFGene->Name, GeneName, MAX_NAME);
        strncpy(GFFGene->DatabaseName, DatabaseName, MAX_NAME);
        Parser->CurrentGene = GFFGene;
        if (!Builder->FirstGFFGene)
        {
            Builder->FirstGFFGene = GFFGene;
        }
        else
        {
            Builder->LastGFFGene->Next = GFFGene;
        }
        Builder->LastGFFGene = GFFGene;
        Parser->PrevExon = NULL;
    }
    else
    {
        // We're continuing along the same GFFGene, so we can link the previous exon to this one.
        // (ASSUMPTION: Exons for the same gene are linked by introns, and come IN ORDER!)
    }
    // Append a new GFFExon to the current GFFGene:
    GFFExon = (GFFExonClass*)calloc(1, sizeof(GFFExonClass));
    GFFExon->Start = Start;
    GFFExon->End = End;
    GFFExon->ReadingFrame = ReadingFrame;
    if (!Parser->CurrentGene->FirstExon)
    {
        Parser->CurrentGene->FirstExon = GFFExon;
    }
    else
    {
        Parser->CurrentGene->LastExon->Next = GFFExon;
    }
    Parser->CurrentGene->LastExon = GFFExon;

    // Construct an MS2Exon:
    Exon = HashExon(Builder, Start, End, ReadingFrame);
    
    // Add a link, if necessary, between this exon and the previous:
    if (Parser->PrevExon)
    {
        // Report an error if the exons overlap:
        if (Parser->PrevExon->End > Exon->Start && Parser->PrevExon->Start < Exon->End)
        {
            REPORT_ERROR_IIII(17, Parser->PrevExon->Start, Parser->PrevExon->End, Exon->Start, Exon->End);
        }
        else if (Parser->PrevExon->Start < Exon->Start)
        {
            // Exons listed from low genome coords to high.  Typical for forward strand.
            LinkExonForward(Parser->PrevExon, Exon);
        }
        else
        {
            // Exons listed from high genome coords to low.  Typical for reverse strand.
            LinkExonForward(Exon, Parser->PrevExon);
        }
        
    }
    Parser->PrevExon = Exon;

    // Report syntax errors:
cleanup:
    if (SyntaxErrorFlag)
    {
        REPORT_ERROR_IS(14, LineNumber, Parser->CurrentFileName);
        return 0;
    }
    return 1;

}

// Iterate over all GFF files.  Parse each one, using the HandleGFFFileLine callback to do the real work.
void ParseGFFFiles(MS2Builder* Builder)
{
    StringNode* GFFNode;
    FILE* GFFFile;
    GFFParser* Parser;

    for (GFFNode = GlobalOptions->FirstGFFFileName; GFFNode; GFFNode = GFFNode->Next)
    {
        GFFFile = fopen(GFFNode->String, "rb");
        if (!GFFFile)
        {
            REPORT_ERROR_S(8, GFFNode->String);
            continue;
        }
        Parser = (GFFParser*)calloc(1, sizeof(GFFParser));
        Parser->Builder = Builder;
        Parser->CurrentFileName = GFFNode->String;
        ParseFileByLines(GFFFile, HandleGFFFileLine, Parser, 0);
        free(Parser);
        fclose(GFFFile);
    }
}

// Remove an exon from the master linked-list:
void DeleteMS2Exon(MS2Builder* Builder, MS2Exon* Exon)
{
    if (Exon == Builder->FirstExon)
    {
        Builder->FirstExon = Builder->FirstExon->Next;
    }
    if (Exon == Builder->LastExon)
    {
        Builder->LastExon = Builder->LastExon->Prev;
    }
    if (Exon->Next)
    {
        Exon->Next->Prev = Exon->Prev;
    }
    if (Exon->Prev)
    {
        Exon->Prev->Next = Exon->Next;
    }
    FreeMS2Exon(Exon);
    Builder->ExonCount--;
}

// When we're making our double-iteration over the exon linked list, it's important to keep 
// track of where we are.  (Note that ExonA itself may be deleted, so we can't simply finish 
// the loop and then go to ExonA->Next!)  The variable NextExonA is what ExonA should be set to 
// for the next A-loop.  And normally NextExonA is simply equal to ExonA->Next.  However, 
// if we delete exons A and B, *and* exon B happens to be ExonA->Next, then on the next time through
// the loop, ExonA must be shifted two positions forward.  Trust me.
#define DELETE_EXON_B()\
{\
    if (NextExonA == ExonB)\
    {\
        NextExonA = ExonB->Next;\
    }\
    DeleteMS2Exon(Builder, ExonB);\
}

// Iterate over all pairs of exons in the builder.  If the exons overlap, split them!
// We handle overlap in many different special cases, each of which is straightforward.
void SplitOverlappingExons(MS2Builder* Builder)
{
    MS2Exon* Exon1;
    MS2Exon* Exon2;
    MS2Exon* Exon3;
    MS2Exon* ExonA;
    MS2Exon* NextExonA;
    MS2Exon* ExonB;
    MS2Exon* NextExonB;
    int ReadingFrame;
    int OverlapFlag = 0;
    //
    ExonA = Builder->FirstExon;
    while (ExonA)
    {
        if (OverlapFlag)
        {
            //DebugPrintMS2Builder(Builder, "A");
            OverlapFlag = 0;
        }

        //printf("ExonA[%d]: %d-%d\n", ExonA->ReadingFrame, ExonA->Start, ExonA->End);
        ReadingFrame = ExonA->ReadingFrame;
        NextExonA = ExonA->Next;
        // Loop B:
        ExonB = ExonA->Next;
        while (ExonB)
        {
            if (OverlapFlag)
            {
                //DebugPrintMS2Builder(Builder, "B");
                OverlapFlag = 0;
            }
            // Compare exon A to exon B:
            if (ExonA->ReadingFrame != ExonB->ReadingFrame)
            {
                ExonB = ExonB->Next;
                continue;
            }
            if (ExonA->End <= ExonB->Start)
            {
                ExonB = ExonB->Next;
                continue;
            }
            if (ExonA->Start >= ExonB->End)
            {
                ExonB = ExonB->Next;
                continue;
            }
            NextExonB = ExonB->Next;
            printf("  %d-%d Overlaps with %d-%d\n", ExonA->Start, ExonA->End, ExonB->Start, ExonB->End);
            OverlapFlag = 1;
            ////////////////////////////////////////////////////////////////////////////////////
            // There's overlap.  Handle each case in turn.
            if (ExonA->Start == ExonB->Start && ExonA->End == ExonB->End)
            {
                ExonInheritBackwardEdges(ExonA, ExonB);
                ExonInheritForwardEdges(ExonA, ExonB);
                DELETE_EXON_B();
                ExonB = NextExonB;
                continue;
            }
            if (ExonA->Start == ExonB->Start)
            {
                if (ExonA->End > ExonB->End)
                {
                    // A-----
                    // B---
                    //     11
                    Exon1 = NewExon(ExonB->End, ExonA->End, ReadingFrame);
                    AddMS2Exon(Builder, Exon1);
                    ExonInheritForwardEdges(Exon1, ExonA);
                    ExonInheritBackwardEdges(ExonB, ExonA);
                    DeleteMS2Exon(Builder, ExonA);
                    break;
                }
                else
                {
                    // A---  
                    // B-----
                    //     11
                    Exon1 = NewExon(ExonA->End, ExonB->End, ReadingFrame);
                    AddMS2Exon(Builder, Exon1);
                    ExonInheritForwardEdges(Exon1, ExonB);
                    ExonInheritBackwardEdges(ExonA, ExonB);
                    DELETE_EXON_B();
                    ExonB = NextExonB;
                    continue;
                }
            } // end: ExonA->Start == ExonB->Start
            if (ExonA->End == ExonB->End)
            {
                if (ExonA->Start < ExonB->Start)
                {
                    // A-----
                    // B  ---
                    //  11
                    Exon1 = NewExon(ExonA->Start, ExonB->Start, ReadingFrame);
                    AddMS2Exon(Builder, Exon1);
                    ExonInheritForwardEdges(ExonB, ExonA);
                    ExonInheritBackwardEdges(Exon1, ExonA);
                    DeleteMS2Exon(Builder, ExonA);
                    break;
                }
                else
                {
                    // A  ---  
                    // B-----
                    //  11
                    Exon1 = NewExon(ExonB->Start, ExonA->Start, ReadingFrame);
                    AddMS2Exon(Builder, Exon1);
                    ExonInheritForwardEdges(ExonA, ExonB);
                    ExonInheritBackwardEdges(Exon1, ExonB);
                    DELETE_EXON_B();
                    ExonB = NextExonB;
                    continue;
                }
            } // end: ExonA->End == ExonB->End
            if (ExonA->Start < ExonB->Start && ExonA->End < ExonB->End)
            {
                // A------
                // B   ------
                //  111222333
                Exon1 = NewExon(ExonA->Start, ExonB->Start, ReadingFrame);
                AddMS2Exon(Builder, Exon1);
                Exon2 = NewExon(ExonB->Start, ExonA->End, ReadingFrame);
                AddMS2Exon(Builder, Exon2);
                Exon3 = NewExon(ExonA->End, ExonB->End, ReadingFrame);
                AddMS2Exon(Builder, Exon3);
                ExonInheritBackwardEdges(Exon1, ExonA);
                ExonInheritBackwardEdges(Exon2, ExonB);
                ExonInheritForwardEdges(Exon2, ExonA);
                ExonInheritForwardEdges(Exon3, ExonB);
                DELETE_EXON_B();
                ExonB = NextExonB;
                DeleteMS2Exon(Builder, ExonA);
                break;
            }
            if (ExonA->Start < ExonB->Start && ExonA->End > ExonB->End)
            {
                // A---------
                // B   ---
                //  111   222
                Exon1 = NewExon(ExonA->Start, ExonB->Start, ReadingFrame);
                AddMS2Exon(Builder, Exon1);
                Exon2 = NewExon(ExonB->End, ExonA->End, ReadingFrame);
                AddMS2Exon(Builder, Exon2);
                ExonInheritBackwardEdges(Exon1, ExonA);
                ExonInheritForwardEdges(Exon2, ExonA);
                DeleteMS2Exon(Builder, ExonA);
                break;
            }
            if (ExonA->Start > ExonB->Start && ExonA->End > ExonB->End)
            {
                // A   ------
                // B------
                //  111222333
                Exon1 = NewExon(ExonB->Start, ExonA->Start, ReadingFrame);
                AddMS2Exon(Builder, Exon1);
                Exon2 = NewExon(ExonA->Start, ExonB->End, ReadingFrame);
                AddMS2Exon(Builder, Exon2);
                Exon3 = NewExon(ExonB->End, ExonA->End, ReadingFrame);
                AddMS2Exon(Builder, Exon3);
                ExonInheritBackwardEdges(Exon1, ExonB);
                ExonInheritBackwardEdges(Exon2, ExonA);
                ExonInheritForwardEdges(Exon3, ExonA);
                ExonInheritForwardEdges(Exon2, ExonB);
                DELETE_EXON_B();
                ExonB = NextExonB;
                DeleteMS2Exon(Builder, ExonA);
                break;
            }
            if (ExonA->Start > ExonB->Start && ExonA->End < ExonB->End)
            {
                // A   ---   
                // B---------
                //  111   222
                Exon1 = NewExon(ExonB->Start, ExonA->Start, ReadingFrame);
                AddMS2Exon(Builder, Exon1);
                Exon2 = NewExon(ExonA->End, ExonB->End, ReadingFrame);
                AddMS2Exon(Builder, Exon2);
                ExonInheritBackwardEdges(Exon1, ExonB);
                ExonInheritForwardEdges(Exon2, ExonB);
                DELETE_EXON_B();
                ExonB = NextExonB;
                continue;
            }
            INSPECT_ASSERT(0); // we'd better not reach this point!
        }
        ExonA = NextExonA;
    }
}

// If two exons are adjacent (one begins just after the other ends) and have 
// compatible reading frames, then add a link between them if necessary.
void AddAdjacentExonLinks(MS2Builder* Builder)
{
    MS2Exon* ExonA;
    MS2Exon* ExonB;
    MS2Edge* TestEdge;
    int LinkFound;
    //
    for (ExonA = Builder->FirstExon; ExonA; ExonA = ExonA->Next)
    {
        //printf("AAEL: %d-%d\n", ExonA->Start, ExonA->End);
        for (ExonB = Builder->FirstExon; ExonB; ExonB = ExonB->Next)
        {
            if (ExonA->End == ExonB->Start && ExonA->ReadingFrame == ExonB->ReadingFrame)
            {
                LinkFound = 0;
                for (TestEdge = ExonA->FirstForward; TestEdge; TestEdge = TestEdge->Next)
                {
                    if (TestEdge->LinkTo == ExonB)
                    {
                        LinkFound = 1;
                        break;
                    }
                }
                if (!LinkFound)
                {
                    LinkExonForward(ExonA, ExonB);
                }
            }
        }
    }
}

// Add an MS2Exon to an MS2Gene.  Also, add to the gene all the exons
// which are (recursively) linked to by MS2Exon.  
void AddExonToGene(MS2Gene* Gene, MS2Exon* Exon)
{
    MS2ExonNode* Node;
    MS2Edge* Edge;
    //
    // We follow edges forward as well as edges back, so we'll re-visit the
    // same exons, which stops the recursion:
    if (Exon->Gene == Gene)
    {
        return;
    }
    //printf("[[Add exon %d-%d R%d to gene %d\n", Exon->Start, Exon->End, Exon->ReadingFrame, Gene->Index);
    INSPECT_ASSERT(!Exon->Gene);
    Exon->Gene = Gene;
    Node = (MS2ExonNode*)calloc(1, sizeof(ExonNode));
    Node->Exon = Exon;
    if (!Gene->FirstExon)
    {
        Gene->FirstExon = Node;
    }
    else
    {
        Gene->LastExon->Next = Node;
        Node->Prev = Gene->LastExon;
    }
    Gene->LastExon = Node;
    // Follow edges:
    for (Edge = Exon->FirstForward; Edge; Edge = Edge->Next)
    {
        AddExonToGene(Gene, Edge->LinkTo);
    }
    for (Edge = Exon->FirstBackward; Edge; Edge = Edge->Next)
    {
        AddExonToGene(Gene, Edge->LinkTo);
    }
}

// Assimilate all MS2Exons from the master list into MS2Genes.  Iteratively:
// Take the first exon that's not in a gene.  Build a new gene, and add this exon, 
// and (recursively) add in everything the exon links to.
void GroupExonsIntoGenes(MS2Builder* Builder)
{
    MS2Exon* Exon;
    MS2Gene* Gene;
    //
    // Iterate over exons:
    for (Exon = Builder->FirstExon; Exon; Exon = Exon->Next)
    {
        if (Exon->Gene)
        {
            continue;
        }
        // This exon doesn't have a gene yet.  Create a gene to contain it:
        Gene = (MS2Gene*)calloc(1, sizeof(MS2Gene));
        Gene->Index = Builder->GeneCount;
        AddExonToGene(Gene, Exon);
        
        if (!Builder->FirstGene)
        {
            Builder->FirstGene = Gene;
        }
        else
        {
            Builder->LastGene->Next = Gene;
        }
        Builder->LastGene = Gene;
        Builder->GeneCount++;
    }
    // All exons are now assigned to genes.
}

// Temp-struct for sorting exons by genome-position
typedef struct MS2SortedExonNode
{
    MS2Exon* Exon;
} MS2SortedExonNode;

// Callback for qsort, to sort exons by genome-position, FORWARD strand
int CompareMS2ExonNodesForward(const MS2SortedExonNode* NodeA, const MS2SortedExonNode* NodeB)
{
    if (NodeA->Exon->Start < NodeB->Exon->Start)
    {
        return -1;
    }
    if (NodeA->Exon->Start > NodeB->Exon->Start)
    {
        return 1;
    }
    if (NodeA->Exon->End < NodeB->Exon->End)
    {
        return -1;
    }
    if (NodeA->Exon->End > NodeB->Exon->End)
    {
        return 1;
    }
    return 0;
}

// Callback for qsort, to sort exons by genome-position, REVERSE strand
int CompareMS2ExonNodesBackward(const MS2SortedExonNode* NodeA, const MS2SortedExonNode* NodeB)
{
    if (NodeA->Exon->Start < NodeB->Exon->Start)
    {
        return 1;
    }
    if (NodeA->Exon->Start > NodeB->Exon->Start)
    {
        return -1;
    }
    if (NodeA->Exon->End < NodeB->Exon->End)
    {
        return 1;
    }
    if (NodeA->Exon->End > NodeB->Exon->End)
    {
        return -1;
    }
    return 0;
}

// Read (and translate) the protein sequence for an exon.
void ReadExonSequence(MS2Builder* Builder, MS2Exon* Exon)
{
    FILE* File;
    int DNALength;
    int DNABufferSize = 0;
    char* DNABuffer = NULL;
    char* RCBuffer = NULL;
    int Modulo;
    char* TranslationStart;
    int AAIndex;
    char* TranslateMe;
    int SuffixPos;
    int AALength;
    int LengthPrefix;
    int LengthBody;
    int LengthSuffix;
    int LengthFull;
    //
    File = GlobalOptions->OutputFile;

    // Allocate a buffer to store the DNA sequence (and reverse complement):
    DNALength = Exon->End - Exon->Start;
    if (DNALength + 1 > DNABufferSize)
    {
        SafeFree(DNABuffer);
        SafeFree(RCBuffer);
        DNABufferSize = max(1024, DNALength + 5);
        DNABuffer = (char*)calloc(DNABufferSize, sizeof(char));
        RCBuffer = (char*)calloc(DNABufferSize, sizeof(char));
    }

    // Retrieve the DNA:
    fseek(Builder->GenomeFile, Exon->Start, 0);
    ReadBinary(DNABuffer, sizeof(char), DNALength, Builder->GenomeFile);
    DNABuffer[DNALength] = '\0';
    if (Builder->ForwardFlag)
    {
        Modulo = Exon->Start % CODON_LENGTH;
        if (Modulo == Exon->ReadingFrame)
        {
            TranslationStart = DNABuffer;
            Exon->Prefix[0] = '\0';
        }
        else if ((Exon->ReadingFrame + 1) % CODON_LENGTH == Modulo % CODON_LENGTH)
        {
            TranslationStart = DNABuffer + 2;
            Exon->Prefix[0] = DNABuffer[0];
            Exon->Prefix[1] = DNABuffer[1];
            Exon->Prefix[2] = '\0';
        }
        else
        {
            TranslationStart = DNABuffer + 1;
            Exon->Prefix[0] = DNABuffer[0];
            Exon->Prefix[1] = '\0';
        }
        TranslateMe = DNABuffer + strlen(Exon->Prefix);
    }
    else
    {
        WriteReverseComplement(DNABuffer, RCBuffer);
        Modulo = (Exon->End - 1) % 3;
        if (Modulo == Exon->ReadingFrame)
        {
            TranslationStart = RCBuffer;
            Exon->Prefix[0] = '\0';
        }
        else if ((Exon->ReadingFrame + 1) % CODON_LENGTH == Modulo % CODON_LENGTH)
        {
            TranslationStart = RCBuffer + 1;
            Exon->Prefix[0] = RCBuffer[0];
            Exon->Prefix[1] = '\0';
        }
        else
        {
            TranslationStart = RCBuffer + 2;
            Exon->Prefix[0] = RCBuffer[0];
            Exon->Prefix[1] = RCBuffer[1];
            Exon->Prefix[2] = '\0';
        }
        TranslateMe = RCBuffer + strlen(Exon->Prefix);
    }
    AALength = (DNALength - strlen(Exon->Prefix)) / CODON_LENGTH;
    Exon->Sequence = (char*)calloc(AALength + 1, sizeof(char));
    for (AAIndex = 0; AAIndex < AALength; AAIndex++)
    {
        Exon->Sequence[AAIndex] = TranslateCodon(TranslateMe);
        TranslateMe += CODON_LENGTH;
    }
    // Set the suffix:
    for (SuffixPos = 0; SuffixPos < CODON_LENGTH; SuffixPos++)
    {
        Exon->Suffix[SuffixPos] = *TranslateMe;
        if (!*TranslateMe)
        {
            break;
        }
        TranslateMe++;
    }
    // Double-check lengths:
    LengthPrefix = strlen(Exon->Prefix);
    LengthSuffix = strlen(Exon->Suffix);
    LengthBody = strlen(Exon->Sequence) * CODON_LENGTH;
    LengthFull = LengthPrefix + LengthSuffix + LengthBody;
    if (LengthFull != Exon->End - Exon->Start)
    {
        printf("** Error: Length %d != genomic span %d\n", LengthFull, Exon->End - Exon->Start);
    }
    SafeFree(DNABuffer);
    SafeFree(RCBuffer);
}

// Output the <Exon> tag for this MS2Exon, along with child tags (edges, and sequence)
void OutputMS2Exon(MS2Builder* Builder, MS2Gene* Gene, MS2Exon* Exon)
{
    FILE* File;
    char SpanningCodon[4];
    int LengthA;
    int LengthB;
    char AA;
    MS2Edge* Edge;
    MS2Exon* LinkExon;
    //
    File = GlobalOptions->OutputFile;
    // Start the exon tag:
    fprintf(File, "  <Exon Index=\"%d\" Start=\"%d\" End=\"%d\"", 
        Exon->Index, Exon->Start, Exon->End);
    fprintf(File, ">\n");

    if (Exon->Sequence)
    {
        fprintf(File, "    <ExonSequence Length=\"%d\">%s</ExonSequence>\n", strlen(Exon->Sequence), Exon->Sequence);
    }
    //fprintf(File, "    <ExonSequence>%s</ExonSequence>\n", Exon->Sequence);
    // Write out all the edges linking back from this exon to lower-numbered exons:
    if (Builder->ForwardFlag)
    {
        Edge = Exon->FirstBackward;
    }
    else
    {
        Edge = Exon->FirstForward;
    }
    for (; Edge; Edge = Edge->Next)
    {
        // Start an <ExtendsExon> or a <LinkFrom> tag:
        LinkExon = Edge->LinkTo;
        if (LinkExon->Start == Exon->End || LinkExon->End == Exon->Start)
        {
            fprintf(File, "    <ExtendsExon");
        }
        else
        {
            fprintf(File, "    <LinkFrom");
        }

        // Indicate the exon index:
        fprintf(File, " Index=\"%d\"", LinkExon->Index);

        // Get the amino acid!
        if (Exon->Flags & MS2_EXON_CUSTOMAA_HEAD)
        {
            AA = Exon->CustomAA;
        }
        else if (Exon->Flags & MS2_EXON_CUSTOMAA)
        {
            AA = '\0';
        }
        else if (LinkExon->Flags & MS2_EXON_CUSTOMAA)
        {
            AA = LinkExon->CustomAA;
        }
        else
        {
            // The spanning codon consists of 1 or 2 bases from this exon,
            // and 2 or 1 bases from the linked exon.
            AA = '\0';
            memset(SpanningCodon, 0, sizeof(char) * 4);
            LengthA = strlen(LinkExon->Suffix);
            LengthB = strlen(Exon->Prefix);
            if (LengthA + LengthB == CODON_LENGTH)
            {
                strcpy(SpanningCodon, LinkExon->Suffix);
                strcat(SpanningCodon, Exon->Prefix);
            }
            else if (LengthA + LengthB != 0)
            {
                // Report an error now, if the exons have incompatible reading frames!
                REPORT_ERROR_IIII(16, Exon->Start, Exon->End, LinkExon->Start, LinkExon->End);
            }

            if (SpanningCodon[0])
            {
                AA = TranslateCodon(SpanningCodon);
            }
        }
        if (AA)
        {
            fprintf(File, " AA=\"%c\"", AA);
        }

        // End the tag:
        fprintf(File, " />\n");
    }
    
    // End the exon tag:
    fprintf(File, "  </Exon>\n");
    
}

// Assign exon indexes for this gene, by first sorting the exons:
void SortMS2GeneExons(MS2Builder* Builder, MS2Gene* Gene)
{
    MS2SortedExonNode* SortedExonBlock;
    int ExonIndex;
    int ExonCount;
    MS2ExonNode* Node;

    //   
    ExonCount = 0;
    for (ExonIndex = 0, Node = Gene->FirstExon; Node; ExonIndex++, Node = Node->Next)
    {
        ExonCount++;
    }

    SortedExonBlock = (MS2SortedExonNode*)calloc(ExonCount, sizeof(MS2SortedExonNode));
    for (ExonIndex = 0, Node = Gene->FirstExon; Node; ExonIndex++, Node = Node->Next)
    {
        SortedExonBlock[ExonIndex].Exon = Node->Exon;
    }
    if (Builder->ForwardFlag)
    {
        qsort(SortedExonBlock, ExonCount, sizeof(MS2SortedExonNode), (QSortCompare)CompareMS2ExonNodesForward);
    }
    else
    {
        qsort(SortedExonBlock, ExonCount, sizeof(MS2SortedExonNode), (QSortCompare)CompareMS2ExonNodesBackward);
    }
    for (ExonIndex = 0; ExonIndex < ExonCount; ExonIndex++)
    {
        SortedExonBlock[ExonIndex].Exon->Index = ExonIndex;
        //ReadExonSequence(Builder, Gene, SortedExonBlock[ExonIndex].Exon);
    }
    DebugPrintMS2Builder(Builder, "Exons sorted");
    SafeFree(SortedExonBlock);
}

// Generate XML for an MS2CrossReference.  We need the GFF gene's database and accession number,
// and we need the list of exon indices.
void OutputMS2CrossReference(MS2Builder* Builder, FILE* File, MS2CrossReference* CR)
{
    IntNode* Node;
    //
    fprintf(File, "  <CrossReference Database=\"%s\" ID=\"%s\">\n", CR->GFFGene->DatabaseName, CR->GFFGene->Name);
    fprintf(File, "    <CRExons Index=\"");
    for (Node = CR->FirstExonID; Node; Node = Node->Next)
    {
        if (Node->Next)
        {
            fprintf(File, "%d, ", Node->Value);
        }
        else
        {
            fprintf(File, "%d", Node->Value);
        }
    }
    fprintf(File, "\"/>\n");
    fprintf(File, "  </CrossReference>\n");
}

// Output the XML for this MS2Gene.
void OutputMS2Gene(MS2Builder* Builder, MS2Gene* Gene)
{
    FILE* File;
    int ExonCount;
    int ExonIndex;
    MS2ExonNode* Node;
    MS2CrossReference* CR;
    
    File = GlobalOptions->OutputFile;

    // Count exons in the gene:
    ExonCount = 0;
    for (ExonIndex = 0, Node = Gene->FirstExon; Node; ExonIndex++, Node = Node->Next)
    {
        ExonCount++;
    }

    // Start the Gene tag:
    fprintf(File, "<Gene ExonCount=\"%d\" Chromosome=\"%s\" ForwardFlag=\"%d\">\n", ExonCount, Builder->ChromosomeName, Builder->ForwardFlag);

    // Loop over exons, and ouptut an <Exon> tag for each one:
    for (ExonIndex = 0; ExonIndex < ExonCount; ExonIndex++)
    {
        for (Node = Gene->FirstExon; Node; Node = Node->Next)
        {
            if (Node->Exon->Index != ExonIndex)
            {
                continue;
            }
            OutputMS2Exon(Builder, Gene, Node->Exon);
            break;
        }
    }

    // Output all cross-references for the gene:
    for (CR = Gene->FirstCrossReference; CR; CR = CR->Next)
    {
        OutputMS2CrossReference(Builder, File, CR);
    }

    // Complete the Gene tag:
    fprintf(File, "</Gene>\n\n");
    fflush(File);
    
}

// Convert a codon into a number from 0 to 63.  (We probably could just TRANSLATE the
// codon and use the amino acid value...)
int GetCodonHashValue(char* EncodeCodon)
{
    int Pos;
    int Multiplier[] = {1, 4, 16};
    int Value = 0;
    for (Pos = 0; Pos < CODON_LENGTH; Pos++)
    {
        switch (EncodeCodon[Pos])
        {
        case 'a':
        case 'A':
            Value += 0 * Multiplier[Pos];
            break;
        case 'c':
        case 'C':
            Value += 1 * Multiplier[Pos];
            break;
        case 'g':
        case 'G':
            Value += 2 * Multiplier[Pos];
            break;
        case 't':
        case 'T':
            Value += 3 * Multiplier[Pos];
            break;
        default:
            //printf("* Error in GetCodonHashValue('%c')\n", EncodeCodon[Pos]);
            REPORT_ERROR_I(24, EncodeCodon[Pos]);
            return 0;
        }
    }
    return Value;
}

// Scenario: What if exon X consists of a single base pair!?
// That's tricky if the base pair is the middle of a codon, because the
// linked exons must look PAST this central exon to get their prefix / suffix.
// We'll produce special degree-1 "customAA" exons, one for each codon, to
// get from each predecessor of exon X to each successor of exon X.  
// We produce one CustomAA exon for each possible codon.  
void RepairPromiscuousSingletonExons(MS2Builder* Builder)
{
    MS2Exon* Exon;
    MS2Exon* NextExon;
    int Modulo;
    MS2Exon* CodonExons[64];
    MS2Edge* BackwardEdge;
    MS2Edge* ForwardEdge;
    MS2Edge* Edge;
    char Codon[4];
    char RCCodon[4];
    char* EncodeCodon;
    int CodonValue;
    int BridgedFlag;
    //
    memset(Codon, 0, sizeof(char) * 4);
    memset(RCCodon, 0, sizeof(char) * 4);
    Exon = Builder->FirstExon;
    while (1)
    {
        if (!Exon)
        {
            break;
        }
        // Skip this exon if its length is > 1:
        if (Exon->End > Exon->Start + 1)
        {
            Exon = Exon->Next;
            continue;
        }
        // Skip this exon, unless its base pair is the middle of a codon:
        if (Builder->ForwardFlag)
        {
            Modulo = Exon->Start % CODON_LENGTH;
            if ((Exon->ReadingFrame + 1) % CODON_LENGTH != Modulo)
            {
                Exon = Exon->Next;
                continue;
            }
        }
        else
        {
            Modulo = (Exon->End - 1) % CODON_LENGTH;
            if ((Exon->ReadingFrame - 1) % CODON_LENGTH != Modulo)
            {
                Exon = Exon->Next;
                continue;
            }
        }
        // Skip customAA exons built by previous passes through this loop:
        if (Exon->CustomAA)
        {
            Exon = Exon->Next;
            continue;
        }
        // This is the tricky case: A length-1 exon in the middle of a codon.
        // Consider every pairing of a predecessor and a successor.  For each distinct codon,
        // build one CustomAA exon:
        memset(CodonExons, 0, sizeof(MS2Exon*) * 64);
        BridgedFlag = 0;
        for (BackwardEdge = Exon->FirstBackward; BackwardEdge; BackwardEdge = BackwardEdge->Next)
        {
            for (ForwardEdge = Exon->FirstForward; ForwardEdge; ForwardEdge = ForwardEdge->Next)
            {
                if (Builder->ForwardFlag)
                {
                    Codon[0] = BackwardEdge->LinkTo->Suffix[0];
                    Codon[1] = Exon->Prefix[0];
                    Codon[2] = ForwardEdge->LinkTo->Prefix[0];
                    EncodeCodon = Codon;
                }
                else
                {
                    Codon[0] = ForwardEdge->LinkTo->Suffix[0];
                    Codon[1] = Exon->Prefix[0];
                    Codon[2] = BackwardEdge->LinkTo->Prefix[0];
                    WriteReverseComplement(Codon, RCCodon);
                    EncodeCodon = RCCodon;
                }
                CodonValue = GetCodonHashValue(EncodeCodon);
                INSPECT_ASSERT(CodonValue >= 0 && CodonValue < 64);
                // CodonExons[CodonValue] will hold the custom-aa exon:
                if (!CodonExons[CodonValue])
                {
                    CodonExons[CodonValue] = (MS2Exon*)calloc(1, sizeof(MS2Exon));
                    CodonExons[CodonValue]->CustomAA = TranslateCodon(EncodeCodon);
                    CodonExons[CodonValue]->Start = Exon->Start;
                    CodonExons[CodonValue]->End = Exon->End;
                    CodonExons[CodonValue]->ReadingFrame = Exon->ReadingFrame;
                    AddMS2Exon(Builder, CodonExons[CodonValue]);
                }
                ExonInheritOneForwardEdge(CodonExons[CodonValue], ForwardEdge);
                ExonInheritOneBackwardEdge(CodonExons[CodonValue], BackwardEdge);
                BridgedFlag = 1;
            }
        }

        // Assign flags to these CustomAA exons.  All this work to handle 
        // the case of a single-base-pair exon with out-degree 1 whose 
        // outgoing edge is to an adjacent exon; in that case, we want the 
        // amino acid to be placed on our incoming edges rather than 
        // on the outgoing edges.
        for (CodonValue = 0; CodonValue < 64; CodonValue++)
        {
            if (CodonExons[CodonValue])
            {
                CodonExons[CodonValue]->Flags = MS2_EXON_CUSTOMAA;
                if (Builder->ForwardFlag)
                {
                    Edge = CodonExons[CodonValue]->FirstForward;
                }
                else
                {
                    Edge = CodonExons[CodonValue]->FirstBackward;
                }
                if (!Edge || !Edge->Next)
                {
                    continue;
                }
                if (Edge->LinkFrom->Start == Edge->LinkTo->End || Edge->LinkFrom->End == Edge->LinkTo->Start)
                {
                    CodonExons[CodonValue]->Flags = MS2_EXON_CUSTOMAA_HEAD;
                }
            }
        }
        NextExon = Exon->Next;
        if (BridgedFlag)
        {
            DeleteMS2Exon(Builder, Exon);
        }
        else
        {
            REPORT_WARNING_I(21, Exon->Start);
        }
        Exon = NextExon;
    }
}

void ExonInheritOneForwardEdge(MS2Exon* Exon, MS2Edge* OldEdge)
{
    MS2Exon* DisplacedExon;
    MS2Exon* LinkedExon;
    //
    DisplacedExon = OldEdge->LinkFrom;
    LinkedExon = OldEdge->LinkTo;
    RemoveBackwardEdge(LinkedExon, DisplacedExon);
    RemoveForwardEdge(DisplacedExon, LinkedExon);
    LinkExonForward(Exon, LinkedExon);
}

void ExonInheritOneBackwardEdge(MS2Exon* Exon, MS2Edge* OldEdge)
{
    MS2Exon* DisplacedExon;
    MS2Exon* LinkedExon;
    //
    DisplacedExon = OldEdge->LinkFrom;
    LinkedExon = OldEdge->LinkTo;
    RemoveBackwardEdge(DisplacedExon, LinkedExon);
    RemoveForwardEdge(LinkedExon, DisplacedExon);
    LinkExonForward(LinkedExon, Exon);
}

// Main entry point for building MS2 database.
void BuildMS2DB()
{
    MS2Builder* Builder;
    int ForwardFlag;
    MS2Gene* Gene;
    MS2Exon* Exon;
    //
    Builder = (MS2Builder*)calloc(1, sizeof(MS2Builder));
    // Builder->VerboseFlag = 1; // spewy!
    // Open the genome file:
    Builder->GenomeFile = fopen(GlobalOptions->GenomeFileName, "rb");
    if (!Builder->GenomeFile)
    {
        REPORT_ERROR_S(8, GlobalOptions->GenomeFileName);
        goto cleanup;
    }
    // At least one GFF file must be specified!
    if (!GlobalOptions->FirstGFFFileName)
    {
        REPORT_ERROR(12);
        goto cleanup;
    }
    fprintf(GlobalOptions->OutputFile, "<Database CreatedBy=\"BuildMS2DB.c\">\n");
    // Loop: First the forward strand, then the reverse strand:
    for (ForwardFlag = 1; ForwardFlag >= 0; ForwardFlag--)
    {
        Builder->ForwardFlag = ForwardFlag;
        Builder->ExonHash = (MS2ExonNode**)calloc(EXON_HASH_SIZE, sizeof(ExonNode*));
        strncpy(Builder->ChromosomeName, GlobalOptions->ChromosomeName, 256);

        // Parse exons from GFF files:
        ParseGFFFiles(Builder);

        // Bail out, if we have no exons at all:
        if (!Builder->FirstExon)
        {
            REPORT_ERROR(15);
            continue;
        }

        printf("Parsed GFF files.  We now have %d exons.\n", Builder->ExonCount);
        DebugPrintMS2Builder(Builder, "After GFF parse");

        // Merge and split any overlapping exons as needed.  Note that if we merge exons, then we
        // can't report a cross-reference ("record FOO covers exons 1, 2, 3, 4, 5").  Therefore, 
        // most exons are NOT permitted to be merged.  Only exons produced from EST alignments
        // should be considerd merge-able.
        SplitOverlappingExons(Builder);
        
        DebugPrintMS2Builder(Builder, "After exon split");

        // Add edges between adjacent exons:
        AddAdjacentExonLinks(Builder);

        // Read all exon sequences.  We *could* read just the exons for one gene at a time.
        for (Exon = Builder->FirstExon; Exon; Exon = Exon->Next)
        {
            ReadExonSequence(Builder, Exon);
        }

        // Ensure that length-1 exons (if any exist!) have at most one back-link.
        RepairPromiscuousSingletonExons(Builder);

        // Group exons into genes:
        GroupExonsIntoGenes(Builder);

        DebugPrintMS2Builder(Builder, "After gene grouping");

        // Sort exons within genes, assigning exons index numbers:
        for (Gene = Builder->FirstGene; Gene; Gene = Gene->Next)
        {
            SortMS2GeneExons(Builder, Gene);
        }

        // Add cross-references to genes:
        BuildGFFCrossReferences(Builder);

        // Output XML:
        for (Gene = Builder->FirstGene; Gene; Gene = Gene->Next)
        {
            OutputMS2Gene(Builder, Gene);
        }
        // Free our exon hash, exon lists, gene lists, etc:
        FreeExonHash(Builder);
        FreeMS2Genes(Builder);
        FreeGFFGenes(Builder);
        FreeMS2Exons(Builder);
    }
    fprintf(GlobalOptions->OutputFile, "\n</Database>\n");
cleanup:
    FreeExonHash(Builder);
    if (Builder->GenomeFile)
    {
        fclose(Builder->GenomeFile);
    }
    SafeFree(Builder);
}


// Handy debugging function: Spew out *all* the exons and genes parsed so far!
void DebugPrintMS2Builder(MS2Builder* Builder, char* Notes)
{
    MS2Exon* Exon;
    MS2ExonNode* Node;
    MS2Gene* Gene;
    int GeneExonCount;
    int GeneStart;
    int GeneEnd;
    int GeneIndex;
    MS2Edge* Edge;
    int ExonIndex;
    int ExonCount = 0;
    int ForwardEdgeCount = 0;
    int BackwardEdgeCount = 0;
    MS2Edge* PrevEdge;
    //
    if (!Builder->VerboseFlag)
    {
        return;
    }
    printf("\n=-=-{O}=-=-{O}=-=-{O}=-=-{O}=-=-{O}=-=-{O}=-=-{O}=-=-{O}=-=-{O}=-=-\n");
    if (Notes)
    {
        printf("*-*-> %s\n", Notes);
    }
    else
    {
        printf("*-*-> MS2Builder state:\n");
    }
    
    for (Exon = Builder->FirstExon, ExonIndex = 0; Exon; Exon = Exon->Next, ExonIndex++)
    {
        printf("  Exon %d: %d-%d R %d", ExonIndex, Exon->Start, Exon->End, Exon->ReadingFrame);
        if (Exon->Gene)
        {
            printf(" Gene %d", Exon->Gene->Index);
        }
        printf("\n");
        if (Exon->Sequence)
        {
            printf("    Prefix '%s' Suffix '%s'\n", Exon->Prefix, Exon->Suffix);
            INSPECT_ASSERT(strlen(Exon->Sequence) * 3 + strlen(Exon->Prefix) + strlen(Exon->Suffix) == (Exon->End - Exon->Start));
        }
        PrevEdge = NULL;
        for (Edge = Exon->FirstForward; Edge; Edge = Edge->Next)
        {
            printf("    >>>Link to %d-%d R%d\n", Edge->LinkTo->Start, Edge->LinkTo->End, Edge->LinkTo->ReadingFrame);
            ForwardEdgeCount++;
            PrevEdge = Edge;
        }
        if (PrevEdge != Exon->LastForward)
        {
            printf("   *** Error: LastForward link is corrupt!\n");
        }
        PrevEdge = NULL;
        for (Edge = Exon->FirstBackward; Edge; Edge = Edge->Next)
        {
            printf("    <<<Link from %d-%d R%d\n", Edge->LinkTo->Start, Edge->LinkTo->End, Edge->LinkTo->ReadingFrame);
            BackwardEdgeCount++;
            PrevEdge = Edge;
        }
        if (PrevEdge != Exon->LastBackward)
        {
            printf("   *** Error: LastForward link is corrupt!\n");
        }

        ExonCount++;
    }
    printf("\n");
    for (Gene = Builder->FirstGene, GeneIndex = 0; Gene; Gene = Gene->Next, GeneIndex++)
    {
        GeneExonCount = 0;
        GeneStart = Gene->FirstExon->Exon->Start;
        GeneEnd = Gene->FirstExon->Exon->End;
        for (Node = Gene->FirstExon; Node; Node = Node->Next)
        {
            GeneExonCount++;
            GeneStart = min(GeneStart, Node->Exon->Start);
            GeneEnd = max(GeneEnd, Node->Exon->End);
        }
        printf("Gene %d/%d (%d...%d) has %d exons\n", GeneIndex, Gene->Index, GeneStart, GeneEnd, GeneExonCount);
        for (Node = Gene->FirstExon; Node; Node = Node->Next)
        {
            if (Node->Exon->Gene != Gene)
            {
                printf("** ERROR: Exon %d-%d doesn't link up!\n", Node->Exon->Start, Node->Exon->End);
            }
        }
    }
    printf("\n...total of %d exons, %d/%d edges\n", ExonCount, ForwardEdgeCount, BackwardEdgeCount);
}
