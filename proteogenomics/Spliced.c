
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
#include "Trie.h"
#include "Spliced.h"
#include "Inspect.h"
#include "ExonGraphAlign.h"
#include "Errors.h"
// The left and right extensions of a tag can have at most this many successes.
#define MAX_SIDE_EXTENSIONS 128

// For keeping a linked list of genes in memory.  In practice, we generally DON'T use this list,
// instead we load and search one gene at a time.  (Genomes are big!)
GeneStruct* FirstGene;
GeneStruct* LastGene;

ExonStruct** g_TagExonArray;
int g_TagExonArrayPos;
ExonEdge** g_TagSpliceArray;
int g_TagSpliceArrayPos;
// Tags often contain residues K/Q and I/L.  We wish to report the true residue found
// in the exon, undoing Q->K and I->L substitutions as necessary.  So, we log the 
// matched tag chars in g_TagBuffer.
char* g_TagBuffer;
char* g_TagBufferSpliced;
int g_TagBufferPos;
int g_TagBufferPosSpliced;

static char* MatchedBases;
char* ExtensionBufferLeft;
char* ExtensionBufferRight;
int* ExtensionLeftDecorations;
int* ExtensionRightDecorations;
int* ExtensionGenomicStart;
int* ExtensionGenomicEnd;
MSSpectrum** ExtensionSpectra;

static int MH_MinMatchMass;
static int MH_MaxMatchMass;
static char* MH_MatchBuffer;
static char* MH_MatchBufferSpliced;
static int* MH_MatchDecoration;
static int MH_MatchCount;
static char* MH_Buffer;
static char* MH_BufferSpliced;
static ExonStruct** MH_MatchExons;
static ExonEdge** MH_MatchEdges;
static int MH_MatchExonPos;
static int MH_MatchEdgePos;
//static int MH_MatchSplicePos;

ExonEdge* GetReciprocalExonEdge(ExonEdge* Edge, int ForwardFlag);

// Free one Gene (and its exons)
void FreeGene(GeneStruct* Gene)
{
    int Index;
    //
    if (!Gene)
    {
        return;
    }
    for (Index = 0; Index < Gene->ExonCount; Index++)
    {
        SafeFree(Gene->Exons[Index].ForwardEdges);
        SafeFree(Gene->Exons[Index].BackwardEdges);
        SafeFree(Gene->Exons[Index].Sequence);
    }
    SafeFree(Gene->Exons);
    SafeFree(Gene);
}

// Free the global list of genes.  (Not used in practice, since we load one at a time)
void FreeGenes()
{
    GeneStruct* Gene;
    GeneStruct* Prev = NULL;
    //
    for (Gene = FirstGene; Gene; Gene = Gene->Next)
    {
        if (Prev)
        {
            FreeGene(Prev);
        }
        Prev = Gene;
    }
    if (Prev)
    {
        FreeGene(Prev);
    }
    FirstGene = NULL;
    LastGene = NULL;
}

// For debugging purposes: Print out a list of exons (with partial sequences) and edges.
// (Mostly for verifying database generation worked)
void DebugPrintGene(GeneStruct* Gene)
{
    int ExonIndex;
    ExonStruct* Exon;
    int EdgeIndex;
    //
    printf("*Gene %s (%s) has %d exons\n", Gene->Name, Gene->SprotName, Gene->ExonCount);
    for (ExonIndex = 0; ExonIndex < Gene->ExonCount; ExonIndex++)
    {
        Exon = Gene->Exons + ExonIndex;
        printf("Exon %d from %d-%d cov %d: \n", Exon->Index, Exon->Start, Exon->End, Exon->Occurrences);
        if (Exon->Sequence)
        {
            printf(Exon->Sequence);
        }
        else
        {
            printf("<none>");
        }
        printf("\n");
        //printf("  Exon from %d-%d coverage %d sequence %s...\n", Exon->Start, Exon->End, Exon->Occurrences, Buffer);
        for (EdgeIndex = 0; EdgeIndex < Exon->ForwardEdgeCount; EdgeIndex++)
        {
            printf("  >> (%d) '%c' to exon #%d %d-%d\n", Exon->ForwardEdges[EdgeIndex].Power, 
                Exon->ForwardEdges[EdgeIndex].AA, Exon->ForwardEdges[EdgeIndex].Exon->Index, 
                Exon->ForwardEdges[EdgeIndex].Exon->Start, Exon->ForwardEdges[EdgeIndex].Exon->End);
        }
        for (EdgeIndex = 0; EdgeIndex < Exon->BackEdgeCount; EdgeIndex++)
        {
            printf("  << (%d) '%c' to exon #%d %d-%d\n", Exon->BackwardEdges[EdgeIndex].Power, 
                Exon->BackwardEdges[EdgeIndex].AA, Exon->BackwardEdges[EdgeIndex].Exon->Index, 
                Exon->BackwardEdges[EdgeIndex].Exon->Start, Exon->BackwardEdges[EdgeIndex].Exon->End);
        }
    }
}

// For debugging: Print *all* our genes, and their exons and edges
void DebugPrintGenes()
{
    GeneStruct* Gene;
    //
    printf("Genes:\n");
    for (Gene = FirstGene; Gene; Gene = Gene->Next)
    {
        printf("\n");
        DebugPrintGene(Gene);
    }
}

// Load one gene from the (binary) gene file.  Does some basic error checking, in case of 
// obsolete or broken file formats.
GeneStruct* LoadGene(FILE* File)
{
    char Buffer[1024];
    int Bytes;
    GeneStruct* Gene;
    int ExonIndex;
    int OtherExonIndex;
    ExonStruct* Exon;
    ExonStruct* OtherExon;
    int Length;
    int EdgeIndex;
    char AA;
    int LinkPower;
    //
    Bytes = ReadBinary(Buffer, sizeof(char), GENE_NAME_LENGTH, File);
    if (!Bytes)
    {
        return NULL; // eof
    }   
    Gene = (GeneStruct*)calloc(1, sizeof(GeneStruct));
    strncpy(Gene->Name, Buffer, GENE_NAME_LENGTH);
    ReadBinary(Gene->SprotName, sizeof(char), GENE_NAME_LENGTH, File);
    ReadBinary(&Gene->ChromosomeNumber, sizeof(int), 1, File);
    if (!Gene->ChromosomeNumber)
    {
        printf("** Warning: No chromosome number for gene '%s'\n", Gene->Name);
    }
    ReadBinary(&Gene->ForwardFlag, sizeof(char), 1, File);
    ReadBinary(&Gene->ExonCount, sizeof(int), 1, File);
    if (Gene->ExonCount < 1 || Gene->ExonCount > MAX_GENE_EXONS)
    {
        printf("** Warning: suspicious exon-count %d encountered in LoadGene().  File position is %ld.\n", Gene->ExonCount, ftell(File));
        return NULL;
    }
    //fread(&GIIDBlock, sizeof(int), 10, File);
    Gene->Exons = (ExonStruct*)calloc(Gene->ExonCount, sizeof(ExonStruct));

    // Read the gene's exons:
    for (ExonIndex = 0; ExonIndex < Gene->ExonCount; ExonIndex++)
    {
        //printf("Filepos %d, now read exon %d of %d\n", ftell(File), ExonIndex, Gene->ExonCount);
        Exon = Gene->Exons + ExonIndex;
        Exon->Gene = Gene;
        Bytes = ReadBinary(&Exon->Start, sizeof(int), 1, File);
        if (!Bytes)
        {
            printf("** Error: EOF encountered while reading exon %d of gene '%s'\n", ExonIndex, Gene->Name);
            break;
        }
        Exon->Index = ExonIndex;
        ReadBinary(&Exon->End, sizeof(int), 1, File);
        ReadBinary(&Length, sizeof(int), 1, File);
        if (Length < 0 || Length > 10000)
        {
            printf("** Error: Bogus sequence length %d encountered while reading exon %d of gene '%s'\n", Length, ExonIndex, Gene->Name);
            break;
        }
        ReadBinary(&Exon->Occurrences, sizeof(int), 1, File);
        Exon->Length = Length;
        if (Length)
        {
            Exon->Sequence = (char*)calloc(Length + 1, sizeof(char));
            ReadBinary(Exon->Sequence, sizeof(char), Length, File);
        }
        else
        {
            Exon->Sequence = NULL;
        }
        //printf("%d '%s'\n", ExonIndex, Exon->Sequence); // 
        ReadBinary(&Exon->Prefix, sizeof(char), 2, File);
        ReadBinary(&Exon->Suffix, sizeof(char), 2, File);
        ReadBinary(&Exon->BackEdgeCount, sizeof(int), 1, File);
        if (Exon->BackEdgeCount < 0 || Exon->BackEdgeCount > 500)
        {
            printf("** zomg broken back edge count in LoadGene() exon %d gene '%s'\n", ExonIndex, Gene->Name);
        }
        ReadBinary(&Exon->ForwardEdgeCount, sizeof(int), 1, File);
        if (Exon->ForwardEdgeCount < 0 || Exon->ForwardEdgeCount > 500)
        {
            printf("** zomg broken forward edge count in LoadGene() exon %d gene '%s'\n", ExonIndex, Gene->Name);
        }
        
        if (Exon->ForwardEdgeCount)
        {
            Exon->ForwardEdges = (ExonEdge*)calloc(Exon->ForwardEdgeCount, sizeof(ExonEdge));
        }
        if (Exon->BackEdgeCount)
        {
            Exon->BackwardEdges = (ExonEdge*)calloc(Exon->BackEdgeCount, sizeof(ExonEdge));
        }
        // Read all the edges for this exon.  (Read all the back-edges, and THEN take care of old forward-edges)
        for (EdgeIndex = 0; EdgeIndex < Exon->BackEdgeCount; EdgeIndex++)
        {
            Bytes = ReadBinary(&OtherExonIndex, sizeof(int), 1, File);
            if (!Bytes)
            {
                printf("** Error: EOF encountered while reading exon %d edge %d of gene '%s'\n", ExonIndex, EdgeIndex, Gene->Name);
                break;
            }
            ReadBinary(&LinkPower, sizeof(int), 1, File);
            ReadBinary(&AA, sizeof(char), 1, File);
            if (OtherExonIndex < 0 || OtherExonIndex >= Gene->ExonCount)
            {
                printf("** Error: Illegal exon back-link %d encountered for exon %d edge %d gene '%s'\n", OtherExonIndex, ExonIndex, EdgeIndex, Gene->Name);
            }
            else
            {
                OtherExon = Gene->Exons + OtherExonIndex;
                Exon->BackwardEdges[EdgeIndex].Exon = OtherExon;
                Exon->BackwardEdges[EdgeIndex].AA = AA;
                Exon->BackwardEdges[EdgeIndex].Power = LinkPower;
                Exon->BackwardEdges[EdgeIndex].Source = Exon;
            }
        }
    } // exon loop
    // We set all the back-links while we're reading the exons in.  Now, let's go through
    // and fix all the forward-links.  
    SetExonForwardEdges(Gene);
    return Gene;
}

// INPUT: A gene where the backward edges are populated, and the exon forward edges are allocated but *not* populated.
// Result: Forward edges are populated.
void SetExonForwardEdges(GeneStruct* Gene)
{
    int ExonIndex;
    ExonStruct* Exon;
    ExonStruct* OtherExon;
    int EdgeIndex;
    int OtherEdgeIndex;
    int ForwardEdgeSet;
    char AA;
    int LinkPower;
    //
    for (ExonIndex = 0; ExonIndex < Gene->ExonCount; ExonIndex++)
    {
        Exon = Gene->Exons + ExonIndex;
        for (EdgeIndex = 0; EdgeIndex < Exon->BackEdgeCount; EdgeIndex++)
        {
            // The first empty slot in the OtherExon forward arrays will now be set:
            OtherExon = Exon->BackwardEdges[EdgeIndex].Exon;
            AA = Exon->BackwardEdges[EdgeIndex].AA;
            LinkPower = Exon->BackwardEdges[EdgeIndex].Power;
            ForwardEdgeSet = 0;
            for (OtherEdgeIndex = 0; OtherEdgeIndex < OtherExon->ForwardEdgeCount; OtherEdgeIndex++)
            {
                if (!OtherExon->ForwardEdges[OtherEdgeIndex].Exon)
                {
                    OtherExon->ForwardEdges[OtherEdgeIndex].Exon = Exon;
                    OtherExon->ForwardEdges[OtherEdgeIndex].AA = AA;
                    OtherExon->ForwardEdges[OtherEdgeIndex].Power = LinkPower;
                    OtherExon->ForwardEdges[OtherEdgeIndex].Source = OtherExon;
                    ForwardEdgeSet = 1;
                    break;
                }
            }
            if (!ForwardEdgeSet)
            {
                REPORT_ERROR_IIS(26, OtherExon->Index, Exon->Index, Gene->Name);
            }
        }
    }
}

// Load genes from a binary file, built by running inspect with splicedb arguments.
// (In practice, we don't call this - we just load ONE gene at a time!)
void LoadGenes(char* FileName)
{
    FILE* File;
    GeneStruct* Gene;
    //
    File = fopen(FileName, "rb");
    if (!File)
    {
        printf("** Error: Unable to open gene file '%s'\n", FileName);
        return;
    }
    while (1)
    {
        Gene = LoadGene(File);
        if (!Gene)
        {
            break;
        }

        // Insert new gene into list:
        if (LastGene)
        {
            LastGene->Next = Gene;
            Gene->Prev = LastGene;
        }
        else
        {
            FirstGene = Gene;
        }
        LastGene = Gene;
    }
    fclose(File);
}

// Static structures used in splice-tolerant search.  Tag extension builds up an array
// of extension-matches for the left and for the right, then tries each combination of
// a legal right-extension and a legal left-extension.
char* SLeftMatchBuffer = NULL; // The AAs of the extension
char* SLeftMatchBufferSpliced = NULL; // The AAs of the extension, with splice boundaries
int* SLeftMatchDecoration = NULL; // The decoration to be attached over the extension
int* SLeftGenomicPosition = NULL;
int* SRightGenomicPosition = NULL;
char* SLeftPrefixes = NULL; // The AA just *beyond* the extension
char* SRightMatchBuffer = NULL;
char* SRightMatchBufferSpliced = NULL;
int* SRightMatchDecoration = NULL;
char* SRightSuffixes = NULL;
ExonStruct** SLeftExon = NULL; // The exons reached by prefix extension.  SLeftExon[MatchNumber*16 + ExonIndex]
ExonEdge** SLeftEdge = NULL; // The splice boundaries crossed by prefix extension.
//int* SLeftSpliceScore = NULL; // The scores of splice boundaries used in prefix extension.
int* SLeftExonCount = NULL;
int* SLeftSpliceCount = NULL;
ExonStruct** SRightExon = NULL; // The exons reached by suffix extension.  SLeftExon[MatchNumber*16 + ExonIndex]
ExonEdge** SRightEdge = NULL; // The splice boundaries crossed by suffix extension.
//int* SRightSpliceScore = NULL;
int* SRightExonCount = NULL;
int* SRightSpliceCount = NULL;

void AllocSpliceStructures()
{
    if (SLeftMatchBuffer)
    {
        return; // It seems we've already allocated them.
    }
    
    // The Spliced buffer is made extra-long so that we can afford to add two symbol chars
    // per amino acid.
    SLeftMatchBuffer = (char*)calloc(MAX_EXTENSION_LENGTH * MAX_SIDE_EXTENSIONS, sizeof(char));
    SLeftMatchBufferSpliced = (char*)calloc(MAX_SEXTENSION_LENGTH * MAX_SIDE_EXTENSIONS, sizeof(char));
    SLeftMatchDecoration = (int*)calloc(MAX_SIDE_EXTENSIONS + 1, sizeof(int));
    SLeftGenomicPosition = (int*)calloc(MAX_SIDE_EXTENSIONS + 1, sizeof(int));
    SLeftPrefixes = (char*)calloc(MAX_SIDE_EXTENSIONS, sizeof(char));
    SRightMatchBuffer = (char*)calloc(MAX_EXTENSION_LENGTH * MAX_SIDE_EXTENSIONS, sizeof(char));
    SRightMatchBufferSpliced = (char*)calloc(MAX_SEXTENSION_LENGTH * MAX_SIDE_EXTENSIONS, sizeof(char));
    SRightMatchDecoration = (int*)calloc(MAX_SIDE_EXTENSIONS + 1, sizeof(int));
    SRightGenomicPosition = (int*)calloc(MAX_SIDE_EXTENSIONS + 1, sizeof(int));
    SRightSuffixes = (char*)calloc(MAX_SIDE_EXTENSIONS, sizeof(char));

    SLeftExon = (ExonStruct**)calloc(MAX_SIDE_EXTENSIONS * MAX_EXTENSION_EXONS, sizeof(ExonStruct*));
    SLeftEdge = (ExonEdge**)calloc(MAX_SIDE_EXTENSIONS * MAX_EXTENSION_EXONS, sizeof(ExonEdge*));
    //SLeftSpliceScore = (int*)calloc(MAX_SIDE_EXTENSIONS * MAX_EXTENSION_EXONS, sizeof(int));
    SLeftExonCount = (int*)calloc(MAX_SIDE_EXTENSIONS, sizeof(int));
    SLeftSpliceCount = (int*)calloc(MAX_SIDE_EXTENSIONS, sizeof(int));

    SRightExon = (ExonStruct**)calloc(sizeof(ExonStruct*), MAX_SIDE_EXTENSIONS * MAX_EXTENSION_EXONS);
    SRightEdge = (ExonEdge**)calloc(MAX_SIDE_EXTENSIONS * MAX_EXTENSION_EXONS, sizeof(ExonEdge*));
    //SRightSpliceScore = (int*)calloc(sizeof(int), MAX_SIDE_EXTENSIONS * MAX_EXTENSION_EXONS);
    SRightExonCount = (int*)calloc(sizeof(int), MAX_SIDE_EXTENSIONS);
    SRightSpliceCount = (int*)calloc(sizeof(int), MAX_SIDE_EXTENSIONS);

    g_TagExonArray = (ExonStruct**)calloc(sizeof(ExonStruct*), 16);
    g_TagSpliceArray = (ExonEdge**)calloc(sizeof(ExonEdge*), 16);
    g_TagBuffer = (char*)calloc(sizeof(int), 10);
    g_TagBufferSpliced = (char*)calloc(sizeof(int), MAX_EXTENSION_EXONS);

    MH_Buffer = (char*)calloc(sizeof(char), 128);
    MH_BufferSpliced = (char*)calloc(sizeof(char), 128);
    MH_MatchExons = (ExonStruct**)calloc(sizeof(ExonStruct*), MAX_EXTENSION_EXONS);
    MH_MatchEdges = (ExonEdge**)calloc(sizeof(ExonEdge*), MAX_EXTENSION_EXONS);
    //MH_MatchSplices = (int*)calloc(sizeof(int), MAX_EXTENSION_EXONS);

    MatchedBases = (char*)calloc(sizeof(char), 512);
    ExtensionBufferLeft = (char*)calloc(sizeof(char), MAX_EXTENSION_LENGTH*512);
    ExtensionBufferRight = (char*)calloc(sizeof(char), MAX_EXTENSION_LENGTH*512);
    ExtensionLeftDecorations = (int*)calloc(sizeof(int), 512);
    ExtensionRightDecorations = (int*)calloc(sizeof(int), 512);
    ExtensionGenomicStart = (int*)calloc(sizeof(int), 512);
    ExtensionGenomicEnd = (int*)calloc(sizeof(int), 512);

    ExtensionSpectra = (MSSpectrum**)calloc(sizeof(MSSpectrum*), 512);
}


// Helper function: We've successfully extended a tag either forward (Direction=1) or backward (Direction=-1)
// along the peptide.  Set the genomic endpoint, and the flanking (prefix or suffix) amino acid character.
// What makes the job tricky is that we may have finished at the edge of an exon, either by using up an
// incoming edge (if Pos==-1) or by using up the entire exon (if Pos+Direction falls off the edge).  
// If we used up the full exon and there's an edge, report the AA for the first edge.  (TODO: Sort edges, maybe).  
// If we used up a full exon and there's nothing to link to, report char '-'
void SetMatchPrefixSuffix(ExonStruct* Exon, int Pos, int Direction)
{
    int Length = 0;
    char AA;
    ExonEdge* Edge;
    if (Direction > 0)
    {
        // Direction is 1, so set RightGenomicPosition:
        if (Exon->Start < 0)
        {
            // This exon has no knwon genomic position:
            SRightGenomicPosition[MH_MatchCount] = -1;
        }
        else if (Exon->Sequence)
        {
            if (Exon->Gene->ForwardFlag)
            {
                if (Pos > -1)
                {
                    SRightGenomicPosition[MH_MatchCount] = Exon->Start + (Pos+1)*3 + strlen(Exon->Prefix);
                }
                else
                {
                    SRightGenomicPosition[MH_MatchCount] = Exon->Start + strlen(Exon->Prefix);
                }
            }
            else
            {
                if (Pos > -1)
                {
                    // Yes, still add prefix length here:
                    SRightGenomicPosition[MH_MatchCount] = Exon->End - (Pos+1)*3 - strlen(Exon->Prefix);
                }
                else
                {
                    SRightGenomicPosition[MH_MatchCount] = Exon->End - strlen(Exon->Prefix);
                }
            }
        }
        else
        {
            if (Exon->Gene->ForwardFlag)
            {
                SRightGenomicPosition[MH_MatchCount] = Exon->Start + strlen(Exon->Prefix); 
            }
            else
            {
                SRightGenomicPosition[MH_MatchCount] = Exon->End - strlen(Exon->Prefix);
            }
        }
        Length = Exon->Length;
        if (Pos + Direction < Length)
        {
            SRightSuffixes[MH_MatchCount] = Exon->Sequence[Pos + Direction];
            return;
        }
        // If we have a forward-edge, use the (most common) forward edge aa
        if (Exon->ForwardEdgeCount)
        {
            Edge = Exon->ForwardEdges;
            AA = Edge->AA;
            if (!AA && Edge->Exon->Sequence)
            {
                AA = Edge->Exon->Sequence[0];
            }
            SRightSuffixes[MH_MatchCount] = AA;
            return;
        }
        SRightSuffixes[MH_MatchCount] = '-';
        return;
    }
    else
    {
        // Direction is -1.  Set LeftGenomicPosition:
        if (Exon->Start < 0)
        {
            SLeftGenomicPosition[MH_MatchCount] = -1;
        }
        else if (Exon->Sequence)
        {
            if (Exon->Gene->ForwardFlag)
            {
                if (Pos >= 0)
                {
                    SLeftGenomicPosition[MH_MatchCount] = Exon->Start + Pos*3 + strlen(Exon->Prefix);
                }
                else
                {
                    // We never used any sequence from the exon proper, we used
                    // an incoming aa-edge:
                    SLeftGenomicPosition[MH_MatchCount] = Exon->End - strlen(Exon->Suffix); 
                }
            }
            else
            {
                if (Pos >= 0)
                {
                    // Yes, still add prefix length here:
                    SLeftGenomicPosition[MH_MatchCount] = Exon->End - Pos*3 - strlen(Exon->Prefix);
                }
                else
                {
                    SLeftGenomicPosition[MH_MatchCount] = Exon->Start + strlen(Exon->Suffix);
                }
            }
        }
        else
        {
            if (Exon->Gene->ForwardFlag)
            {
                SLeftGenomicPosition[MH_MatchCount] = Exon->End - strlen(Exon->Suffix); 
            }
            else
            {
                SLeftGenomicPosition[MH_MatchCount] = Exon->Start + strlen(Exon->Suffix); 
            }
        }

        if (Pos + Direction >= 0 && Exon->Sequence)
        {
            SLeftPrefixes[MH_MatchCount] = Exon->Sequence[Pos + Direction];
            return;
        }
        else if (Exon->BackEdgeCount)
        {
            Edge = Exon->BackwardEdges;
            AA = Edge->AA;
            if (!AA && Edge->Exon->Sequence)
            {
                Length = strlen(Edge->Exon->Sequence);
                AA = Edge->Exon->Sequence[Length-1];
            }
            SLeftPrefixes[MH_MatchCount] = AA;
            return;
        }
        SLeftPrefixes[MH_MatchCount] = '-';
        return;
    }
 }

// Copy the exon list from MH_MatchExons into either SLeftExon or SRightExon,
// initializing the left-over entries to NULL
void MatchHelperSetExons(int Direction)
{
    ExonStruct** MatchExons;
    ExonEdge** MatchEdges;
    int Index;
    if (Direction < 0)
    {
        MatchExons = SLeftExon;
        MatchEdges = SLeftEdge;
        //MatchSplices = SLeftSpliceScore;
    }
    else
    {
        MatchExons = SRightExon;
        MatchEdges = SRightEdge;
        //MatchSplices = SRightSpliceScore;
    }
    for (Index = 0; Index < MAX_EXTENSION_EXONS; Index++)
    {
        if (Index >= MH_MatchExonPos)
        {
            MatchExons[MH_MatchCount * MAX_EXTENSION_EXONS + Index] = NULL;    
            //MatchSplices[MH_MatchCount*MAX_EXTENSION_EXONS + Index] = -1;    
        }
        else
        {
            MatchExons[MH_MatchCount * MAX_EXTENSION_EXONS + Index] = MH_MatchExons[Index];
            //MatchSplices[MH_MatchCount*MAX_EXTENSION_EXONS + Index] = MH_MatchSplices[Index];
        }
        if (Index >= MH_MatchEdgePos)
        {
            MatchEdges[MH_MatchCount * MAX_EXTENSION_EXONS + Index] = NULL; 
        }
        else
        {
            MatchEdges[MH_MatchCount * MAX_EXTENSION_EXONS + Index] = MH_MatchEdges[Index]; 
        }
    }
}

// Recursion counter - tracks calls to MatchFlankingMassSpliceHelper, so that
// we can bail out if we take absurdly long.  (If we have a 3000Da+ flanking mass
// and a lot of SNPs, then the search time becomes unacceptable)
int g_SpliceHelperRecursionCount = 0;
// The largest count we ever saw before the limit was added: 1862528167
// Second-largest 816910, 99.99% were <30000
#define MAX_HELPER_RECURSION_COUNT 30000
// Recursable function for matching MatchMass.  We start out with decoration DecorationMatchIndex, and we try
// smaller decorations (with smaller index number) as we go.  We start out with FlankingMass = 0 on the first
// call; it's nonzero if we hit a splice junction and recurse.  
int MatchFlankingMassSpliceHelper(MSSpectrum* Spectrum, TrieTag* Tag, ExonStruct* Exon, 
    int StartPos, int Direction, int MatchMass, int ModsRemaining, 
    int DecorationMassIndex, int FlankingMass, int BufferPos, int BufferPosSpliced)
{
    int Pos;
    int AAMass;
    int Diff;
    int AbsDiff;
    int MandatoryDecorationChange = 0;
    int BridgeBufferPos;
    int BridgeDMI;
    ExonStruct* BridgeExon;
    int BridgeMass;
    //char* EdgeAA;
    int EdgeCount;
    //ExonStruct** EdgeExon;
    //int* EdgePower;
    ExonEdge* Edges;
    int EdgeIndex;
    int Length;
    int VerboseFlag = 0;
    int OldMatchExonPos;
    //int OldMatchSplicePos;
    int OldMatchEdgePos;
    int BridgeBufferPosSpliced;

    g_SpliceHelperRecursionCount++;

    //////////////////////////
    // StartPos < 0 if we're starting at the edge of the exon and working inward.
    if (StartPos < 0) 
    {
        if (Direction>0)
        {
            Pos = 0;
        }
        else
        {
            if (Exon->Sequence)
            {
                Pos = Exon->Length - 1;
            }
            else
            {
                Pos = -1;
            }
        }
        MH_MatchExons[MH_MatchExonPos] = Exon;
        MH_MatchExonPos++;
        if (MH_MatchExonPos >= MAX_EXTENSION_EXONS)
        {
            // Bail out!  We extended across too many exons!
            return 0;
        }
    }
    else
    {
        // The tag includes the character at StartPos, so move to the next character:
        Pos = StartPos + Direction;
    }
    Length = Exon->Length;
    
    // First, we'll extend out as far as possible WITHOUT bridging:
    while (1)
    {
        if (Pos < 0 || Pos >= Length)
        {
            break;
        }
        if (DecorationMassIndex < 0)
        {
            break;
        }
        AAMass = PeptideMass[Exon->Sequence[Pos]];
        if (!AAMass)
        {
            // We've reached a stop codon.
            DecorationMassIndex = -1;
            break;
        }
        FlankingMass += AAMass;
        MH_Buffer[BufferPos++] = Exon->Sequence[Pos];
        MH_BufferSpliced[BufferPosSpliced++] = Exon->Sequence[Pos];
        Diff = MatchMass  - (FlankingMass + AllDecorations[DecorationMassIndex].Mass);
        AbsDiff = abs(Diff);
        if (AbsDiff < GlobalOptions->FlankingMassEpsilon)
        {
            // Aha!  This is *probably* a match.  Check to be sure we have the bases we need:
            if (CheckForPTAttachmentPoints(DecorationMassIndex, MH_Buffer, 0, BufferPos - 1, 1))
            {
                if (VerboseFlag)
                {
                    printf("Side is match!  Dec-index %d, flank %.2f.\n", DecorationMassIndex, FlankingMass / (float)MASS_SCALE);
                    printf("Copy to buffer.  Match count is %d, bufferpos is %d\n", MH_MatchCount, BufferPos);
                }
                strncpy(MH_MatchBuffer + MAX_EXTENSION_LENGTH * MH_MatchCount, MH_Buffer, BufferPos);
                MH_MatchBuffer[MAX_EXTENSION_LENGTH * MH_MatchCount + BufferPos] = '\0';
                strncpy(MH_MatchBufferSpliced + MAX_SEXTENSION_LENGTH * MH_MatchCount, MH_BufferSpliced, BufferPosSpliced);
                MH_MatchBufferSpliced[MAX_SEXTENSION_LENGTH * MH_MatchCount + BufferPosSpliced] = '\0';
                // Set prefix or suffix for this extension:
                SetMatchPrefixSuffix(Exon, Pos, Direction);
                MH_MatchDecoration[MH_MatchCount] = DecorationMassIndex;
                MatchHelperSetExons(Direction);
                MH_MatchCount++;
                if (MH_MatchCount >= MAX_SIDE_EXTENSIONS)
                {
                    return MH_MatchCount;
                }
            }
        }
        // Move the DecorationMassIndex, if needed.
        while (MandatoryDecorationChange || FlankingMass + AllDecorations[DecorationMassIndex].Mass > MH_MinMatchMass)
        {
            // The flanking sequence's mass is significantly bigger than our (decorated) target mass.
            // Move to a smaller decoration:
            MandatoryDecorationChange = 0;
            DecorationMassIndex--;
            if (DecorationMassIndex<0)
            {
                break;
            }
            // Skip any decorations that include phosphorylation, if we're not on phospho mode.
            if (!GlobalOptions->PhosphorylationFlag && g_PhosphorylationMod > -1 && AllDecorations[DecorationMassIndex].Mods[g_PhosphorylationMod])
            {
                MandatoryDecorationChange = 1;
                continue;
            }
            if (AllDecorations[DecorationMassIndex].TotalMods > ModsRemaining)
            {
                continue;
            }
            // This decoration is acceptable.  Check for a match:
            Diff = MatchMass  - (FlankingMass + AllDecorations[DecorationMassIndex].Mass);
            AbsDiff = abs(Diff);
            if (AbsDiff < GlobalOptions->FlankingMassEpsilon) 
            {
                // Aha!  This is *probably* a match.  Check to be sure we have the bases we need:
                if (CheckForPTAttachmentPoints(DecorationMassIndex, MH_Buffer, 0, BufferPos-1, 1))
                {
                    if (VerboseFlag)
                    {
                        printf("Left is match!  Dec-index %d, flank %.2f.\n", DecorationMassIndex, FlankingMass / (float)MASS_SCALE);
                    }
                    strncpy(MH_MatchBuffer + MAX_EXTENSION_LENGTH * MH_MatchCount, MH_Buffer, BufferPos);
                    MH_MatchBuffer[MAX_EXTENSION_LENGTH * MH_MatchCount + BufferPos] = '\0';
                    strncpy(MH_MatchBufferSpliced + MAX_SEXTENSION_LENGTH * MH_MatchCount, MH_BufferSpliced, BufferPosSpliced);
                    MH_MatchBufferSpliced[MAX_SEXTENSION_LENGTH * MH_MatchCount + BufferPosSpliced] = '\0';

                    // Set prefix or suffix for this extension:
                    SetMatchPrefixSuffix(Exon, Pos, Direction);
                    MH_MatchDecoration[MH_MatchCount] = DecorationMassIndex;
                    MatchHelperSetExons(Direction);
                    MH_MatchCount++;
                    if (MH_MatchCount >= MAX_SIDE_EXTENSIONS)
                    {
                        return MH_MatchCount;
                    }
                    MandatoryDecorationChange = 1;
                }
            }
        }
        Pos += Direction;
    }

    // If DMI < 0, then our flanking mass became too large or we hit a stop codon:
    if (DecorationMassIndex<0)
    {
        return MH_MatchCount;
    }

    // Now: We reached the end of the exon, so next we'll try each edge:
    if (Direction > 0)
    {
        Edges = Exon->ForwardEdges;
        EdgeCount = Exon->ForwardEdgeCount;
    }
    else
    {
        Edges = Exon->BackwardEdges;
        EdgeCount = Exon->BackEdgeCount;
    }
    // Save our current state (FlankingMass, BufferPos, and DecorationMassIndex).  After trying each edge,
    // we return to this state.
    BridgeMass = FlankingMass;
    BridgeBufferPos = BufferPos;
    BridgeBufferPosSpliced = BufferPosSpliced;
    BridgeDMI = DecorationMassIndex;
    OldMatchExonPos = MH_MatchExonPos;
    OldMatchEdgePos = MH_MatchEdgePos;

    for (EdgeIndex = 0; EdgeIndex < EdgeCount; EdgeIndex++)
    {
        FlankingMass = BridgeMass;
        BufferPos = BridgeBufferPos;
        BufferPosSpliced = BridgeBufferPosSpliced;
        DecorationMassIndex = BridgeDMI;
        MH_MatchExonPos = OldMatchExonPos;
        MH_MatchEdgePos = OldMatchEdgePos;
        
        BridgeExon = Edges[EdgeIndex].Exon;  
        MH_MatchEdges[MH_MatchEdgePos] = Edges + EdgeIndex;
        MH_MatchEdgePos++;
        // Extend with the edge amino acid:
        if (Edges[EdgeIndex].AA)
        {
            AAMass = PeptideMass[Edges[EdgeIndex].AA];
            if (!AAMass)
            {
                continue; // terminator
            }
            FlankingMass += AAMass;
            // If this is a "true edge" (not an adjacent-edge), then note the splicing:
            if (Edges[EdgeIndex].Power)
            {
                MH_BufferSpliced[BufferPosSpliced++] = ';';
                MH_BufferSpliced[BufferPosSpliced++] = Edges[EdgeIndex].AA;
                MH_BufferSpliced[BufferPosSpliced++] = ';';
            }
            else
            {
                MH_BufferSpliced[BufferPosSpliced++] = Edges[EdgeIndex].AA;
            }
            MH_MatchExons[MH_MatchExonPos] = BridgeExon;
            MH_MatchExonPos++;
            if (MH_MatchExonPos >= MAX_EXTENSION_EXONS)
            {
                // Bail out!  We extended across too many exons!
                MH_MatchExonPos--;
                continue;
            }

            MH_Buffer[BufferPos++] = Edges[EdgeIndex].AA; //EdgeAA[EdgeIndex];
            Diff = MatchMass  - (FlankingMass + AllDecorations[DecorationMassIndex].Mass);
            AbsDiff = abs(Diff);
            if (AbsDiff < GlobalOptions->FlankingMassEpsilon)
            {
                // Aha!  This is *probably* a match.  Check to be sure we have the bases we need:
                if (CheckForPTAttachmentPoints(DecorationMassIndex, MH_Buffer, 0, BufferPos-1, 1))
                {
                    if (VerboseFlag)
                    {
                        printf("Side is match!  Dec-index %d, flank %.2f.\n", DecorationMassIndex, FlankingMass / (float)MASS_SCALE);
                    }
                    strncpy(MH_MatchBuffer + MAX_EXTENSION_LENGTH * MH_MatchCount, MH_Buffer, BufferPos);
                    MH_MatchBuffer[MAX_EXTENSION_LENGTH * MH_MatchCount + BufferPos] = '\0';
                    strncpy(MH_MatchBufferSpliced + MAX_SEXTENSION_LENGTH * MH_MatchCount, MH_BufferSpliced, BufferPosSpliced);
                    MH_MatchBufferSpliced[MAX_SEXTENSION_LENGTH * MH_MatchCount + BufferPosSpliced] = '\0';
                    // Set prefix or suffix for this extension:
                    if (Direction > 0)
                    {
                        // Direction > 0: set suffix!
                        if (BridgeExon->Sequence)
                        {
                            SRightSuffixes[MH_MatchCount] = BridgeExon->Sequence[0];
                        }
                        else
                        {
                            SRightSuffixes[MH_MatchCount] = '-';
                        }
                        if (Exon->Start < 0)
                        {
                            SRightGenomicPosition[MH_MatchCount] = -1;
                        }
                        else if (Exon->Gene->ForwardFlag)
                        {
                            SRightGenomicPosition[MH_MatchCount] = BridgeExon->Start + strlen(BridgeExon->Prefix);
                        }
                        else
                        {
                            SRightGenomicPosition[MH_MatchCount] = BridgeExon->End - strlen(BridgeExon->Prefix);
                        }

                    }
                    else
                    {
                        // Direction < 0: set prefix!
                        if (BridgeExon->Sequence)
                        {
                            if (strlen(BridgeExon->Sequence) < 1)
                            {
                                MH_MatchCount = MH_MatchCount;
                            }
                            SLeftPrefixes[MH_MatchCount] = BridgeExon->Sequence[strlen(BridgeExon->Sequence)-1];
                        }
                        else
                        {
                            SLeftPrefixes[MH_MatchCount] = '-';
                        }
                        if (Exon->Start < 0)
                        {
                            SLeftGenomicPosition[MH_MatchCount] = -1;
                        }
                        else if (Exon->Gene->ForwardFlag)
                        {
                            SLeftGenomicPosition[MH_MatchCount] = BridgeExon->End - strlen(BridgeExon->Suffix);
                        }
                        else
                        {
                            SLeftGenomicPosition[MH_MatchCount] = BridgeExon->Start + strlen(BridgeExon->Suffix);
                        }
                    }
                    MH_MatchDecoration[MH_MatchCount] = DecorationMassIndex;
                    MatchHelperSetExons(Direction);
                    MH_MatchCount++;
                    if (MH_MatchCount >= MAX_SIDE_EXTENSIONS)
                    {
                        return MH_MatchCount;
                    }
                }
            }
            // Move the DecorationMassIndex, if needed.
            while (MandatoryDecorationChange || FlankingMass + AllDecorations[DecorationMassIndex].Mass > MH_MinMatchMass)
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
                    if (CheckForPTAttachmentPoints(DecorationMassIndex, MH_Buffer, 0, BufferPos-1, 1))
                    {
                        if (VerboseFlag)
                        {
                            printf("Left is match!  Dec-index %d, flank %.2f.\n", DecorationMassIndex, FlankingMass / (float)MASS_SCALE);
                        }
                        strncpy(MH_MatchBuffer + MAX_EXTENSION_LENGTH*MH_MatchCount, MH_Buffer, BufferPos);
                        MH_MatchBuffer[MAX_EXTENSION_LENGTH * MH_MatchCount + BufferPos] = '\0';
                        strncpy(MH_MatchBufferSpliced + MAX_SEXTENSION_LENGTH * MH_MatchCount, MH_BufferSpliced, BufferPosSpliced);
                        MH_MatchBufferSpliced[MAX_SEXTENSION_LENGTH * MH_MatchCount + BufferPosSpliced] = '\0';
                        // Set prefix or suffix for this extension:
                        if (Direction > 0)
                        {
                            if (BridgeExon->Sequence)
                            {
                                SRightSuffixes[MH_MatchCount] = BridgeExon->Sequence[0];
                            }
                            else
                            {
                                SRightSuffixes[MH_MatchCount] = '-';
                            }
                            if (BridgeExon->Start < 0)
                            {
                                SRightGenomicPosition[MH_MatchCount] = -1;
                            }
                            else if (Exon->Gene->ForwardFlag)
                            {
                                SRightGenomicPosition[MH_MatchCount] = BridgeExon->Start + strlen(BridgeExon->Prefix);
                            }
                            else
                            {
                                SRightGenomicPosition[MH_MatchCount] = BridgeExon->End - strlen(BridgeExon->Prefix);
                            }
                        }
                        else
                        {
                            if (BridgeExon->Sequence)
                            {
                                SLeftPrefixes[MH_MatchCount] = BridgeExon->Sequence[strlen(BridgeExon->Sequence)-1];
                            }
                            else
                            {
                                SLeftPrefixes[MH_MatchCount] = '-';
                            }
                            if (BridgeExon->Start < 0)
                            {
                                SLeftGenomicPosition[MH_MatchCount] = -1;
                            }
                            else if (Exon->Gene->ForwardFlag)
                            {
                                SLeftGenomicPosition[MH_MatchCount] = BridgeExon->End - strlen(BridgeExon->Suffix);
                            }
                            else
                            {
                                SLeftGenomicPosition[MH_MatchCount] = BridgeExon->Start + strlen(BridgeExon->Suffix);
                            }
                        }

                        MH_MatchDecoration[MH_MatchCount] = DecorationMassIndex;
                        MatchHelperSetExons(Direction);
                        MH_MatchCount++;
                        if (MH_MatchCount >= MAX_SIDE_EXTENSIONS)
                        {
                            return MH_MatchCount;
                        }
                        MandatoryDecorationChange = 1;
                    }
                }
            }
            MH_MatchExonPos--;
        } // If the edge has an AA
        else
        {
            if (Edges[EdgeIndex].Power)
            {
                MH_BufferSpliced[BufferPosSpliced++] = ':';
            }
        }

        // Recurse!  Call MatchFlankingMassSpliceHelper again:
        MatchFlankingMassSpliceHelper(Spectrum,  Tag, BridgeExon, -1, Direction, MatchMass, ModsRemaining,
            DecorationMassIndex, FlankingMass, BufferPos, BufferPosSpliced);
        if (MH_MatchCount >= MAX_SIDE_EXTENSIONS)
        {
            return MH_MatchCount;
        }
        if (g_SpliceHelperRecursionCount >= MAX_HELPER_RECURSION_COUNT)
        {
            return MH_MatchCount;
        }
    }  // Iteration over edges
    return MH_MatchCount;
}

int MatchFlankingMassSpliced(MSSpectrum* Spectrum, TrieTag* Tag, ExonStruct* Exon, int StartPos, int Direction, 
    int MatchMass, int ModsRemaining)
{
    static int DecorationMassIndex;
    static int AAMass;
    //
    /////////////////////////////////////////////////////////
    // If prefix mass is zero, that qualifies as a match always.
    MH_MatchCount = 0;
    if (MatchMass < GlobalOptions->FlankingMassEpsilon) 
    {  
        if (Direction < 0)
        {
            SLeftMatchDecoration[0] = PlainOldDecorationIndex;
            SLeftExon[0] = NULL;
            SLeftEdge[0] = NULL;
            SLeftMatchBuffer[0] = '\0';
            SLeftMatchBufferSpliced[0] = '\0';
            SetMatchPrefixSuffix(Exon, StartPos, Direction);
        }
        else
        {
            SRightMatchDecoration[0] = PlainOldDecorationIndex;
            SRightExon[0] = NULL;
            SRightEdge[0] = NULL;
            SRightMatchBuffer[0] = '\0';
            SRightMatchBufferSpliced[0] = '\0';
            SetMatchPrefixSuffix(Exon, StartPos, Direction);
        }
        return 1;
    }

    MH_MinMatchMass = MatchMass - GlobalOptions->FlankingMassEpsilon;
    MH_MaxMatchMass = MatchMass + GlobalOptions->FlankingMassEpsilon;
    if (Direction < 0)
    {
        MH_MatchBuffer = SLeftMatchBuffer;
        MH_MatchBufferSpliced = SLeftMatchBufferSpliced;
        MH_MatchDecoration = SLeftMatchDecoration;
        MH_MatchExonPos = 0;
        //MH_MatchSplices = SLeftSpliceScore;
        //MH_MatchSplicePos = 0;
        MH_MatchEdgePos = 0;
    }
    else
    {
        MH_MatchBuffer = SRightMatchBuffer;
        MH_MatchBufferSpliced = SRightMatchBufferSpliced;
        MH_MatchDecoration = SRightMatchDecoration;
        MH_MatchExonPos = 0;
        //MH_MatchSplices = SRightSpliceScore;
        //MH_MatchSplicePos = 0;
        MH_MatchEdgePos = 0;
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

    MH_MatchExonPos = 0;
    MH_MatchEdgePos = 0;
    // Perform tag extension, following edges as needed:
    g_SpliceHelperRecursionCount = 0;
    MatchFlankingMassSpliceHelper(Spectrum, Tag, Exon, StartPos, Direction,
        MatchMass, ModsRemaining, DecorationMassIndex, 0, 0, 0);
    return MH_MatchCount;

}

// Copies a string to a destination, in reverse character order.
void ReverseStringCopy(char* Target, char* Source)
{
    int Length;
    char* SourceChar;
    //
    Length = strlen(Source);
    for (SourceChar = Source + Length - 1; SourceChar >= Source; SourceChar--)
    {
        *Target = *SourceChar;
        Target++;
    }

}

#define MINIMUM_EXON_LENGTH 4

// A tag has been matched.  Its left edge lies at LeftExonPos in LeftExon (or at -1 if its leftmost character
// comes from an edge).  Its right edge lies at RightExonPos in RightExon (or at -1 if its rightmost character
// comes from an edge).  Try to extend out to a prefix/suffix mass match.  Analogous to the GetMatches() function
// in standard trie search.  The difference is that our extension can follow an exon edge.
void GetSplicedMatches(SearchInfo* Info, TrieNode* Node, ExonStruct* LeftExon, int LeftExonPos, 
    ExonStruct* RightExon, int RightExonPos)
{
    int LeftMatchCount;
    int RightMatchCount;
    int LeftMatchIndex;
    int RightMatchIndex;
    int ModIndex;
    int Length;
    int ModsRemaining;
    int Pos;
    Peptide* Match;
    int VerboseFlag = 0;
    int ForwardFlag;
    MSSpectrum* Spectrum;
    static int PTMLimit[MAX_PT_MODTYPE];
    TrieTagHanger* TagNode;
    int ExtensionIndex;
    int ExtensionCount = 0;
    int ExtensionFound;
    int UsedTooMany;
    int ExIndex;
    int ExonCount;
    int SpliceScoreCount;
    ExonStruct* TempExon;
    ExonStruct* AllExons[256];
    ExonEdge* AllEdges[256];
    char SplicedBases[256];
    int AllEdgeCount;
    int AllExonCount;
    PeptideSpliceNode* SpliceTail;
    PeptideSpliceNode* SpliceNode;
    //int GenomicLocation;
    PeptideSpliceNode* PrevSpliceNode;
    int EdgeIndex;
    ExonEdge* TempEdge;
    int GenomicStart;
    int GenomicEnd;
    char* ShortExonCheck;
    int DistanceFromLastJunction;
    int InvalidExonFlag;
    //////////////
    //printf("GetSplicedMatches() called for tag %s\n", Node->FirstTag->Tag->Tag); 
    if (!Node->FirstTag)
    {
        return;
    }
    ForwardFlag = LeftExon->Gene->ForwardFlag;
    for (TagNode = Node->FirstTag; TagNode; TagNode = TagNode->Next)
    {
        if (VerboseFlag)
        {
            printf("Matched tag '%s' (pre %.2f post %.2f).\n  Left exon %d pos %d, right exon %d pos %d\n",
                TagNode->Tag->Tag, TagNode->Tag->PrefixMass / (float)MASS_SCALE, TagNode->Tag->SuffixMass / (float)MASS_SCALE,
                LeftExon->Start, LeftExonPos, RightExon->Start, RightExonPos);
        }
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
        LeftMatchCount = MatchFlankingMassSpliced(Spectrum, TagNode->Tag, LeftExon, LeftExonPos, -1, TagNode->Tag->PrefixMass, ModsRemaining);
        if (LeftMatchCount == 0)
        {
            continue;
        }
        RightMatchCount = MatchFlankingMassSpliced(Spectrum, TagNode->Tag, RightExon, RightExonPos, 1, TagNode->Tag->SuffixMass, ModsRemaining);
        if (RightMatchCount == 0)
        {
            continue;
        }

        // Consider each combination of left-decoration and right-decoration:
        for (LeftMatchIndex = 0; LeftMatchIndex < LeftMatchCount; LeftMatchIndex++)
        {
            for (RightMatchIndex = 0; RightMatchIndex < RightMatchCount; RightMatchIndex++)
            {
                if (VerboseFlag)
                {
                    printf("LMI %d RMI %d Count %d\n", LeftMatchIndex, RightMatchIndex, ExtensionCount);
                }
                UsedTooMany = 0;
                for (ModIndex = 0; ModIndex < AllPTModCount; ModIndex++)
                {

                    if (AllDecorations[SLeftMatchDecoration[LeftMatchIndex]].Mods[ModIndex] + 
                        AllDecorations[SRightMatchDecoration[RightMatchIndex]].Mods[ModIndex] > PTMLimit[ModIndex])
                    {
                        UsedTooMany = 1;
                        break;
                    }
                }
                if (UsedTooMany)
                {
                    continue;
                }
                if (AllDecorations[SLeftMatchDecoration[LeftMatchIndex]].TotalMods + 
                    AllDecorations[SRightMatchDecoration[RightMatchIndex]].TotalMods > ModsRemaining)
                {
                    continue;
                }
                if (GlobalOptions->MandatoryModIndex > -1 && 
                    !TagNode->Tag->MandatoryModUsed &&
                    AllDecorations[SLeftMatchDecoration[LeftMatchIndex]].Mods[GlobalOptions->MandatoryModIndex] == 0 &&
                    AllDecorations[SRightMatchDecoration[RightMatchIndex]].Mods[GlobalOptions->MandatoryModIndex] == 0)
                {
                    continue; // We don't have our mandatory PTM (biotin, or whatever)
                }
                if (LeftExon->Gene->ForwardFlag)
                {
                    GenomicStart = SLeftGenomicPosition[LeftMatchIndex];
                    GenomicEnd = SRightGenomicPosition[RightMatchIndex];
                }
                else
                {
                    GenomicStart = SRightGenomicPosition[RightMatchIndex];
                    GenomicEnd = SLeftGenomicPosition[LeftMatchIndex];
                }               
                // Don't produce the same extension multiple times:
                ExtensionFound = 0;
                for (ExtensionIndex = 0; ExtensionIndex < ExtensionCount; ExtensionIndex++)
                {
                    if (!strcmp(ExtensionBufferLeft + ExtensionIndex*MAX_EXTENSION_LENGTH, SLeftMatchBuffer + LeftMatchIndex*MAX_EXTENSION_LENGTH)
                    && !strcmp(ExtensionBufferRight + ExtensionIndex*MAX_EXTENSION_LENGTH, SRightMatchBuffer + RightMatchIndex*MAX_EXTENSION_LENGTH)
                    && ExtensionLeftDecorations[ExtensionIndex] == SLeftMatchDecoration[LeftMatchIndex]
                    && ExtensionRightDecorations[ExtensionIndex] == SRightMatchDecoration[RightMatchIndex]
                    && ExtensionSpectra[ExtensionIndex] == TagNode->Tag->PSpectrum)
                    {
                        // Gosh, looks like we found the same peptide again (probably by starting with
                        // another valid tag).  Let's check whether the genomic endpoints are the
                        // same as well:
                        if (GenomicStart == ExtensionGenomicStart[ExtensionIndex] && GenomicEnd == ExtensionGenomicEnd[ExtensionIndex])
                        {
                            ExtensionFound = 1;
                            break;
                        }
                    }
                }
                if (ExtensionFound)
                {
                    continue;
                }
                ExtensionLeftDecorations[ExtensionCount] = SLeftMatchDecoration[LeftMatchIndex];
                ExtensionRightDecorations[ExtensionCount] = SRightMatchDecoration[RightMatchIndex];
                strcpy(ExtensionBufferLeft + ExtensionCount * MAX_EXTENSION_LENGTH, SLeftMatchBuffer + LeftMatchIndex * MAX_EXTENSION_LENGTH);
                strcpy(ExtensionBufferRight + ExtensionCount * MAX_EXTENSION_LENGTH, SRightMatchBuffer + RightMatchIndex * MAX_EXTENSION_LENGTH);
                ExtensionSpectra[ExtensionCount] = TagNode->Tag->PSpectrum;

                // MatchedBases is concatenated together from five sources:
                // prefixAA
                //    |        TAG
                //    A BBBBBB CCC DDDDDDD E
                //      left       right   |
                //       ext.       ext.   suffix
                Pos = strlen(SLeftMatchBuffer + MAX_EXTENSION_LENGTH * LeftMatchIndex);
                MatchedBases[0] = SLeftPrefixes[LeftMatchIndex];
                ReverseStringCopy(MatchedBases + 1, SLeftMatchBuffer + MAX_EXTENSION_LENGTH*LeftMatchIndex);
                g_TagBuffer[g_TagBufferPos] = '\0';
                strcpy(MatchedBases + 1 + Pos, g_TagBuffer);
                //strcpy(MatchedBases + 1 + Pos, TagNode->Tag->Tag);
                strcpy(MatchedBases + 1 + Pos + strlen(TagNode->Tag->Tag), SRightMatchBuffer + MAX_EXTENSION_LENGTH*RightMatchIndex);
                Length = strlen(MatchedBases+1);
                MatchedBases[strlen(MatchedBases+1)+2] = '\0';
                MatchedBases[strlen(MatchedBases+1)+1] = SRightSuffixes[RightMatchIndex];

                // Set SplicedBases, and check for unacceptably short exons:
                Pos = strlen(SLeftMatchBufferSpliced + MAX_SEXTENSION_LENGTH * LeftMatchIndex);
                ReverseStringCopy(SplicedBases, SLeftMatchBufferSpliced + MAX_SEXTENSION_LENGTH * LeftMatchIndex);
                g_TagBufferSpliced[g_TagBufferPosSpliced] = '\0';
                strcpy(SplicedBases + Pos, g_TagBufferSpliced);
                strcpy(SplicedBases + Pos + strlen(g_TagBufferSpliced), SRightMatchBufferSpliced + MAX_SEXTENSION_LENGTH * RightMatchIndex);
                DistanceFromLastJunction = 999;
                InvalidExonFlag = 0;
                for (ShortExonCheck = SplicedBases; *ShortExonCheck; ShortExonCheck++)
                {
                    switch (*ShortExonCheck)
                    {
                    case ';':
                        if (DistanceFromLastJunction < MINIMUM_EXON_LENGTH)
                        {
                            InvalidExonFlag = 1;
                        }
                        else
                        {
                            ShortExonCheck += 2; // We're at the start of ;x;, skip over aa and other ;
                            DistanceFromLastJunction = 0;
                        }
                        break;
                    case ':':
                        if (DistanceFromLastJunction < MINIMUM_EXON_LENGTH)
                        {
                            InvalidExonFlag = 1;
                        }
                        DistanceFromLastJunction = 0;
                        break;
                    default:
                        DistanceFromLastJunction++;
                        break;
                    }
                    if (InvalidExonFlag)
                    {
                        break;
                    }
                }
                // Reject, if unacceptably short exons were used:
                if (InvalidExonFlag)
                {
                    continue;
                }

                ExtensionGenomicStart[ExtensionCount] = GenomicStart;
                ExtensionGenomicEnd[ExtensionCount] = GenomicEnd;
                Match = AddNewMatch(Info, -1, TagNode->Tag, 
                    MatchedBases + 1, Length,
                    strlen(SLeftMatchBuffer + MAX_EXTENSION_LENGTH * LeftMatchIndex),
                    SLeftMatchDecoration[LeftMatchIndex], SRightMatchDecoration[RightMatchIndex],
                    GenomicStart, GenomicEnd);

                if (Match)
                {
                    // We might have some splice nodes stored here.  If so, free them:
                    if (Match->SpliceHead)
                    {
                        PrevSpliceNode = NULL;
                        for (SpliceNode = Match->SpliceHead; SpliceNode; SpliceNode = SpliceNode->Next)
                        {
                            SafeFree(PrevSpliceNode);
                            PrevSpliceNode = SpliceNode;
                        }
                        SafeFree(PrevSpliceNode);
                        Match->SpliceHead = NULL;
                    }
                    //Match->GenomicLocation = GenomicLocation;

                    Match->ChromosomeNumber = LeftExon->Gene->ChromosomeNumber;
                    Match->ChromosomeForwardFlag = LeftExon->Gene->ForwardFlag;
                    Match->RecordNumber = Info->RecordNumber;
                    // Copy in the list of exons and the splice scores:
                    ExonCount = 0;
                    SpliceScoreCount = 0;
                    AllExonCount = 0;
                    AllEdgeCount = 0;
                    // Read exons from the prefix:
                    for (ExIndex = 0; ExIndex < MAX_EXTENSION_EXONS; ExIndex++)
                    {
                        TempExon = SLeftExon[LeftMatchIndex * MAX_EXTENSION_EXONS + ExIndex];
                        if (!TempExon)
                        {
                            ExIndex--;
                            break;
                        }
                    }
                    while (ExIndex > -1)
                    {
                        AllExons[AllExonCount] = SLeftExon[LeftMatchIndex * MAX_EXTENSION_EXONS + ExIndex];
                        AllExonCount++;
                        ExIndex--;
                    }
                    // Read edges from the prefix:
                    for (ExIndex = 0; ExIndex < MAX_EXTENSION_EXONS; ExIndex++)
                    {
                        TempEdge = SLeftEdge[LeftMatchIndex * MAX_EXTENSION_EXONS + ExIndex];
                        if (!TempEdge)
                        {
                            ExIndex--;
                            break;
                        }
                    }
                    while (ExIndex > -1)
                    {
                        AllEdges[AllEdgeCount] = GetReciprocalExonEdge(SLeftEdge[LeftMatchIndex * MAX_EXTENSION_EXONS + ExIndex], 0);
                        AllEdgeCount++;
                        ExIndex--;
                    }
                    // Read exons from the tag:
                    for (ExIndex = 0; ExIndex < g_TagExonArrayPos; ExIndex++)
                    {
                        if (AllExonCount && (AllExons[AllExonCount-1] == g_TagExonArray[ExIndex]))
                        {
                            continue;
                        }
                        AllExons[AllExonCount] = g_TagExonArray[ExIndex];
                        AllExonCount++;
                    }
                    // Read edges from the tag:
                    for (ExIndex = 0; ExIndex < g_TagSpliceArrayPos; ExIndex++)
                    {
                        AllEdges[AllEdgeCount] = g_TagSpliceArray[ExIndex];
                        AllEdgeCount++;
                    }
                    // Read exons from the suffix:
                    for (ExIndex = 0; ExIndex < MAX_EXTENSION_EXONS; ExIndex++)
                    {
                        TempExon = SRightExon[RightMatchIndex * MAX_EXTENSION_EXONS + ExIndex];
                        if (TempExon)
                        {
                            if (AllExonCount && (AllExons[AllExonCount-1] == TempExon))
                            {
                                continue;
                            }
                            AllExons[AllExonCount] = TempExon;
                            AllExonCount++;
                        }
                        else
                        {
                            break; // After the first null exon comes undefined rubbish data
                        }
                    }
                    // Read edges from the suffix:
                    for (ExIndex = 0; ExIndex < MAX_EXTENSION_EXONS; ExIndex++)
                    {
                        TempEdge = SRightEdge[RightMatchIndex * MAX_EXTENSION_EXONS + ExIndex];
                        if (TempEdge)
                        {
                            AllEdges[AllEdgeCount] = TempEdge;
                            AllEdgeCount++;
                        }
                        else
                        {
                            break; // After the first null exon comes undefined rubbish data
                        }
                    }
                    // Store the sequence, with splice boundaries indicated:
                    SafeFree(Match->SplicedBases);
                    Match->SplicedBases = (char*)calloc(sizeof(char), 256);
                    strncpy(Match->SplicedBases, SplicedBases, 256);
                    
                    // We know the exons, now we'll store all the splicing info for the match:
                    SpliceTail = NULL;
                    for (EdgeIndex = 0; EdgeIndex < AllEdgeCount; EdgeIndex++)
                    {
                        if (AllEdges[EdgeIndex]->Power)
                        {
                            SpliceNode = (PeptideSpliceNode*)calloc(sizeof(PeptideSpliceNode), 1);
                            if (ForwardFlag)
                            {
                                SpliceNode->DonorPos = AllEdges[EdgeIndex]->Source->End;
                                SpliceNode->AcceptorPos = AllEdges[EdgeIndex]->Exon->Start;
                            }
                            else
                            {
                                SpliceNode->DonorPos = AllEdges[EdgeIndex]->Source->Start;
                                SpliceNode->AcceptorPos = AllEdges[EdgeIndex]->Exon->End;
                            }
                            SpliceNode->ChromosomeNumber = LeftExon->Gene->ChromosomeNumber;
                            if (SpliceTail)
                            {
                                SpliceTail->Next = SpliceNode;
                            }
                            else
                            {
                                Match->SpliceHead = SpliceNode;
                            }
                            SpliceTail = SpliceNode;
                        }
                    }

                    //// %%% SANITY CHECK SPLICING %%%
                    //if (Match->SpliceHead && (!strstr(Match->SplicedBases, ";") && !strstr(Match->SplicedBases, ":")))
                    //{
                    //    printf("Warning: Match found with no true splicing, but splice junction stored!\n");
                    //    printf("%s %s\n", Match->Bases, Match->SplicedBases);
                    //    printf("SpliceNode: %d-%d\n", Match->SpliceHead->DonorPos, Match->SpliceHead->AcceptorPos);
                    //    DebugPrintGene(LeftExon->Gene);
                    //}
                    //if (!Match->SpliceHead && (strstr(Match->SplicedBases, ";") || strstr(Match->SplicedBases, ":")))
                    //{
                    //    printf("Warning: Match found with true splicing, but splice junction not stored!\n");
                    //    printf("%s %s\n", Match->Bases, Match->SplicedBases);
                    //    //printf("SpliceNode: %d-%d\n", Match->SpliceHead->DonorPos, Match->SpliceHead->AcceptorPos);
                    //    DebugPrintGene(LeftExon->Gene);
                    //}
                } // if match
                ExtensionCount = min(511, ExtensionCount + 1);
            } // RightMatchIndex
        } // LeftMatchIndex
    } // Tag loop
    return;
}

// Integrity checking of a gene.  (For debugging use only)
void CheckGene(GeneStruct* Gene)
{
    int ExonIndex;
    int EdgeIndex;
    ExonStruct* Exon;
    ExonStruct* Exon2;
    //
    if (!Gene)
    {
        return;
    }
    for (ExonIndex = 0; ExonIndex < Gene->ExonCount; ExonIndex++)
    {
        Exon = Gene->Exons + ExonIndex;
        if (Exon->Start < 0 || Exon->End < 0 || Exon->Start >= Exon->End)
        {
            printf("*ERROR\n");
        }
        for (EdgeIndex = 0; EdgeIndex < Exon->BackEdgeCount; EdgeIndex++)
        {
            Exon2 = Exon->BackwardEdges[EdgeIndex].Exon;
            if (!Exon2)
            {
                printf("*ERROR!\n");
            }
        }
        for (EdgeIndex = 0; EdgeIndex < Exon->ForwardEdgeCount; EdgeIndex++)
        {
            Exon2 = Exon->ForwardEdges[EdgeIndex].Exon;
            if (!Exon2)
            {
                printf("*ERROR!\n");
            }
        }        
    }
}

// Given an exon, search it for tag matches.  If AnchoredFlag is true, then we've already
// matched part of the tag (and Root isn't the root of the entire trie)
void GetSplicedTagMatches(SearchInfo* Info, ExonStruct* LeftExon, int LeftExonPos, ExonStruct* Exon, 
    TrieNode* Root, int AnchoredFlag)
{
    int AnchorMax;
    int AnchorPos;
    int OldExonPos;
    int OldSplicePos;
    char AA;
    ExonStruct* BridgedExon;
    TrieNode* CurrentNode;
    TrieNode* SubNode;
    int EdgeIndex;
    int SequenceLength;
    int SequencePos;
    int OldTagBufferPos;
    int OldTagBufferPos2;
    int OldTagBufferPosSpliced;
    int OldTagBufferPosSpliced2;
    //
    OldExonPos = g_TagExonArrayPos;
    OldSplicePos = g_TagSpliceArrayPos;
    OldTagBufferPos = g_TagBufferPos;
    OldTagBufferPosSpliced = g_TagBufferPosSpliced;
    SequenceLength = Exon->Length;
    if (AnchoredFlag)
    {
        AnchorMax = min(1, SequenceLength); // it's possible that sequencelength is 0!
    }
    else
    {
        AnchorMax = SequenceLength;
    }
    //printf("AnchoredFlag %d sequencelen %d anchormax %d\n", AnchoredFlag, SequenceLength, AnchorMax);
    for (AnchorPos = 0; AnchorPos < AnchorMax; AnchorPos++)
    {
        CurrentNode = Root->Children[Exon->Sequence[AnchorPos] - 'A'];
        if (!CurrentNode)
        {
            continue;
        }

        SequencePos = AnchorPos;
        g_TagBufferPos = OldTagBufferPos;
        g_TagBufferPosSpliced = OldTagBufferPosSpliced;
        g_TagBuffer[g_TagBufferPos++] = Exon->Sequence[AnchorPos];
        g_TagBufferSpliced[g_TagBufferPosSpliced++] = Exon->Sequence[AnchorPos];

        // If we're performing a tagless search, then our tag may have length 1, 
        // so we could get matches right now:
        if (CurrentNode->FirstTag)
        {
            if (AnchoredFlag)
            {
                GetSplicedMatches(Info, CurrentNode, LeftExon, LeftExonPos, Exon, AnchorPos);
            }
            else
            {
                GetSplicedMatches(Info, CurrentNode, Exon, AnchorPos, Exon, AnchorPos);
            }
        }

        while (1)
        {
            SequencePos++;
            //printf("Exon %d anchor %d sequence pos %d\n", Exon->Index, AnchorPos, SequencePos);
            if (SequencePos >= SequenceLength)
            {
                // Try to follow any edges forward
                OldTagBufferPos2 = g_TagBufferPos;
                OldTagBufferPosSpliced2 = g_TagBufferPosSpliced;
                for (EdgeIndex = 0; EdgeIndex < Exon->ForwardEdgeCount; EdgeIndex++)
                {
                    g_TagBufferPos = OldTagBufferPos2;
                    g_TagBufferPosSpliced = OldTagBufferPosSpliced2;
                    AA = Exon->ForwardEdges[EdgeIndex].AA;
                    if (AA)
                    {
                        SubNode = CurrentNode->Children[AA-'A'];
                        g_TagBuffer[g_TagBufferPos++] = AA;
                        if (Exon->ForwardEdges[EdgeIndex].Power)
                        {
                            g_TagBufferSpliced[g_TagBufferPosSpliced++] = ';';
                            g_TagBufferSpliced[g_TagBufferPosSpliced++] = AA;
                            g_TagBufferSpliced[g_TagBufferPosSpliced++] = ';';
                        }
                        else
                        {
                            g_TagBufferSpliced[g_TagBufferPosSpliced++] = AA;
                        }
                    }
                    else
                    {
                        SubNode = CurrentNode;
                        if (Exon->ForwardEdges[EdgeIndex].Power > 0)
                        {
                            g_TagBufferSpliced[g_TagBufferPosSpliced++] = ':';
                        }

                    }
                    if (!SubNode)
                    {
                        continue;
                    }
                    BridgedExon = Exon->ForwardEdges[EdgeIndex].Exon;
                    if (AA)
                    {
                        g_TagExonArray[g_TagExonArrayPos++] = BridgedExon;
                        g_TagSpliceArray[g_TagSpliceArrayPos++] = Exon->ForwardEdges + EdgeIndex;  //Exon->ForwardEdgePower[EdgeIndex];
                    }
                    if (SubNode->FirstTag && AA)
                    {
                        if (AnchoredFlag)
                        {
                            GetSplicedMatches(Info, SubNode, LeftExon, LeftExonPos, BridgedExon, -1);
                        }
                        else
                        {
                            GetSplicedMatches(Info, SubNode, Exon, AnchorPos, BridgedExon, -1);
                        }
                    }
                    if (!AA)
                    {
                        g_TagExonArray[g_TagExonArrayPos++] = BridgedExon;
                        g_TagSpliceArray[g_TagSpliceArrayPos++] = Exon->ForwardEdges + EdgeIndex; //Exon->ForwardEdgePower[EdgeIndex];
                    }

                    // We've now spanned an edge with our tag. 
                    if (AnchoredFlag)
                    {
                        GetSplicedTagMatches(Info, LeftExon, LeftExonPos, BridgedExon, SubNode, 1);
                    }
                    else
                    {
                        GetSplicedTagMatches(Info, Exon, AnchorPos, BridgedExon, SubNode, 1);
                    }
                    g_TagExonArrayPos = OldExonPos;
                    g_TagSpliceArrayPos = OldSplicePos;
                }
                break;
            } // following an edge forward
            else
            {
                CurrentNode = CurrentNode->Children[Exon->Sequence[SequencePos] - 'A'];
                g_TagBuffer[g_TagBufferPos++] = Exon->Sequence[SequencePos];
                g_TagBufferSpliced[g_TagBufferPosSpliced++] = Exon->Sequence[SequencePos];
                if (!CurrentNode)
                {
                    break;
                }
                if (CurrentNode->FirstTag)
                {
                    if (AnchoredFlag)
                    {
                        GetSplicedMatches(Info, CurrentNode, LeftExon, LeftExonPos, Exon, SequencePos);
                    }
                    else
                    {
                        GetSplicedMatches(Info, CurrentNode, Exon, AnchorPos, Exon, SequencePos);
                    }
                }
            }
        } // sequencepos iteration
    } // anchorpos
}

// Given an edge record for one exon, get the corresponding edge struct for the linked exon.
// If ForwardFlag is 1, then the edge passed is a forward edge (and the reciprocal edge is a backward edge)
// If ForwardFlag is 0, then the edge passed is a backward edge (and the reciprocal edge is a forward edge)
ExonEdge* GetReciprocalExonEdge(ExonEdge* Edge, int ForwardFlag)
{
    int EdgeIndex;
    ExonEdge* OtherEdge;
    if (ForwardFlag)
    {
        for (EdgeIndex = 0; EdgeIndex < Edge->Exon->BackEdgeCount; EdgeIndex++)
        {
            OtherEdge = Edge->Exon->BackwardEdges + EdgeIndex;
            if (OtherEdge->Exon == Edge->Source && OtherEdge->AA == Edge->AA)
            {
                return Edge->Exon->BackwardEdges + EdgeIndex;
            }
        }
    }
    else
    {
        for (EdgeIndex = 0; EdgeIndex < Edge->Exon->ForwardEdgeCount; EdgeIndex++)
        {
            OtherEdge = Edge->Exon->ForwardEdges + EdgeIndex;
            if (OtherEdge->Exon == Edge->Source && OtherEdge->AA == Edge->AA)
            {
                return Edge->Exon->ForwardEdges + EdgeIndex;
            }
        }
    }
    INSPECT_ASSERT(0);
    return NULL;
}

void SearchSplicableGene(SearchInfo* Info, GeneStruct* Gene)
{
    TrieNode* CurrentNode;
    int ExonIndex;
    int SequencePos;
    int SequenceLength;
    int EdgeIndex;
    int VerboseFlag = 0;
    ExonStruct* Exon;
    ExonStruct* ActiveExon;
    int TotalSequenceLength = 0;
    int TotalExonCount = 0;
    ExonEdge* Edge;

    // CheckGene(Gene);
    if (VerboseFlag)
    {
        printf("Gene %d: '%s' (%d exons)\n", Info->RecordNumber, Gene->Name, Gene->ExonCount);
    }
    for (ExonIndex = 0; ExonIndex < Gene->ExonCount; ExonIndex++)
    {
        SequencePos = 0;
        Exon = Gene->Exons + ExonIndex;
        ActiveExon = Exon;
        if (VerboseFlag)
        {
            printf("Search exon %d: '%s'\n", ExonIndex, Exon->Sequence);
        }

        SequenceLength = Exon->Length;
        TotalExonCount++;
        TotalSequenceLength += SequenceLength;

        ////////////////////////////////////////////////////////////
        // Try starting with each edge INTO the exon.  XXX-T-AGXX
        for (EdgeIndex = 0; EdgeIndex < Exon->BackEdgeCount; EdgeIndex++)
        {
            Edge = Exon->BackwardEdges + EdgeIndex;
            if (!Edge->AA)
            {
                continue;
            }
            if (Edge->Source->Start == 31887068 && Edge->AA == 'D' && Edge->Source->Index == 463)
            {
                Edge = Edge;
            }
            g_TagExonArray[0] = Edge->Exon;
            g_TagExonArray[1] = Exon;
            g_TagExonArrayPos = 2;
            g_TagSpliceArray[0] = GetReciprocalExonEdge(Edge, 0);
            g_TagSpliceArrayPos = 1;
            g_TagBuffer[0] = Edge->AA;
            g_TagBufferPos = 1;
            if (Edge->Power)
            {
                g_TagBufferSpliced[0] = ';';
                g_TagBufferSpliced[1] = Edge->AA;
                g_TagBufferSpliced[2] = ';';
                g_TagBufferPosSpliced = 3;
            }
            else
            {
                g_TagBufferSpliced[0] = Edge->AA;
                g_TagBufferPosSpliced = 1;
            }
            CurrentNode = Info->Root->Children[Edge->AA - 'A'];
            if (CurrentNode)
            {
                GetSplicedTagMatches(Info, Edge->Exon, -1, Exon, CurrentNode, 1);
                // Special for tagless search:
                if (CurrentNode->FirstTag)
                {
                    GetSplicedMatches(Info, CurrentNode, Edge->Exon, -1, Exon, -1);
                }
            }
        }
        g_TagExonArray[0] = Exon;
        g_TagExonArrayPos = 1;
        g_TagSpliceArrayPos = 0;
        g_TagBufferPos = 0;
        g_TagBufferPosSpliced = 0;
        GetSplicedTagMatches(Info, Exon, 0, Exon, Info->Root, 0);
    } // loop over exons
}

// Main method for Spliced.c: Given a collection of tags (Root) and a binary splicedb (FileName), search
// for matches to the current list of spectra (list head GlobalOptions->FirstSpectrum, but we get to them
// via back-links from tags).  Score matches with Scorer.  If GeneNames is not null, then it's an array
// of gene names to be searched, and we skip any gene whose name isn't on the list.
void SearchSplicableGenes(SearchInfo* Info)
{
    FILE* File;
    GeneStruct* Gene;
    int VerboseFlag = 0;
    int RecordNumber = -1;
    int TotalSequenceLength = 0;
    int TotalExonCount = 0;

    //
    AllocSpliceStructures();

    File = Info->DB->DBFile;
    if (!File)
    {
        printf("** Erorr: Unable to open gene database '%s'.  No search performed.\n", Info->DB->FileName);
        return;
    }
    fseek(File, 0, 0);
    while (1)
    {
        RecordNumber++;
        //StartTime = clock();
        Gene = LoadGene(File);
        if (!Gene)
        {
            break;
        }
        SearchSplicableGene(Info, Gene);

        FreeGene(Gene);

    } // while genes
    printf("Searched %d genes, %d exons, %d bases\n", RecordNumber, TotalExonCount, TotalSequenceLength);
}

// For debugging: Exercise the splice functionality.
void TestSplicing()
{
    if (!SLeftMatchBuffer)
    {
        AllocSpliceStructures();
    }
    LoadGenes("C:\\source\\Inspect\\Inspect\\SpliceDB\\Ch1.dat");
    DebugPrintGenes();
}


// inspect test splicedb <DBPath> <Start> <End> [<DesiredProtein>]
// Print out all genes which overlap the region of interest.  Useful for asking
// why a particular known protein isn't present in its entirety.  If a protein
// sequence is provided as well, we align that sequence against the exon graph,
// determining how much of it is present in database.
void TestSpliceDB(int argc, char** argv)
{
    int StartPos;
    int EndPos;
    FILE* DBFile;
    GeneStruct* Gene;
    int GeneStart;
    int GeneEnd;
    int ExonIndex;
    char* DesiredProtein = NULL;
    int DesiredProteinLength;
    FILE* ProteinFile;
    int BytesRead;
    //
    if (argc < 4)
    {
        printf("** Not enough args - bailing out\n");
        return;
    }
    DBFile = fopen(argv[3], "rb");
    if (!DBFile)
    {
        printf("** Error: Can't open splicedb at '%s'\n", argv[3]);
        return;
    }
    StartPos = 0;
    EndPos = -1;
    if (argc > 4)
    {
        StartPos = atoi(argv[4]);
    }
    if (argc > 5)
    {
        EndPos = atoi(argv[5]);
    }
    if (argc > 6)
    {
        // Read protein sequence:
        ProteinFile = fopen(argv[6], "rb");
        if (!ProteinFile)
        {
            printf("** Error: Can't read target protein sequence from '%s'\n", argv[6]);
            return;
        }
        fseek(ProteinFile, 0, 2);
        DesiredProteinLength = ftell(ProteinFile);
        DesiredProtein = (char*)calloc(DesiredProteinLength, sizeof(char));
        fseek(ProteinFile, 0, 0);
        BytesRead = ReadBinary(DesiredProtein, sizeof(char), DesiredProteinLength, ProteinFile);
        DesiredProtein[BytesRead] = '\0';
        fclose(ProteinFile);
    }
    ////////////////////////
    while (1)
    {
        Gene = LoadGene(DBFile);
        if (!Gene)
        { 
            break;
        }
        // Decide whether to print the gene:
        GeneStart = -1;
        GeneEnd = -1;
        for (ExonIndex = 0; ExonIndex < Gene->ExonCount; ExonIndex++)
        {
            if (GeneStart < 0 || GeneStart > Gene->Exons[ExonIndex].Start)
            {
                GeneStart = Gene->Exons[ExonIndex].Start;
            }
            if (GeneEnd < 0 || GeneEnd < Gene->Exons[ExonIndex].End)
            {
                GeneEnd = Gene->Exons[ExonIndex].End;
            }
        }
        if (GeneEnd >= StartPos && (EndPos < 0 || GeneStart < EndPos))
        {
            //DebugPrintGene(Gene);
            if (DesiredProtein)
            {
                AlignSequenceAgainstExonGraph(Gene, DesiredProtein, NULL, -10, -1);
            }
            else
            {
                DebugPrintGene(Gene);
            }
        }
        FreeGene(Gene);
    }

}
