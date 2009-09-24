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
#include <stdio.h>
#include <memory.h>
//#include <malloc.h>
#include "Utils.h"
#include "Spliced.h"

// This file implements an alignment algorithm between a sequence and 
// an exon graph, or between two exon graphs. 
// Exon-graph alignment is very similar to the smith-waterman alignment
// algorithm.  The main difference is that in the recurrence relation,
// we may have several "previous" cells to move to, due to the makeup
// of the exon graph.

#define AA_COUNT 26

// IntNode: For handling lists of successors to each row / column
typedef struct IntNode
{
    int Value;
    struct IntNode* Next;
} IntNode;


// Forward declarations:
int AlignExonGraphAgainstExonGraph(GeneStruct* GeneA, GeneStruct* GeneB,
    char* ScoringMatrixFileName, int StartGapPenalty, int ExtendGapPenalty);
void FreePrevCellTable(IntNode** PrevCell, int Size);

// Default distance matrix for aligning protein sequences: Hamming distance,
// a bonus for each match and a penalty for each mismatch.  This is adequate
// for most purposes.
int* GenerateHammingDistanceMatrix()
{
    int X;
    int Y;
    int* Matrix;
    //
    Matrix = (int*)calloc(AA_COUNT*AA_COUNT, sizeof(int));
    for (X = 0; X < AA_COUNT; X++)
    {
        for (Y = 0; Y < AA_COUNT; Y++)
        {
            if (X == Y)
            {   
                Matrix[Y*AA_COUNT + X] = 10;
            }
            else
            {
                Matrix[Y*AA_COUNT + X] = -10;
            }
        }
    }
    return Matrix;
}

int* LoadScoringMatrix(char* ScoringMatrixFileName)
{
    //FILE* ScoringMatrixFile;
    printf("** Scoring matrix support not implemented yet - use Hamming disatnce for now\n");
    return NULL;
}


void XALinkRowToRow(IntNode** PrevY, int Y, int TargetY)
{
    IntNode* INode;
    IntNode* NewINode;
    NewINode = (IntNode*)calloc(1, sizeof(IntNode));
    NewINode->Value = TargetY;
    //printf("Back-link from row %d to target %d\n", Y, TargetY);
    if (PrevY[Y])
    {
        for (INode = PrevY[Y]; INode->Next; INode = INode->Next)
        {
            ;
        }
        INode->Next = NewINode;
    }
    else
    {
        PrevY[Y] = NewINode;
    }
}

// This function links the specified row to the specified exon.
// The catch: If the exon has length 0, then we link to the specified exon's predecessors.
void XALinkRowToExon(GeneStruct* Gene, IntNode** PrevY, int Y, int ExonIndex, int* ExonOffsets, int* ExonEdgeOffsets)
{
    ExonStruct* Exon;
    int EdgeIndex;
    int AAEdgeCount;
    ExonEdge* Edge;
    //
    Exon = Gene->Exons + ExonIndex;
    // Standard case: The exon is non-empty, so we link back to it.
    if (Exon->Length)
    {
        XALinkRowToRow(PrevY, Y, ExonOffsets[Exon->Index] + Exon->Length - 1);
        return;
    }
    // Special case: The exon is empty, so we link to the exon's predecessors.
    AAEdgeCount = 0;
    for (EdgeIndex = 0; EdgeIndex < Exon->BackEdgeCount; EdgeIndex++)
    {
        Edge = Exon->BackwardEdges + EdgeIndex;
        if (Edge->AA)
        {
            XALinkRowToRow(PrevY, Y, ExonEdgeOffsets[ExonIndex] + AAEdgeCount);
            AAEdgeCount++;
        }
        else
        {
            XALinkRowToExon(Gene, PrevY, Y, Edge->Exon->Index, ExonOffsets, ExonEdgeOffsets);
        }
    }
}

// Test scaffolding for graph-based alignment
void TestExonGraphAlignment(int argc, char* argv[])
{
    GeneStruct* GeneA;
    GeneStruct* GeneB;
    FILE* GeneFile;
    //
    GeneFile = fopen(argv[1], "rb");
    if (!GeneFile)
    {
        printf("** Error: Can't open gene file '%s'.\n", argv[1]);
        return;
    }
    GeneA = LoadGene(GeneFile);
    fclose(GeneFile);
    GeneFile = fopen(argv[2], "rb");
    if (!GeneFile)
    {
        printf("** Error: Can't open gene file '%s'.\n", argv[1]);
        return;
    }
    GeneB = LoadGene(GeneFile);
    AlignExonGraphAgainstExonGraph(GeneA, GeneB, NULL, -10, -3);
    //AlignSequenceAgainstExonGraph(Gene, Sequence, NULL, -10, -3);
    printf("\n\nAlignment complete.\n");
}

#define Z_STANDARD 0
#define Z_GAP_IN_X 1
#define Z_GAP_IN_Y 2

// Count the links and the total amino acids in this exon graph.  Used
// in determining the table size in alignment.
void GetExonGraphSize(GeneStruct* Gene, int* pLinkCount, int* pSize)
{
    int ExonIndex;
    int LinkIndex;
    ExonStruct* Exon;
    ExonEdge* Edge;
    //
    *pSize = 0;
    for (ExonIndex = 0; ExonIndex < Gene->ExonCount; ExonIndex++)
    {
        Exon = Gene->Exons + ExonIndex;
        *pSize += Exon->Length;
        
        for (LinkIndex = 0; LinkIndex < Exon->BackEdgeCount; LinkIndex++)
        {
            Edge = Exon->BackwardEdges + LinkIndex;
            if (Edge->AA)
            {
                (*pLinkCount)++;
                (*pSize)++;
            }
        }
    }
}

int CompareExonsForward(const ExonStruct* ExonA, const ExonStruct* ExonB)
{
    if (ExonA->Start < ExonB->Start)
    {
        return -1;
    }
    if (ExonA->Start > ExonB->Start)
    {
        return 1;
    }
    // ExonA->Start == ExonB->Start    
    if (ExonA->End < ExonB->End)
    {
        return -1;
    }
    if (ExonA->End > ExonB->End)
    {
        return 1;
    }
    // Same coordinates?  Arbitrary sort:
    if (ExonA < ExonB)
    {
        return -1;
    }
    else 
    {
        return 1;
    }
}
int CompareExonsReverse(const ExonStruct* ExonA, const ExonStruct* ExonB)
{
    if (ExonA->End > ExonB->End)
    {
        return -1;
    }
    if (ExonA->End < ExonB->End)
    {
        return 1;
    }
    // ExonA->End == ExonB->End
    if (ExonA->Start > ExonB->Start)
    {
        return -1;
    }
    if (ExonA->Start > ExonB->Start)
    {
        return 1;
    }
    // Same coordinates?  Arbitrary sort:
    if (ExonA > ExonB)
    {
        return 1;
    }
    else 
    {
        return -1;
    }
}


// Build the necessary arrays for a d.p. alignment against an exon graph:
// - ExonEdgeOffsets[n] is the row for the first back-link-with-aa of exon n
// - ExonOffsets[n] is the row for the first aa in exon n
// - YSequence[n] is the nth character in the flattened graph
// - PrevY[n] is a linked list of predecessors for row n.  If n is within an exon,
//  there's just one entry, n-1.  If n is the start of an exon, there may be several
//  entries.  If n comes from an edge, there'll be exactly one entry, for the earlier exon.
void FlattenExonsForAlignment(GeneStruct* Gene, int* ExonOffsets, int* ExonEdgeOffsets, 
    char* YSequence, IntNode** PrevY, char** YRowInfo)
{
    int Y;
    int ExonIndex;
    int StartExonIndex;
    int ExonIterateDir;
    int ExonCount = Gene->ExonCount;
    ExonStruct* Exon;
    int AALinkCount;
    int LinkIndex;
    IntNode* NewINode;
    int Pos;
    ExonEdge* Edge;
    //
    StartExonIndex = 0;
    ExonIterateDir = 1;
    Y = 0;
    for (ExonIndex = StartExonIndex; ExonIndex >= 0 && ExonIndex < ExonCount; ExonIndex += ExonIterateDir)
    {
        Exon = Gene->Exons + ExonIndex;
        // Add a row for each edge with AA:
        AALinkCount = 0;
        ExonEdgeOffsets[ExonIndex] = -1;
        for (LinkIndex = 0; LinkIndex < Exon->BackEdgeCount; LinkIndex++)
        {
            // For each exon with an associated amino acid:
            // - Add a row for the edge, with one back-link.
            // - Remember the row number, so that the exon can be linked back to this row
            Edge = Exon->BackwardEdges + LinkIndex;
            if (Edge->AA)
            {
                if (ExonEdgeOffsets[ExonIndex] < 0)
                {
                    ExonEdgeOffsets[ExonIndex] = Y;
                }
                //EdgeOffsets[AALinkCount] = Y;
                YSequence[Y] = Edge->AA;
                sprintf(YRowInfo[Y], "X%d backlink%d to %d", ExonIndex, LinkIndex, Edge->Exon->Index);
                XALinkRowToExon(Gene, PrevY, Y, Edge->Exon->Index, ExonOffsets, ExonEdgeOffsets);
                //NewINode = (IntNode*)calloc(1, sizeof(IntNode));
                //NewINode->Value = ExonOffsets[Exon->BackExon[LinkIndex]->Index] + Exon->BackExon[LinkIndex]->Length;
                //PrevY[Y] = NewINode;
                AALinkCount++;
                if (PrevY[Y])
                {
                    printf("Y %d: exon %d link %d, back to exon %d y %d\n", Y, ExonIndex, LinkIndex, Edge->Exon->Index, PrevY[Y]->Value);
                }
                Y++;
            }
        }
        if (!Exon->Length)
        {
            continue;
        }
        printf("Y %d: start exon %d body.\n", Y, ExonIndex);
        // Add back-links for the first AA in the exon:
        ExonOffsets[ExonIndex] = Y;
        AALinkCount = 0;
        for (LinkIndex = 0; LinkIndex < Exon->BackEdgeCount; LinkIndex++)
        {
            Edge = Exon->BackwardEdges + LinkIndex;
            if (Edge->AA)
            {
                XALinkRowToRow(PrevY, Y, ExonEdgeOffsets[ExonIndex] + AALinkCount);
                printf("ExonEdge offset %d.  Link %d goes to AARow %d\n", ExonEdgeOffsets[ExonIndex], LinkIndex, PrevY[Y]->Value);
                AALinkCount++;
            }
            else
            {
                XALinkRowToExon(Gene, PrevY, Y, Edge->Exon->Index, ExonOffsets, ExonEdgeOffsets);
                printf("Y %d: start exon %d.  Link %d goes to exon row %d\n", Y, ExonIndex, LinkIndex, PrevY[Y]->Value);
            }
        }
        // Add one row for each AA in the exon proper:
        for (Pos = 0; Pos < Exon->Length; Pos++)
        {
            sprintf(YRowInfo[Y], "X%d pos %d/%d", ExonIndex, Pos, Exon->Length);
            YSequence[Y] = Exon->Sequence[Pos];
            if (Pos)
            {
                NewINode = (IntNode*)calloc(1, sizeof(IntNode));
                NewINode->Value = Y-1;
                PrevY[Y] = NewINode;
            }
            //printf("%d %c pos %d in exon #%d\n", Y, YSequence[Y], Pos, Exon->Index);
            Y++;
        }
    }
}

void SortGeneExons(GeneStruct* Gene)
{
    return;
}

// ExonGraphAlign extends the Smith-Waterman alignment algorithm to 
// handle local alignment of a sequence with an exon graph.  The scoring
// matrix (such as BLOSUM) can be specified, or Hamming distance can be
// used.  Gap penalties should also be specified.  The function will 
// return the score of the best alignment, and (optionally) set 'verbose' 
// strings specifying the alignment itself, like this:
//   EAM--APK
//   ***  * *
//   EAMCGARK
// The data structure we use in implementing this alignment is a grid
// (stored as an array), where each row has a linked list of zero or more
// nodes specifying the legal predecessor rows.
int AlignSequenceAgainstExonGraph(GeneStruct* Gene, char* Sequence, 
    char* ScoringMatrixFileName, int StartGapPenalty, int ExtendGapPenalty)
{
    int* ScoreTable;
    int* NextX;
    int* NextY;
    int* NextZ;
    int X;
    int Y;
    int Z;
    IntNode** PrevY;
    int TableIndex;
    int PrevTableIndex;
    int TableSize;
    int* ScoringMatrix = NULL;
    int* ExonOffsets;
    int* ExonEdgeOffsets;
    int ExonCount = 0;
    int LinkCount = 0;
    int SequenceLength;
    int MaxY;
    int AlignScore;
    int Score;
    char* YSequence;
    IntNode* INode;
    int XBlockSize;
    int YBlockSize;
    char ResidueA;
    char ResidueB;
    int BestX = 0;
    int BestY = 0;
    int BestZ = 0;
    int AlignStringLength;
    char* AlignStringA;
    char* AlignStringB;
    char* AlignStringC;
    int BestScore;
    int StartExonIndex;
    int NearY;
    int ExonIterateDir;
    char** YRowInfo;
    // Ensure gap penalties are NEGATIVE numbers.  Negative is bad.
    if (StartGapPenalty > 0)
    {
        StartGapPenalty = -StartGapPenalty;
    }
    if (ExtendGapPenalty > 0)
    {
        ExtendGapPenalty = -ExtendGapPenalty;
    }

    // Load the scoring matrix (or use default hamming distance)
    if (ScoringMatrixFileName)
    {
        ScoringMatrix = LoadScoringMatrix(ScoringMatrixFileName);
    }
    if (!ScoringMatrix)
    {
        ScoringMatrix = GenerateHammingDistanceMatrix();
    }
    SequenceLength = strlen(Sequence);
    ////////////////////////////////////////////////////////////
    // Count the exons and edges (with aa):
    ExonCount = Gene->ExonCount;
    
    if (Gene->ForwardFlag)
    {
        ExonIterateDir = 1;
        StartExonIndex = 0;
    }
    else
    {
        ExonIterateDir = -1;
        StartExonIndex = Gene->ExonCount - 1;
    }
    GetExonGraphSize(Gene, &LinkCount, &MaxY);
    ////////////////////////////////////////////////////////////
    // Allocate arrays:
    TableSize = MaxY * SequenceLength * 3;
    ScoreTable = (int*)calloc(TableSize, sizeof(int));
    NextX = (int*)calloc(TableSize, sizeof(int));
    NextY = (int*)calloc(TableSize, sizeof(int));
    NextZ = (int*)calloc(TableSize, sizeof(int));
    PrevY = (IntNode**)calloc(MaxY, sizeof(IntNode*));
    YSequence = (char*)calloc(MaxY, sizeof(char));
    ExonOffsets = (int*)calloc(ExonCount, sizeof(int));
    ExonEdgeOffsets = (int*)calloc(ExonCount, sizeof(int));
    YRowInfo = (char**)calloc(MaxY, sizeof(char*));
    for (Y = 0; Y < MaxY; Y++)
    {
        YRowInfo[Y] = (char*)calloc(64, sizeof(char));
    }
    ////////////////////////////////////////////////////////////
    // Initialize the linked lists giving predecessors at each point.
    SortGeneExons(Gene);
    //DebugPrintGene(Gene); 
    FlattenExonsForAlignment(Gene, ExonOffsets, ExonEdgeOffsets, YSequence, PrevY, YRowInfo);
    //////////////////////////////////////////////////////////////////////
    // Debug print:
    for (Y = 0; Y < MaxY; Y++)
    {
        printf("%d ", Y);
        for (NearY = Y - 5; NearY < Y + 6; NearY++)
        {
            if (NearY >= 0 && NearY < MaxY)
            {
                printf("%c", YSequence[NearY]);
            }
        }
        for (INode = PrevY[Y]; INode; INode = INode->Next)
        {
            printf(" ->%d (", INode->Value);
            for (NearY = INode->Value - 3; NearY < INode->Value; NearY++)
            {
                if (NearY >= 0 && NearY < MaxY)
                {
                    printf("%c", YSequence[NearY]);
                }
            }
            printf(" %c ", YSequence[INode->Value]);
            for (NearY = INode->Value + 1; NearY < INode->Value + 4; NearY++)
            {
                if (NearY >= 0 && NearY < MaxY)
                {
                    printf("%c", YSequence[NearY]);
                }
            }
            printf(")");
        }
        printf("\n");
    }
    ////////////////////////////////////////////////////////////
    // Carry out dynamic programming:
    XBlockSize = 3;
    YBlockSize = XBlockSize * SequenceLength;
    for (Y = 0; Y < MaxY; Y++)
    {
        ResidueB = YSequence[Y] - 'A';
        if (ResidueB < 0 || ResidueB > 26)
        {
            ResidueB = 23; //'X';
        }
        for (X = 0; X < SequenceLength; X++)
        {
            ResidueA = Sequence[X] - 'A';
            if (ResidueA < 0 || ResidueA > 26)
            {
                ResidueA = 23; //'X';
            }
            ////////////////////////////
            // Z == 0, the alignment table:
            TableIndex = Y*YBlockSize + X*XBlockSize + Z_STANDARD;
            // Default: Jump in
            BestScore = 0;
            ScoreTable[TableIndex] = BestScore;
            NextX[TableIndex] = -1;
            NextY[TableIndex] = -1;
            NextZ[TableIndex] = -1;
            // Consider aligning:
            AlignScore = ScoringMatrix[ResidueA * AA_COUNT + ResidueB];
            // Aligning at the edges of the world is allowed:
            if (!X || !PrevY[Y])
            {
                if (AlignScore > BestScore)
                {
                    ScoreTable[TableIndex] = AlignScore;
                    BestScore = AlignScore;
                }
            }
            else
            {
                // Consider each predecessor row:
                for (INode = PrevY[Y]; INode; INode = INode->Next)
                {
                    PrevTableIndex = INode->Value * YBlockSize + (X-1)*XBlockSize + 0;
                    Score = AlignScore + ScoreTable[PrevTableIndex];
                    if (Score > BestScore)
                    {
                        BestScore = Score;
                        ScoreTable[TableIndex] = BestScore;
                        NextX[TableIndex] = X - 1;
                        NextY[TableIndex] = INode->Value;
                        NextZ[TableIndex] = 0;
                    }
                }
            }
            // Consider gapping in x:
            if (X)
            {
                PrevTableIndex = Y * YBlockSize + (X-1) * XBlockSize + Z_GAP_IN_X;
                Score = StartGapPenalty + ScoreTable[PrevTableIndex];
                if (Score > BestScore)
                {
                    BestScore = Score;
                    ScoreTable[TableIndex] = BestScore;
                    NextX[TableIndex] = X - 1;
                    NextY[TableIndex] = Y;
                    NextZ[TableIndex] = Z_GAP_IN_X;
                }
            }
            // Consider gapping in y:
            for (INode = PrevY[Y]; INode; INode = INode->Next)
            {
                PrevTableIndex = INode->Value * YBlockSize + X * XBlockSize + Z_GAP_IN_Y;
                Score = StartGapPenalty + ScoreTable[PrevTableIndex];
                if (Score > BestScore)
                {
                    BestScore = Score;
                    ScoreTable[TableIndex] = BestScore;
                    NextX[TableIndex] = X;
                    NextY[TableIndex] = INode->Value;
                    NextZ[TableIndex] = Z_GAP_IN_Y;
                }
            }
            //printf("At %d, %d, 0: Score %d, prev %d, %d, %d\n", X, Y, ScoreTable[TableIndex],
            //    NextX[TableIndex], NextY[TableIndex], NextZ[TableIndex]);
            ////////////////////////////
            // Z=1, gapping in x:
            // By default, close the gap...but also consider extending it (unless x == 0)
            TableIndex = Y*YBlockSize + X*XBlockSize + Z_GAP_IN_X;
            PrevTableIndex = Y*YBlockSize + X*XBlockSize + Z_STANDARD;
            BestScore = ScoreTable[PrevTableIndex];
            ScoreTable[TableIndex] = BestScore;
            NextX[TableIndex] = X;
            NextY[TableIndex] = Y;
            NextZ[TableIndex] = Z_STANDARD;
            if (X)
            {
                Score = ExtendGapPenalty + ScoreTable[Y*YBlockSize + (X-1)*XBlockSize + Z_GAP_IN_X];
                if (Score > BestScore)
                {
                    ScoreTable[TableIndex] = BestScore;
                    NextX[TableIndex] = X - 1;
                    NextY[TableIndex] = Y;
                    NextZ[TableIndex] = Z_GAP_IN_X;
                }
            }
            ////////////////////////////
            // Z=2, gapping in y:
            // By default, close the gap...but also consider extending it
            TableIndex = Y*YBlockSize + X*XBlockSize + Z_GAP_IN_Y;
            PrevTableIndex = Y*YBlockSize + X*XBlockSize + Z_STANDARD;
            BestScore = ScoreTable[PrevTableIndex];
            ScoreTable[TableIndex] = BestScore;
            NextX[TableIndex] = X;
            NextY[TableIndex] = Y;
            NextZ[TableIndex] = Z_STANDARD;
            for (INode = PrevY[Y]; INode; INode = INode->Next)
            {
                Score = ExtendGapPenalty + ScoreTable[INode->Value*YBlockSize + X*XBlockSize + Z_GAP_IN_Y];
                if (Score > BestScore)
                {
                    ScoreTable[TableIndex] = BestScore;
                    NextX[TableIndex] = X;
                    NextY[TableIndex] = INode->Value;
                    NextZ[TableIndex] = Z_GAP_IN_Y;
                }
            }
        }
    }
    ////////////////////////////////////////////////////////////
    // Find where the best alignment ends:
    BestScore = -9999;
    for (X = 0; X < SequenceLength; X++)
    {
        for (Y = 0; Y < MaxY; Y++)
        {
            for (Z = 0; Z < 3; Z++)
            {
                Score = ScoreTable[Y*YBlockSize + X*XBlockSize + Z];
                if (Score > BestScore)
                {
                    BestScore = Score;
                    BestX = X;
                    BestY = Y;
                    BestZ = Z;
                }
            }
        }
    }
    ////////////////////////////////////////////////////////////
    // Produce strings for the optimal alignment:
    X = BestX;
    Y = BestY;
    Z = BestZ;
    AlignStringLength = 0;
    while (X >= 0)
    {
        TableIndex = Y*YBlockSize + X*XBlockSize + Z;
        // Each step we take will add to the string...except closing a gap.
        if (!Z || NextZ[TableIndex])
        {
            AlignStringLength++;
        }
        X = NextX[TableIndex];
        Y = NextY[TableIndex];
        Z = NextZ[TableIndex];
    }

    AlignStringA = (char*)calloc(AlignStringLength + 1, sizeof(char));
    AlignStringB = (char*)calloc(AlignStringLength + 1, sizeof(char));
    AlignStringC = (char*)calloc(AlignStringLength + 1, sizeof(char));
    X = BestX;
    Y = BestY;
    Z = BestZ;
    while (X >= 0)
    {
        AlignStringLength--;
        TableIndex = Y*YBlockSize + X*XBlockSize + Z;
        switch (Z)
        {
        case Z_STANDARD:
            switch (NextZ[TableIndex])
            {
            case Z_STANDARD:
            default:
                ResidueA = Sequence[X];
                ResidueB = YSequence[Y];
                AlignStringA[AlignStringLength] = Sequence[X];
                AlignStringC[AlignStringLength] = YSequence[Y];
                if (ResidueA == ResidueB)
                {
                    AlignStringB[AlignStringLength] = '*';
                }
                else
                {
                    AlignStringB[AlignStringLength] = ' ';
                }
                printf("X %d (%c) Y %d (%c) %s\n", X, ResidueA, Y, ResidueB, YRowInfo[Y]);
                break;
            case Z_GAP_IN_X:
                AlignStringA[AlignStringLength] = Sequence[X];
                AlignStringB[AlignStringLength] = ' ';
                AlignStringC[AlignStringLength] = '-';
                break;
            case Z_GAP_IN_Y:
                AlignStringA[AlignStringLength] = '-';
                AlignStringB[AlignStringLength] = ' ';
                AlignStringC[AlignStringLength] = YSequence[Y];
                break;
            }
            break;
        case Z_GAP_IN_X:
            if (NextZ[TableIndex])
            {
                AlignStringA[AlignStringLength] = Sequence[X];
                AlignStringB[AlignStringLength] = ' ';
                AlignStringC[AlignStringLength] = '-';
            }
            break;            
        case Z_GAP_IN_Y:
            if (NextZ[TableIndex])
            {
                AlignStringA[AlignStringLength] = '-';
                AlignStringB[AlignStringLength] = ' ';
                AlignStringC[AlignStringLength] = YSequence[Y];
            }
            break;
        }
        
        // Each step we take will add to the string...except closing a gap.
        if (Z && !NextZ[TableIndex])
        {
            AlignStringLength++;
        }
        X = NextX[TableIndex];
        Y = NextY[TableIndex];
        Z = NextZ[TableIndex];
    }
    printf("Alignment score %d.  Alignment follows:\n", BestScore);
    printf("%s\n", AlignStringA);
    printf("%s\n", AlignStringB);
    printf("%s\n", AlignStringC);

    ////////////////////////////////////////////////////////////
    // cleanup:
    SafeFree(ScoringMatrix);
    SafeFree(ScoreTable);
    SafeFree(ExonOffsets);
    SafeFree(NextX);
    SafeFree(NextY);
    SafeFree(NextZ);
    SafeFree(YSequence);
    SafeFree(ExonEdgeOffsets);
    SafeFree(AlignStringA);
    SafeFree(AlignStringB);
    SafeFree(AlignStringC);
    FreePrevCellTable(PrevY, MaxY);
    for (Y = 0; Y < MaxY; Y++)
    {
        SafeFree(YRowInfo[Y]);
    }
    SafeFree(YRowInfo);

    return BestScore;
}

// Free an array of linked-lists providing predecessor cells
void FreePrevCellTable(IntNode** PrevCell, int Size)
{
    int Index;
    IntNode* Node;
    IntNode* Prev;
    if (!PrevCell)
    {
        return;
    }
    for (Index = 0; Index < Size; Index++)
    {
        Prev = NULL;
        for (Node = PrevCell[Index]; Node; Node = Node->Next)
        {
            SafeFree(Prev);
            Prev = Node;
        }
        SafeFree(Prev);
    }
    SafeFree(PrevCell);
}


// Align an exon graph against another exon graph.
int AlignExonGraphAgainstExonGraph(GeneStruct* GeneA, GeneStruct* GeneB,
    char* ScoringMatrixFileName, int StartGapPenalty, int ExtendGapPenalty)
{
    int MaxX;
    int MaxY;
    int LinkCountA;
    int LinkCountB;
    int* ScoringMatrix = NULL;
    int ExonCountA;
    int ExonCountB;
    int* ScoreTable;
    int* NextX;
    int* NextY;
    int* NextZ;
    IntNode** PrevX;
    IntNode** PrevY;
    char* XSequence;
    char* YSequence;
    int* ExonOffsetsA;
    int* ExonOffsetsB;
    int* ExonEdgeOffsetsA;
    int* ExonEdgeOffsetsB;
    int XBlockSize;
    int YBlockSize;
    IntNode* PrevNodeX;
    IntNode* PrevNodeY;
    int TableSize;
    int X;
    int Y;
    int Z;
    int BestX = 0;
    int BestY = 0;
    int BestZ = 0;
    char ResidueA;
    char ResidueB;
    int Score;
    int AlignScore;
    int BestScore;
    char* AlignStringA;
    char* AlignStringB;
    char* AlignStringC;
    int TableIndex;
    int PrevTableIndex;
    int AlignStringLength;
    char** RowInfoA;
    char** RowInfoB;
    // Ensure gap penalties are NEGATIVE numbers.  Negative is bad.
    if (StartGapPenalty > 0)
    {
        StartGapPenalty = -StartGapPenalty;
    }
    if (ExtendGapPenalty > 0)
    {
        ExtendGapPenalty = -ExtendGapPenalty;
    }

    // Load the scoring matrix (or use default hamming distance)

    if (ScoringMatrixFileName)
    {
        ScoringMatrix = LoadScoringMatrix(ScoringMatrixFileName);
    }
    if (!ScoringMatrix)
    {
        ScoringMatrix = GenerateHammingDistanceMatrix();
    }
    printf("\n\nGene A:\n");
    //DebugPrintGene(GeneA); 
    printf("\n\nGene B:\n");
    //DebugPrintGene(GeneB); 

    ////////////////////////////////////////////////////////////
    // Count the exons and edges (with aa):
    ExonCountA = GeneA->ExonCount;
    ExonCountB = GeneB->ExonCount;
    
    GetExonGraphSize(GeneA, &LinkCountA, &MaxX);
    GetExonGraphSize(GeneB, &LinkCountB, &MaxY);
    ////////////////////////////////////////////////////////////
    // Allocate arrays:
    TableSize = MaxY * MaxX * 3;
    ScoreTable = (int*)calloc(TableSize, sizeof(int));
    NextX = (int*)calloc(TableSize, sizeof(int));
    NextY = (int*)calloc(TableSize, sizeof(int));
    NextZ = (int*)calloc(TableSize, sizeof(int));
    XSequence = (char*)calloc(MaxX + 1, sizeof(char));
    PrevX = (IntNode**)calloc(MaxX, sizeof(IntNode*));
    YSequence = (char*)calloc(MaxY + 1, sizeof(char));
    PrevY = (IntNode**)calloc(MaxY, sizeof(IntNode*));
    ExonOffsetsA = (int*)calloc(ExonCountA, sizeof(int));
    ExonEdgeOffsetsA = (int*)calloc(ExonCountA, sizeof(int));
    ExonOffsetsB = (int*)calloc(ExonCountB, sizeof(int));
    ExonEdgeOffsetsB = (int*)calloc(ExonCountB, sizeof(int));
    RowInfoA = (char**)calloc(MaxX, sizeof(char*));
    RowInfoB = (char**)calloc(MaxY, sizeof(char*));
    for (X = 0; X < MaxX; X++)
    {
        RowInfoA[X] = (char*)calloc(64, sizeof(char));
    }
    for (Y = 0; Y < MaxY; Y++)
    {
        RowInfoB[Y] = (char*)calloc(64, sizeof(char));
    }

    ////////////////////////////////////////////////////////////
    // Initialize the linked lists giving predecessors at each point.
    SortGeneExons(GeneA);
    SortGeneExons(GeneB);
    FlattenExonsForAlignment(GeneA, ExonOffsetsA, ExonEdgeOffsetsA, XSequence, PrevX, RowInfoA);
    FlattenExonsForAlignment(GeneB, ExonOffsetsB, ExonEdgeOffsetsB, YSequence, PrevY, RowInfoB);
    ////////////////////////////////////////////////////////////
    // Carry out dynamic programming:
    XBlockSize = 3;
    YBlockSize = XBlockSize * MaxX;
    for (Y = 0; Y < MaxY; Y++)
    {
        ResidueB = YSequence[Y] - 'A';
        if (ResidueB < 0 || ResidueB > 26)
        {
            ResidueB = 23; //'X';
        }
        for (X = 0; X < MaxX; X++)
        {
            ResidueA = XSequence[X] - 'A';
            if (ResidueA < 0 || ResidueA > 26)
            {
                ResidueA = 23; //'X';
            }
            ////////////////////////////
            // Z == 0, the alignment table:
            TableIndex = Y*YBlockSize + X*XBlockSize + Z_STANDARD;
            // Default: Jump in
            BestScore = 0;
            ScoreTable[TableIndex] = BestScore;
            NextX[TableIndex] = -1;
            NextY[TableIndex] = -1;
            NextZ[TableIndex] = -1;
            // Consider aligning:
            AlignScore = ScoringMatrix[ResidueA * AA_COUNT + ResidueB];
            // Aligning at the edges of the world is allowed:
            if (!PrevX[X] || !PrevY[Y])
            {
                if (AlignScore > BestScore)
                {
                    ScoreTable[TableIndex] = AlignScore;
                    BestScore = AlignScore;
                }
            }
            else
            {
                // Consider each predecessor cell (x, y):
                for (PrevNodeX = PrevX[X]; PrevNodeX; PrevNodeX = PrevNodeX->Next)
                {
                    for (PrevNodeY = PrevY[Y]; PrevNodeY; PrevNodeY = PrevNodeY->Next)
                    {
                        PrevTableIndex = PrevNodeY->Value * YBlockSize + PrevNodeX->Value * XBlockSize + 0;
                        Score = AlignScore + ScoreTable[PrevTableIndex];
                        if (Score > BestScore)
                        {
                            BestScore = Score;
                            ScoreTable[TableIndex] = BestScore;
                            NextX[TableIndex] = PrevNodeX->Value;
                            NextY[TableIndex] = PrevNodeY->Value;
                            NextZ[TableIndex] = 0;
                        }
                    }
                }
            }
            // Consider gapping in x:
            for (PrevNodeX = PrevX[X]; PrevNodeX; PrevNodeX = PrevNodeX->Next)
            {
                PrevTableIndex = Y * YBlockSize + PrevNodeX->Value * XBlockSize + Z_GAP_IN_X;
                Score = StartGapPenalty + ScoreTable[PrevTableIndex];
                if (Score > BestScore)
                {
                    BestScore = Score;
                    ScoreTable[TableIndex] = BestScore;
                    NextX[TableIndex] = PrevNodeX->Value;
                    NextY[TableIndex] = Y;
                    NextZ[TableIndex] = Z_GAP_IN_X;
                }
            }
            // Consider gapping in y:
            for (PrevNodeY = PrevY[Y]; PrevNodeY; PrevNodeY = PrevNodeY->Next)
            {
                PrevTableIndex = PrevNodeY->Value * YBlockSize + X * XBlockSize + Z_GAP_IN_Y;
                Score = StartGapPenalty + ScoreTable[PrevTableIndex];
                if (Score > BestScore)
                {
                    BestScore = Score;
                    ScoreTable[TableIndex] = BestScore;
                    NextX[TableIndex] = X;
                    NextY[TableIndex] = PrevNodeY->Value;
                    NextZ[TableIndex] = Z_GAP_IN_Y;
                }
            }
            //printf("At %d, %d, 0: Score %d, prev %d, %d, %d\n", X, Y, ScoreTable[TableIndex],
            //    NextX[TableIndex], NextY[TableIndex], NextZ[TableIndex]);
            ////////////////////////////
            // Z=1, gapping in x:
            // By default, close the gap...but also consider extending it (unless x == 0)
            TableIndex = Y*YBlockSize + X*XBlockSize + Z_GAP_IN_X;
            PrevTableIndex = Y*YBlockSize + X*XBlockSize + Z_STANDARD;
            BestScore = ScoreTable[PrevTableIndex];
            ScoreTable[TableIndex] = BestScore;
            NextX[TableIndex] = X;
            NextY[TableIndex] = Y;
            NextZ[TableIndex] = Z_STANDARD;
            for (PrevNodeX = PrevX[X]; PrevNodeX; PrevNodeX = PrevNodeX->Next)
            {
                Score = ExtendGapPenalty + ScoreTable[Y*YBlockSize + PrevNodeX->Value * XBlockSize + Z_GAP_IN_X];
                if (Score > BestScore)
                {
                    ScoreTable[TableIndex] = BestScore;
                    NextX[TableIndex] = PrevNodeX->Value;
                    NextY[TableIndex] = Y;
                    NextZ[TableIndex] = Z_GAP_IN_X;
                }
            }
            ////////////////////////////
            // Z=2, gapping in y:
            // By default, close the gap...but also consider extending it
            TableIndex = Y*YBlockSize + X*XBlockSize + Z_GAP_IN_Y;
            PrevTableIndex = Y*YBlockSize + X*XBlockSize + Z_STANDARD;
            BestScore = ScoreTable[PrevTableIndex];
            ScoreTable[TableIndex] = BestScore;
            NextX[TableIndex] = X;
            NextY[TableIndex] = Y;
            NextZ[TableIndex] = Z_STANDARD;
            for (PrevNodeY = PrevY[Y]; PrevNodeY; PrevNodeY = PrevNodeY->Next)
            {
                Score = ExtendGapPenalty + ScoreTable[PrevNodeY->Value*YBlockSize + X*XBlockSize + Z_GAP_IN_Y];
                if (Score > BestScore)
                {
                    ScoreTable[TableIndex] = BestScore;
                    NextX[TableIndex] = X;
                    NextY[TableIndex] = PrevNodeY->Value;
                    NextZ[TableIndex] = Z_GAP_IN_Y;
                }
            }
        }
    }
    ////////////////////////////////////////////////////////////
    // Find where the best alignment ends:
    BestScore = -9999;
    for (X = 0; X < MaxX; X++)
    {
        for (Y = 0; Y < MaxY; Y++)
        {
            for (Z = 0; Z < 3; Z++)
            {
                Score = ScoreTable[Y*YBlockSize + X*XBlockSize + Z];
                if (Score > BestScore)
                {
                    BestScore = Score;
                    BestX = X;
                    BestY = Y;
                    BestZ = Z;
                }
            }
        }
    }
    ////////////////////////////////////////////////////////////
    // Produce strings for the optimal alignment:
    X = BestX;
    Y = BestY;
    Z = BestZ;
    AlignStringLength = 0;
    while (X >= 0)
    {
        TableIndex = Y*YBlockSize + X*XBlockSize + Z;
        // Each step we take will add to the string...except closing a gap.
        if (!Z || NextZ[TableIndex])
        {
            AlignStringLength++;
        }
        X = NextX[TableIndex];
        Y = NextY[TableIndex];
        Z = NextZ[TableIndex];
    }

    AlignStringA = (char*)calloc(AlignStringLength + 1, sizeof(char));
    AlignStringB = (char*)calloc(AlignStringLength + 1, sizeof(char));
    AlignStringC = (char*)calloc(AlignStringLength + 1, sizeof(char));
    X = BestX;
    Y = BestY;
    Z = BestZ;
    while (X >= 0)
    {
        AlignStringLength--;
        TableIndex = Y*YBlockSize + X*XBlockSize + Z;
        switch (Z)
        {
        case Z_STANDARD:
            switch (NextZ[TableIndex])
            {
            case Z_STANDARD:
            default:
                ResidueA = XSequence[X];
                ResidueB = YSequence[Y];
                AlignStringA[AlignStringLength] = ResidueA;
                AlignStringC[AlignStringLength] = ResidueB;
                if (ResidueA == ResidueB)
                {
                    AlignStringB[AlignStringLength] = '*';
                }
                else
                {
                    AlignStringB[AlignStringLength] = ' ';
                }
                break;
            case Z_GAP_IN_X:
                AlignStringA[AlignStringLength] = XSequence[X];
                AlignStringB[AlignStringLength] = ' ';
                AlignStringC[AlignStringLength] = '-';
                break;
            case Z_GAP_IN_Y:
                AlignStringA[AlignStringLength] = '-';
                AlignStringB[AlignStringLength] = ' ';
                AlignStringC[AlignStringLength] = YSequence[Y];
                break;
            }
            break;
        case Z_GAP_IN_X:
            if (NextZ[TableIndex])
            {
                AlignStringA[AlignStringLength] = XSequence[X];
                AlignStringB[AlignStringLength] = ' ';
                AlignStringC[AlignStringLength] = '-';
            }
            break;            
        case Z_GAP_IN_Y:
            if (NextZ[TableIndex])
            {
                AlignStringA[AlignStringLength] = '-';
                AlignStringB[AlignStringLength] = ' ';
                AlignStringC[AlignStringLength] = YSequence[Y];
            }
            break;
        }
        
        // Each step we take will add to the string...except closing a gap.
        if (Z && !NextZ[TableIndex])
        {
            AlignStringLength++;
        }
        X = NextX[TableIndex];
        Y = NextY[TableIndex];
        Z = NextZ[TableIndex];
    }
    printf("Alignment score %d.  Alignment follows:\n", BestScore);
    printf("%s\n", AlignStringA);
    printf("%s\n", AlignStringB);
    printf("%s\n", AlignStringC);
    
    ////////////////////////////////////////////////////////////
    // cleanup:
    SafeFree(ScoringMatrix);
    SafeFree(ScoreTable);
    SafeFree(ExonOffsetsA);
    SafeFree(ExonOffsetsB);
    SafeFree(NextX);
    SafeFree(NextY);
    SafeFree(NextZ);
    SafeFree(YSequence);
    SafeFree(XSequence);
    SafeFree(ExonEdgeOffsetsA);
    SafeFree(ExonEdgeOffsetsB);
    if (AlignStringA)
    {
        SafeFree(AlignStringA);
        SafeFree(AlignStringB);
        SafeFree(AlignStringC);
    }
    FreePrevCellTable(PrevY, MaxY);
    FreePrevCellTable(PrevX, MaxX);
    for (Y = 0; Y < MaxY; Y++)
    {
        SafeFree(RowInfoB[Y]);
    }
    SafeFree(RowInfoB);
    for (X = 0; X < MaxX; X++)
    {
        SafeFree(RowInfoA[X]);
    }
    SafeFree(RowInfoA);
    return BestScore;
}


