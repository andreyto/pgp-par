
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


// SpliceDB.c constructs a splice-tolerant database, starting from a collection of INTERVALS with LINKS.
// Translated from the original Python script, CollectExons.py, for efficiency
#include "CMemLeak.h"
#include "Utils.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "Trie.h"
#include "Inspect.h"
#include "Spliced.h"
#include "SpliceDB.h"
#include "SNP.h"

#define MAX_INTERVALS_PER_GENE 2000
#define MAX_INTERVAL_LENGTH 100000

// Reject any splice junctions from an EST which don't have this good of 
// a splice consensus score...unless we see more than one EST supporting
// the junction.
#define SPLICE_SIGNAL_SCORE_CUTOFF -15
#define DEFAULT_MINIMUM_ORF_LENGTH 50
#define IFLAG_FRAME_0 1
#define IFLAG_FRAME_1 2
#define IFLAG_FRAME_2 4
#define IFLAG_ALL_READING_FRAMES 7

// We might be parsing many, many intervals from two different sources...which means
// they'll be read out-of-order.  We maintain an index g_IntervalIndex such that
// g_IntervalIndex[n] is the first interval whose start is at least n*10000.
// When it comes time to insert a new interval, we check g_IntervalIndex.  If 
// the entry is NULL, the interval goes at the end of the global list.  If the
// interval isn't null, the interval goes NEAR that interval (maybe a little earlier, 
// maybe a little later, but the scan is cheap and that's the key idea)
IntervalNode** g_IntervalIndex = NULL;

// Linked list of all intervals in a chromosome (+ orientation):
IntervalNode* g_FirstInterval;
IntervalNode* g_LastInterval;

// Linked list of GeneNode structs for the current gene. 
GeneNode* g_GeneFirst;
GeneNode* g_GeneLast;
int GeneNodeCount;

// Int variables for reporting statistics on the database generation algorithms:
int g_StatsIncompleteGeneCount = 0;
int g_StatsLargestGeneSize = 0;
int g_StatsLargestGeneRecordNumber = 0;
int g_StatsIntervalsBeforeMerge = 0;
int g_StatsEdgesBeforeMerge = 0;
int g_StatsIntervalsAfterMerge = 0;
int g_StatsEdgesAfterMerge = 0;
int g_StatsIntervalsAfterIntersect = 0;
int g_StatsEdgesAfterIntersect = 0;
int g_StatsTotalExonsWritten = 0;
int g_StatsTotalEdgesWritten = 0;

// Forward declarations:
IntervalNode* InsertIntervalBefore(IntervalNode* Interval, IntervalNode* Before);
IntervalNode* InsertIntervalAfter(IntervalNode* Interval, IntervalNode* After);
int BuildAndWriteExons(FILE* GenomicFile, FILE* OutputFile, int ReverseFlag, char* GeneName, int ChromosomeNumber, int MinORFLength);
void MaskBrokenSequence(char* Protein, int MinORFLength);
void IntegrityCheckGene();
void PruneShortORFs(int ReverseFlag, int MinimumORFLength);
void DeleteExonLink(ExonNode* Exon, ExonLink* Link, int ForwardFlag);
void PurgeNonCodingExonChunks();
void GenomeDAGLinkBack(GenomeDAGNode* DAGNode, GenomeDAGNode* BackDAGNode, int Count);
void FreeIntervalExons(IntervalNode* Interval);

// Print the current GeneNode list to stdout, for debugging.
void DebugPrintBuiltGene()
{
    GeneNode* GNode;
    ExonNode* Exon;
    ExonNode* OtherExon;
    ExonLink* Link;
    char Buffer[512];
    int IntervalLen;
    int TrueLen;
    GenomeDAGNode* DAGNode;
    GenomeDAGLink* DAGLink;
    int DAGNodeIndex;
    //
    printf("\n--== Current gene ==--\n");
    for (GNode = g_GeneFirst; GNode; GNode = GNode->Next)
    {
        printf("  Interval from %d to %d flag %d\n", GNode->Interval->Start, GNode->Interval->End, GNode->Interval->Flags);
        for (DAGNodeIndex = 0; DAGNodeIndex < GNode->Interval->DAGNodeCount; DAGNodeIndex++)
        {
            DAGNode = GNode->Interval->DAGNodes + DAGNodeIndex;
            // Skip any extra allocation (null DAG nodes)
            if (!DAGNode->Sequence)
            {
                continue;
            }
            printf("    DAG node from %d to %d (%d bases)\n", DAGNode->Start, DAGNode->End , DAGNode->End - DAGNode->Start);
            for (DAGLink = DAGNode->FirstBack; DAGLink; DAGLink = DAGLink->Next)
            {
                printf("    << Link back %d to DAG node %d-%d\n", DAGLink->Count, DAGLink->Node->Start, DAGLink->Node->End);
            }
            for (DAGLink = DAGNode->FirstForward; DAGLink; DAGLink = DAGLink->Next)
            {
                printf("    >> Link forw %d to DAG node %d-%d\n", DAGLink->Count, DAGLink->Node->Start, DAGLink->Node->End);
            }
        }
        for (Exon = GNode->Interval->FirstExon; Exon; Exon = Exon->Next)
        {
            printf("    Exon from %d to %d (%dAA)\n", Exon->Start, Exon->End, Exon->Length);
            IntervalLen = Exon->End - Exon->Start;
            if (IntervalLen <= 0)
            {
                printf("** WARNING: Exon is * * E M P T Y * * \n");
            }
            TrueLen = strlen(Exon->Prefix) + Exon->Length*3 + strlen(Exon->Suffix);
            if (IntervalLen != TrueLen)
            {
                printf("** Warning: %d-%d is length %d, but true exon length is %d+%d+%d\n", 
                    Exon->Start, Exon->End, IntervalLen, 
                    strlen(Exon->Prefix), Exon->Length*3, strlen(Exon->Suffix));
            }
            if (Exon->Sequence)
            {
                strncpy(Buffer, Exon->Sequence, 512);
                Buffer[511] = '\0';
                printf("    Sequence(partial): %s\n", Buffer);
            }
            for (Link = Exon->FirstBack; Link; Link = Link->Next)
            {
                OtherExon = Link->Exon;
                printf("      Link back %d to an exon from %d to %d (%dAA)\n", Link->Power, 
                    OtherExon->Start, OtherExon->End, OtherExon->Length);
                if (OtherExon->Sequence)
                {
                    strncpy(Buffer, OtherExon->Sequence, 50);
                    Buffer[50] = '\0';
                    printf("      Ls: %s\n", Buffer);
                }
                if ((OtherExon->Start != Exon->End) && (OtherExon->Start != Exon->End + 1) && 
                    (OtherExon->End != Exon->Start) && (OtherExon->End != Exon->Start - 1))
                {
                    if (!Link->Power)
                    {
                        printf("** Warning: Link with no power!\n");
                    }
                }
            }
            
            for (Link = Exon->FirstForward; Link; Link = Link->Next)
            {
                OtherExon = Link->Exon;
                printf("      Link forward %d to an exon from %d to %d (%dAA)\n", Link->Power, OtherExon->Start, 
                    OtherExon->End, OtherExon->Length);
                if (OtherExon->Sequence)
                {
                    strncpy(Buffer, OtherExon->Sequence, 50);
                    Buffer[50] = '\0';
                    printf("      Ls: %s\n", Buffer);
                }
                if ((OtherExon->Start != Exon->End) && (OtherExon->Start != Exon->End + 1) && 
                    (OtherExon->End != Exon->Start) && (OtherExon->End != Exon->Start - 1))
                {
                    if (!Link->Power)
                    {
                        printf("** Warning: Link with no power!\n");
                    }
                }

            }
        }
    }
}

// Print all intervals to stdout.  (VERY verbose, if done for the whole chromosome!)
void DebugPrintIntervals(int IncludeLinks, int CountingFlag,
    int CoverageStart, int CoverageEnd)
{
    IntervalNode* Interval;
    EdgeNode* Edge;
    int IntervalCount = 0;
    int ForwardCount = 0;
    int BackwardCount = 0;
    int IForwardCount = 0;
    int IBackwardCount = 0;

    printf("\n\n=-=-=-=-=- Intervals =-=-=-=-=-\n");
    for (Interval = g_FirstInterval; Interval; Interval = Interval->Next)
    {
        IntervalCount++;
        // Skip output of intervals not in the range CoverageStart...CoverageEnd
        if (CoverageStart >= 0 && Interval->End < CoverageStart)
        {
            continue;
        }
        if (CoverageEnd >= 0 && Interval->Start > CoverageEnd)
        {
            continue;
        }
        if (IncludeLinks >= 0)
        {
            IForwardCount = 0;
            IBackwardCount = 0;
            for (Edge = Interval->FirstForward; Edge; Edge = Edge->Next)
            {
                IForwardCount++;
            }
            for (Edge = Interval->FirstBack; Edge; Edge = Edge->Next)
            {
                IBackwardCount++;
            }
            printf("%d-%d %d <%d >%d\n", Interval->Start, Interval->End, Interval->Occurrences,
                IBackwardCount, IForwardCount);
        }
        
        for (Edge = Interval->FirstForward; Edge; Edge = Edge->Next)
        {
            if (IncludeLinks > 0)
            {
                printf("  -> %d-%d (%d)\n", Edge->Interval->Start, Edge->Interval->End, Edge->Count);
            }
            if (Edge->Interval->Start < Interval->Start)
            {
                printf("** Corruption: Forward link goes to an interval EARLIER along the chrom\n");
                printf("** Start %d-%d, edge to %d-%d\n", Interval->Start, Interval->End, 
                    Edge->Interval->Start, Edge->Interval->End);
            }
            ForwardCount++;
        }
        for (Edge = Interval->FirstBack; Edge; Edge = Edge->Next)
        {
            if (IncludeLinks > 0)
            {
                printf("  <- %d-%d (%d)\n", Edge->Interval->Start, Edge->Interval->End, Edge->Count);
            }
            if (Edge->Interval->Start > Interval->Start)
            {
                printf("** Corruption: Forward link goes to an interval EARLIER along the chrom\n");
                printf("** Start %d-%d, edge to %d-%d\n", Interval->Start, Interval->End, 
                    Edge->Interval->Start, Edge->Interval->End);
            }
            BackwardCount++;
        }
    }
    printf("Total: %d intervals, %d forward links, %d backward links.\n", IntervalCount, ForwardCount, BackwardCount);
    switch (CountingFlag)
    {
    case 1:
        g_StatsIntervalsBeforeMerge = IntervalCount;
        g_StatsEdgesBeforeMerge = ForwardCount;
        break;
    case 2:
        g_StatsIntervalsAfterMerge = IntervalCount;
        g_StatsEdgesAfterMerge = ForwardCount;
        break;
    case 3:
        g_StatsIntervalsAfterIntersect = IntervalCount;
        g_StatsEdgesAfterIntersect = ForwardCount;
        break;        
    default:
        break;
    }
}

// Add a new interval to the master list.  Or, if that interval has already been
// seen, increment its count.  We use g_IntervalIndex to jump to *approximately* the right place 
// in the global list of intervals, then scan forward or backward to find exactly the right spot.
IntervalNode* AddInterval(int Start, int End, int Flags)
{
    IntervalNode* OldInterval;
    IntervalNode* NewInterval;
    int Bin;
    int IterateBin;
    //
    if (!g_IntervalIndex)
    {
        // Somewhat hacky: Hard-coded size of 25000, large enough to cover human chromosome #1
        g_IntervalIndex = (IntervalNode**)calloc(25000, sizeof(IntervalNode*));
    }
    Bin = Start / 10000;
    OldInterval = g_IntervalIndex[Bin];
    if (!OldInterval)
    {
        // This interval's start position is larger than any seen before! 
        // Insert the interval at the end of the global list:
        NewInterval = (IntervalNode*)calloc(sizeof(IntervalNode), 1);
        NewInterval->Occurrences = 1;
        NewInterval->Start = Start;
        NewInterval->End = End;
        NewInterval->Flags = Flags;
        if (g_LastInterval)
        {
            g_LastInterval->Next = NewInterval;
        }
        NewInterval->Prev = g_LastInterval;
        if (!g_FirstInterval)
        {
            g_FirstInterval = NewInterval;
        }
        g_LastInterval = NewInterval;
        // Update the index:
        for (IterateBin = Bin; IterateBin >= 0; IterateBin--)
        {
            if (g_IntervalIndex[IterateBin])
            {
                break;
            }
            g_IntervalIndex[IterateBin] = NewInterval;
        }
        return NewInterval;
    }
    // Next case: OldInterval is exactly right:
    if (Start == OldInterval->Start && End == OldInterval->End)
    {
        OldInterval->Occurrences++;
        OldInterval->Flags |= Flags;
        return OldInterval;
    }
    // Next case: OldInterval precedes this interval.
    if (Start > OldInterval->Start || (Start == OldInterval->Start && End > OldInterval->End))
    {
        // Iterate forward until OldInterval is NULL or OldInterval comes AFTER the new interval:
        for (; OldInterval; OldInterval = OldInterval->Next)
        {
            if (OldInterval->Start > Start)
            {
                break;
            }
            if (OldInterval->Start == Start && OldInterval->End > End)
            {
                break;
            }
            if (OldInterval->Start == Start && OldInterval->End == End)
            {
                OldInterval->Occurrences++;
                OldInterval->Flags |= Flags;
                return OldInterval;
            }
        }
        NewInterval = (IntervalNode*)calloc(sizeof(IntervalNode), 1);
        NewInterval->Occurrences = 1;
        NewInterval->Start = Start;
        NewInterval->End = End;
        NewInterval->Flags = Flags;
        if (!OldInterval)
        {
            // The new interval comes at the END of the list:
            if (g_LastInterval)
            {
                g_LastInterval->Next = NewInterval;
            }
            NewInterval->Prev = g_LastInterval;
            g_LastInterval = NewInterval;
        }
        else
        {
            // Insert new interval just before OldInterval:
            if (OldInterval->Prev)
            {
                OldInterval->Prev->Next = NewInterval;
            }
            NewInterval->Prev = OldInterval->Prev;
            NewInterval->Next = OldInterval;
            OldInterval->Prev = NewInterval;
        }
        return NewInterval;
    }
    else
    {
        // Last case: The new interval immediately precedes OldInterval.
        NewInterval = (IntervalNode*)calloc(sizeof(IntervalNode), 1);
        NewInterval->Occurrences = 1;
        NewInterval->Start = Start;
        NewInterval->End = End;
        NewInterval->Flags = Flags;
        if (OldInterval->Prev)
        {
            OldInterval->Prev->Next = NewInterval;
        }
        if (g_FirstInterval == OldInterval)
        {
            g_FirstInterval = NewInterval;
        }
        NewInterval->Prev = OldInterval->Prev;
        NewInterval->Next = OldInterval;
        OldInterval->Prev = NewInterval;
        for (IterateBin = Bin; IterateBin >= 0; IterateBin--)
        {
            if (g_IntervalIndex[IterateBin] && (g_IntervalIndex[IterateBin]->Start < Start || (g_IntervalIndex[IterateBin]->Start == Start && g_IntervalIndex[IterateBin]->End < End)))
            {
                break;
            }
            g_IntervalIndex[IterateBin] = NewInterval;
        }
        return NewInterval;
    }
    //if (!FirstInterval)
    //{
    //    NewInterval = (IntervalNode*)calloc(sizeof(IntervalNode), 1);
    //    NewInterval->Occurrences = 1;
    //    NewInterval->Start = StartPos;
    //    NewInterval->End = EndPos;
    //    NewInterval->Flags = Flags;
    //    FirstInterval = NewInterval;
    //    LastInterval = NewInterval;
    //    return NewInterval;
    //}
    //// After this loop, Interval is the last one before the new guy (or NULL, if the new guy
    //// belongs at the start of the list), and NextInterval is the first one after the
    //// new guy (or NULL, if the new guy belongs at the end of the list).
    //for (Interval = LastInterval; Interval; Interval = Interval->Prev)
    //{
    //    if (Interval->Start == StartPos)
    //    {
    //        if (Interval->End > EndPos)
    //        {
    //            NextInterval = Interval;
    //            continue;
    //        }
    //        if (Interval->End == EndPos)
    //        {
    //            Interval->Occurrences++;
    //            Interval->Flags |= Flags;
    //            return Interval;
    //        }
    //        break;
    //    }
    //    if (Interval->Start < StartPos)
    //    {
    //        break;
    //    }
    //    NextInterval = Interval;
    //}

    //NewInterval = (IntervalNode*)calloc(sizeof(IntervalNode), 1);
    //NewInterval->Occurrences = 1;
    //NewInterval->Start = StartPos;
    //NewInterval->End = EndPos;
    //NewInterval->Flags = Flags;

    //if (!Interval)
    //{
    //    FirstInterval->Prev = NewInterval;
    //    NewInterval->Next = FirstInterval;
    //    FirstInterval = NewInterval;
    //    return NewInterval;
    //}
    //Interval->Next = NewInterval;
    //NewInterval->Prev = Interval;
    //NewInterval->Next = NextInterval;
    //if (NextInterval)
    //{
    //    NextInterval->Prev = NewInterval;
    //}
    //else
    //{
    //    LastInterval = NewInterval;
    //}
    //return NewInterval;
    
}

// Link forward from interval A to interval B.
void LinkIntervals(IntervalNode* A, IntervalNode* B, int Count, float Score)
{
    EdgeNode* OldEdge;
    EdgeNode* NewEdge;
    int Linked;
    //
    Linked = 0;
    for (OldEdge = A->FirstForward; OldEdge; OldEdge = OldEdge->Next)
    {
        if (OldEdge->Interval == B)
        {
            OldEdge->Count += Count;
            Linked = 1;
        }
    }
    if (!Linked)
    {
        NewEdge = (EdgeNode*)calloc(sizeof(EdgeNode), 1);
        NewEdge->Count = Count;
        NewEdge->Score = Score;
        NewEdge->Interval = B;
        if (!A->FirstForward)
        {
            A->FirstForward = NewEdge;
        }
        else
        {
            A->LastForward->Next = NewEdge;
            NewEdge->Prev = A->LastForward;
        }
        A->LastForward = NewEdge;
    }
    Linked = 0;
    for (OldEdge = B->FirstBack; OldEdge; OldEdge = OldEdge->Next)
    {
        if (OldEdge->Interval == A)
        {
            OldEdge->Count += Count;
            Linked = 1;
        }
    }
    if (!Linked)
    {
        NewEdge = (EdgeNode*)calloc(sizeof(EdgeNode), 1);
        NewEdge->Count = Count;
        NewEdge->Score = Score;
        NewEdge->Interval = A;
        if (!B->FirstBack)
        {
            B->FirstBack = NewEdge;
        }
        else
        {
            B->LastBack->Next = NewEdge;
            NewEdge->Prev = B->LastBack;
        }
        B->LastBack = NewEdge;
    }
}

// Copied-and-modified from ParseIntervalsESTBinaryFile.
// Parse intervals from a binary file, with reading-frame flags attached.
void ParseIntervalsGeneFindBinaryFile(char* FileName)
{
    FILE* File;
    IntervalNode* Interval;
    IntervalNode* BackInterval;
    int Start;
    int End;
    int BytesRead;
    int Score;
    int Flags;
    int FilePos = 0;
    int TotalIntervals = 0;
    int JunctionCount = 0;
    int JunctionIndex;
    int JunctionStart;
    float JunctionScore;
    int TotalJunctions = 0;
    int BackIntervalFound;
    //
    File = fopen(FileName, "rb");
    if (!File)
    {
        printf("** Error in ParseIntervalsBinaryFile: Can't open '%s'\n", FileName);
        return;
    }
    while (1)
    {

        BytesRead = ReadBinary(&Start, sizeof(int), 1, File);
        if (!BytesRead)
        {
            break;
        }
        FilePos += BytesRead;
        BytesRead += ReadBinary(&End, sizeof(int), 1, File);
        // Sanity check:
        if (Start<0 || End<0 || End<=Start)
        {
            printf("** BARF: Gene finder output reports interval from %d to %d!\n", Start, End);
        }
        BytesRead += ReadBinary(&Flags, sizeof(int), 1, File);
        BytesRead += ReadBinary(&Score, sizeof(int), 1, File);
        
        Interval = AddInterval(Start, End, Flags);
        TotalIntervals++;
        //FilePos += ReadBinary(&Interval->Occurrences, sizeof(int), 1, File);
        FilePos += ReadBinary(&JunctionCount, sizeof(int), 1, File);
        // Read a list of junctions that END at this interval.
        for (JunctionIndex = 0; JunctionIndex < JunctionCount; JunctionIndex++)
        {
            FilePos += ReadBinary(&JunctionStart, sizeof(int), 1, File);
            //FilePos += ReadBinary(&JunctionOccurrences, sizeof(int), 1, File);
            FilePos += ReadBinary(&JunctionScore, sizeof(float), 1, File);
            // Right here is where we filter crummy splice junctions:
            //if (JunctionOccurrences < 2 && JunctionScore < SPLICE_SIGNAL_SCORE_CUTOFF)
            //{
            //    continue;
            //}
            TotalJunctions++;
            // Find an interval which ends at the junction's splice point:
            BackIntervalFound = 0;
            for (BackInterval = Interval->Prev; BackInterval; BackInterval = BackInterval->Prev)
            {
                if (BackInterval->End == JunctionStart)
                {
                    BackIntervalFound = 1;
                    LinkIntervals(BackInterval, Interval, 1, JunctionScore);
                    break;
                }
            }
            if (!BackIntervalFound)
            {
                printf("** Warning: Found a junction with no back-interval!\n");
                printf("  Junction %d %f\n", JunctionStart,  JunctionScore);
                printf("  Interval %d-%d\n", Interval->Start, Interval->End);
            }
        }
    }
    fclose(File);
    printf("Read %d intervals, %d junctions.\n", TotalIntervals, TotalJunctions);
}

// Parse intervals from a binary file containing Interval records. Each Interval
// record may contain a list of Junction records.
// Interval record: IntervalStart, IntervalEnd, IntervalCount, JunctionCount
// Junction record: Start, Count, Score
// The list contains all junctions that END at the START of this interval.  (That way,
// we can look *back* through the list to find the splice donor)
// We can filter any junctions that don't have a good occurrence-Count or a good 
// consensus splice signal Score.
// Note that EST intervals have no particular reading frame specified.
void ParseIntervalsESTBinaryFile(char* FileName)
{
    FILE* File;
    IntervalNode* Interval;
    IntervalNode* BackInterval;
    int Start;
    int End;
    int BytesRead;
    int FilePos = 0;
    int TotalIntervals = 0;
    int JunctionCount;
    int JunctionIndex;
    int JunctionStart;
    int JunctionOccurrences;
    float JunctionScore;
    int BackIntervalFound;
    int TotalJunctions = 0;
    //
    File = fopen(FileName, "rb");
    if (!File)
    {
        printf("** Error in ParseIntervalsBinaryFile: Can't open '%s'\n", FileName);
        return;
    }
    while (1)
    {

        BytesRead = ReadBinary(&Start, sizeof(int), 1, File);
        if (!BytesRead)
        {
            break;
        }
        FilePos += BytesRead;
        BytesRead += ReadBinary(&End, sizeof(int), 1, File);
        Interval = AddInterval(Start, End, IFLAG_ALL_READING_FRAMES);
        TotalIntervals++;
        FilePos += ReadBinary(&Interval->Occurrences, sizeof(int), 1, File);
        FilePos += ReadBinary(&JunctionCount, sizeof(int), 1, File);
        for (JunctionIndex = 0; JunctionIndex < JunctionCount; JunctionIndex++)
        {
            FilePos += ReadBinary(&JunctionStart, sizeof(int), 1, File);
            FilePos += ReadBinary(&JunctionOccurrences, sizeof(int), 1, File);
            FilePos += ReadBinary(&JunctionScore, sizeof(float), 1, File);
            // Right here is where we filter crummy splice junctions:
            if (JunctionOccurrences < 2 && JunctionScore < SPLICE_SIGNAL_SCORE_CUTOFF)
            {
                continue;
            }
            TotalJunctions++;
            // Find an interval which ends at the junction's splice point:
            BackIntervalFound = 0;
            for (BackInterval = Interval->Prev; BackInterval; BackInterval = BackInterval->Prev)
            {
                if (BackInterval->End == JunctionStart)
                {
                    BackIntervalFound = 1;
                    LinkIntervals(BackInterval, Interval, JunctionOccurrences, JunctionScore);
                    break;
                }
            }
            if (!BackIntervalFound)
            {
                printf("** Warning: Found a junction with no back-interval!\n");
                printf("  Junction %d %d %f\n", JunctionStart, JunctionOccurrences, JunctionScore);
                printf("  Interval %d-%d\n", Interval->Start, Interval->End);
            }
        }
    }
    fclose(File);
    printf("Read %d intervals, %d junctions.\n", TotalIntervals, TotalJunctions);
}

// B inherits all backward links from A.
// before:           after:
//  C<->A           C<-\   A
//      B               \->B
void AssimilateLinksBack(IntervalNode* A, IntervalNode* B)
{
    EdgeNode* NodeA;
    EdgeNode* PrevA;
    EdgeNode* NodeB;
    EdgeNode* NodeC;
    EdgeNode* Next;
    IntervalNode* C;
    int Found;
    int ACStrength = 0;
    int BCStrength = 0;
    //
    
    for (NodeA = A->FirstBack; NodeA; NodeA = NodeA->Next)
    {
        ACStrength = NodeA->Count;
        BCStrength = 0;
        // Ensure that B has a link to this target:
        Found = 0;
        for (NodeB = B->FirstBack; NodeB; NodeB = NodeB->Next)
        {
            if (NodeB->Interval == NodeA->Interval)
            {
                BCStrength = NodeB->Count;
                //NodeB->Count += ACStrength; // Counts are already full
                Found = 1;
                break;
            }
        }
        // If B didn't already link back to the target, add an EdgeNode to B's list:
        if (!Found)
        {
            NodeB = (EdgeNode*)calloc(sizeof(EdgeNode), 1);
            NodeB->Interval = NodeA->Interval;
            NodeB->Count = ACStrength;
            if (B->LastBack)
            {
                NodeB->Prev = B->LastBack;
                B->LastBack->Next = NodeB;
            }
            else
            {
                B->FirstBack = NodeB;
            }
            B->LastBack = NodeB;
        }
        // Switch C to point to B.  It's possible that C has a pointer to B already, in which case
        // we'll free the old one.
        C = NodeA->Interval;
        for (NodeC = C->FirstForward; NodeC; NodeC = NodeC->Next)
        {
            if (NodeC->Interval == B)
            {
                //FoundCount += NodeC->Count;
                Next = NodeC->Next;
                if (Next)
                {
                    Next->Prev = NodeC->Prev;
                }
                if (NodeC->Prev)
                {
                    NodeC->Prev->Next = Next;
                }
                if (C->FirstForward == NodeC)
                {
                    C->FirstForward = NodeC->Next;
                }
                if (C->LastForward == NodeC)
                {
                    C->LastForward = NodeC->Prev;
                }
                SafeFree(NodeC);
                NodeC = Next;
                if (!NodeC)
                {
                    break;
                }
            }
        }
        Found = 0;
        for (NodeC = C->FirstForward; NodeC; NodeC= NodeC->Next)
        {
            if (NodeC->Interval == A)
            {
                //NodeC->Count += FoundCount;
                //NodeC->Count = ACStrength + BCStrength; // Counts are already full
                NodeC->Count = ACStrength;
                NodeC->Interval = B;
                Found = 1;
            }
        }
        if (!Found)
        {
            printf("*** Corruption!  %d-%d should link forward to %d-%d\n", C->Start, C->End, A->Start, A->End);
        }
    }
    
    // Free the old nodes:
    PrevA = NULL;
    for (NodeA = A->FirstBack; NodeA; NodeA = NodeA->Next)
    {
        SafeFree(PrevA);
        PrevA = NodeA;
    }
    SafeFree(PrevA);
    A->FirstBack = NULL;
    A->LastBack = NULL;

}

// B inherits all the forward links from A.
void AssimilateLinksForward(IntervalNode* A, IntervalNode* B)
{
    EdgeNode* NodeA;
    EdgeNode* PrevA;
    EdgeNode* NodeB;
    EdgeNode* NodeC;
    EdgeNode* Next;
    IntervalNode* C;
    int Found;
    int ACStrength = 0;
    int BCStrength = 0;
    //
    for (NodeA = A->FirstForward; NodeA; NodeA = NodeA->Next)
    {
        ACStrength = NodeA->Count;
        BCStrength = 0;
        // Ensure that B has a link to this target:
        Found = 0;
        for (NodeB = B->FirstForward; NodeB; NodeB = NodeB->Next)
        {
            if (NodeB->Interval == NodeA->Interval)
            {
                BCStrength = NodeB->Count;
                //NodeB->Count += ACStrength; // Counts are already full
                Found = 1;
                break;
            }
        }
        if (!Found)
        {
            NodeB = (EdgeNode*)calloc(sizeof(EdgeNode), 1);
            NodeB->Interval = NodeA->Interval;
            NodeB->Count = ACStrength;
            if (B->LastForward)
            {
                NodeB->Prev = B->LastForward;
                B->LastForward->Next = NodeB;
            }
            else
            {
                B->FirstForward = NodeB;
            }
            B->LastForward = NodeB;
        }
        // Switch C to point to B.  It's possible that C has a pointer to B already, in which case
        // we'll free the old one.
        C = NodeA->Interval;
        for (NodeC = C->FirstBack; NodeC; NodeC= NodeC->Next)
        {
            if (NodeC->Interval == B)
            {
                //FoundCount += NodeC->Count;
                Next = NodeC->Next;
                if (Next)
                {
                    Next->Prev = NodeC->Prev;
                }
                if (NodeC->Prev)
                {
                    NodeC->Prev->Next = Next;
                }
                if (C->FirstBack == NodeC)
                {
                    C->FirstBack = NodeC->Next;
                }
                if (C->LastBack == NodeC)
                {
                    C->LastBack = NodeC->Prev;
                }
                SafeFree(NodeC);
                NodeC = Next;
                if (!NodeC)
                {
                    break;
                }
            }
        }
        Found = 0;
        for (NodeC = C->FirstBack; NodeC; NodeC = NodeC->Next)
        {
            if (NodeC->Interval == A)
            {
                //NodeC->Count = ACStrength + BCStrength;// Counts are already full
                NodeC->Count = ACStrength;
                NodeC->Interval = B;
                Found = 1;
            }
        }
        if (!Found)
        {
            printf("*** Corruption!  %d-%d should link backward to %d-%d\n", C->Start, C->End, A->Start, A->End);
        }
    }
    
    // Free the old nodes:
    PrevA = NULL;
    for (NodeA = A->FirstForward; NodeA; NodeA = NodeA->Next)
    {
        SafeFree(PrevA);
        PrevA = NodeA;
    }
    SafeFree(PrevA);
    A->FirstForward = NULL;
    A->LastForward = NULL;
}
void FreeIntervalDAG(IntervalNode* Interval)
{
    int DAGNodeIndex;
    GenomeDAGNode* DAGNode;
    GenomeDAGLink* Link;
    GenomeDAGLink* PrevLink;
    //
    if (!Interval || !Interval->DAGNodes)
    {
        return;
    }
    for (DAGNodeIndex = 0; DAGNodeIndex < Interval->DAGNodeCount; DAGNodeIndex++)
    {
        DAGNode = Interval->DAGNodes + DAGNodeIndex;
        SafeFree(DAGNode->Sequence);
        SafeFree(DAGNode->Exons);
        // Free links back:
        PrevLink = NULL;
        for (Link = DAGNode->FirstBack; Link; Link = Link->Next)
        {
            SafeFree(PrevLink);
            PrevLink = Link;
        }
        SafeFree(PrevLink);
        // Free links forward:
        PrevLink = NULL;
        for (Link = DAGNode->FirstForward; Link; Link = Link->Next)
        {
            SafeFree(PrevLink);
            PrevLink = Link;
        }
        SafeFree(PrevLink);
    }
    SafeFree(Interval->DAGNodes);
    Interval->DAGNodes = NULL;
    Interval->DAGNodeCount = 0;
}

void FreeInterval(IntervalNode* Interval)
{
    //
    FreeIntervalDAG(Interval);
    FreeIntervalExons(Interval);

    Interval->FirstForward = NULL;
    Interval->FirstBack = NULL;
    Interval->LastForward = NULL;
    Interval->LastBack = NULL;
    Interval->Start = -1;
    Interval->End = -1;
    SafeFree(Interval);
}

// Remove an interval from the master list.  And USUALLY, free the
// interval and its edges.  If DontFree is true, then don't free
// any memory yet.
void RemoveInterval(IntervalNode* Interval, int DontFree)
{
    EdgeNode* Prev;
    EdgeNode* Edge;
    EdgeNode* NeighborEdge;
    if (Interval == g_FirstInterval)
    {
        g_FirstInterval = Interval->Next;
    }
    if (Interval == g_LastInterval)
    {
        g_LastInterval = Interval->Prev;
    }
    if (Interval->Prev)
    {
        Interval->Prev->Next = Interval->Next;
    }
    if (Interval->Next)
    {
        Interval->Next->Prev = Interval->Prev;
    }
    if (!DontFree)
    {
        Prev = NULL;
        Edge = Interval->FirstForward;
        while (Edge)
        {
            SafeFree(Prev);
            Prev = Edge;
            // If someone points at us, free their pointer (to avoid corruption!)
            for (NeighborEdge = Edge->Interval->FirstBack; NeighborEdge; NeighborEdge = NeighborEdge->Next)
            {
                if (NeighborEdge->Interval == Interval)
                {
                    if (Edge->Interval->FirstBack == NeighborEdge)
                    {
                        Edge->Interval->FirstBack = NeighborEdge->Next;
                    }
                    if (Edge->Interval->LastBack == NeighborEdge)
                    {
                        Edge->Interval->LastBack = NeighborEdge->Prev;
                    }
                    if (NeighborEdge->Prev)
                    {
                        NeighborEdge->Prev->Next = NeighborEdge->Next;
                    }
                    if (NeighborEdge->Next)
                    {
                        NeighborEdge->Next->Prev = NeighborEdge->Prev;
                    }
                    SafeFree(NeighborEdge);
                    break;
                }
            }

            Edge = Edge->Next;
        }
        SafeFree(Prev);
        //
        Prev = NULL;
        Edge = Interval->FirstBack;
        while (Edge)
        {
            SafeFree(Prev);
            Prev = Edge;
            // If someone points at us, free their pointer (to avoid corruption!)
            for (NeighborEdge = Edge->Interval->FirstForward; NeighborEdge; NeighborEdge = NeighborEdge->Next)
            {
                if (NeighborEdge->Interval == Interval)
                {
                    if (Edge->Interval->FirstForward == NeighborEdge)
                    {
                        Edge->Interval->FirstForward = NeighborEdge->Next;
                    }
                    if (Edge->Interval->LastForward == NeighborEdge)
                    {
                        Edge->Interval->LastForward = NeighborEdge->Prev;
                    }
                    if (NeighborEdge->Prev)
                    {
                        NeighborEdge->Prev->Next = NeighborEdge->Next;
                    }
                    if (NeighborEdge->Next)
                    {
                        NeighborEdge->Next->Prev = NeighborEdge->Prev;
                    }
                    SafeFree(NeighborEdge);
                    break;
                }
            }
            Edge = Edge->Next;
        }
        SafeFree(Prev);
        FreeInterval(Interval);
    }
}

// Merge all redundant intervals.  Intervals which overlap, and
// have no incompatible edges, can be merged into one large(r) interval.
// The merged interval inherits all reading frames of the subintervals.
// (This could add some redundancy, but not much, especially if we later
// prune short ORFs)
void MergeIntervals()
{
    IntervalNode* MergeA;
    IntervalNode* NextMergeA;
    IntervalNode* MergeB;
    int MergePerformed = 0;
    int TotalMergesPerformed = 0;
    //
    
    NextMergeA = g_FirstInterval;
    while (1)
    {
        
        MergeA = NextMergeA;
        if (!MergeA)
        {
            break;
        }
        MergePerformed = 0;
        MergeB = MergeA->Next;
        while (1)
        {
            if (MergePerformed)
            {
                TotalMergesPerformed++;
                //DebugPrintIntervals(-1. 0);
                break;
            }
            if (!MergeB || MergeB->Start > MergeA->End)
            {
                NextMergeA = MergeA->Next;
                break;
            }
            // Case 0: Identical!
            if (MergeA->Start == MergeB->Start && MergeA->End == MergeB->End)
            {
                //printf("%% [0] Merge two identical intervals %d-%d\n", MergeA->Start, MergeA->End);
                AssimilateLinksBack(MergeB, MergeA);
                AssimilateLinksForward(MergeB, MergeA);
                MergeA->Occurrences += MergeB->Occurrences;
                MergeA->Flags |= MergeB->Flags;
                RemoveInterval(MergeB, 0);
                MergePerformed = 1;
            }
            // Case 1: Same starting point, A doesn't link forward:
            else if (MergeA->Start == MergeB->Start && !MergeA->FirstForward)
            {
                //printf("%% [1] Same starting point: %d-%d, %d-%d\n", MergeA->Start, MergeA->End, MergeB->Start, MergeB->End);
                AssimilateLinksBack(MergeA, MergeB);
                MergeB->Occurrences += MergeA->Occurrences;
                MergeB->Flags |= MergeA->Flags;
                NextMergeA = MergeA->Next;
                RemoveInterval(MergeA, 0);
                MergePerformed = 1;
            }
            // Case 2: Same ending point, B doesn't link backward:
            else if (MergeA->End == MergeB->End && !MergeB->FirstBack)
            {
                //printf("%% [2] Same ending point: %d-%d, %d-%d\n", MergeA->Start, MergeA->End, MergeB->Start, MergeB->End);
                AssimilateLinksForward(MergeB, MergeA);
                MergeA->Occurrences += MergeB->Occurrences;
                MergeA->Flags |= MergeB->Flags;
                NextMergeA = MergeA;
                RemoveInterval(MergeB, 0);
                MergePerformed = 1;
            }
            // Case 3: Full overlap, no links in B:
            else if (MergeA->Start < MergeB->Start && MergeB->End < MergeA->End && !MergeB->FirstBack && !MergeB->FirstForward)
            {
                //printf("%% [3] full overlap: %d-%d, %d-%d\n", MergeA->Start, MergeA->End, MergeB->Start, MergeB->End);
                MergeA->Occurrences += MergeB->Occurrences;
                MergeA->Flags |= MergeB->Flags;
                NextMergeA = MergeA;
                RemoveInterval(MergeB, 0);
                MergePerformed = 1;
            }
            // Case 4: 'proper' overlap, A no forward, B no backward:
            else if (MergeB->Start > MergeA->Start && MergeB->End > MergeA->End && !MergeA->FirstForward && !MergeB->FirstBack)
            {
                //printf("%% [4] Proper overlap: %d-%d, %d-%d\n", MergeA->Start, MergeA->End, MergeB->Start, MergeB->End);
                MergeA->End = MergeB->End;
                AssimilateLinksForward(MergeB, MergeA);
                MergeA->Occurrences += MergeB->Occurrences;
                MergeA->Flags |= MergeB->Flags;
                NextMergeA = MergeA;
                RemoveInterval(MergeB, 0);
                MergePerformed = 1;
            }
            else
            {
                // Default case: Non-mergable
                MergeB = MergeB->Next;
            }
        }
    } // Iterate MergeA
    printf("Performed a total of %d merges.\n", TotalMergesPerformed);
}

// If two intervals intersect, then we don't want to store separate exons for each one.  That would be a lot of 
// redundant sequence data!  So, after calling MergeIntervals, we call IntersectIntervals().  
// The routine IntersectIntervals will produce a (near-)minimal disjoint set of intervals covering all the 
// ESTs and splice boundaries we've ever seen.  The intersection interval inherits all reading frames from its
// parents.
void IntersectIntervals()
{
    IntervalNode* A;
    IntervalNode* NextA;
    IntervalNode* B;
    IntervalNode* C;
    IntervalNode* D;
    int IntersectPerformed = 0;
    //
    
    NextA = g_FirstInterval;
    while (1)
    {
        A = NextA;
        if (!A)
        {
            break;
        }
        if (IntersectPerformed)
        {
            //DebugPrintIntervals(1, 0);
        }
        IntersectPerformed = 0;
        B = A->Next;
        if (!B)
        {
            break;
        }
        // Easy case: B starts after A ends.  No intersection required:
        if (B->Start >= A->End)
        {
            NextA = A->Next;
            continue;
        }
        if (B->Start == A->Start && B->End == A->End)
        {
            //printf("%d-%d is identical to %d-%d\n", A->Start, A->End, B->Start, B->End);
            AssimilateLinksForward(B, A);
            AssimilateLinksBack(B, A);
            A->Occurrences += B->Occurrences;
            A->Flags |= B->Flags;
            NextA = A;
            RemoveInterval(B, 0);
            IntersectPerformed = 1;
            continue;
        }
        if (B->Start == A->Start)
        {
            // |----| A    
            // |-----------| B
            //
            // |----|------|
            //   A     C
            //printf("%d-%d has same START as %d-%d\n", A->Start, A->End, B->Start, B->End);
            C = (IntervalNode*)calloc(sizeof(IntervalNode), 1);
            C->Start = A->End;
            C->End = B->End;
            C->Occurrences = B->Occurrences;
            C->Flags = B->Flags;
            A->Flags |= B->Flags;
            AssimilateLinksForward(B, C);
            AssimilateLinksBack(B, A);
            LinkIntervals(A, C, 0, 0);
            RemoveInterval(B, 0);
            C = InsertIntervalAfter(C, A);
            NextA = A;
            IntersectPerformed = 1;
            continue;
        }
        if (B->End == A->End)
        {
            // |-----------| A    
            //       |-----| B
            //
            // |----|------|
            //   C     B

            //printf("%d-%d has same END as %d-%d\n", A->Start, A->End, B->Start, B->End);
            NextA = A->Prev;
            C = (IntervalNode*)calloc(sizeof(IntervalNode), 1);
            C->Start = A->Start;
            C->End = B->Start;
            C->Occurrences = A->Occurrences;
            C->Flags = A->Flags;
            B->Flags |= A->Flags;
            AssimilateLinksForward(A, B);
            AssimilateLinksBack(A, C);
            LinkIntervals(C, B, 0, 0);
            RemoveInterval(A, 0);
            C = InsertIntervalBefore(C, B);
            if (!NextA)
            {
                NextA = g_FirstInterval;
            }
            IntersectPerformed = 1;
            continue;
        }
        // |---------------|
        //       |---|
        //
        // |-----|---|-----|
        //    C    B    D
        if (B->End < A->End)
        {
            //printf("%d-%d CONTAINS %d-%d\n", A->Start, A->End, B->Start, B->End);
            C = (IntervalNode*)calloc(sizeof(IntervalNode), 1);
            C->Start = A->Start;
            C->End = B->Start;
            C->Flags = A->Flags;
            D = (IntervalNode*)calloc(sizeof(IntervalNode), 1);
            D->Start = B->End;
            D->End = A->End;
            D->Flags = A->Flags;
            B->Flags |= A->Flags;
            AssimilateLinksBack(A, C);
            AssimilateLinksForward(A, D);
            LinkIntervals(C, B, 0, 0);
            LinkIntervals(B, D, 0, 0);
            C = InsertIntervalBefore(C, B);
            D = InsertIntervalAfter(D, B);
            RemoveInterval(A, 0);
            NextA = C;
            IntersectPerformed = 1;
            continue;
        }
        // |-------------| A
        //          |---------| B
        //
        // |--------|----|----|
        //     C      B     D
        //printf("%d-%d has PROPER OVERLAP with %d-%d\n", A->Start, A->End, B->Start, B->End);
        C = (IntervalNode*)calloc(sizeof(IntervalNode), 1);
        C->Start = A->Start;
        C->End = B->Start;
        C->Occurrences = A->Occurrences;
        C->Flags = A->Flags;
        D = (IntervalNode*)calloc(sizeof(IntervalNode), 1);
        D->Start = A->End;
        D->End = B->End;
        D->Occurrences = B->Occurrences;
        D->Flags = B->Flags;
        B->Flags |= A->Flags;
        //B2 = (IntervalNode*)calloc(sizeof(IntervalNode), 1);
        AssimilateLinksBack(A, C);
        AssimilateLinksForward(B, D);
        AssimilateLinksForward(A, B);
        LinkIntervals(C, B, 0, 0);
        LinkIntervals(B, D, 0, 0);
        C = InsertIntervalBefore(C, B);
        D = InsertIntervalAfter(D, B);
        RemoveInterval(B, 1);
        B->End = A->End;
        RemoveInterval(A, 0);
        B = InsertIntervalBefore(B, D);
        NextA = C;
        IntersectPerformed = 1;
        continue;
    }
}

// Insert Interval into the master list.  It comes after After.
IntervalNode* InsertIntervalAfter(IntervalNode* Interval, IntervalNode* After)
{
    IntervalNode* Node;
    //
    Node = After;
    if (!Node)
    {
        Node = g_FirstInterval;
    }
    while (Node)
    {
        if (Node->Start > Interval->Start)
        {
            break;
        }
        if (Node->Start == Interval->Start)
        {
            if (Node->End == Interval->End)
            {
                AssimilateLinksForward(Interval, Node);
                AssimilateLinksBack(Interval, Node);
                Node->Occurrences += Interval->Occurrences;
                SafeFree(Interval);
                return Node;
            }
            if (Node->End > Interval->End)
            {
                break;
            }
        }
        Node = Node->Next;
    }
    // At this point, Node is the guy that Interval will be inserted before:
    if (!Node)
    {
        g_LastInterval->Next = Interval;
        Interval->Prev = g_LastInterval;
        g_LastInterval = Interval;
    }
    else
    {
        if (Node->Prev)
        {
            Node->Prev->Next = Interval;
        }
        Interval->Prev = Node->Prev;
        Interval->Next = Node;
        Node->Prev = Interval;
    }
    return Interval;
}

// Insert Interval into the master list.  It comes before Before.
IntervalNode* InsertIntervalBefore(IntervalNode* Interval, IntervalNode* Before)
{
    IntervalNode* Node;
    //
    Node = Before;
    if (!Node)
    {
        Node = g_LastInterval;
    }
    while (Node)
    {
        if (Node->Start < Interval->Start)
        {
            break;
        }
        if (Node->Start == Interval->Start)
        {
            if (Node->End == Interval->End)
            {
                AssimilateLinksForward(Interval, Node);
                AssimilateLinksBack(Interval, Node);
                Node->Occurrences += Interval->Occurrences;
                SafeFree(Interval);
                return Node;
            }
            if (Node->End < Interval->End)
            {
                break;
            }
        }
        Node = Node->Prev;
    }
    // At this point, Node is the guy that Interval will be inserted after:
    if (!Node)
    {
        g_FirstInterval->Prev = Interval;
        Interval->Next = g_FirstInterval;
        g_FirstInterval = Interval;
    }
    else
    {
        if (Node->Next)
        {
            Node->Next->Prev = Interval;
        }
        Interval->Next = Node->Next;
        Interval->Prev = Node;
        Node->Next = Interval;
    }
    return Interval;
}

// Add Interval to the current gene sometime after Start
GeneNode* AddIntervalToGeneAfter(GeneNode* Start, IntervalNode* Interval)
{
    GeneNode* Node;
    GeneNode* NewNode;
    //
    for (Node = Start; Node; Node = Node->Next)
    {
        if (Node->Interval->Start == Interval->Start)
        {
            // Already on list, good.
            return Node;
        }
        if (Node->Interval->Start > Interval->Start)
        {
            NewNode = (GeneNode*)calloc(sizeof(GeneNode), 1);
            NewNode->Interval = Interval;
            Interval->GNode = NewNode;
            if (Node->Prev)
            {
                Node->Prev->Next = NewNode;
                NewNode->Prev = Node->Prev;
            }
            NewNode->Next = Node;
            Node->Prev = NewNode;
            GeneNodeCount++;
            return NewNode;
        }
    }
    // We ran off the edge of the list without seeing something that comes after the new interval.
    // So, the new interval becomes the last one of the gene:
    NewNode = (GeneNode*)calloc(sizeof(GeneNode), 1);
    NewNode->Interval = Interval;
    Interval->GNode = NewNode;
    NewNode->Prev = g_GeneLast;
    g_GeneLast->Next = NewNode;
    g_GeneLast = NewNode;
    GeneNodeCount++;
    return NewNode;
}

// Add Interval to the current gene sometime before Start
GeneNode* AddIntervalToGeneBefore(GeneNode* Start, IntervalNode* Interval)
{
    GeneNode* Node;
    GeneNode* NewNode;
    //
    for (Node = Start; Node; Node = Node->Prev)
    {
        if (Node->Interval->Start == Interval->Start)
        {
            // Already on list, good.
            return Node;
        }
        if (Node->Interval->Start < Interval->Start)
        {
            NewNode = (GeneNode*)calloc(sizeof(GeneNode), 1);
            NewNode->Interval = Interval;
            Interval->GNode = NewNode;
            if (Node->Next)
            {
                Node->Next->Prev = NewNode;
                NewNode->Next = Node->Next;
            }
            NewNode->Prev = Node;
            Node->Next = NewNode;
            GeneNodeCount++;
            return NewNode;
        }
    }
    // We ran off the edge of the list without seeing something that comes after the new interval.
    // So, the new interval becomes the first one of the gene:
    NewNode = (GeneNode*)calloc(sizeof(GeneNode), 1);
    NewNode->Interval = Interval;
    Interval->GNode = NewNode;
    NewNode->Next = g_GeneFirst;
    g_GeneFirst->Prev = NewNode;
    g_GeneFirst = NewNode;
    GeneNodeCount++;
    return NewNode;
}

// Add new GeneNodes to handle any peptides that start in Node->Interval and extend forward
// GNode is the bookmark where we started the satisfaction effort, so when (if) we insert new nodes,
// we'll insert
int SatisfyIntervalForward(GeneNode* GNode, int CharsSoFar)
{
    EdgeNode* Edge;
    int Chars;
    GeneNode* SubGNode;
    int RX;
    int MinRX;
    //
    // If this node has already been satisfied, then return immediately:
    if (GNode->RX + CharsSoFar > 60)
    {
        return GNode->RX;
    }
    MinRX = 9999;
    // Iterate over all 'forward intervals' that this interval links to:
    for (Edge = GNode->Interval->FirstForward; Edge; Edge = Edge->Next)
    {
        // Find (or create) the GeneNode for the forward interval:
        SubGNode = Edge->Interval->GNode;
        if (!SubGNode)
        {
            SubGNode = AddIntervalToGeneAfter(GNode, Edge->Interval);
        }
        RX = Edge->Interval->End - Edge->Interval->Start;
        Chars = CharsSoFar + (Edge->Interval->End - Edge->Interval->Start);
        if (Chars < 60)
        {
            // We're not yet satisfied along this path, so continue adding intervals:
            RX += SatisfyIntervalForward(SubGNode, Chars);
        }
        MinRX = min(MinRX, RX);
    }
    // Sanity check: RX cannot decrease when you add more intervals, it can only improve.
    if (MinRX < GNode->RX)
    {
        printf("%d < %d???\n", MinRX, GNode->RX);
    }
    GNode->RX = MinRX;
    return MinRX;
}

// Add new GeneNodes to handle any peptides that start in Node->Interval and extend backward
int SatisfyIntervalBack(GeneNode* GNode, int CharsSoFar)
{
    EdgeNode* Edge;
    int Chars;
    GeneNode* SubGNode;
    int LX;
    int MinLX;
    //
    if (GNode->LX + CharsSoFar > 60)
    {
        return GNode->LX;
    }
    MinLX = 9999;
    for (Edge = GNode->Interval->FirstBack; Edge; Edge = Edge->Next)
    {
        SubGNode = Edge->Interval->GNode;
        if (!SubGNode)
        {
            SubGNode = AddIntervalToGeneBefore(GNode, Edge->Interval);
        }
        LX = Edge->Interval->End - Edge->Interval->Start;
        Chars = CharsSoFar + (Edge->Interval->End - Edge->Interval->Start);
        if (Chars < 60)
        {
            LX += SatisfyIntervalBack(SubGNode, Chars);
        }
        MinLX = min(MinLX, LX);
    }
    if (MinLX < GNode->LX)
    {
        printf("%d < %d???\n", MinLX, GNode->LX);
    }

    GNode->LX = MinLX;
    return MinLX;

}

// Free all the GeneNode instances in the global list.
void FreeGeneNodes()
{
    GeneNode* Prev;
    GeneNode* Node;
    // Free all the gene nodes:
    Prev = NULL;
    for (Node = g_GeneFirst; Node; Node = Node->Next)
    {
        Node->Interval->GNode = NULL;
        SafeFree(Prev);
        Prev = Node;
    }
    SafeFree(Prev);
    g_GeneFirst = NULL;
    g_GeneLast = NULL;
    GeneNodeCount = 0;
}

// Take this interval, and 'satisfy' it by adding linked intervals until (a) there are 
// no more links to follow, or (b) we have extended a considerable distance (in amino acids!).
// Take the resulting pool of intervals, build exons for them, and write out one 
// "gene" record.  
int SatisfyIntervalAndWriteGene(IntervalNode* NextUnsatisfied, FILE* SequenceFile, FILE* OutputFile, 
    int RecordNumber, int ChromosomeNumber, int ReverseFlag, int MinORFLength)
{
    GeneNode* Node;
    GeneNode* Dissatisfied;
    char GeneName[256];
    int AllSatisfied;
    char DirectionChar;
    int ValidGeneFlag;
    //

    // First gene node wraps NextUnsatisfied:
    g_GeneFirst = (GeneNode*)calloc(sizeof(GeneNode), 1);
    g_GeneFirst->Interval = NextUnsatisfied;
    NextUnsatisfied->GNode = g_GeneFirst;
    g_GeneLast = g_GeneFirst;
    GeneNodeCount = 1;

    // Add the necessary gene nodes to satisfy:
    Node = g_GeneFirst;
    SatisfyIntervalForward(Node, 0);
    SatisfyIntervalBack(Node, 0);
    Node->Interval->Satisfied = 1;
    // Iterate: If there are unsatisfied intervals in the gene, and the gene isn't too large, satisfy some more.
    AllSatisfied = 0;
    while (GeneNodeCount < MAX_INTERVALS_PER_GENE)
    {
        Dissatisfied = NULL;
        // Find the first not-yet-satisfied interval:
        for (Node = g_GeneFirst; Node; Node = Node->Next)
        {
            if (!Node->Interval->Satisfied)
            {
                Dissatisfied = Node;
                break;
            }
        }
        if (!Dissatisfied)
        {
            AllSatisfied = 1;
            break; // Done!
        }
        SatisfyIntervalForward(Dissatisfied, 0);
        SatisfyIntervalBack(Dissatisfied, 0);
        Dissatisfied->Interval->Satisfied = 1;
    }
    if (!AllSatisfied)
    {
        g_StatsIncompleteGeneCount++;
        IntegrityCheckGene();
    }

    // Write this gene out:
    if (ReverseFlag)
    {
        DirectionChar = '-';
    }
    else
    {
        DirectionChar = '+';
    }
    sprintf(GeneName, "%d%c Gene %d, %d-%d", ChromosomeNumber, DirectionChar, RecordNumber, g_GeneFirst->Interval->Start, g_GeneLast->Interval->End);
    // *************************
    ValidGeneFlag = BuildAndWriteExons(SequenceFile, OutputFile, ReverseFlag, GeneName, ChromosomeNumber, MinORFLength);

    FreeGeneNodes();
    return ValidGeneFlag;
}


// Once the master interval list has been prepared, we can write out genes.
// The procedure works like this:
// We'll build a linked list of GeneNode structs, from g_GeneFirst to g_GeneLast, with size GeneNodeCount.  The intervals
//   contained in this list of GeneNodes are what we'll write out as a gene record.
// Iterate:
// - Take the first interval not yet satisfied, A.  
// - Find all neighbors necessary in order to satisfy A, and add them to the gene.
// - Iterate:
// -- If the gene contains too many intervals, stop.
// -- If every interval in the gene has now been satisfied, stop.
// -- Otherwise, take the first interval in the gene which has not yet been flagged satisfied, and add the necessary 
//    intervals to satisfy it.  (It's possible that we already have the necessary intervals, and just need to 
//    discoverify that fact) 
// - For each interval in the active range: Construct exons
// - Using the exon structs, write out a gene record
// - Free the exon structs (they're bloaty, containing sequence strings) and the GeneNode list
//
// If the interval graph is well-behaved, and consists of small connected components, then we write one gene for each
// connected component.  If the interval graph is messy, then our iterative procedure will cover the graph in
// a reasonably efficient way.  (We're guaranteed to satisfy one interval with each gene, and we're likely to satisfy
// several)
int WriteGenesForIntervals(char* SequenceFileName, char* OutputFileName, int ChromosomeNumber, int ReverseFlag, int MinORFLength)
{
    FILE* SequenceFile;
    FILE* OutputFile;
    IntervalNode* NextUnsatisfied;
    int RecordNumber = 0;
    int ValidGeneFlag;
    //
    SequenceFile = fopen(SequenceFileName, "rb");
    if (!SequenceFile)
    {
        printf("** ERROR: Unable to open chromosome database '%s'\n", SequenceFileName);
        return 0;
    }
    OutputFile = fopen(OutputFileName, "wb");
    if (!OutputFile)
    {
        printf("** ERROR: Unable to open output file '%s'\n", OutputFileName);
        return 0;
    }
    // Iterate over intervals, skipping over intervals that have already been satisfied.
    NextUnsatisfied = g_FirstInterval;
    while (1)
    {
        if (!NextUnsatisfied)
        {
            break;
        }
        if (NextUnsatisfied->Satisfied)
        {
            NextUnsatisfied = NextUnsatisfied->Next;
            continue;
        }
        
        //printf("\n  - - - Satisfy the next interval: %d-%d\n", NextUnsatisfied->Start, NextUnsatisfied->End);
        // if SatisfyIntervalAndWriteGene returns 0, then there's no real gene here (short ORFs were pruned).
        ValidGeneFlag = SatisfyIntervalAndWriteGene(NextUnsatisfied, SequenceFile, OutputFile, RecordNumber, ChromosomeNumber, ReverseFlag, MinORFLength);
        if (ValidGeneFlag)
        {
            RecordNumber++;
            //printf("Wrote gene record %d\n", RecordNumber);
        }
    } // main loop for writing out genes.
    printf("Genes have been written out.  Statistics:\n");
    printf("%d\t", ChromosomeNumber);
    printf("%d\t", ReverseFlag);
    printf("%d\t", RecordNumber);
    printf("%d\t", g_StatsIncompleteGeneCount);
    printf("%d\t", g_StatsLargestGeneSize);
    printf("%d\t", g_StatsLargestGeneRecordNumber);
    printf("%d\t", g_StatsIntervalsBeforeMerge);
    printf("%d\t", g_StatsEdgesBeforeMerge);
    printf("%d\t", g_StatsIntervalsAfterMerge);
    printf("%d\t", g_StatsEdgesAfterMerge);
    printf("%d\t", g_StatsIntervalsAfterIntersect);
    printf("%d\t", g_StatsEdgesAfterIntersect);
    printf("%d\t", g_StatsTotalExonsWritten);
    printf("%d\t", g_StatsTotalEdgesWritten);
    printf("\n");
    return RecordNumber;
}

// Free exons for an interval
void FreeIntervalExons(IntervalNode* Interval)
{
    ExonNode* PrevExon;
    ExonNode* Exon;
    ExonLink* PrevLink;
    ExonLink* Link;

    PrevExon = NULL;
    for (Exon = Interval->FirstExon; Exon; Exon = Exon->Next)
    {
        // Free forward links:
        PrevLink = NULL;
        for (Link = Exon->FirstForward; Link; Link = Link->Next)
        {
            SafeFree(PrevLink);
            PrevLink = Link;
        }
        SafeFree(PrevLink);

        // Free backward links:
        PrevLink = NULL;
        for (Link = Exon->FirstBack; Link; Link = Link->Next)
        {
            SafeFree(PrevLink);
            PrevLink = Link;
        }
        SafeFree(PrevLink);

        // Free exon itself:
        SafeFree(PrevExon->Sequence);
        SafeFree(PrevExon);
        PrevExon = Exon;
    }
    SafeFree(PrevExon->Sequence);
    SafeFree(PrevExon);
    Interval->FirstExon = NULL;
    Interval->LastExon = NULL;
}

void AddExonToInterval(IntervalNode* Interval, ExonNode* Exon)
{
    if (Interval->LastExon)
    {
        Interval->LastExon->Next = Exon;
    }
    else
    {
        Interval->FirstExon = Exon;
    }
    Interval->LastExon = Exon;
    Exon->Interval = Interval;
}



// Given an exon and its dna sequence, translate into amino acids.  
// The exon's prefix has already been set, but we'll set the suffix (with the 'leftovers')
// If MinORFLength>0, then call MaskBrokenSequence to MASK OUT the interval between two stop 
// codons separated by less than MinORFLength.
void GetExonSequence(ExonNode* Exon, char* DNA, int MinORFLength)
{
    char ProteinBuffer[MAX_INTERVAL_LENGTH];
    int Pos;
    int Length;
    char* Peptide;
    int SuffixPos;

    if (!DNA || !*DNA)
    {
        Exon->Suffix[0] = '\0';
        Exon->Sequence = NULL;
        Exon->Length = 0;
        return;
    }
    Length = strlen(DNA);
    Pos = 0;
    Peptide = ProteinBuffer;
    while (Pos + 2 < Length)
    {
        *Peptide = TranslateCodon(DNA + Pos);
        Peptide++;
        Pos += 3;
    }
    *Peptide = '\0';
    MaskBrokenSequence(ProteinBuffer, MinORFLength);
    //Exon->Length = strlen(ProteinBuffer);
    Exon->Length = strlen(ProteinBuffer);
    if (Exon->Length)
    {
        Exon->Sequence = (char*)calloc(sizeof(char), Exon->Length + 1);
        strcpy(Exon->Sequence, ProteinBuffer);
    }
    SuffixPos = 0;
    while (Pos < Length)
    {
        Exon->Suffix[SuffixPos] = DNA[Pos];
        SuffixPos++;
        Pos++;
    }
    Exon->Suffix[SuffixPos] = '\0';
}

// If an exon's protein sequence has two stop codons, with fewer than 50 residues in between,
// then CUT OUT that section.  (Because we interpret genomic intervals in multiple reading 
// frames, we often get reads of very short length; we don't believe that such short peptides
// are reasonable!)
// Todo: Try encoding some stop codons as selenocysteines (U).
// Update: We no longer CUT the sequence, because that would ruin our genomic coordinates.
// Rather, we MASK super-short reading frames!
void MaskBrokenSequence(char* Protein, int MinORFLength)
{
    int AnchorPos = -1;
    int Pos;
    char AA;
    int MaskPos;
    //
    // if MinORFLength <= 0, then don't filter.
    if (MinORFLength <= 0)
    {
        return;
    }
    Pos = 0;
    while (1)
    {
        AA = Protein[Pos];
        if (!AA)
        {
            break;
        }
        if (AA == 'X')
        {
            if (AnchorPos == -1 || (Pos - AnchorPos >= MinORFLength))
            {
                AnchorPos = Pos;
            }
            else
            {
                for (MaskPos = AnchorPos + 1; MaskPos < Pos; MaskPos++)
                {
                    Protein[MaskPos] = 'X';
                }
                AnchorPos = Pos;
                Pos++;
                continue;
            }
        }
        Pos++;
    }
}

// Add a forward link from exon A to exon B
void AddExonLink(ExonNode* A, ExonNode* B, char AA, int Power)
{
    ExonLink* Link;
    //
    Link = (ExonLink*)calloc(sizeof(ExonLink), 1);
    Link->Exon = B;
    Link->AA = AA;
    Link->Power = Power;
    if (A->LastForward)
    {
        A->LastForward->Next = Link;
    }
    else
    {
        A->FirstForward = Link;
    }
    A->LastForward = Link;
    //
    Link = (ExonLink*)calloc(sizeof(ExonLink), 1);
    Link->Exon = A;
    Link->AA = AA;
    Link->Power = Power;
    if (B->LastBack)
    {
        B->LastBack->Next = Link;
    }
    else
    {
        B->FirstBack = Link;
    }
    B->LastBack = Link;
}

// Link up the "edge" DAG nodes for intervals, if their parent intervals are linked.
void LinkDAGAcrossIntervals(IntervalNode* Interval, EdgeNode* Edge, int ReverseFlag)
{
    IntervalNode* OtherInterval;
    GenomeDAGNode* DAGNode;
    GenomeDAGNode* OtherDAGNode;
    int DAGNodeIndex;
    int OtherDAGNodeIndex;
    //
    OtherInterval = Edge->Interval;
    for (DAGNodeIndex = 0; DAGNodeIndex < Interval->DAGNodeCount; DAGNodeIndex++)
    {
        DAGNode = Interval->DAGNodes + DAGNodeIndex;
        if (!DAGNode->Sequence)
        {
            continue;
        }
        // We link forward only from dag nodes that touch the edge:
        if (DAGNode->End < Interval->End)
        {
            continue;
        }
        for (OtherDAGNodeIndex = 0; OtherDAGNodeIndex < OtherInterval->DAGNodeCount; OtherDAGNodeIndex++)
        {
            OtherDAGNode = OtherInterval->DAGNodes + OtherDAGNodeIndex;
            if (!OtherDAGNode->Sequence)
            {
                continue;
            }
            if (OtherDAGNode->Start > OtherInterval->Start)
            {
                continue;
            }
            GenomeDAGLinkBack(OtherDAGNode, DAGNode, Edge->Count); // link!
        }
    }
}

// Add links between edges, according to how the 'parent' DAG is linked.
// Each DAG can have three exons: prefix-length 0, prefix 1, prefix 2.
// If you're length 1, then you get one with no prefix (and length-1 suffix), one with prefix (and no suffix).
// If you're longer than 1, then the exons you get depend on your reading frame flags: EST-dervied exons get 
// three exons, most gene-finding-derived exons get a single exon for the single plausible reading frame.
void LinkIntervalExons(IntervalNode* Interval, int ReverseFlag)
{
    ExonNode* Exon;
    ExonNode* NextExon = NULL;
    GenomeDAGNode* DAGNode;
    GenomeDAGLink* Edge;
    GenomeDAGLink* NextEdge;
    int DAGNodeIndex;
    int ExonIndex;
    int SuffixLength;
    char DNA[4];
    char AA = 0;
    int Power = 0;
    //
    DNA[3] = '\0';
    for (DAGNodeIndex = 0; DAGNodeIndex < Interval->DAGNodeCount; DAGNodeIndex++)
    {
        DAGNode = Interval->DAGNodes + DAGNodeIndex;
        if (ReverseFlag)
        {
            Edge = DAGNode->FirstBack;
        }
        else
        {
            Edge = DAGNode->FirstForward;
        }
        while (Edge)
        {
            // This DAG has one, two, or three exons to join.
            for (ExonIndex = 0; ExonIndex < 3; ExonIndex++)
            {
                Exon = DAGNode->Exons[ExonIndex];
                if (!Exon)
                {
                    continue;
                }
                SuffixLength = strlen(Exon->Suffix);
                switch (SuffixLength)
                {
                case 0:
                    NextExon = Edge->Node->Exons[0];
                    if (!NextExon)
                    {
                        continue;
                    }
                    AA = '\0';
                    Power = Edge->Count;
                    break;
                case 1:
                    // A length-1 suffix.  We link to a length-2 prefix, if available:
                    NextExon = Edge->Node->Exons[2];
                    if (NextExon)
                    {
                        // Combine our length-1 suffix with the length-2 prefix:
                        DNA[0] = Exon->Suffix[0];
                        DNA[1] = NextExon->Prefix[0];
                        DNA[2] = NextExon->Prefix[1];
                        AA = TranslateCodon(DNA);
                        Power = Edge->Count;
                        break;
                    }
                    else
                    {
                        // There's no length-2 prefix available.  If that's because the next
                        // exon has length <1, then we "leapfrog" through it:
                        if ((Edge->Node->End == Edge->Node->Start + 1) && Edge->Node->Exons[0])
                        {
                            // The ugly case.  Take our suffix char, add the forward interval's base, 
                            // and then add one base from the "far interval".
                            DNA[0] = Exon->Suffix[0];
                            DNA[1] = Edge->Node->Exons[0]->Suffix[0];
                            if (ReverseFlag)
                            {
                                NextEdge = Edge->Node->FirstBack;
                            }
                            else
                            {
                                NextEdge = Edge->Node->FirstForward;
                            }
                            for (; NextEdge; NextEdge = NextEdge->Next)
                            {
                                NextExon = NextEdge->Node->Exons[1];
                                if (NextExon)
                                {
                                    DNA[2] = NextExon->Prefix[0];
                                    AA = TranslateCodon(DNA);
                                    Power = max(Edge->Count, NextEdge->Count);
                                    AddExonLink(Exon, NextExon, AA, Power);
                                }
                            }
                        }
                        continue;
                    }
                case 2:
                    NextExon = Edge->Node->Exons[1];
                    if (!NextExon)
                    {
                        continue;
                    }
                    DNA[0] = Exon->Suffix[0];
                    DNA[1] = Exon->Suffix[1];
                    DNA[2] = NextExon->Prefix[0];
                    AA = TranslateCodon(DNA);
                    Power = Edge->Count;
                    break;
                }
                AddExonLink(Exon, NextExon, AA, Power);

            } // exon loop
            Edge = Edge->Next;
        } // edge loop
    } // DAG node loop
}

// Write out one exon record in binary format. 
void WriteExonRecord(ExonNode* Exon, FILE* OutputFile, int ReverseFlag)
{
    int Length;
    int BackLinkCount;
    int ForwardLinkCount;
    ExonLink* Link;
    //
    WriteBinary(&Exon->Start, sizeof(int), 1, OutputFile);
    WriteBinary(&Exon->End, sizeof(int), 1, OutputFile);
    Length = Exon->Length;
    g_StatsTotalExonsWritten++;

    WriteBinary(&Length, sizeof(int), 1, OutputFile);
    WriteBinary(&Exon->Interval->Occurrences, sizeof(int), 1, OutputFile);
    if (Length)
    {
        WriteBinary(Exon->Sequence, sizeof(char), Length, OutputFile);
    }
    BackLinkCount = 0;
    for (Link = Exon->FirstBack; Link; Link = Link->Next)
    {
        BackLinkCount++;
    }
    ForwardLinkCount = 0;
    for (Link = Exon->FirstForward; Link; Link = Link->Next)
    {
        ForwardLinkCount++;
    }

    if (0) //ReverseFlag)
    {
        WriteBinary(Exon->Prefix, sizeof(char), 2, OutputFile);
        WriteBinary(Exon->Suffix, sizeof(char), 2, OutputFile);

        WriteBinary(&ForwardLinkCount, sizeof(int), 1, OutputFile);
        WriteBinary(&BackLinkCount, sizeof(int), 1, OutputFile);
        for (Link = Exon->FirstForward; Link; Link = Link->Next)
        {
            g_StatsTotalEdgesWritten++;
            WriteBinary(&Link->Exon->Index, sizeof(int), 1, OutputFile);
            WriteBinary(&Link->Power, sizeof(int), 1, OutputFile);
            WriteBinary(&Link->AA, sizeof(char), 1, OutputFile);
        }

    }
    else
    {
        WriteBinary(Exon->Prefix, sizeof(char), 2, OutputFile);
        WriteBinary(Exon->Suffix, sizeof(char), 2, OutputFile);

        WriteBinary(&BackLinkCount, sizeof(int), 1, OutputFile);
        WriteBinary(&ForwardLinkCount, sizeof(int), 1, OutputFile);
        for (Link = Exon->FirstBack; Link; Link = Link->Next)
        {
            g_StatsTotalEdgesWritten++;
            WriteBinary(&Link->Exon->Index, sizeof(int), 1, OutputFile);
            WriteBinary(&Link->Power, sizeof(int), 1, OutputFile);
            WriteBinary(&Link->AA, sizeof(char), 1, OutputFile);
        }
    }
}

// Given a range of intervals (from First to Last), with exons built, write
// out the binary gene record to the splice-tolerant database file.
int WriteGeneRecord(int ChromosomeNumber, char* GeneName, int ReverseFlag, FILE* OutputFile)
{
    int ExonCount = 0;
    IntervalNode* Interval;
    ExonNode* Exon;
    GeneNode* Node;
    char ForwardFlag;
    for (Node = g_GeneFirst; Node; Node = Node->Next)
    {
        Interval = Node->Interval;
        for (Exon = Interval->FirstExon; Exon; Exon = Exon->Next)
        {
            ExonCount++;
        }
    }
    if (!ExonCount)
    {
        // No exons!?  That can happen if we pruned them all away due to short ORFs.
        // No need to write anything at all:
        return 0;
    }

    WriteBinary(GeneName, sizeof(char), 256, OutputFile);
    WriteBinary(GeneName, sizeof(char), 256, OutputFile);
    WriteBinary(&ChromosomeNumber, sizeof(int), 1, OutputFile);
    if (ReverseFlag)
    {
        ForwardFlag = 0;
    }
    else
    {
        ForwardFlag = 1;
    }
    WriteBinary(&ForwardFlag, sizeof(char), 1, OutputFile);
    WriteBinary(&ExonCount, sizeof(int), 1, OutputFile);
    if (ReverseFlag)
    {
        // Re-index all the exons:
        ExonCount = 0;
        for (Node = g_GeneLast; Node; Node = Node->Prev)
        {
            for (Exon = Node->Interval->FirstExon; Exon; Exon = Exon->Next)
            {
                Exon->Index = ExonCount;
                ExonCount++;
            }
        }
        for (Node = g_GeneLast; Node; Node = Node->Prev)
        {
            Interval = Node->Interval;
            for (Exon = Interval->FirstExon; Exon; Exon = Exon->Next)
            {
                WriteExonRecord(Exon, OutputFile, ReverseFlag);
            }
        }
    }
    else
    {
        // Re-index all the exons:
        ExonCount = 0;
        for (Node = g_GeneFirst; Node; Node = Node->Next)
        {
            for (Exon = Node->Interval->FirstExon; Exon; Exon = Exon->Next)
            {
                Exon->Index = ExonCount;
                ExonCount++;
            }
        }

        for (Node = g_GeneFirst; Node; Node = Node->Next)
        {
            Interval = Node->Interval;
            for (Exon = Interval->FirstExon; Exon; Exon = Exon->Next)
            {
                WriteExonRecord(Exon, OutputFile, ReverseFlag);
            }
        }
    }
    return 1;
}

// Verify that our exon construction is valid.  The number of exon forward-links and
// backward-links should match.  And no exons should go outside the active range...
void IntegrityCheckGene()
{
    int ForwardCount = 0;
    int BackwardCount = 0;
    int ExonCount = 0;
    int IntervalCount = 0;
    ExonNode* Exon;
    ExonLink* Link;
    ExonLink* RecipLink;
    int FoundFlag;
    IntervalNode* Interval;
    GeneNode* Node;
    EdgeNode* Edge;
    int Count;
    //
    printf("\n===Integrity check: Intervals from %d to %d\n", g_GeneFirst->Interval->Start, g_GeneLast->Interval->End);

    for (Node = g_GeneFirst; Node; Node = Node->Next)
    {
        Interval = Node->Interval;
        if (Interval->Satisfied)
        {
            printf("%d - %d SATISFIED ", Interval->Start, Interval->End);
        }
        else
        {
            printf("%d - %d unsatisfied ", Interval->Start, Interval->End);
        }
        Count = 0;
        for (Edge = Interval->FirstForward; Edge; Edge = Edge->Next)
        {
            Count++;
        }
        switch (Count)
        {
        case 0:
            break;
        case 1:
            printf("to %d", Interval->FirstForward->Interval->Start);
            break;
        case 2:
            printf("to %d, %d", Interval->FirstForward->Interval->Start, Interval->FirstForward->Next->Interval->Start);
            break;
        default:
            printf("to %d, %d, +%d", Interval->FirstForward->Interval->Start, 
                Interval->FirstForward->Next->Interval->Start, Count-2);
            break;

        }
        printf("\n");
        IntervalCount++;
        for (Exon = Interval->FirstExon; Exon; Exon = Exon->Next)
        {
            ExonCount++;
            for (Link = Exon->FirstForward; Link; Link = Link->Next)
            {
                ForwardCount++;
                if (!Link->Exon->Interval->GNode)
                {
                    printf("** Warning: Exon %d links forward out of this world.\n", Exon->Index);
                }
                // Check for a reciprocal link, too:
                FoundFlag = 0;
                for (RecipLink = Link->Exon->FirstBack; RecipLink; RecipLink = RecipLink->Next)
                {
                    if (RecipLink->Exon == Exon)
                    {
                        FoundFlag = 1;
                        if (Link->Power != RecipLink->Power)
                        {
                            printf("** Warning: Exon link %d to %d has inconsistent strength %d, %d\n", 
                                Exon->Index, Link->Exon->Index, Link->Power, RecipLink->Power);
                        }
                        break;
                    }
                }
                if (!FoundFlag)
                {
                    printf("** Warning: Exon %d has a non-reciprocated forward link.\n", Exon->Index);
                }
            }
            for (Link = Exon->FirstBack; Link; Link = Link->Next)
            {
                BackwardCount++;
                if (!Link->Exon->Interval->GNode)
                {
                    printf("** Warning: Exon %d links backward out of this world.\n", Exon->Index);
                }
                // Check for a reciprocal link, too:
                FoundFlag = 0;
                for (RecipLink = Link->Exon->FirstForward; RecipLink; RecipLink = RecipLink->Next)
                {
                    if (RecipLink->Exon == Exon)
                    {
                        FoundFlag = 1;
                        break;
                    }
                }
                if (!FoundFlag)
                {
                    printf("** Warning: Exon %d has a non-reciprocated backward link.\n", Exon->Index);
                }
            }
        } // exon loop
    } // interval loop
    printf("Saw %d intervals, %d exons, %d links.\n", IntervalCount, ExonCount, ForwardCount);
    if (ForwardCount != BackwardCount)
    {
        printf("** Warning: Total forward links is %d != backward links %d\n", ForwardCount, BackwardCount);
    }
}

int g_UFTotalExons = 0;
int g_UFTotalAA = 0;
int g_UFTotalEdges = 0;
int g_TotalExons = 0;
int g_TotalAA = 0;
int g_TotalEdges = 0;
int g_TotalTrueExons = 0;
int g_TotalTrueEdges = 0;

typedef struct ExonSortNode
{
    ExonNode* Exon;
} ExonSortNode;

int CompareExonNodesForward(const ExonSortNode* NodeA, const ExonSortNode* NodeB)
{
    if (NodeA->Exon->Start < NodeB->Exon->Start)
    {
        return -1;
    }
    if (NodeA->Exon->Start > NodeB->Exon->Start)
    {
        return 1;
    }
    // arbitrary:
    if (NodeA->Exon < NodeB->Exon)
    {
        return -1;
    }
    else
    {
        return 1;
    }
}
int CompareExonNodesBackward(const ExonSortNode* NodeA, const ExonSortNode* NodeB)
{
    if (NodeA->Exon->Start < NodeB->Exon->Start)
    {
        return 1;
    }
    if (NodeA->Exon->Start > NodeB->Exon->Start)
    {
        return -1;
    }
    // arbitrary:
    if (NodeA->Exon < NodeB->Exon)
    {
        return 1;
    }
    else
    {
        return -1;
    }
}

// It is desirable for an exon's back-links to always hit exons with LOWER index numbers.
void SortExons(int ReverseFlag)
{
    ExonSortNode* ExonNodes;
    GeneNode* GNode;
    int ExonCount;
    int ExonIndex;
    ExonNode* Exon;
    //
    for (GNode = g_GeneFirst; GNode; GNode = GNode->Next)
    {
        ExonCount = 0;
        for (Exon = GNode->Interval->FirstExon; Exon; Exon = Exon->Next)
        {
            ExonCount++;
        }
        if (ExonCount)
        {
            ExonNodes = (ExonSortNode*)calloc(ExonCount, sizeof(ExonSortNode));
            ExonIndex = 0;
            for (Exon = GNode->Interval->FirstExon; Exon; Exon = Exon->Next)
            {
                ExonNodes[ExonIndex].Exon = Exon;
                ExonIndex++;
            }
            if (ReverseFlag)
            {
                qsort(ExonNodes, ExonCount, sizeof(ExonSortNode), (QSortCompare)CompareExonNodesBackward);
            }
            else
            {
                qsort(ExonNodes, ExonCount, sizeof(ExonSortNode), (QSortCompare)CompareExonNodesForward);
            }
            GNode->Interval->FirstExon = ExonNodes[0].Exon;
            GNode->Interval->LastExon = ExonNodes[ExonCount - 1].Exon;
            for (ExonIndex = 0; ExonIndex < ExonCount; ExonIndex++)
            {
                if (ExonIndex < ExonCount - 1)
                {
                    ExonNodes[ExonIndex].Exon->Next = ExonNodes[ExonIndex + 1].Exon;
                }
                else
                {
                    ExonNodes[ExonIndex].Exon->Next = NULL;
                }
            }
            SafeFree(ExonNodes);
        }
    }
}

// For reporting purposes, count how many exons and edges and amino acids are in our db:
void CountExons(int PreFilterFlag)
{
    GeneNode* GNode;
    ExonNode* Exon;
    ExonLink* Link;
    int Pos;
    int TrueExonFlag;
    //
    for (GNode = g_GeneFirst; GNode; GNode = GNode->Next)
    {
        for (Exon = GNode->Interval->FirstExon; Exon; Exon = Exon->Next)
        {
            if (PreFilterFlag)
            {
                g_UFTotalExons++;
                for (Pos = 0; Pos < Exon->Length; Pos++)
                {
                    if (Exon->Sequence[Pos]!='X')
                    {
                        g_UFTotalAA++;
                    }
                }
                
            }
            else
            {
                g_TotalExons++;
                TrueExonFlag = 1;
                for (Link = Exon->FirstBack; Link; Link = Link->Next)
                {
                    if (Link->Exon->End == Exon->Start || Link->Exon->Start == Exon->End)
                    {
                        TrueExonFlag = 0;
                    }
                    else
                    {
                        g_TotalTrueEdges++;
                    }
                }
                g_TotalTrueExons += TrueExonFlag;
                for (Pos = 0; Pos < Exon->Length; Pos++)
                {
                    if (Exon->Sequence[Pos]!='X')
                    {
                        g_TotalAA++;
                    }
                }
            }
            for (Link = Exon->FirstBack; Link; Link = Link->Next)
            {
                if (PreFilterFlag)
                {
                    g_UFTotalEdges++;
                    if (Link->AA)
                    {
                        g_UFTotalAA++;
                    }
                }
                else
                {
                    g_TotalEdges++;
                    if (Link->AA)
                    {
                        g_TotalAA++;
                    }
                }
            }
            for (Link = Exon->FirstForward; Link; Link = Link->Next)
            {
                if (PreFilterFlag)
                {
                    g_UFTotalEdges++;
                    if (Link->AA)
                    {
                        g_UFTotalAA++;
                    }
                }
                else
                {
                    g_TotalEdges++;
                    if (Link->AA)
                    {
                        g_TotalAA++;
                    }
                }
            } // iterate forward links
        } // iterate exons
    } // iterate GNodes
}

// Create a link between two genome DAG nodes.
void GenomeDAGLinkBack(GenomeDAGNode* DAGNode, GenomeDAGNode* BackDAGNode, int Count)
{
    GenomeDAGLink* NewLink;
    GenomeDAGLink* Link;
    //
    // Add the back-link:
    NewLink = (GenomeDAGLink*)calloc(1, sizeof(GenomeDAGLink));
    NewLink->Node = BackDAGNode;
    NewLink->Count = Count;
    if (!DAGNode->FirstBack)
    {
        DAGNode->FirstBack = NewLink;
    }
    else
    {
        for (Link = DAGNode->FirstBack; Link->Next; Link = Link->Next)
        {
            ;
        }
        Link->Next = NewLink;
    }
    // Add the forward-link:
    NewLink = (GenomeDAGLink*)calloc(1, sizeof(GenomeDAGLink));
    NewLink->Node = DAGNode;
    NewLink->Count = Count;
    if (!BackDAGNode->FirstForward)
    {
        BackDAGNode->FirstForward = NewLink;
    }
    else
    {
        for (Link = BackDAGNode->FirstForward; Link->Next; Link = Link->Next)
        {
            ;
        }
        Link->Next = NewLink;
    }
}

int GetReadingFrameFlag(int Start, int End, int Offset, int ReverseFlag)
{
    int ReadingFrameFlag;
    //
    if (ReverseFlag)
    {
        switch ((End - 1 - Offset) % 3)
        {
        case 0:
            ReadingFrameFlag = IFLAG_FRAME_0;
            break;
        case 1:
            ReadingFrameFlag = IFLAG_FRAME_1;
            break;
        case 2:
            ReadingFrameFlag = IFLAG_FRAME_2;
            break;
        }
    }
    else
    {
        switch ((Start + Offset) % 3)
        {
        case 0:
            ReadingFrameFlag = IFLAG_FRAME_0;
            break;
        case 1:
            ReadingFrameFlag = IFLAG_FRAME_1;
            break;
        case 2:
            ReadingFrameFlag = IFLAG_FRAME_2;
            break;
        }
    }
    return ReadingFrameFlag;
}

// Build a DAG for this genomic interval.  Then we'll generate three exons (two, for length-1 nodes)
// for each node of the DAG.  The DAG is generally just one node, corresponding to genomic DNA.
// But the DAG may have extra nodes and edges if there are SNPs that fall within the interval.
void BuildDAGForInterval(IntervalNode* Interval, FILE* GenomicFile, int MinORFLength, int ReverseFlag)
{
    int DAGNodeCount;
    int DAGNodeIndex;
    int NextDAGStart;
    GenomeDAGNode* DAGNode;
    int Length;
    ExonNode* Exon;
    char RCDNABuffer[MAX_INTERVAL_LENGTH + 1];
    int PolyIndex;
    int FirstPolyIndex;
    Polymorphism* Poly;
    int PrevNodesStart = -1;
    int PrevNodesEnd = -1;
    int PrevNodeIndex;
    int NewPrevNodesStart;
    int SNPIndex;
    int ReadingFrameFlag;
    //

    FirstPolyIndex = FindPolyInInterval(Interval->Start, Interval->End);

    ////////////////////////////////////////////////////////////
    // How many nodes in our DAG?  Assume that no polymorphisms overlap, and so we will need 1 dag node
    // plus three per polymorphism (or two for a polymorphism directly after another polymorphism, but
    // just alloc three anyway)
    DAGNodeCount = 1;
    if (FirstPolyIndex >= 0)
    {
        for (PolyIndex = FirstPolyIndex; PolyIndex < g_PolymorphismCount; PolyIndex++)
        {
            if (g_Polymorphisms[PolyIndex].Pos >= Interval->End)
            {
                break;
            }
            // A polymorphism means one node for each SNP allele, and (USUALLY) a node
            // for the continued 'core track'
            //       a
            //  XXXX b xxxx
            //       c 
            DAGNodeCount++;
            DAGNodeCount += strlen(g_Polymorphisms[PolyIndex].SNP);
        }
    }
    Interval->DAGNodeCount = DAGNodeCount;
    Interval->DAGNodes = (GenomeDAGNode*)calloc(DAGNodeCount, sizeof(GenomeDAGNode));
    ////////////////////////////////////////////////////////////
    // Initialize all the DAG nodes:
    NextDAGStart = Interval->Start;
    DAGNodeIndex = 0;
    PolyIndex = FirstPolyIndex;
    while (1)
    {
        if (PolyIndex < 0 || PolyIndex >= g_PolymorphismCount || g_Polymorphisms[PolyIndex].Pos >= Interval->End)
        {
            // There are no more polymorphisms.  
            if (NextDAGStart < Interval->End)
            {
                // Generate an interval that extends to the end:
                DAGNode = Interval->DAGNodes + DAGNodeIndex;
                DAGNode->Start = NextDAGStart;
                DAGNode->End = Interval->End; 
                Length = DAGNode->End - DAGNode->Start;
                DAGNode->Sequence = (char*)calloc(Length + 1, sizeof(char)); // +1 for null terminator
                fseek(GenomicFile, DAGNode->Start, 0);
                ReadBinary(DAGNode->Sequence, sizeof(char), Length, GenomicFile);
                // If we have some nodes already, that's because there's a polymorphism.  Link to the 
                // previous two nodes:
                if (DAGNodeIndex)
                {
                    for (PrevNodeIndex = PrevNodesStart; PrevNodeIndex <= PrevNodesEnd; PrevNodeIndex++)
                    {
                        GenomeDAGLinkBack(DAGNode, Interval->DAGNodes + PrevNodeIndex, 0);
                    }
                }
            }
            // And we're done!
            break;
        }
        // There is another polymorphism.  
        Poly = g_Polymorphisms + PolyIndex;
        if (NextDAGStart < Poly->Pos)
        {
            // If there's non-polymorphic sequence before the next poly, then
            // generate a DAG node for it.
            DAGNode = Interval->DAGNodes + DAGNodeIndex;
            DAGNode->Start = NextDAGStart;
            DAGNode->End = Poly->Pos; 
            Length = DAGNode->End - DAGNode->Start;
            DAGNode->Sequence = (char*)calloc(Length + 1, sizeof(char)); // +1 for null terminator
            fseek(GenomicFile, DAGNode->Start, 0);
            ReadBinary(DAGNode->Sequence, sizeof(char), Length, GenomicFile);
            // If we have some nodes already, that's because there's a polymorphism.  Link to the 
            // previous nodes:
            if (DAGNodeIndex)
            {
                for (PrevNodeIndex = PrevNodesStart; PrevNodeIndex <= PrevNodesEnd; PrevNodeIndex++)
                {
                    GenomeDAGLinkBack(DAGNode, Interval->DAGNodes + PrevNodeIndex, 0);
                }
            }
            PrevNodesStart = DAGNodeIndex;
            PrevNodesEnd = DAGNodeIndex;
            DAGNodeIndex++;
        } 
        // Nodes for the two (or more) alleles:
        NewPrevNodesStart = DAGNodeIndex;
        for (SNPIndex = 0; SNPIndex < 4; SNPIndex++)
        {
            if (!Poly->SNP[SNPIndex])
            {
                break;
            }
            DAGNode = Interval->DAGNodes + DAGNodeIndex;
            DAGNode->Start = Poly->Pos;
            DAGNode->End = Poly->Pos + 1;
            DAGNode->Sequence = (char*)calloc(2, sizeof(char));
            DAGNode->Sequence[0] = Poly->SNP[SNPIndex];
            if (PrevNodesStart > -1)
            {
                for (PrevNodeIndex = PrevNodesStart; PrevNodeIndex <= PrevNodesEnd; PrevNodeIndex++)
                {
                    GenomeDAGLinkBack(DAGNode, Interval->DAGNodes + PrevNodeIndex, 0);
                }
            }
            DAGNodeIndex++;
        }
        PrevNodesStart = NewPrevNodesStart;
        PrevNodesEnd = DAGNodeIndex - 1;        
        NextDAGStart = Poly->Pos + 1;
        PolyIndex++;
    }
    ////////////////////////////////////////////////////////////////////////////////
    // The DAG for the interval has now been constructed.  Build exons for all DAG nodes.
    for (DAGNodeIndex = 0; DAGNodeIndex < Interval->DAGNodeCount; DAGNodeIndex++)
    {
        DAGNode = Interval->DAGNodes + DAGNodeIndex;
        if (!DAGNode->Sequence)
        {
            continue; // not a real DAG node.
        }
        if (DAGNode->End <= DAGNode->Start)
        {
            DAGNode = DAGNode;
        }
        // Reverse-complement the DAG's sequence, if necessary:
        Length = DAGNode->End - DAGNode->Start;
        if (Length >= MAX_INTERVAL_LENGTH)
        {
            printf("** Warning: Genomic interval from %d to %d is MUCH too long to process; truncating!", DAGNode->Start, DAGNode->End);
            DAGNode->Sequence[MAX_INTERVAL_LENGTH] = '\0';
        }
        if (ReverseFlag)
        {
            strcpy(RCDNABuffer, DAGNode->Sequence);
            WriteReverseComplement(RCDNABuffer, DAGNode->Sequence);
        }
        DAGNode->Exons = (ExonNode**)calloc(3, sizeof(ExonNode*));
        // Add two or three exons for this DAG node.
        // Check the reading frame of the interval to decide where codons are supposed to begin.
        ReadingFrameFlag = GetReadingFrameFlag(DAGNode->Start, DAGNode->End, 0, ReverseFlag);
        if (ReadingFrameFlag & Interval->Flags)
        {
            // Reading Frame 0:
            Exon = (ExonNode*)calloc(1, sizeof(ExonNode));
            Exon->Prefix[0] = '\0';
            GetExonSequence(Exon, DAGNode->Sequence, MinORFLength);
            Exon->Start = DAGNode->Start;
            Exon->End = DAGNode->End;
            Exon->DAGNode = DAGNode;
            DAGNode->Exons[0] = Exon;
            AddExonToInterval(Interval, Exon);
        }
        ReadingFrameFlag = GetReadingFrameFlag(DAGNode->Start, DAGNode->End, 1, ReverseFlag);
        if (ReadingFrameFlag & Interval->Flags)
        {
            // Reading frame 1:
            Exon = (ExonNode*)calloc(sizeof(ExonNode), 1);
            Exon->Prefix[0] = DAGNode->Sequence[0];
            Exon->Prefix[1] = '\0';
            GetExonSequence(Exon, DAGNode->Sequence + 1, MinORFLength);
            Exon->Start = DAGNode->Start;
            Exon->End = DAGNode->End;
            Exon->DAGNode = DAGNode;
            DAGNode->Exons[1] = Exon;
            AddExonToInterval(Interval, Exon);
        }
        // Reading frame 2:
        if (Length > 1)
        {
            ReadingFrameFlag = GetReadingFrameFlag(DAGNode->Start, DAGNode->End, 2, ReverseFlag);
            if (ReadingFrameFlag & Interval->Flags)
            {
                Exon = (ExonNode*)calloc(sizeof(ExonNode), 1);
                Exon->Prefix[0] = DAGNode->Sequence[0];
                Exon->Prefix[1] = DAGNode->Sequence[1];
                Exon->Prefix[2] = '\0';
                GetExonSequence(Exon, DAGNode->Sequence + 2, MinORFLength);
                Exon->Start = DAGNode->Start;
                Exon->End = DAGNode->End;
                Exon->DAGNode = DAGNode;
                DAGNode->Exons[2] = Exon;
                AddExonToInterval(Interval, Exon);
            }
        }
    }
}

// Every interval gives rise to three exons (two, if it's only one base long).
// If interval A links to interval B, then A's exons each link to a compatible
// exon in B.  Exception: If an exon with suffix length 1 links to an interval
// of length 1, then we must go to the NEXT-next interval to complete a codon.
//
// GenomicFile is the file containing the genomic sequence.
int BuildAndWriteExons(FILE* GenomicFile, FILE* OutputFile, int ReverseFlag, 
    char* GeneName, int ChromosomeNumber, int MinORFLength)
{
    IntervalNode* Interval;
    EdgeNode* Edge;
    GeneNode* Node;
    int IntervalCount = 0;
    int ValidGeneFlag;
    int VerboseFlag = 0;
    // 
    // Construct 1-3 exons for each interval within the gene:
    for (Node = g_GeneFirst; Node; Node = Node->Next)
    {
        BuildDAGForInterval(Node->Interval, GenomicFile, MinORFLength, ReverseFlag);
        IntervalCount++;
    }
    if (IntervalCount > g_StatsLargestGeneSize)
    {
        g_StatsLargestGeneSize = IntervalCount;
    }
    //DebugPrintBuiltGene();
    // Link up the DAG graphs for all the intervals:
    for (Node = g_GeneFirst; Node; Node = Node->Next)
    {
        Interval = Node->Interval;

        for (Edge = Interval->FirstForward; Edge; Edge = Edge->Next)
        {
            // Ignore edges that extend out of this gene (we'll get them in overlap)
            if (!Edge->Interval->GNode)
            {
                continue;
            }
            LinkDAGAcrossIntervals(Interval, Edge, ReverseFlag);
        }
    }
    //printf("\nLinked DAG across intervals:\n");
    //DebugPrintBuiltGene();
    // Link up the exons, in accordance with the DAG graph linkage:
    for (Node = g_GeneFirst; Node; Node = Node->Next)
    {
        Interval = Node->Interval;
        LinkIntervalExons(Interval, ReverseFlag);
    }

    CountExons(1);

    // If an exon isn't part of any long reading frame, it can be dropped.  And if it contains a stop
    // codon, or the stop codon's prefix (or suffix), and/or some edges, can be dropped.  Perform
    // that filtering now:
    //printf("\nLinked interval exons:\n");
    if (VerboseFlag)
    {
        DebugPrintBuiltGene();
    }
    PruneShortORFs(ReverseFlag, MinORFLength);
    // Exons may include 'masked-out' sequence blocks between stop codons.  These 
    // sequences aren't needed for search, so delete them, splitting the exons if necessary.
    //printf("\nPruned short ORFs:\n");
    //DebugPrintBuiltGene();
    
    PurgeNonCodingExonChunks(ReverseFlag);
    //printf("\nPruned non-coding chunks:\n");
    //DebugPrintBuiltGene();
    CountExons(0);
    SortExons(ReverseFlag);
    // Write out a gene record:
    if (VerboseFlag)
    {
        DebugPrintBuiltGene();
    }
    ValidGeneFlag = WriteGeneRecord(ChromosomeNumber, GeneName, ReverseFlag, OutputFile);

    // Go back and free all the exon records:
    for (Node = g_GeneFirst; Node; Node = Node->Next)
    {
        Interval = Node->Interval;
        FreeIntervalExons(Interval);
        FreeIntervalDAG(Interval);
    }
    return ValidGeneFlag;
}

// Delete an exon entirely!
void DeleteExon(ExonNode* Exon)
{
    ExonNode* OtherExon;
    ExonLink* Link;
    ExonLink* NextLink;
    
    // First, fix the pointers from the parent interval:
    if (Exon->Interval->FirstExon == Exon)
    {
        Exon->Interval->FirstExon = Exon->Next;
        if (Exon->Interval->LastExon == Exon)
        {
            Exon->Interval->LastExon = NULL;
        }
    }
    else
    {
        for (OtherExon = Exon->Interval->FirstExon; OtherExon; OtherExon = OtherExon->Next)
        {
            if (OtherExon->Next == Exon)
            {
                OtherExon->Next = Exon->Next;
                if (Exon->Interval->LastExon == Exon)
                {
                    Exon->Interval->LastExon = OtherExon;
                }
                break;
            }
        }
    }
    // Now, free all the edges (and reciprocal edges):
    Link = Exon->FirstBack;
    while (Link)
    {
        NextLink = Link->Next;
        DeleteExonLink(Exon, Link, 0);
        Link = NextLink;
    }
    Link = Exon->FirstForward;
    while (Link)
    {
        NextLink = Link->Next;
        DeleteExonLink(Exon, Link, 1);
        Link = NextLink;
    }

    // Now, free the exon itself:
    SafeFree(Exon->Sequence);
    SafeFree(Exon);
}

// Delete the specified Link from this Exon.  ForwardFlag indicates
// whether it's a forward link.
void DeleteExonLink(ExonNode* Exon, ExonLink* Link, int ForwardFlag)
{
    ExonLink* OtherLink;
    ExonNode* OtherExon;
    ExonLink* Prev;
    //
    if (ForwardFlag)
    {
        // Update the exon's linked list of edges, removing Link:
        for (OtherLink = Exon->FirstForward; OtherLink; OtherLink = OtherLink->Next)
        {
            if (OtherLink->Next == Link)
            {
                OtherLink->Next = Link->Next;
                if (Exon->LastForward == Link)
                {
                    Exon->LastForward = OtherLink;
                }
                break;
            }
        }
        if (Exon->FirstForward == Link)
        {
            Exon->FirstForward = Link->Next;
        }
        if (Exon->LastForward == Link)
        {
            Exon->LastForward = NULL;
        }

        // Remove the link from the other exon:
        OtherExon = Link->Exon;
        Prev = NULL;
        for (OtherLink = OtherExon->FirstBack; OtherLink; OtherLink = OtherLink->Next)
        {
            if (OtherLink->Exon == Exon && OtherLink->AA == Link->AA)
            {
                if (OtherExon->LastBack == OtherLink)
                {
                    OtherExon->LastBack = Prev;
                }
                if (Prev)
                {
                    Prev->Next = OtherLink->Next;
                }
                else
                {
                    OtherExon->FirstBack = OtherLink->Next;
                }
                SafeFree(OtherLink);
                break;
            }
            Prev = OtherLink;
        }
        SafeFree(Link);
    } // forward link
    else
    {
        // Update the exon's linked list of edges, removing Link:
        for (OtherLink = Exon->FirstBack; OtherLink; OtherLink = OtherLink->Next)
        {
            if (OtherLink->Next == Link)
            {
                OtherLink->Next = Link->Next;
                if (Exon->LastBack == Link)
                {
                    Exon->LastBack = OtherLink;
                }
                break;
            }
        }
        if (Exon->FirstBack == Link)
        {
            Exon->FirstBack = Link->Next;
        }
        if (Exon->LastBack == Link)
        {
            Exon->LastBack = NULL;
        }

        // Remove the link from the other exon:
        OtherExon = Link->Exon;
        Prev = NULL;
        for (OtherLink = OtherExon->FirstForward; OtherLink; OtherLink = OtherLink->Next)
        {
            if (OtherLink->Exon == Exon && OtherLink->AA == Link->AA)
            {
                if (OtherExon->LastForward == OtherLink)
                {
                    OtherExon->LastForward = Prev;
                }
                if (Prev)
                {
                    Prev->Next = OtherLink->Next;
                }
                else
                {
                    OtherExon->FirstForward = OtherLink->Next;
                }
                SafeFree(OtherLink);
                break;
            }
            Prev = OtherLink;
        }
        SafeFree(Link);
    } // backward link
}

// if Link is set, we've counted the exon itself and we're on this link.  
// if Link is null, we're entering the exon:
int GeneFindLongestExtension(ExonNode* OldExon, ExonLink* Link, int LongEnough, int ForwardFlag)
{
    int Length;
    ExonNode* Exon;
    ExonLink* OtherLink;
    int Extension;
    int BestExtension;
    int Pos;
    //
    if (Link && Link->AA)
    {
        Length = 1;
        if (Length >= LongEnough)
        {
            return Length;
        }
    }
    else
    {
        Length = 0;
    }
    Exon = Link->Exon;

    // Iterate over bases in the exon, and add to our length:
    if (ForwardFlag)
    {
        for (Pos = 0; Pos < Exon->Length; Pos++)
        {
            if (Exon->Sequence[Pos] == 'X')
            {
                return Length;
            }
            Length++;
        }
    }
    else
    {
        for (Pos = Exon->Length - 1; Pos >= 0; Pos--)
        {
            if (Exon->Sequence[Pos] == 'X')
            {
                return Length;
            }
            Length++;
        }
    }
    if (Length >= LongEnough)
    {
        return Length;
    }

    // Continue following edges:
    if (ForwardFlag)
    {
        OtherLink = Exon->FirstForward;
    }
    else
    {
        OtherLink = Exon->FirstBack;
    }
    BestExtension = 0;
    while (OtherLink)
    {
        Extension = GeneFindLongestExtension(Exon, OtherLink, LongEnough - Length, ForwardFlag);
        if (Length + Extension >= LongEnough)
        {
            return (Length + Extension);
        }
        BestExtension = max(BestExtension, Extension);
        OtherLink = OtherLink->Next;
    }
    return (Length + BestExtension);
}

// For each exon:
// if its length is zero, we're done
// if it's all stop codons, delete the exon
// if it's all 'real' residues, we're done
// if it starts with one or more stop codons, delete them (and back-edges)
// otherwise, it starts with one or more real residues.  Split them to a separate exon.
void PurgeNonCodingExonChunks(int ReverseFlag)
{
    GeneNode* GNode;
    ExonNode* Exon;
    ExonNode* NextExon = NULL;
    ExonNode* NewExon;
    int FirstReal;
    int FirstStop;
    char AA = 0;
    ExonLink* Link;
    ExonLink* OtherLink;
    ExonLink* NextLink;
    char* NewSequence;
    int AminoIndex;
    int AAEdgeBack;
    int AAEdgeForward;
    //
    for (GNode = g_GeneFirst; GNode; GNode = GNode->Next)
    {
        Exon = GNode->Interval->FirstExon;
        while (Exon)
        {
            // An EMPTY exon can be nuked:
            if (Exon->Start == Exon->End)
            {
                NextExon = Exon->Next;
                DeleteExon(Exon);
                Exon = NextExon;
                continue;
            }
            // This exon has no amino acids, but it has edges.  Keep it:
            if (Exon->Length == 0)
            {
                Exon = Exon->Next;
                continue;
            }

            // Loop over residues to find the first stop codon (-1 if none) and the first real codon (-1 if none)
            FirstReal = -1;
            FirstStop = -1;
            for (AminoIndex = 0; AminoIndex < Exon->Length; AminoIndex++)
            {
                AA = Exon->Sequence[AminoIndex];
                if (AA == 'X')
                {
                    if (FirstStop < 0)
                    {
                        FirstStop = AminoIndex;
                    }
                }
                else
                {
                    if (FirstReal < 0)
                    {
                        FirstReal = AminoIndex;
                    }
                }
            }
            // Count the number of AAEdges (edges with an amino acid char attached) back and forward
            AAEdgeBack = 0;
            AAEdgeForward = 0;
            for (Link = Exon->FirstBack; Link; Link = Link->Next)
            {
                if (Link->AA)
                {
                    AAEdgeBack++;
                }
            }
            for (Link = Exon->FirstForward; Link; Link = Link->Next)
            {
                if (Link->AA)
                {
                    AAEdgeForward++;
                }
            }
            //printf("FirstReal %d FirstStop %d EdgeBack %d EdgeForward %d\n", FirstReal, FirstStop, AAEdgeBack, AAEdgeForward);
            if (FirstReal == -1)
            {
                // This exon contains nothing but stop codons!  
                // We probably mustn't simply delete this exon, or else the aa-edges leading into and out of 
                // it will be broken.  But let's truncate its sequence, and delete any
                // non-AA edges.
                if (AAEdgeBack)
                {
                    NewExon = (ExonNode*)calloc(1, sizeof(ExonNode));
                    NewExon->Interval = Exon->Interval;
                    NewExon->Length = 0;
                    if (ReverseFlag)
                    {
                        NewExon->Start = Exon->End - strlen(Exon->Prefix);
                        NewExon->End = Exon->End;                        
                    }
                    else
                    {
                        NewExon->Start = Exon->Start;
                        NewExon->End = Exon->Start + strlen(Exon->Prefix);
                    }
                    strcpy(NewExon->Prefix, Exon->Prefix);
                    // Assimilate all edges with AA into the new exon:
                    Link = Exon->FirstBack;
                    while (Link)
                    {
                        if (!Link->AA)
                        {
                            NextLink = Link->Next;
                            DeleteExonLink(Exon, Link, 0);
                            Link = NextLink;
                            continue;
                        }
                        // Fix the reciprocal links to point to the new exon:
                        for (OtherLink = Link->Exon->FirstForward; OtherLink; OtherLink = OtherLink->Next)
                        {
                            if (OtherLink->Exon == Exon && OtherLink->AA == Link->AA)
                            {
                                OtherLink->Exon = NewExon;
                            }
                        }
                        if (NewExon->FirstBack)
                        {
                            NewExon->LastBack->Next = Link;
                        }
                        else
                        {
                            NewExon->FirstBack = Link;
                        }
                        NewExon->LastBack = Link;
                        Link = Link->Next;
                    }
                    if (NewExon->LastBack)
                    {
                        NewExon->LastBack->Next = NULL;
                    }
                    Exon->FirstBack = NULL;
                    Exon->LastBack = NULL;
                    GNode->Interval->LastExon->Next = NewExon;
                    GNode->Interval->LastExon = NewExon;
                }
                if (AAEdgeForward)
                {
                    NewExon = (ExonNode*)calloc(1, sizeof(ExonNode));
                    NewExon->Interval = Exon->Interval;
                    NewExon->Length = 0;
                    if (ReverseFlag)
                    {
                        NewExon->Start = Exon->Start;
                        NewExon->End = Exon->Start + strlen(Exon->Suffix);
                        strcpy(NewExon->Suffix, Exon->Suffix);
                    }
                    else
                    {
                        NewExon->Start = Exon->End - strlen(Exon->Suffix);
                        NewExon->End = Exon->End;
                        strcpy(NewExon->Suffix, Exon->Suffix);
                    }
                    // Assimilate all edges with AA into the new exon:
                    Link = Exon->FirstForward;
                    while (Link)
                    {
                        if (!Link->AA)
                        {
                            NextLink = Link->Next;
                            DeleteExonLink(Exon, Link, 1);
                            Link = NextLink;
                            continue;
                        }
                        // Fix the reciprocal links to point to the new exon:
                        for (OtherLink = Link->Exon->FirstBack; OtherLink; OtherLink = OtherLink->Next)
                        {
                            if (OtherLink->Exon == Exon && OtherLink->AA == Link->AA)
                            {
                                OtherLink->Exon = NewExon;
                            }
                        }
                        if (NewExon->FirstForward)
                        {
                            NewExon->LastForward->Next = Link;
                        }
                        else
                        {
                            NewExon->FirstForward = Link;
                        }
                        NewExon->LastForward = Link;
                        Link = Link->Next;
                    }
                    if (NewExon->LastForward)
                    {
                        NewExon->LastForward->Next = NULL;
                    }
                    Exon->FirstForward = NULL;
                    Exon->LastForward = NULL;
                    GNode->Interval->LastExon->Next = NewExon;
                    GNode->Interval->LastExon = NewExon;
                }
                NextExon = Exon->Next;
                DeleteExon(Exon);
                Exon = NextExon;
                continue;
            } // if exon contains only stop codons
            if (FirstStop == -1)
            {
                // This exon contains no stop codons.  Leave it alone, move on.
                Exon = Exon->Next;
                continue;
            }
            if (FirstStop == 0)
            {
                // This exon begins with one or more stop codons.  Delete all back edges except
                // those containing an amino acid:
                if (AAEdgeBack)
                {
                    NewExon = (ExonNode*)calloc(1, sizeof(ExonNode));
                    NewExon->Interval = Exon->Interval;
                    NewExon->Length = 0;
                    if (ReverseFlag)
                    {
                        NewExon->Start = Exon->End - strlen(Exon->Prefix);
                        NewExon->End = Exon->End;
                    }
                    else
                    {
                        NewExon->Start = Exon->Start;
                        NewExon->End = Exon->Start + strlen(Exon->Prefix);
                    }
                    strcpy(NewExon->Prefix, Exon->Prefix);
                    // Assimilate all edges with AA into the new exon:
                    Link = Exon->FirstBack;
                    while (Link)
                    {
                        if (!Link->AA)
                        {
                            NextLink = Link->Next;
                            DeleteExonLink(Exon, Link, 0);
                            Link = NextLink;
                            continue;
                        }
                        // Fix the reciprocal links to point to the new exon:
                        for (OtherLink = Link->Exon->FirstForward; OtherLink; OtherLink = OtherLink->Next)
                        {
                            if (OtherLink->Exon == Exon)
                            {
                                OtherLink->Exon = NewExon;
                            }
                        }
                        if (NewExon->FirstBack)
                        {
                            NewExon->LastBack->Next = Link;
                        }
                        else
                        {
                            NewExon->FirstBack = Link;
                        }
                        NewExon->LastBack = Link;
                        Link = Link->Next;
                    }
                    if (NewExon->LastBack)
                    {
                        NewExon->LastBack->Next = NULL;
                    }
                    Exon->FirstBack = NULL;
                    Exon->LastBack = NULL;
                    GNode->Interval->LastExon->Next = NewExon;
                    GNode->Interval->LastExon = NewExon;
                }
                else
                {
                    Link = Exon->FirstBack;
                    while (Link)
                    {
                        NextLink = Link->Next;
                        DeleteExonLink(Exon, Link, 0);
                        Link = NextLink;
                    }
                }
                // Delete the exon's prefix, and move its start position up:
                if (ReverseFlag)
                {
                    Exon->End -= strlen(Exon->Prefix) + 3 * FirstReal;
                }
                else
                {
                    Exon->Start += strlen(Exon->Prefix) + 3 * FirstReal;
                }
                Exon->Prefix[0] = '\0';
                NewSequence = (char*)calloc(Exon->Length - FirstReal + 1, sizeof(char)); // 1 byte for null
                strcpy(NewSequence, Exon->Sequence + FirstReal);
                SafeFree(Exon->Sequence);
                Exon->Sequence = NewSequence;
                Exon->Length = Exon->Length - FirstReal;
                // We'll revisit this exon in the next loop pass, in case it has more stop codons later on.
                continue;
            }  // if sequence starts with stop codon

            // This exon contains a stop codon, preceded by 'real' AAs.  Build a new exon 
            // to hold our suffix.  The old exon gets truncated and gets its genomic 
            // pos changed.
            NewExon = (ExonNode*)calloc(1, sizeof(ExonNode));
            NewExon->Interval = Exon->Interval;
            NewExon->Length = Exon->Length - FirstStop - 1;
            NewExon->Sequence = (char*)calloc(NewExon->Length + 1, sizeof(char));
            strcpy(NewExon->Sequence, Exon->Sequence + FirstStop + 1);
            if (ReverseFlag)
            {
                NewExon->Start = Exon->Start;
                NewExon->End = Exon->End - strlen(Exon->Prefix) - (FirstStop + 1) * 3;
                Exon->Start = NewExon->End + 3;
            }
            else
            {
                NewExon->End = Exon->End;
                NewExon->Start = Exon->Start + strlen(Exon->Prefix) + (FirstStop + 1) * 3;
                Exon->End = NewExon->Start - 3;
            }

            NewExon->FirstForward = Exon->FirstForward;
            NewExon->LastForward = Exon->LastForward;
            for (Link = Exon->FirstForward; Link; Link = Link->Next)
            {
                for (OtherLink = Link->Exon->FirstBack; OtherLink; OtherLink = OtherLink->Next)
                {
                    if (OtherLink->Exon == Exon)
                    {
                        OtherLink->Exon = NewExon;
                    }
                }
            }
            Exon->FirstForward = NULL;
            Exon->LastForward = NULL;
            Exon->Sequence[FirstStop] = '\0';
            Exon->Length = FirstStop;
            strcpy(NewExon->Suffix, Exon->Suffix);
            Exon->Suffix[0] = '\0';
            GNode->Interval->LastExon->Next = NewExon;
            GNode->Interval->LastExon = NewExon;
            Exon = Exon->Next;
            continue;
        } // Iteration over exons
    } // Iteration over GNodes (intervals)
} // PurgeNonCodingExonChunks

int SetExonLinkExtensionLengthsBack(ExonNode* Exon, int MinimumORFLength, int IncludeBody)
{
    ExonLink* Link;
    int Length;
    int Pos;
    //

    // First case: We're starting INTO the exon sequence.  We may stop partway.
    if (IncludeBody && Exon->Sequence)
    {
        Length = 0;
        for (Pos = Exon->Length - 1; Pos >= 0; Pos--)
        {
            if (Exon->Sequence[Pos] == 'X')
            {
                return Length;
                break;
            }
            Length++;
            if (Length >= MinimumORFLength)
            {
                // That's long enough already!
                return Length;
            }
        }
    }
    else
    {
        Length = 0;
    }

    if (Exon->MaxBackOverall != -1)
    {
        return Length + Exon->MaxBackOverall;
    }

    // Set Link->MaxLength for each link back.  We always do this when we're called with
    // IncludeBody == 0, we MAY do it for IncludeBody == 1 (as necessary)
    Exon->MaxBackOverall = 0;
    for (Link = Exon->FirstBack; Link; Link = Link->Next)
    {
        // If the extension for the link is already known, note it and continue
        if (Link->MaxLength != -1)
        {
            Exon->MaxBackOverall = max(Exon->MaxBackOverall, Link->MaxLength);
            continue;
        }
        if (Link->AA)
        {
            Link->MaxLength = 1;
        }
        else
        {
            Link->MaxLength = 0;
        }
        Link->MaxLength += SetExonLinkExtensionLengthsBack(Link->Exon, MinimumORFLength, 1);
        Exon->MaxBackOverall = max(Exon->MaxBackOverall, Link->MaxLength);
    }
    return Length + Exon->MaxBackOverall;
}

int SetExonLinkExtensionLengthsForward(ExonNode* Exon, int MinimumORFLength, int IncludeBody)
{
    ExonLink* Link;
    int Length;
    int Pos;
    //

    // First case: We're starting INTO the exon sequence.  We may stop partway.
    if (IncludeBody && Exon->Sequence)
    {
        Length = 0;
        for (Pos = 0; Pos < Exon->Length; Pos++)
        {
            if (Exon->Sequence[Pos] == 'X')
            {
                return Length;
                break;
            }
            Length++;
            if (Length >= MinimumORFLength)
            {
                // That's long enough already!
                return Length;
            }
        }
    }
    else
    {
        Length = 0;
    }

    if (Exon->MaxForwardOverall != -1)
    {
        return Length + Exon->MaxForwardOverall;
    }

    // Set Link->MaxLength for each link back.  We always do this when we're called with
    // IncludeBody==0, we MAY do it for IncludeBody==1 (as necessary)
    Exon->MaxForwardOverall = 0;
    for (Link = Exon->FirstForward; Link; Link = Link->Next)
    {
        // If the extension for the link is already known, note it and continue
        if (Link->MaxLength != -1)
        {
            Exon->MaxForwardOverall = max(Exon->MaxForwardOverall, Link->MaxLength);
            continue;
        }
        if (Link->AA)
        {
            Link->MaxLength = 1;
        }
        else
        {
            Link->MaxLength = 0;
        }
        Link->MaxLength += SetExonLinkExtensionLengthsForward(Link->Exon, MinimumORFLength, 1);
        Exon->MaxForwardOverall = max(Exon->MaxForwardOverall, Link->MaxLength);
    }
    return Length + Exon->MaxForwardOverall;
}

// Short open reading frame pruning:
// - Determine the maximum length of each reading frame attainable by linking forward along the graph F1...Fm
// - Determine the maximum length of each reading frame attainable by linking backward along the graph B1...Bn
// - If there are no stop codons: 
//   If len + max(F1...Fm) + max(B1...Bn) < N, kill the exon and all links
//   Else:
//     If len + max(B1...Bn) + Fi < N, kill forward link i
//     If len + max(F1...Fm) + Bi < N, kill backward link i
// - Let S be the length of the suffix (all AAs after the last stop codon)
//   If S + max(F1...Fm) < N, mask S and remove all forward links
//   Else if S + Fi < N, remove Fi
// - Let P be the length of the prefix (all AAs up to the first stop codon)
//   If P + max(B1...Bn) < N, mask P and remove all backward links
//   Else if P + Bi < N, remove Bi
void PruneShortORFs(int ReverseFlag, int MinimumORFLength)
{
    GeneNode* GNode;
    ExonNode* Exon;
    ExonNode* NextExon = NULL;
    ExonLink* Link;
    ExonLink* NextLink;
    int LinkIndex;
    int PrefixLength;
    int SuffixLength;
    int Pos;
    char* NewSequence;
    int CutsPerformed;
    int Flag;
    int GeneNodeIndex = 0;
    int ExonIndex;
    //
    // if MinOrfLength <= 0, then don't filter.
    if (MinimumORFLength <= 0)
    {
        return;
    }
    // Iterate over all exons in all intervals.  Init the link lengths.
    for (GNode = g_GeneFirst; GNode; GNode = GNode->Next)
    {
        for (Exon = GNode->Interval->FirstExon; Exon; Exon = Exon->Next)
        {
            for (Link = Exon->FirstForward; Link; Link = Link->Next)
            {
                Link->MaxLength = -1;
            }
            for (Link = Exon->FirstBack; Link; Link = Link->Next)
            {
                Link->MaxLength = -1;
            }
            // Exon->MaxBackOverall is set to -1 to indicate that it hasn't been
            // processed yet.
            Exon->MaxBackOverall = -1;
            Exon->MaxForwardOverall = -1;
        }
    }
    // Iterate over all exons in all intervals.  Set the max lengths of all
    // their links.
    for (GNode = g_GeneFirst; GNode; GNode = GNode->Next)
    {
        for (Exon = GNode->Interval->FirstExon; Exon; Exon = Exon->Next)
        {
            if (Exon->MaxBackOverall == -1)
            {
                SetExonLinkExtensionLengthsBack(Exon, MinimumORFLength, 0);
            }
            if (Exon->MaxForwardOverall == -1)
            {
                SetExonLinkExtensionLengthsForward(Exon, MinimumORFLength, 0);
            }
        }
    }

    for (GNode = g_GeneFirst; GNode; GNode = GNode->Next)
    {
        //printf("Start gene node #%d\n", GeneNodeIndex);
        Exon = GNode->Interval->FirstExon;
        ExonIndex = 0;
        while (Exon)
        {
            //printf("Start GeneNode#%d exon#%d\n", GeneNodeIndex, ExonIndex);
            ExonIndex++;
            
            // Measure the exon, its prefix, and its suffix.
            PrefixLength = 0;
            SuffixLength = 0;
            if (Exon->Sequence)
            {
                for (Pos = 0; Pos < Exon->Length; Pos++)
                {
                    if (Exon->Sequence[Pos] == 'X')
                    {
                        break;
                    }
                    PrefixLength++;
                }
                for (Pos = Exon->Length - 1; Pos >= 0; Pos--)
                {
                    if (Exon->Sequence[Pos] == 'X')
                    {
                        break;
                    }
                    SuffixLength++;
                }
            }
            // Consider removing the exon entirely:
            if (Exon->Length + Exon->MaxBackOverall + Exon->MaxForwardOverall < MinimumORFLength)
            {
                //printf("*Delete the exon entirely!\n");
                // Zap!  Free the exon, and its links.
                NextExon = Exon->Next;
                DeleteExon(Exon);
                Exon = NextExon;
                continue;
            }

            if (PrefixLength == Exon->Length)
            {
                // This exon contains no stop codons.  And we cannot delete it entirely,
                // but we can perhaps still prune a few links.
                // Try removing forward links:
                Link = Exon->FirstForward;
                LinkIndex = 0;
                while (Link)
                {
                    if (Exon->MaxBackOverall + Exon->Length + Link->MaxLength < MinimumORFLength)
                    {
                        // This link can't be part of a full-length ORF.  So, let's remove the link:
                        NextLink = Link->Next;
                        //printf("*Delete a forward link\n");
                        DeleteExonLink(Exon, Link, 1);
                        Link = NextLink;
                        LinkIndex++;
                        continue;
                    }
                    LinkIndex++;
                    Link = Link->Next;
                }
                // Try removing backward links:
                Link = Exon->FirstBack;
                LinkIndex = 0;
                while (Link)
                {
                    if (Exon->MaxForwardOverall + Exon->Length + Link->MaxLength < MinimumORFLength)
                    {
                        // This link can't be part of a full-length ORF.  So, let's remove the link:
                        NextLink = Link->Next;
                        //printf("*Delete a backward link\n");
                        DeleteExonLink(Exon, Link, 0);
                        Link = NextLink;
                        LinkIndex++;
                        continue;
                    }
                    LinkIndex++;
                    Link = Link->Next;
                }
            }
            else
            {
                // This exon contains at least one stop codon.  We'll consider pruning the 
                // prefix (everything up to and including the stop), the suffix (the stop
                // and everything beyond it).
                CutsPerformed = 0;
                // First, consider removing the prefix (or some incoming links):
                if (PrefixLength + Exon->MaxBackOverall < MinimumORFLength)
                {
                    // We can cut the prefix!  First delete all backward links:
                    Link = Exon->FirstBack;
                    while (Link)
                    {
                        NextLink = Link->Next;
                        DeleteExonLink(Exon, Link, 0);
                        Link = NextLink;
                    }
                    // If there's at least one non-stop character after the prefix,
                    // then the exon still has a sequence:
                    if (Exon->Length > PrefixLength + 1)
                    {
                        NewSequence = (char*)calloc(Exon->Length - PrefixLength, sizeof(char)); 
                        strcpy(NewSequence, Exon->Sequence + PrefixLength + 1);
                        Exon->Length = Exon->Length - PrefixLength - 1;
                        SafeFree(Exon->Sequence);
                        Exon->Sequence = NewSequence;
                        // Fix the genomic start coordinate:
                        if (ReverseFlag)
                        {
                            Exon->End -= strlen(Exon->Prefix) + (PrefixLength + 1)*3;
                        }
                        else
                        {
                            Exon->Start += strlen(Exon->Prefix) + (PrefixLength + 1)*3;
                        }
                    }
                    else
                    {
                        // The exon's body is gone!
                        if (ReverseFlag)
                        {
                            Exon->End -= strlen(Exon->Prefix) + Exon->Length * 3;
                        }
                        else
                        {
                            Exon->Start += strlen(Exon->Prefix) + Exon->Length * 3;
                        }
                        Exon->Length = 0;
                        SafeFree(Exon->Sequence);
                        Exon->Sequence = NULL;
                    }
                    CutsPerformed++;
                    Exon->Prefix[0] = '\0';
                }
                // Consider removing the suffix (or some outgoing links):
                if (SuffixLength + Exon->MaxForwardOverall < MinimumORFLength)
                {
                    // Delete all forward links:
                    //printf("*Delete all forward links\n");
                    Link = Exon->FirstForward;
                    while (Link)
                    {
                        NextLink = Link->Next;
                        DeleteExonLink(Exon, Link, 1);
                        Link = NextLink;
                    }
                    if (Exon->Length > SuffixLength + 1)
                    {
                        NewSequence = (char*)calloc(Exon->Length - SuffixLength, sizeof(char));
                        strncpy(NewSequence, Exon->Sequence, Exon->Length - SuffixLength - 1);
                        Exon->Length = Exon->Length - SuffixLength - 1;
                        NewSequence[Exon->Length] = '\0';
                        SafeFree(Exon->Sequence);
                        Exon->Sequence = NewSequence;
                        if (ReverseFlag)
                        {
                            Exon->Start += strlen(Exon->Suffix) + (SuffixLength + 1) * 3;
                        }
                        else
                        {
                            Exon->End -= strlen(Exon->Suffix) + (SuffixLength + 1) * 3;
                        }
                    }
                    else
                    {
                        // The exon's body is gone!
                        if (ReverseFlag)
                        {
                            Exon->Start += strlen(Exon->Suffix) + Exon->Length * 3;
                        }
                        else
                        {
                            Exon->End -= strlen(Exon->Suffix) + Exon->Length * 3;
                        }
                        Exon->Length = 0;
                        SafeFree(Exon->Sequence);
                        Exon->Sequence = NULL;

                    }
                    Exon->Suffix[0] = '\0';
                    CutsPerformed++;
                }
                // If we removed the exon body, and there's no prefix or suffix left, then cut the exon itself:
                if (Exon->Start == Exon->End)
                {
                    NextExon = Exon->Next;
                    DeleteExon(Exon);
                    Exon = NextExon;
                    continue;
                }
                // If we cut a prefix and cut a suffix, then we have no links, and it's entirely possible 
                // that no real sequence remains. Check to see that we have a non-stop-codon:
                if (CutsPerformed == 2) 
                {
                    Flag = 0;
                    for (Pos = 0; Pos < Exon->Length; Pos++)
                    {
                        if (Exon->Sequence[Pos] != 'X')
                        {
                            Flag = 1;
                            break;
                        }
                    }
                    if (!Flag)
                    {
                        NextExon = Exon->Next;
                        //printf("*Delete the exon itself\n");
                        DeleteExon(Exon);
                        Exon = NextExon;
                        continue;
                    }
                }
                // Even if we didn't cut the prefix (or suffix) outright, we may be able to remove
                // some incoming and outgoing links:
                Link = Exon->FirstBack;
                LinkIndex = 0;
                while (Link)
                {
                    if (Link->MaxLength + PrefixLength < MinimumORFLength)
                    {
                        NextLink = Link->Next;
                        //printf("*Delete ONE backward link\n");
                        DeleteExonLink(Exon, Link, 0);
                        Link = NextLink;
                        LinkIndex++;
                        continue;
                    }
                    Link = Link->Next;
                    LinkIndex++;
                }
                Link = Exon->FirstForward;
                LinkIndex = 0;
                while (Link)
                {
                    if (Link->MaxLength + SuffixLength < MinimumORFLength)
                    {
                        NextLink = Link->Next;
                        //printf("*Delete ONE forward link\n");
                        DeleteExonLink(Exon, Link, 1);
                        Link = NextLink;
                        LinkIndex++;
                        continue;
                    }
                    Link = Link->Next;
                    LinkIndex++;
                }
            } // if the exon contains a stop codon
            Exon = Exon->Next;
        } // exon loop
        GeneNodeIndex++;
    } // gene node loop
}

// Prepare a splice-db for a particular target gene.  We'll read all the intervals for the given
// chromosome number, perform the merge and interstect algorithms, and then construct a set 
// of genomic intervals which 'satisfy' the target.
void PrepareOneGeneSpliceDB(int ChromosomeNumber, int ReverseFlag, int IntervalStart, 
    int IntervalEnd, char* CustomFileName, char* GeneName, int MinORFLength)
{
    char FileName[1024];
    char GenomeFileName[1024];
    FILE* GenomicFile;
    FILE* CustomFile;
    GeneNode* GNode;
    IntervalNode* Interval;
    int SatisfiedOne;
    char ReverseChar;
    //
    if (ReverseFlag)
    {
        ReverseChar = '-';
    }
    else
    {
        ReverseChar = '+';
    }
    sprintf(FileName, "NewSpliceDB\\%d%c.filtered", ChromosomeNumber, ReverseChar);
    //ParseIntervalsFlatFile(FileName, -1);
    ParseIntervalsESTBinaryFile(FileName);
    
    /////////////////////////////////////////////////////////////////
    //// For purposes of debugging, we can trim the list of intervals a bit.  (Debug printout of a whole
    //// chromosome is unwieldy!)  In production, we MUST NOT trim, because one of the intervals in the
    //// master-interval may be linked to a far-away interval.
    //PruneEdge = IntervalStart - 5000;
    //while (g_FirstInterval->End < PruneEdge)
    //{
    //    RemoveInterval(g_FirstInterval, 0);
    //}
    //PruneEdge = IntervalEnd + 5000;
    //while (g_LastInterval->Start > PruneEdge)
    //{
    //    RemoveInterval(g_LastInterval, 0);
    //}

    printf("BEFORE merge:\n");
    DebugPrintIntervals(1, 1, -1, -1);
    MergeIntervals();
    printf("AFTER merge:\n");
    DebugPrintIntervals(1, 2, -1, -1);
    IntersectIntervals();
    printf("AFTER intersect:\n");
    DebugPrintIntervals(1, 3, -1, -1);

    //sprintf(GenomeFileName, "C:\\source\\Bafna\\Splice\\chromFa\\chr%d.trie", ChromosomeNumber);
    sprintf(GenomeFileName, "e:\\Chromosome\\chr%d.trie", ChromosomeNumber);
    GenomicFile = fopen(GenomeFileName, "rb");
    CustomFile = fopen(CustomFileName, "wb");
    // Create the gene node list.  First, add every interval that overlaps
    // the requested 'master interval':
    for (Interval = g_FirstInterval; Interval; Interval = Interval->Next)
    {
        if (Interval->Start > IntervalEnd)
        {
            break;
        }
        if (Interval->End < IntervalStart)
        {
            continue;
        }
        GNode = (GeneNode*)calloc(1, sizeof(GeneNode));
        GNode->Interval = Interval;
        Interval->GNode = GNode;
        if (g_GeneFirst)
        {
            g_GeneLast->Next = GNode;
            GNode->Prev = g_GeneLast;
        }
        else
        {
            g_GeneFirst = GNode;
        }
        g_GeneLast = GNode;
    }
    // Iterate: Find the first interval overlapping the master which is not satisfied.  Then, satisfy it.
    // Break after every interval overlapping the master has been satisfied.
    while (1)
    {
        SatisfiedOne = 0;
        for (GNode = g_GeneFirst; GNode; GNode = GNode->Next)
        {
            if (GNode->Interval->End < IntervalStart || GNode->Interval->Start > IntervalEnd)
            {
                continue; // We needn't satisfy this one, since it's not in the master-interval.
            }
            if (!GNode->Interval->Satisfied)
            {
                SatisfyIntervalForward(GNode, 0);
                SatisfyIntervalBack(GNode, 0);
                GNode->Interval->Satisfied = 1;
                SatisfiedOne = 1;
                break;
            }
        }
        if (!SatisfiedOne)
        {
            // Everyone's happy, so stop now.
            break;
        }
    }
   
    BuildAndWriteExons(GenomicFile, CustomFile, ReverseFlag, GeneName, ChromosomeNumber, MinORFLength);
    fclose(GenomicFile);
    fclose(CustomFile);

    FreeGeneNodes();
    // Free the interval list!
    while (g_FirstInterval)
    {
        RemoveInterval(g_FirstInterval, 0);
    }

}

// Parse a binary file listiing genomic intervals, with links between
// them.  Convert it into an exon graph and write it out.
void PrepareSpliceDB(int ChromosomeNumber, int ReverseFlag, int MinORFLength)
{
    FILE* StatsFile;
    char ReverseChar;
    char FileName[1024];
    char ChromosomeFileName[1024];
    char OutputFileName[1024];
    int GeneCount;

    if (ReverseFlag)
    {
        ReverseChar = '-';
    }
    else
    {
        ReverseChar = '+';
    }
    
    ////////////////////////////////////////////////////////////////////////////
    // HUMAN data sources:
    // We can parse ESTs, or gene finder output, or BOTH.  (Both may be slow)
    // ESTREF, if reference sequences are included, or EST, if only ESTs are included:
    //sprintf(FileName, "ESTREF\\%d%c.filtered", ChromosomeNumber, ReverseChar); // %%% hard-coded path
    sprintf(FileName, "EST\\%d%c.filtered", ChromosomeNumber, ReverseChar); // %%% hard-coded path
    ParseIntervalsESTBinaryFile(FileName);
    sprintf(FileName, "GeneFindDB\\%d%c.dat", ChromosomeNumber, ReverseChar); // %%% hard-coded path
    ParseIntervalsGeneFindBinaryFile(FileName);
    sprintf(ChromosomeFileName, "e:\\Chromosome\\chr%d.trie", ChromosomeNumber);
    sprintf(OutputFileName, "ESTSpliceDB\\%d%c.dat", ChromosomeNumber, ReverseChar);

    printf("BEFORE merge:\n");
    DebugPrintIntervals(-1, 1, -1, -1); 
    MergeIntervals();
    printf("AFTER merge:\n");
    DebugPrintIntervals(-1, 2, -1, -1); 
    IntersectIntervals();
    printf("AFTER intersect:\n");
    DebugPrintIntervals(-1, 3, -1, -1); 
    
    GeneCount = WriteGenesForIntervals(ChromosomeFileName, OutputFileName, ChromosomeNumber, ReverseFlag, MinORFLength);
    //WriteGenesForIntervals("C:\\source\\Bafna\\Splice\\chromFa\\chr11.trie", "ESTSpliceDB\\11-.dat", 11, 1);
    StatsFile = fopen("SplicePrepStats.txt", "a+");
    //fprintf(StatsFile, "Genes have been written out.  Statistics:\n");
    fprintf(StatsFile, "%d\t", ChromosomeNumber);
    fprintf(StatsFile, "%d\t", ReverseFlag);
    fprintf(StatsFile, "%d\t", GeneCount);
    fprintf(StatsFile, "%d\t", g_StatsIncompleteGeneCount);
    fprintf(StatsFile, "%d\t", g_StatsLargestGeneSize);
    fprintf(StatsFile, "%d\t", g_StatsLargestGeneRecordNumber);
    fprintf(StatsFile, "%d\t", g_StatsIntervalsBeforeMerge);
    fprintf(StatsFile, "%d\t", g_StatsEdgesBeforeMerge);
    fprintf(StatsFile, "%d\t", g_StatsIntervalsAfterMerge);
    fprintf(StatsFile, "%d\t", g_StatsEdgesAfterMerge);
    fprintf(StatsFile, "%d\t", g_StatsIntervalsAfterIntersect);
    fprintf(StatsFile, "%d\t", g_StatsEdgesAfterIntersect);
    fprintf(StatsFile, "%d\t", g_StatsTotalExonsWritten);
    fprintf(StatsFile, "%d\t", g_StatsTotalEdgesWritten);
    printf("Exon counts:\n");
    fprintf(StatsFile, "\t%d\t", g_PolymorphismCount);
    printf("%d\t%d\t%d\t\t%d\t%d\t%d\t", g_UFTotalExons, g_UFTotalEdges, g_UFTotalAA, g_TotalExons, g_TotalEdges, g_TotalAA);
    fprintf(StatsFile, "\t%d\t%d\t%d\t\t%d\t%d\t%d\t", g_UFTotalExons, g_UFTotalEdges, g_UFTotalAA, g_TotalExons, g_TotalEdges, g_TotalAA);
    // How many exons are there...if you count adjacent exons as a single 'real' exon?
    fprintf(StatsFile, "\t%d\t%d\t", g_TotalTrueExons, g_TotalTrueEdges);
    fprintf(StatsFile, "\n");
    fclose(StatsFile);
    
    // Free the interval list!
    while (g_FirstInterval)
    {
        RemoveInterval(g_FirstInterval, 0);
    }
}
