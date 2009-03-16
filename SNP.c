
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
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include "Trie.h"
#include "Utils.h"
#include "Run.h"
#include "Tagger.h"
#include "Score.h"
#include "FreeMod.h"
#include "Spliced.h"
#include "BN.h"
#include "SVM.h"
#include "Scorpion.h"
#include "ChargeState.h"
#include "SNP.h"

// Code to support the inclusion of POLYMORPHISMS, particularly SNPs, in a protein database.
// The motivation: If a protein has several polymorphic sites, we can include all of them in 
// a string-table proteomic database only by including multiple isoforms.  Since we're already
// using a DAG data-structure to capture alternative splicing, let's include SNPs in the DAG,
// and capture polymorphic variability as well as splicing variability.  

// We use PolyNodes during database construction (not needed during search) to keep track of all
// polymorphisms.  The polynodes are read from a binary file (currently written out by 
// ParseSNPDatabase.py while parsing snp.txt from the ucsc genome browser).  They're ordered 
// by genomic position.

PolyNode* g_FirstPolyNode = NULL;
PolyNode* g_LastPolyNode = NULL;

Polymorphism* g_Polymorphisms = NULL; // array
int g_PolymorphismCount;

// Search for the FIRST polymorphism that lies within the given interval.
// Return its index.  Return -1 if there's no polymorphism in that interval.
// Simple binary search.
int FindPolyInInterval(int Start, int End)
{
    int Low;
    int High;
    int Mid;
    int Pos;
    //
    if (!g_PolymorphismCount)
    {
        return -1;
    }
    Low = 0;
    High = g_PolymorphismCount - 1;
    while (1)
    {
        // If we're down to a minimally-sized poly-interval, check it and return:
        if (Low + 1 >= High)
        {
            Pos = g_Polymorphisms[Low].Pos;
            if (Pos >= Start && Pos < End)
            {
                return Low;
            }
            Pos = g_Polymorphisms[High].Pos;
            if (Pos >= Start && Pos < End)
            {
                return High;
            }
            return -1;
        }

        Mid = (Low + High) / 2;
        Pos = g_Polymorphisms[Mid].Pos;
        if (Pos < Start)
        {
            Low = Mid;
            continue;
        }
        if (Pos >= End)
        {
            High = Mid;
            continue;
        }
        // We found one!  Make sure we have the FIRST one.
        for (Low = Mid; Low >= 0; Low--)
        {
            if (g_Polymorphisms[Low].Pos < Start)
            {
                return (Low + 1);
            }
        }
        return 0;
    }
}

// Free the full linked list of poly nodes. 
void FreePolyNodes()
{
    PolyNode* Node;
    PolyNode* Prev = NULL;
    //
    if (!g_FirstPolyNode)
    {
        return;
    }
    for (Node = g_FirstPolyNode; Node; Node = Node->Next)
    {
        SafeFree(Prev);
        Prev = Node;
    }
    SafeFree(Prev);
    g_FirstPolyNode = NULL;
    g_LastPolyNode = NULL;
}

// Parse polymorphism nodes for the current chromosome.
void ParsePolyNodes(char* FileName)
{
    FILE* File;
    PolyNode* Node;
    int BytesRead;
    int GenomicPosition;
    int RecordNumber;
    //
    File = fopen(FileName, "rb");
    if (!File)
    {
        printf("** Error: Unable to open polymorphism database '%s'\n", FileName);
        return;
    }
    RecordNumber = 0;
    while (1)
    {
        BytesRead = ReadBinary(&GenomicPosition, sizeof(int), 1, File);
        if (!BytesRead)
        {
            break;
        }
        Node = (PolyNode*)calloc(sizeof(PolyNode), 1);
        Node->Pos = GenomicPosition;
        ReadBinary(&Node->Type, sizeof(char), 1, File);
        switch (Node->Type)
        {
        case 0:
            ReadBinary(Node->SNP, sizeof(char), 2, File);
            break;
        case 1:
            ReadBinary(Node->SNP, sizeof(char), 3, File);
            break;
        case 2:
            ReadBinary(Node->SNP, sizeof(char), 4, File);
            break;
        default:
            printf("** Error: Unable to parse polymorphism node %d type '%d'\n", RecordNumber, Node->Type);
            break;
        }
        if (g_LastPolyNode)
        {
            g_LastPolyNode->Next = Node;
            // Sanity check: These nodes MUST come in order.
            if (g_LastPolyNode->Pos >= Node->Pos)
            {
                printf("** Error parsing polymorphism data: Record %d is out of order!  (Start %d vs %d)\n", RecordNumber, Node->Pos, g_LastPolyNode->Pos);
            }
        }
        else
        {
            g_FirstPolyNode = Node;
        }
        g_LastPolyNode = Node;
        RecordNumber++;
    }
    fclose(File);
    ////////////////////////////////////////////////////////////////////////
    // Now, put all those nodes into an array:
    g_PolymorphismCount = RecordNumber;
    g_Polymorphisms = (Polymorphism*)calloc(g_PolymorphismCount, sizeof(Polymorphism));
    RecordNumber = 0;
    for (Node = g_FirstPolyNode; Node; Node = Node->Next)
    {
        g_Polymorphisms[RecordNumber].Pos = Node->Pos;
        memcpy(g_Polymorphisms[RecordNumber].SNP, Node->SNP, sizeof(char) * 4);
        RecordNumber++;
    }
    FreePolyNodes();
}

// For debugging: Print out all the polymorphism nodes.
void DebugPrintPolyNodes(int FirstRecord, int LastRecord)
{
    PolyNode* Node;
    int RecordNumber;
    //
    RecordNumber = 0;
    for (Node = g_FirstPolyNode; Node; Node = Node->Next)
    {
        if (FirstRecord >= 0 && RecordNumber < FirstRecord)
        {
            continue;
        }
        if (LastRecord >= 0 && RecordNumber > LastRecord)
        {
            continue;
        }
        printf("SNP record %d: Pos %d can be %c or %c\n", RecordNumber, Node->Pos, Node->SNP[0], Node->SNP[1]);
        RecordNumber++;
    }
}

void SNPTestMain()
{
    ParsePolyNodes("SNP\\1.snp");
    DebugPrintPolyNodes(-1, -1);
}
