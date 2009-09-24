
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
#include "Inspect.h"
#include "Utils.h"
#include "Mods.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

// AllKnownPTMods is initialized once and only once.  AllowedPTMods is a sub-array, 
// set before doing a search or batch of searches.  AllPTModCount is the size
// of array AllKnownPTMods, and AllowedPTModCount is the size of array AllowedPTMods.
PTMod AllKnownPTMods[MAX_PT_MODTYPE];
int AllPTModCount = 0;

int g_PhosphorylationMod = -1;

// PTMLimit[n] is a limit on how many modifications of type n can be placed
// on a peptide.  For each n, PTMLimit[n] <= GlobalOptions->MaxPTMods.
int g_PTMLimit[MAX_PT_MODTYPE];

int PlainOldDecorationIndex = 0;

Decoration* AllDecorations = NULL;
int AllDecorationCount = 0;
int AllDecorationAllocation = 0;

int CompareDecorations(const Decoration* A, const Decoration* B)
{
    if (A->Mass < B->Mass)
    {
        return -1;
    }
    if (A->Mass > B->Mass)
    {
        return 1;
    }
    return 0;
}

void ExpandDecorationList(int SourceDecorationIndex, int MinPTModIndex, int* PTMRemaining, int ModsLeft)
{
    int PTModIndex;
    int Decor;
    //
    if (ModsLeft <= 0)
    {
        return;
    }
    for (PTModIndex = MinPTModIndex; PTModIndex < AllPTModCount; PTModIndex++)
    {
        if (PTMRemaining[PTModIndex] <= 0)
        {
            continue;
        }
        // If we have a lot of decorations, expand the memory available for them:
        if (AllDecorationCount == AllDecorationAllocation-1)
        {
            AllDecorationAllocation *= 2;
            AllDecorations = (Decoration*)realloc(AllDecorations, sizeof(Decoration) * AllDecorationAllocation);
        }
        Decor = AllDecorationCount;
        AllDecorationCount++;
        
        memcpy(AllDecorations[Decor].Mods, AllDecorations[SourceDecorationIndex].Mods, sizeof(int) * MAX_PT_MODTYPE);
        AllDecorations[Decor].Mods[PTModIndex]++;
        AllDecorations[Decor].TotalMods = AllDecorations[SourceDecorationIndex].TotalMods + 1;
        AllDecorations[Decor].Mass = AllDecorations[SourceDecorationIndex].Mass + MassDeltaByIndex[AMINO_ACIDS*MAX_PT_MODTYPE + PTModIndex]->RealDelta; 

        //printf("Added decoration %d (%.2f) ", Decor, AllDecorations[Decor].Mass);
        //for (ModIndex = 0; ModIndex < AllPTModCount; ModIndex++)
        //{
        //    printf("%d ", AllDecorations[Decor].Mods[ModIndex]);
        //}
        //printf("\n");

        PTMRemaining[PTModIndex] -= 1;
        ExpandDecorationList(Decor, PTModIndex, PTMRemaining, ModsLeft - 1);
        PTMRemaining[PTModIndex] += 1;
    }
}

// Entries of form IsSubDecoration[DecorIndex][OtherDecorIndex]
int** IsSubDecoration = NULL;


// After reading the definitions of all the post-translational modifications, we construct 
// a list of decorations.
// Special case:  If GlobalOptions->MandatoryModName is set, then we set MandatoryModIndex, and
// we only allow decorations that *do* contain that mod.
void BuildDecorations()
{
    int DecorIndex;
    int OtherDecorIndex;
    int ModIndex;
    int ValidSubDecoration;
    int PTMRemaining[MAX_PT_MODTYPE];
    int TotalPTMsPermitted;
    //

    // Free the old IsSubDecoration array, if allocated:
    if (IsSubDecoration)
    {
        for (DecorIndex = 0; DecorIndex < AllDecorationCount; DecorIndex++)
        {
            SafeFree(IsSubDecoration[DecorIndex]);
        }
        SafeFree(IsSubDecoration);
        IsSubDecoration = NULL;
    }
    AllDecorationAllocation = 100;
    SafeFree(AllDecorations); // Remove old ones!
    AllDecorations = NULL;
    AllDecorations = (Decoration*)calloc(AllDecorationAllocation, sizeof(Decoration));
    // AllDecorations[0] is now prepared.  (Mass 0, no mods)
    AllDecorationCount = 1;
    memcpy(PTMRemaining, g_PTMLimit, sizeof(int) * MAX_PT_MODTYPE);
    TotalPTMsPermitted = GlobalOptions->MaxPTMods;
    ExpandDecorationList(0, 0, PTMRemaining, TotalPTMsPermitted);
    qsort(AllDecorations, AllDecorationCount, sizeof(Decoration), (QSortCompare)CompareDecorations);
    // Locate the index of the unmodified null-decoration.  (Usually it's #0, because
    // it has mass 0, but it's possible for PTMs to have a *negative* mass)
    for (DecorIndex = 0; DecorIndex < AllDecorationCount; DecorIndex++)
    {
        if (AllDecorations[DecorIndex].TotalMods == 0)
        {
            PlainOldDecorationIndex = DecorIndex;
            break;
        }
    }
    for (ModIndex = 0; ModIndex < AllPTModCount; ModIndex++)
    {
        if (!CompareStrings(GlobalOptions->MandatoryModName, MassDeltaByIndex[AMINO_ACIDS*MAX_PT_MODTYPE + ModIndex]->Name))
        {
            GlobalOptions->MandatoryModIndex = ModIndex;
        }
    }

    IsSubDecoration = (int**)calloc(AllDecorationCount, sizeof(int*));
    for (DecorIndex = 0; DecorIndex < AllDecorationCount; DecorIndex++)
    {
        IsSubDecoration[DecorIndex] = (int*)calloc(AllDecorationCount, sizeof(int));
        for (OtherDecorIndex = 0; OtherDecorIndex < AllDecorationCount; OtherDecorIndex++)
        {
            ValidSubDecoration = 1; // default
            for (ModIndex = 0; ModIndex < AllPTModCount; ModIndex++)
            {
                if (AllDecorations[OtherDecorIndex].Mods[ModIndex] < AllDecorations[DecorIndex].Mods[ModIndex])
                {
                    ValidSubDecoration = 0;
                    break;
                }
            }
            if (ValidSubDecoration)
            {
                IsSubDecoration[DecorIndex][OtherDecorIndex] = 1;
            }
        }
    }
}

void FreeIsSubDecoration()
{
    int ModIndex;
    for (ModIndex = 0; ModIndex < AllDecorationCount; ModIndex++)
    {
        SafeFree(IsSubDecoration[ModIndex]);
        IsSubDecoration[ModIndex] = NULL;
    }
    SafeFree(IsSubDecoration);
    IsSubDecoration = NULL;
}

// Returns a PTM with this name.  Returns NULL if no match found.
// Case-insensitive (pHoSpHoRyLaTiOn is ok).
MassDelta* FindPTModByName(char Amino, char* Name)
{
    int ModIndex;
    int AminoIndex = Amino - 'A';
    for (ModIndex = 0; ModIndex < GlobalOptions->DeltasPerAA; ModIndex++)
    {
        if (!MassDeltas[AminoIndex][ModIndex].Flags)
        {
            break;
        }
        if (!CompareStrings(MassDeltas[AminoIndex][ModIndex].Name, Name))
        {
            return &MassDeltas[AminoIndex][ModIndex];
        }
    }
    return NULL;
}
