
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

#ifndef MODS_H
#define MODS_H

// Structs to support search with post-translational modifications.

#include "Inspect.h"

#define DELTA_FLAG_VALID 1
#define DELTA_FLAG_PHOSPHORYLATION 2
#define DELTA_FLAG_C_TERMINAL 4
#define DELTA_FLAG_N_TERMINAL 8

// A MassDeltaNode is part of a linked list, each node of which wraps
// a MassDelta.  Given a modification that affects multiple amino acids (e.g.
// oxidation of M or W), we build one MassDelta struct...but there's one
// MassDeltaNode for MassDeltas[M] and one for MassDeltas[W].
typedef struct MassDeltaNode
{
    struct MassDelta* Delta;
    struct MassDeltaNode* Next;
} MassDeltaNode;

typedef struct MassDelta
{ 
    float Score; 
    int Delta;  // in bin-units
    int RealDelta; // in actual mass-units.  RealDelta = Delta * 10
    char Name[21]; 
    char Amino; // if this is a mutation to an amino acid
    int Flags; // Used for noting which is phosphorylation.  If flags == 0, this record is an end-of-array marker!

    // Index of type.  For instance, all phosphorylations have same index.  Offset into AllKnownPTMods.
    int Index; 
} MassDelta;

extern MassDelta** MassDeltas;
extern MassDelta** MassDeltaByIndex;
// A decoration is a collection of post-translational modification.  This includes the
// 'empty decoration', with no modifications, and mass 0.  When we examine the flanking
// regions of a tag match to see whether the masses are valid, we consider each possible
// decoration.  (For instance, if the prefix mass is too low by 80 but phosphorylation
// is available - and there's a phosphorylatable residue in the prefix - then we have a match
// via the decoration of mass 80)

typedef struct PTMod
{
    char Name[40];
    int Mass;
    // How many of this modification can be attached to a base?  (Generally zero or one!)
    int Allowed[TRIE_CHILD_COUNT]; 
    int Flags;
} PTMod;

typedef struct Decoration
{
    int Mass;
    int TotalMods;
    int Mods[MAX_PT_MODTYPE]; // Decoration->Mods[n] is how many MassDeltas with Index of n are in this decoration.
} Decoration;

extern Decoration* AllDecorations;
extern int AllDecorationCount;
extern int PlainOldDecorationIndex;

// AllKnownPTMods - initialized at startup
extern PTMod AllKnownPTMods[MAX_PT_MODTYPE];
extern int AllPTModCount;
extern int g_PTMLimit[MAX_PT_MODTYPE];
extern int g_PhosphorylationMod; // index of the phosphorylation PTM

// Returns index of the PTM with this name.  Returns -1 if no match found.
// Case-insensitive (pHoSpHoRyLaTiOn is ok).
MassDelta* FindPTModByName(char Amino, char* Name);

void BuildDecorations();
void FreeIsSubDecoration();
extern int** IsSubDecoration;
#endif // MODS_H
