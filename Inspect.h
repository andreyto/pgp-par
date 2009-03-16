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


#ifndef INSPECT_H
#define INSPECT_H

#define INSPECT_VERSION_NUMBER "20080404"

#include <stdio.h>
#include "Utils.h"

//////////////////////////////////////////////////////////////////////////////
// General-purpose #definitions.

// Multiplier for scaling floating-point masses up to integers.
// Represent 123.456Da as the integer 123456:
#define MASS_SCALE 1000
#define DALTON 1000
#define HALF_DALTON 500
#define DECI_DALTON 100

// Mass (in amu) of hydrogen.  (Used for, e.g., finding the PRM of a b peak)
//#define HYDROGEN_MASS (float)1.0078
#define HYDROGEN_MASS (int)1008
#define TWO_HYDROGEN_MASS (int)2016
#define GLYCINE_MASS 57000
#define WATER_MASS 18000
#define CAM_MASS 57000
#define PHOSPHATE_MASS 79966
#define PHOSPHATE_WATER_MASS 97966

// The parent mass boost is equal to the difference in mass between a precursor ion and
// the parent *residue* mass (sum of amino acid masses).
#define PARENT_MASS_BOOST (int)19000

// Maximum length of a peptide tag that can be indexed.  (Char-arrays are limited to this size)
// (Somewhat absurdly large, because we may want to use the trie to search for peptides 
// outside the MS/MS context)
#define MAX_TAG_LENGTH 50

// How far can peaks be from theoretical prediction?
#define DEFAULT_EPSILON 500

#define DEFAULT_PARENT_MASS_EPSILON 2500
#define DEFAULT_PARENT_MASS_PPM 2000
#define DEFAULT_FLANKING_MASS_EPSILON 3000

// How large are the tags we generate, by default?  (Overridable by -x, -y)
#define DEFAULT_TAG_LENGTH 3

// Maximum number of post-translational modification types:
#define MAX_PT_MODTYPE 42

// Maximum number of post-tranlational mods that can EVER be allowed:
#define MAX_PT_MODS 8

// How many entries in the match hash-table. (If we have many more matches than this,
// performance will be slowed a bit)
#define MATCH_HASH_SIZE 1000

// Trie node's child array has one entry for each letter (some slots are wasted, 
// since there are only 20 peptides, but it makes for fast searching)
// The index into the array is the peptide char minus 'A'.  (alanine 0, cysteine 2, etc)
#define TRIE_CHILD_COUNT 26 

#define FILENAME_AMINO_ACID_MASSES "AminoAcidMasses.txt"
#define FILENAME_PTMS "PTMods.txt"
#define FILENAME_MASTER_TAGGING_MODEL "PRM.dat"
#define FILENAME_MASS_DELTAS "Mutations.dat"
#define FILENAME_PVALUE "PValue.dat"
#define FILENAME_PVALUE_TRYPTIC "PValueTryptic.dat"
#define FILENAME_SCORING_MODEL "ScoringModel.dat"
#define FILENAME_ISOTOPE_PATTERNS "IsotopePatterns.txt"
#define FILENAME_INTENSITY_RANK_ODDS "IntensityRankIonOdds.txt"
#define FILENAME_WITNESS_SCORES "IonWitnessScores.dat"
#define FILENAME_PRM_MODEL "PRMModel.dat"

#define TWEAK_COUNT 6

//used as switches for fragmentation models
#define FRAGMENTATION_NORMAL 0
#define FRAGMENTATION_PHOSPHO 1

// We may try two or three different charge/parent-mass combinations for one 
// spectrum.  We use SVMs to determine parent mass and charge state, but in
// borderline cases, we try both.  
typedef struct SpectrumTweak
{
    int ParentMass;
    int Charge;
    // Intensities(S, L) is the frequency of intensity level L in sector S
    float Intensities[12]; // SECTOR_COUNT
    int* PRMScores;
    int PRMScoreMax;
} SpectrumTweak;

#define SPECTRUM_FORMAT_INVALID -1
#define SPECTRUM_FORMAT_DTA 0
#define SPECTRUM_FORMAT_PKL 1
#define SPECTRUM_FORMAT_MS2 2
#define SPECTRUM_FORMAT_MGF 3
#define SPECTRUM_FORMAT_MS2_COLONS 4
#define SPECTRUM_FORMAT_MZXML 5
#define SPECTRUM_FORMAT_MZDATA 6

// Create one InputFileNode for each file being searched.  
// If the input file is a standard .dta file, then we create one child SpectrumNode.
// If the input file is a .ms2 file, then we create many child SpectrumNodes.
typedef struct InputFileNode
{
    char FileName[MAX_FILENAME_LEN];
    int SpectrumCount;
    int Format; // 0 dta, 1 pkl, 2 ms2, 3 mgf
    struct InputFileNode* Prev;
    struct InputFileNode* Next;
} InputFileNode;

typedef struct SpectrumNode
{
    struct MSSpectrum* Spectrum;
    struct SpectrumNode* Next;
    SpectrumTweak Tweaks[TWEAK_COUNT];
    int PMCFlag; // Set to 1 after PMC is done and our tweak-array is populated.
    int FilePosition; // seek to here before parsing 
    int ScanNumber;
    int MatchCount;
    struct Peptide* FirstMatch;
    struct Peptide* LastMatch;
    InputFileNode* InputFile; // the file name (and file type)
} SpectrumNode;

// The Stats object is for keeping track of cumulative info (tags generated,
// bytes read, spectra scored, that sort of thing)
typedef struct InspectStats
{
    // Tags generated - raw count of all tripeptide paths through the PRM graph
    long long TagsGenerated;
    // Tag hits in the database (How many tripeptide tag matches were extended?)
    long long TagMatches;
    // Number of candidate peptides that were scored against the source spectrum
    long long CandidatesScored;
    long long TagGraphNodes;
    long long TagGraphEdges;
} InspectStats;

extern InspectStats* GlobalStats;

typedef void (*TrainingCallback)(SpectrumNode*, int, int, struct Peptide*);
void TrainOnOracleFile(char* OracleFileName, char* SpectrumDir, TrainingCallback Callback);
void AddSpectrumToList(InputFileNode* InputFile, int FilePos, int ScanNumber);

#endif // INSPECT_H
