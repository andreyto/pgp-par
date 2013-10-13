
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

#ifndef FREE_MOD_H
#define FREE_MOD_H
#include "Tagger.h"

// FreeMod.h includes code and classes for handling mutations and large sets
// of modifications.  (Most references to "mods" can be taken to mean "mutations and 
// post-translational modifications")  This is a more powerful (but much slower) way to search 
// spectra, and is most appropriate for second-pass searching.  (In multipass
// searching, the database contains only the proteins identified with high
// confidence during a restrictive first-pass search of a large database, like an IPI
// species database or swiss-prot)  

// DELTA_BIN_COUNT is the number of mass bins in the range [MIN_DELTA_AMU, MAX_DELTA_AMU], 400*10 = 4000
// This bin count is the size of each MassDeltaByMass[AA] array.
//#define DELTA_BIN_COUNT 4000

// MASS_TO_BIN and BIN_TO_MASS convert between masses and mass-bins
#define MASS_TO_BIN(mass) (int)((mass + 50) / 100)
#define BIN_TO_MASS(bin) (int)((bin) * 100)

#define MDBI_ALL_MODS 26

// Scaling factor, compensating for the different score ranges of 
// quick PRM-based scoring and final match-scoring.
#define DELTA_SCORE_SCALER 200
#define DELTA_SCORE_SCALER_FINAL 0.5

// Search a database, using *no* tag-based filtering at all.  This is much slower than searching
// with tag-based filters, but also more sensitive, particularly since tagging is harder in the presence of mods.  
void SearchDatabaseTagless(SearchInfo* Info, int MaxMods, int VerboseFlag, SpectrumTweak* Tweak);

// Set Spectrum->PRMScores, using the PRM scoring model.  When extending in blind mode,
// we use the scores of these PRMs as an initial score for our peptides.
void SetPRMScores(MSSpectrum* Spectrum);

// Read, from the binary file Mutations.dat, the definitions of all mass modifications we will consider.
void LoadMassDeltas(char* FileName, int ReadFlag);

// Initialize the hash MassDeltaByMass.  The table entry MassDeltaByMass[AA][Delta] points to a linked list
// of mass deltas for amino acid AA matching Delta.  
void InitMassDeltaByMass();

// Re-score spectral matches.  The matches in the list Spectrum->FirstMatch have been 
// quick-scored, but we can sort them better if we score them more meticulously.  
// Let's do so, and re-sort the list based on the new scores.
void MQScoreSpectralMatches(SpectrumNode* Spectrum);

// Print out a list of matches for the spectrum (Spectrum->FirstMatch through Spectrum->LastMatch).
void DebugPrintMatchList(SpectrumNode* Spectrum);

// Attach edges moving back by one, two, or three amino acid masses to nodes in the TagGraph 
void TagGraphPopulateBackEdges(TagGraph* Graph);

void FreeMassDeltaByMass();
void FreeMassDeltas();
void AddBlindMods();
void AllocMassDeltaByIndex();

#endif // FREE_MOD_H
