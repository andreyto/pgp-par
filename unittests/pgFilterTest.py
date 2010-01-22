#!/usr/bin/env python

"""
Class FilterTest: this serves to test some of the basic functionality of
the filters in PGORFFilters.  We want to have a test set up for each and 
every filter, and for each definable part of a filter.  This will make 
our lives much easier

- Make a method to test a 'unit' of code.  Name the method intuitively
and comment it's purpose
"""

import unittest

import PGORFFilters
import GenomicLocations

class Test(unittest.TestCase):

    def testSequenceComplexity_exclusivelyGA(self):
        """Name: testSequenceComplexity_exclusivelyGA
        Description: This method is meant to test whether we remove peptides
        from ORF objects is they are exclusively GA in content. That's it.
        
        NOTE: currently set up with the old GenomicLocations objects
        and not the new ones that Eli is writing. but we need tests now.
        """
        Filter = PGORFFilters.SequenceComplexityFilter()
        
        P1 = GenomicLocations.GenomicLocationForPeptide()
        P1.Aminos = "GGGGGGGGGGG"
        P2 = GenomicLocations.GenomicLocationForPeptide()
        P2.Aminos = "GGAGAGAGGGGAGAG"
        P3 = GenomicLocations.GenomicLocationForPeptide()
        P3.Aminos = "AGGAGQWEIUPIUFGBVA"
        P4 = GenomicLocations.GenomicLocationForPeptide()
        P4.Aminos = "AGGAGGQHGGGAG"
        
        #now make the encapsulating ORF
        FastaLine = ">Protein0.Chr:NC_001263.Frame1.StartNuc1.Strand+" #needed for the constructor
        AALen = 55 # again, a dummy value just for the constructor
        ORF = GenomicLocations.GenomicLocationForORF(FastaLine, AALen)
        ListOfPeptides = [P1, P2, P3, P4] #brackets make a list Sam, parens make a tuple
        ORF.PeptideLocationList = ListOfPeptides
        
        #now run the filter
        Filter.RemoveExclusivelyLowMWPeptides(ORF)
        #after running this filter, the one thing we care about for this unittest
        #is that P1 and P2 are deleted.  So let's test the length
        self.assertEqual(len(ORF.PeptideLocationList), 2)
        #now because it's easy, we want to test the identity of the peptides left, just
        #to make sure
        AcceptablePeptides = ["AGGAGQWEIUPIUFGBVA", "AGGAGGQHGGGAG"]
        for Peptide in ORF.PeptideLocationList:
            self.assert_(Peptide.Aminos in AcceptablePeptides)
        

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test6Frame']
    unittest.main()
