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
    
    def setUp(self):
        """This is some generic code that I don't feel like writing over again for each test
        """
        #these first four are appropriate for a sequence complexity filter
        self.P1 = GenomicLocations.GenomicLocationForPeptide()
        self.P1.Aminos = "GGGGGGGGGGG"
        self.P1.Unique = 1 #isn't that amazing.  it's unique
        self.P2 = GenomicLocations.GenomicLocationForPeptide()
        self.P2.Aminos = "GGAGAGAGGGGAGAG"
        self.P3 = GenomicLocations.GenomicLocationForPeptide()
        self.P3.Aminos = "AGGAGQWEIUPIUFGBVA"
        self.P3.Unique = 1 # also unique
        self.P4 = GenomicLocations.GenomicLocationForPeptide()
        self.P4.Aminos = "AGGAGGQHGGGAG"
        
        #these next are all tryptic at the c-term
        self.P5 = GenomicLocations.GenomicLocationForPeptide()
        self.P5.Aminos = "MSTAQWSTR"
        

    def testSequenceComplexity_exclusivelyGA(self):
        """Name: testSequenceComplexity_exclusivelyGA
        Description: This method is meant to test whether we remove peptides
        from ORF objects is they are exclusively GA in content. That's it.
        
        NOTE: currently set up with the old GenomicLocations objects
        and not the new ones that Eli is writing. but we need tests now.
        """
        Filter = PGORFFilters.SequenceComplexityFilter()
                
        #now make the encapsulating ORF
        FastaLine = ">Protein0.Chr:NC_001263.Frame1.StartNuc1.Strand+" #needed for the constructor
        AALen = 55 # again, a dummy value just for the constructor
        ORF = GenomicLocations.GenomicLocationForORF(FastaLine, AALen)
        ListOfPeptides = [self.P1, self.P2, self.P3, self.P4] #brackets make a list Sam, parens make a tuple
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
            
    def testMinPeptide(self):
        """Name: testMinPeptide
        Description: this is just to test that our code surrounding the 
        MinPeptide filter works.  You need at least X peptides to pass
        """
        FilterRequire2 = PGORFFilters.MinPeptideFilter(2) #min is two
        FilterRequire1 = PGORFFilters.MinPeptideFilter(1)
        FilterRequire5 = PGORFFilters.MinPeptideFilter(5)
        #make ORFs
        FastaLine = ">Protein0.Chr:NC_001263.Frame1.StartNuc1.Strand+" #needed for the constructor
        AALen = 55 # again, a dummy value just for the constructor
        ListOf1 = [self.P5]
        ListOf2 = [self.P1, self.P4]
        ListOf5 = [self.P1, self.P2, self.P3, self.P4, self.P5] #brackets make a list Sam, parens make a tuple
        ORF0 = GenomicLocations.GenomicLocationForORF(FastaLine, AALen)
        ORF1 = GenomicLocations.GenomicLocationForORF(FastaLine, AALen)
        ORF2 = GenomicLocations.GenomicLocationForORF(FastaLine, AALen)
        ORF5 = GenomicLocations.GenomicLocationForORF(FastaLine, AALen)
        ORF1.PeptideLocationList = ListOf1
        ORF2.PeptideLocationList = ListOf2
        ORF5.PeptideLocationList = ListOf5
        
        #now we run the tests
        self.assertEqual(FilterRequire1.apply(ORF0), 1) # 1 says DELETEME NOW
        self.assertEqual(FilterRequire1.apply(ORF1), 0) 
        self.assertEqual(FilterRequire5.apply(ORF2), 1) 
        self.assertEqual(FilterRequire2.apply(ORF2), 0) 
        self.assertEqual(FilterRequire5.apply(ORF5), 0) 


    def testUnique(self):
        """Name: testUnique
        Description: just tests whether our uniqueness filter is working
        in that it should require a unique peptide.
        """
        Filter = PGORFFilters.UniquenessFilter()
        
        FastaLine = ">Protein0.Chr:NC_001263.Frame1.StartNuc1.Strand+" #needed for the constructor
        AALen = 55 # again, a dummy value just for the constructor
        ListOf5 = [self.P1, self.P2, self.P3, self.P4, self.P5] #brackets make a list Sam, parens make a tuple
        #now this is tricky, but important, if I don't declare something as unique, I assume that it is not
        #so only P1 and P3 are unique.
        ListOfBad = [self.P2, self.P4]
        ORFHasUnique = GenomicLocations.GenomicLocationForORF(FastaLine, AALen)
        ORFHasUnique.PeptideLocationList = ListOf5
        ORFLacksUnique = GenomicLocations.GenomicLocationForORF(FastaLine, AALen)
        ORFLacksUnique.PeptideLocationList = ListOfBad
        
        self.assertEqual(Filter.apply(ORFHasUnique), 0)
        self.assertEqual(Filter.apply(ORFLacksUnique), 1)
        

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test6Frame']
    unittest.main()
