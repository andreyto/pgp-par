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
import PGPeptide

class Test(unittest.TestCase):
    
    def setUp(self):
        """This is some generic code that I don't feel like writing over again for each test
        """
        #these first four are appropriate for a sequence complexity filter

        self.P1 = PGPeptide.LocatedPeptide()
        self.P1.aminos = "GGGGGGGGGGG"
        self.P1.isUnique = 1 #isn't that amazing.  it's unique
        self.P2 = PGPeptide.LocatedPeptide()
        self.P2.aminos = "GGAGAGAGGGGAGAG"
        self.P3 = PGPeptide.LocatedPeptide()
        self.P3.aminos = "AGGAGQWEIUPIUFGBVA"
        self.P3.isUnique = 1 # also unique
        self.P4 = PGPeptide.LocatedPeptide()
        self.P4.aminos = "AGGAGGQHGGGAG"
        
        #these next are all tryptic at the c-term
        self.P5 = PGPeptide.LocatedPeptide()
        self.P5.aminos = "MSTAQWSTR"
        

    def testSequenceComplexity_exclusivelyGA(self):
        """Name: testSequenceComplexity_exclusivelyGA
        Description: This method is meant to test whether we remove peptides
        from ORF objects is they are exclusively GA in content. That's it.
        
        NOTE: currently set up with the old GenomicLocations objects
        and not the new ones that Eli is writing. but we need tests now.
        """
        Filter = PGORFFilters.SequenceComplexityFilter()
                
        #now make the encapsulating ORF
        ORF = PGPeptide.OpenReadingFrame()
        ORF.addLocatedPeptides([self.P1, self.P2, self.P3, self.P4])
        
        #now run the filter
        ORF.filterPeptides( Filter.lowComplexFilter )
        #after running this filter, the one thing we care about for this unittest
        #is that P1 and P2 are deleted.  So let's test the length
        self.assertEqual(ORF.numPeptides(), 2)
        #now because it's easy, we want to test the identity of the peptides left, just
        #to make sure
        AcceptablePeptides = ["AGGAGQWEIUPIUFGBVA", "AGGAGGQHGGGAG"]
        for Peptide in ORF.peptideIter():
            self.assert_(Peptide.aminos in AcceptablePeptides)
            
    def testMinPeptide(self):
        """Name: testMinPeptide
        Description: this is just to test that our code surrounding the 
        MinPeptide filter works.  You need at least X peptides to pass
        """
        FilterRequire2 = PGORFFilters.MinPeptideFilter(2) #min is two
        FilterRequire1 = PGORFFilters.MinPeptideFilter(1)
        FilterRequire5 = PGORFFilters.MinPeptideFilter(5)
        #make ORFs
        ListOf1 = [self.P5]
        ListOf2 = [self.P1, self.P4]
        ListOf5 = [self.P1, self.P2, self.P3, self.P4, self.P5] #brackets make a list Sam, parens make a tuple
        ORF0 = PGPeptide.OpenReadingFrame() 
        ORF1 = PGPeptide.OpenReadingFrame()
        ORF2 = PGPeptide.OpenReadingFrame()
        ORF5 = PGPeptide.OpenReadingFrame()
        ORF1.addLocatedPeptides( ListOf1 )
        ORF2.addLocatedPeptides( ListOf2 )
        ORF5.addLocatedPeptides( ListOf5 )
        
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
        
        ListOf5 = [self.P1, self.P2, self.P3, self.P4, self.P5] #brackets make a list Sam, parens make a tuple
        #now this is tricky, but important, if I don't declare something as unique, I assume that it is not
        #so only P1 and P3 are unique.
        ListOfBad = [self.P2, self.P4]
        ORFHasUnique = PGPeptide.OpenReadingFrame()
        ORFHasUnique.addLocatedPeptides( ListOf5 )
        ORFLacksUnique = PGPeptide.OpenReadingFrame()
        ORFLacksUnique.addLocatedPeptides( ListOfBad )
        
        self.assertEqual(Filter.apply(ORFHasUnique), 0)
        self.assertEqual(Filter.apply(ORFLacksUnique), 1)
        

if __name__ == "__main__":
    unittest.main()
