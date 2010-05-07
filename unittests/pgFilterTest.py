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

        self.fastaLine = ">Protein0.Chr:NC_001263.Frame1.StartNuc1.Strand+" #needed for the construct or 

        self.loc = PGPeptide.GenomicLocation(0,0,'+','NC_001263') #total hack because in this unit test we are NOT testing location.
        self.orfSeq = "THISISADUMMYSEQUENCEHOPEITWORKSOK"
        self.P1 = PGPeptide.LocatedPeptide("GGGGGGGGGGG",self.loc)
        self.P1.isUnique = 1 #isn't that amazing.  it's unique
        self.P2 = PGPeptide.LocatedPeptide("GGAGAGAGGGGAGAG",self.loc)
        self.P3 = PGPeptide.LocatedPeptide("AGGAGQWEIUPIUFGBVA",self.loc)
        self.P3.isUnique = 1 # also unique
        self.P4 = PGPeptide.LocatedPeptide("AGGAGGQHGGGAG",self.loc)

        #these next are all tryptic at the c-term
        self.P5 = PGPeptide.LocatedPeptide("MSTAQWSTR",self.loc)
        self.P6 = PGPeptide.LocatedPeptide("WEROTYSDSWH", self.loc) #NOT tryptic
        self.P5.SetTryptic("K") #making this fully tryptic
        self.P6.SetTryptic("M") #making this not tryptic

    def testSequenceComplexity_exclusivelyGA(self):
        """Name: testSequenceComplexity_exclusivelyGA
        Description: This method is meant to test whether we remove peptides
        from ORF objects is they are exclusively GA in content. That's it.
        
        NOTE: currently set up with the old GenomicLocations objects
        and not the new ones that Eli is writing. but we need tests now.
        """
        Filter = PGORFFilters.SequenceComplexityFilter()
                
        #now make the encapsulating ORF
        ORF = PGPeptide.OpenReadingFrame(self.fastaLine, self.orfSeq)
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
        ORF0 = PGPeptide.OpenReadingFrame(self.fastaLine, self.orfSeq) 
        ORF1 = PGPeptide.OpenReadingFrame(self.fastaLine, self.orfSeq)
        ORF2 = PGPeptide.OpenReadingFrame(self.fastaLine, self.orfSeq)
        ORF5 = PGPeptide.OpenReadingFrame(self.fastaLine, self.orfSeq)
        ORF1.addLocatedPeptides( ListOf1 )
        ORF2.addLocatedPeptides( ListOf2 )
        ORF5.addLocatedPeptides( ListOf5 )
        
        #now we run the tests
        self.assertEqual(FilterRequire1.apply(ORF0), 1) # 1 says DELETEME NOW
        self.assertEqual(FilterRequire1.apply(ORF1), 0) 
        self.assertEqual(FilterRequire5.apply(ORF2), 1) 
        self.assertEqual(FilterRequire2.apply(ORF2), 0) 
        self.assertEqual(FilterRequire5.apply(ORF5), 0) 

        # Also test the FilterList
        filters = PGORFFilters.FilterList([FilterRequire5])
        keep = filters.ApplyAllFilters( { 'ORF5': ORF5, 'ORF2': ORF2 })
        self.assertEqual( 1, len(keep) )
        self.assertEqual( 'ORF5', keep.keys()[0] )

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
        ORFHasUnique = PGPeptide.OpenReadingFrame(self.fastaLine, self.orfSeq)
        ORFHasUnique.addLocatedPeptides( ListOf5 )
        ORFLacksUnique = PGPeptide.OpenReadingFrame(self.fastaLine, self.orfSeq)
        ORFLacksUnique.addLocatedPeptides( ListOfBad )
        
        self.assertEqual(Filter.apply(ORFHasUnique), 0)
        self.assertEqual(Filter.apply(ORFLacksUnique), 1)

    def testTrypticFilter(self):
        "Name: test the tryptic filters of open reading frames"
        Filter = PGORFFilters.TrypticFilter()
        #set tryptic on peptides done in setUp()
        ListHasTryptic = [self.P1, self.P5]
        ListLacksTryptic = [self.P1, self.P6]
        ListUnset = [self.P1, self.P2]
        ORFHasTryptic = PGPeptide.OpenReadingFrame(self.fastaLine, self.orfSeq)
        ORFHasTryptic.addLocatedPeptides( ListHasTryptic)
        ORFLacksTryptic = PGPeptide.OpenReadingFrame(self.fastaLine, self.orfSeq)
        ORFLacksTryptic.addLocatedPeptides( ListLacksTryptic )
        ORFUnset = PGPeptide.OpenReadingFrame(self.fastaLine, self.orfSeq)
        ORFUnset.addLocatedPeptides( ListUnset )
       
        self.assertEqual(Filter.apply(ORFHasTryptic), 0)
        self.assertEqual(Filter.apply(ORFLacksTryptic), 1)
        self.assertEqual(Filter.apply(ORFUnset), 1)
        
        
        

if __name__ == "__main__":
    unittest.main()
