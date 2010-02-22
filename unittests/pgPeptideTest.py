#!/usr/bin/env python

'''
Created on Oct 8, 2009

@author: eli

Some unit tests for the Peptide location code
'''
import unittest
import filecmp
import os

import PGPeptide

class Test(unittest.TestCase):

    def testGenomicLocationOverlap(self):
        "GenomicLocation.overlap() correctly identifies several overlaps."
        p1 = PGPeptide.GenomicLocation(2,11,'+')
        p2 = PGPeptide.GenomicLocation(3,9,'+')
        p3 = PGPeptide.GenomicLocation(5,17,'-')
        p4 = PGPeptide.GenomicLocation(21,367,'-')
        self.assertEqual( 2, p1.start )
        self.assertEqual( 11, p1.stop )

        res = p1.overlap( p2 )
        self.assertEqual( 3, res[0] )
        self.assertEqual( 9, res[1] )

        res = p1.overlap( p3 )
        self.assertEqual( 5, res[0] )
        self.assertEqual( 11, res[1] )
        
        res = p3.overlap( p1 )
        self.assertEqual( 5, res[0] )
        self.assertEqual( 11, res[1] )

        res = p4.overlap( p3 )
        self.assertEqual( None, res )
        
    def testFrameAccessors(self):
        'Verify correct behaviour of frame accessors'
        loc = PGPeptide.GenomicLocation(1,10,'+')
        loc.frame = 2
        self.assertEqual( 2, loc.frame )
        def raiseerror():
            loc.frame = 5
        self.assertRaises(ValueError, raiseerror)

    def testReadGFF(self):
        'Reading ORF peptides from a GFF file.'
        GFF = 'membrane.NC_004837.gff'
        gffReader = PGPeptide.GFF( GFF )
        orf = PGPeptide.OpenReadingFrame(
                'Protein256370.Chr:Chr1.Frame1.StartNuc831372.Strand-',
                'ABOGUSTESTORF')
        gffReader.populateORF( orf )
        i = 0
        for peptide in orf.peptideIter():
            i += 1
            if i == 1:
                self.assertEqual( peptide.aminos, 'PNFEANSITIALPHYVDLPGR' )
                self.assertEqual( peptide.GetStart(), 5745 )
                self.assertEqual( peptide.GetStop(), 5807 )
                self.assertEqual( peptide.Strand(), '-' )
            elif i == 104:
                self.assertEqual( peptide.aminos, 'YYGTVINAGY' )
                self.assertEqual( peptide.GetStart(), 7394 )
                self.assertEqual( peptide.GetStop(), 7423 )
                self.assertEqual( peptide.Strand(), '+' )

        self.assertEqual( 104, i )
        self.assertEqual( orf.numPeptides(), i )

if __name__ == "__main__":
    unittest.main()
