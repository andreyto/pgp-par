#!/usr/bin/env python

'''
Created on Oct 8, 2009

@author: eli

Some unit tests for the database setup code.
'''
import unittest
import filecmp
import os

import PGPeptide

class Test(unittest.TestCase):

    def testGenomicLocation(self):
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

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test6Frame']
    unittest.main()
