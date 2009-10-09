#!/usr/bin/env python

'''
Created on Oct 8, 2009

@author: eli

Some unit tests for the database setup code.
'''
import unittest
import filecmp
import os

import SixFrameFasta
import bioseq

class Test(unittest.TestCase):


    def test6Frame(self):
        translate = SixFrameFasta.AbacusClass()
        translate.ParseCommandLine(["-r","NC_004837.fa","-w","test1.fa","-c","Plasmid1"])
        translate.Main()       
        self.assert_(filecmp.cmp("test1.fa","NC_004837.6frame.fa"))
        os.remove("test1.fa")
        os.remove("Temp.RT.fasta")

    def testFastaReader(self):
        reader = bioseq.FastaReader('NC_004837.fa')
        for seq in reader:
            self.assertEqual('gi|31795327|ref|NC_004837.1|',seq.id)
            self.assertEqual('TGTAACGAACGGTG',seq.seq[0:14])
            self.assertEqual('TACCCCGACCCCTG',seq.seq[-14::])
            self.assertEqual('Yersinia pestis KIM plasmid pPCP1, complete sequence',seq.desc)
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test6Frame']
    unittest.main()
