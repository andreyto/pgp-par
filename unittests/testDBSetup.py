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

    IN = 'NC_004837.fa'
    OUT = 'test.fa'
    SIX = 'NC_004837.6frame.fa'

    def test6Frame(self):
        translate = SixFrameFasta.AbacusClass()
        translate.ParseCommandLine(["-r",Test.IN,"-w",Test.OUT,"-c","Plasmid1"])
        translate.Main()       
        self.assert_(filecmp.cmp(Test.OUT,Test.SIX))
        os.remove(Test.OUT)
        os.remove("Temp.RT.fasta")

    def testFastaReader(self):
        reader = bioseq.FastaReader(Test.IN)
        for seq in reader:
            self.assertEqual('gi|31795327|ref|NC_004837.1|',seq.acc)
            self.assertEqual('TGTAACGAACGGTG',seq.seq[0:14])
            self.assertEqual('TACCCCGACCCCTG',seq.seq[-14:])
            self.assertEqual('Yersinia pestis KIM plasmid pPCP1, complete sequence',seq.desc)

    def testFastaMultiReader(self):
        reader = bioseq.FastaReader(Test.SIX)
        i = 0
        for seq in reader:
            i += 1

        self.assertEqual( 596, i )

    def testFastaOut(self):
        reader = bioseq.FastaReader(Test.IN)
        output = bioseq.FastaOut(Test.OUT)
        output.linesize = 70
        for seq in reader:
            output.write(seq)
            
        output.close()
        self.assert_(filecmp.cmp(Test.OUT,Test.IN))
        os.remove(Test.OUT)
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test6Frame']
    unittest.main()