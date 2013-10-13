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
        "The 6 frame translation works the same as the old run by hand."
        translate = SixFrameFasta.AbacusClass()
        translate.ParseCommandLine(["-r",Test.IN,"-w",Test.OUT,"-c","NC_004837"])
        translate.Main()
        self.assert_(filecmp.cmp(Test.OUT,Test.SIX))
        os.remove(Test.OUT)
        os.remove("Temp.RT.fasta")

if __name__ == "__main__":
    unittest.main()
