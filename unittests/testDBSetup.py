'''
Created on Oct 8, 2009

@author: eli
'''
import unittest
import filecmp
import os

import SixFrameFasta

class Test(unittest.TestCase):


    def test6Frame(self):
        translate = SixFrameFasta.AbacusClass()
        translate.ParseCommandLine(["-r","NC_004837.fa","-w","test1.fa","-c","Plasmid1"])
        translate.Main()       
        self.assert_(filecmp.cmp("test1.fa","NC_004837.6frame.fa"))
        os.remove("test1.fa")
        os.remove("Temp.RT.fasta")
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test6Frame']
    unittest.main()