#!/usr/bin/env python

"""
Class Test: this serves to test some of the basic functionality of
the mapping peptides onto DNA. Testing that coordinates on the + and
- strand of DNA are correct.  Using my old stuff, the coordinates 
are as follows

loading /export/software/proteogenomics/unittests/NC_004837.6frame.trie
ResultsParser:/export/software/proteogenomics/unittests/inspect2.txt.bz2
(0/1) /export/software/proteogenomics/unittests/inspect2.txt.bz2
I found 31 peptides from 52 spectra
GenomicLocationForPeptide object, NSGDSVSIGGDAAGISNK
  found in protein Protein229 (unique=1)
  Chr:Plasmid1, Start 7511, Stop 7564


GenomicLocationForPeptide object, ITVETGQEK
  found in protein Protein453 (unique=1)
  Chr:Plasmid1, Start 5583, Stop 5609


GenomicLocationForPeptide object, AGITAGYQETR
  found in protein Protein229 (unique=1)
  Chr:Plasmid1, Start 7106, Stop 7138


GenomicLocationForPeptide object, IEGLFNDANIGLR
  found in protein Protein453 (unique=1)
  Chr:Plasmid1, Start 5040, Stop 5078


GenomicLocationForPeptide object, TEGTVSYEQK
  found in protein Protein453 (unique=1)
  Chr:Plasmid1, Start 5610, Stop 5639


GenomicLocationForPeptide object, GTVSYEQK
  found in protein Protein453 (unique=1)
  Chr:Plasmid1, Start 5610, Stop 5633


GenomicLocationForPeptide object, TVETGQEK
  found in protein Protein453 (unique=1)
  Chr:Plasmid1, Start 5583, Stop 5606
"""

import unittest
import PeptideMapper


class Test(unittest.TestCase):
    

    def testMappingStartStop(self):
        """Name: testMappingStartStop
        Description:Call out to the PeptideMapper methods to make sure that 
        they are working correctly 

        """
        #copied from above
        AminoList = ["NSGDSVSIGGDAAGISNK", "ITVETGQEK", "AGITAGYQETR", "IEGLFNDANIGLR", "TEGTVSYEQK", "GTVSYEQK", "TVETGQEK"]
        StartList = [7511, 5583, 7106, 5040, 5610, 5610, 5583]
        StopList =  [7564, 5609, 7138, 5078, 5639, 5633, 5606]
        Databases = ["NC_004837.6frame.trie",]
        Mapper = PeptideMapper.PeptideMappingClass()
        Mapper.LoadDatabases(Databases)
        for Index in range(len(AminoList)):
            Aminos = AminoList[Index]
            Start = StartList[Index]
            Stop = StopList[Index]
            (Peptide, ) = Mapper.MapMe(Aminos, 0.01) #hack because they are all unique, it should be a one member list
            self.assertEqual(Peptide.location.start, Start)
            self.assertEqual(Peptide.location.stop, Stop)
        
        

if __name__ == "__main__":
    unittest.main()
