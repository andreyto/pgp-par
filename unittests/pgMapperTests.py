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
import PGPeptide



class Test(unittest.TestCase):
    

    def SetUpPeptides(self):
        """ reusable code to get peptides and coords
        """
        self.AminoList = ["NSGDSVSIGGDAAGISNK", "ITVETGQEK", "AGITAGYQETR", "IEGLFNDANIGLR", "TEGTVSYEQK", "GTVSYEQK", "TVETGQEK"]
        self.StartList = [7511, 5583, 7106, 5040, 5610, 5610, 5583]
        self.StopList =  [7564, 5609, 7138, 5078, 5639, 5633, 5606]
        self.Databases = ["NC_004837.6frame.trie",]
        self.PeptidesInORF = ["Protein229", "Protein453", "Protein229", "Protein453", "Protein453", "Protein453", "Protein453"]

    def testMappingStartStop(self):
        """Name: testMappingStartStop
        Description:Call out to the PeptideMapper methods to make sure that 
        they are working correctly 

        """
        self.SetUpPeptides()
        Mapper = PeptideMapper.PeptideMappingClass()
        Mapper.LoadDatabases(self.Databases)
        for Index in range(len(self.AminoList)):
            Aminos = self.AminoList[Index]
            Start = self.StartList[Index]
            Stop = self.StopList[Index]
            ShouldMapTo = self.PeptidesInORF[Index]
            (Peptide, ) = Mapper.MapPeptide(Aminos, 0.01) #hack because they are all unique, it should be a one member list
            self.assertEqual(Peptide.location.start, Start)
            self.assertEqual(Peptide.location.stop, Stop)
            self.assertEqual(Peptide.ORFName, ShouldMapTo)
        

    def SetUpProteins(self):
        """ put in some proteins so that I don't have to deal with crap"""
        self.Proteins = {}
        self.ProteinStarts = {}
        self.ProteinStops = {}
        self.ProteinORFLocation = {}
        Name1 = "gi|39980761|ref|NP_951039.1| putative replication regulatory protein [Yersinia pestis KIM]"
        Seq1 = "MNKQQQTALNMARFIRSQSLILLEKLDALDADEQAAMCERLHELAEELQNSIQARFEAESETGT"
        self.Proteins[Name1] = Seq1
        self.ProteinStarts[Name1] = 2925
        self.ProteinStops[Name1] = 3119
        self.ProteinORFLocation[Name1] = "Protein87"
        Name2 = "gi|39980762|ref|NP_951040.1| hypothetical protein YPKp06 [Yersinia pestis KIM]"
        Seq2 = "MKFHFCDLNHSYKNQEGKIRSRKTAPGNIRKKQKGDNVSKTKSGRHRLSKTDKRLLAALVVAGYEERTARDLIQKHVYTLTQADLRHLVSEISNGVGQSQAYDAIYQARRIRLARKYLSGKKPEGVEPREGQEREDLP"
        self.Proteins[Name2] = Seq2
        self.ProteinStarts[Name2] = 6006
        self.ProteinStops[Name2] = 6422
        self.ProteinORFLocation[Name2] = "Protein189"
        Name3 = "gi|39980763|ref|NP_951041.1| putative transcriptional regulator [Yersinia pestis KIM]"
        Seq3 = "MRTLDEVIASRSPESQTRIKEMADEMILEVGLQMMREELQLSQKQVAEAMGISQPAVTKLEQRGNDLKLATLKRYVEAMGGKLSLDVELPTGRRVAFHV"
        self.Proteins[Name3] = Seq3
        self.ProteinStarts[Name3] = 7790
        self.ProteinStops[Name3] = 8089
        self.ProteinORFLocation[Name3] = "Protein352"
     
    def testMappingProteins(self):
        """Name: testMappingProteins
        Description:does the Protein mapping code give me good LocatedProtein objects
        """
        self.SetUpProteins()
        ORFDatabases = ["NC_004837.6frame.trie",]
        Mapper = PeptideMapper.PeptideMappingClass()
        Mapper.LoadDatabases(ORFDatabases)
        for (Name, Seq) in self.Proteins.items():
            LocatedProtein = Mapper.MapProtein(Seq, Name)
            self.assertEqual(LocatedProtein.GetStart(), self.ProteinStarts[Name])
            self.assertEqual(LocatedProtein.GetStop(), self.ProteinStops[Name])
            self.assertEqual(LocatedProtein.GetORFName(), self.ProteinORFLocation[Name])

    def testOpenReadingFrame_GetFivePrime(self):
        """Name: testMappingOpenReadingFrame
        Description: does the OpenReadingFrame get created correctly?
        """
        self.SetUpPeptides()
        ORFDatabases = ["NC_004837.6frame.trie",]
        Mapper = PeptideMapper.PeptideMappingClass()
        Mapper.LoadDatabases(ORFDatabases)
        #Make all the Peptides
        AllPeptides = []
        AllORFs = []
        for Aminos in self.AminoList:
            (Peptide, ) = Mapper.MapPeptide(Aminos, 0.01) #hack because they are all unique, it should be a one member list
            AllPeptides.append(Peptide)
            
        for ID in Mapper.ORFDB.ProteinSequences.keys():
            Fasta = Mapper.ORFDB.ProteinNames[ID]
            Seq = Mapper.ORFDB.ProteinSequences[ID]
            Object = PGPeptide.OpenReadingFrame(Fasta, Seq)
            #now go throug the peptides, and add them to this object
            for LocatedPeptide in AllPeptides:
                if LocatedPeptide.ORFName == Object.name:
                    Object.addLocatedPeptide(LocatedPeptide)
                    
            AllORFs.append(Object)
        #now we are done with making our ORFs, lets' get the 5' peptide for each
        # as a crib sheet 
        for ORF in AllORFs:
            Peptide = ORF.GetFivePrimePeptide()
            if ORF.name == "Protein229":
                self.assertEqual(Peptide.aminos, "AGITAGYQETR")
            if ORF.name == "Protein453":
                self.assertEqual(Peptide.aminos, "TEGTVSYEQK")    
        

if __name__ == "__main__":
    unittest.main()
