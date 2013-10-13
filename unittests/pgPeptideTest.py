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

    GFF = 'membrane.NC_004837.gff'
    SixFrame = 'NC_004837.6frame.fa'

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

        genome    = PGPeptide.Genome()
        gffReader = PGPeptide.GFFPeptide( Test.GFF )
        acc = 'NC_004837'
        chrom = genome.makeChromosome(acc)
        gffReader.generateORFs( Test.SixFrame, genome )

        self.assertEqual( 2, len(chrom.pepOnlyOrfs) )
        orf229 = chrom.getOrf('Protein229')
        orf453 = chrom.getOrf('Protein453')
        self.assertEqual( 68, orf229.numPeptides())
        self.assertEqual( 36, orf453.numPeptides())

        i = 0
        for peptide in orf453.peptideIter():
            i += 1
            if i == 1:
                self.assertEqual( peptide.aminos, 'PNFEANSITIALPHYVDLPGR' )
                self.assertEqual( peptide.GetStart(), 5745 )
                self.assertEqual( peptide.GetStop(), 5807 )
                self.assertEqual( peptide.Strand(), '-' )
                self.assertEqual( peptide.bestScore, 0.00564015792442 )
            elif i == 36:
                self.assertEqual( peptide.aminos, 'TVETGQEKDGVK' )
                self.assertEqual( peptide.GetStart(), 5571 )
                self.assertEqual( peptide.GetStop(), 5606 )
                self.assertEqual( peptide.Strand(), '-' )
                self.assertEqual( peptide.bestScore, 0.00465920088483 )

        i = 0
        for peptide in orf229.peptideIter():
            i += 1
            if i == 1:
                self.assertEqual( peptide.aminos, 'YYGTVINAGYYVTPNAK' )
                self.assertEqual( peptide.GetStart(), 7394 )
                self.assertEqual( peptide.GetStop(), 7444 )
                self.assertEqual( peptide.Strand(), '+' )
                self.assertEqual( peptide.bestScore, 0.00465920088483 )
            elif i == 68:
                self.assertEqual( peptide.aminos, 'YYGTVINAGY' )
                self.assertEqual( peptide.GetStart(), 7394 )
                self.assertEqual( peptide.GetStop(), 7423 )
                self.assertEqual( peptide.Strand(), '+' )
                self.assertEqual( peptide.bestScore, 0.016171175805 )

    def testRoundTripGFF(self):
        'Reading then writing ORF peptides from a GFF file.'

        acc = 'NC_004837'
        genome1   = PGPeptide.Genome()
        chrom     = genome1.makeChromosome(acc)
        gffReader = PGPeptide.GFFPeptide( Test.GFF )
        gffReader.generateORFs( Test.SixFrame, genome1 )

        gffOut = 't.gff'
        gffWriter = PGPeptide.GFFPeptide( gffOut, 'w' )
        chrom1Peps = chrom.numORFsWithPeptides()
        for orf in chrom.pepOnlyOrfs.values():
            gffWriter.writeORFPeptides( orf )

        gffWriter.close()
        genome2   = PGPeptide.Genome()
        chrom     = genome2.makeChromosome(acc)
        gffReader = PGPeptide.GFFPeptide( gffOut )
        gffReader.generateORFs( Test.SixFrame, genome2 )
        for name,orf2 in chrom.pepOnlyOrfs.items():
            orf1 = genome1.chromosomes[acc].getOrf( name )
            self.assertEqual( orf1.name, orf2.name )
            self.assertEqual( orf1.chromosome, orf2.chromosome )
            self.assertEqual( orf1.numPeptides(), orf2.numPeptides() )

            for peps in zip(orf1.peptideIter(), orf2.peptideIter()):
                self.assertEqual( peps[0].aminos,     peps[1].aminos)
                self.assertEqual( peps[0].bestScore , peps[1].bestScore)
                self.assertEqual( peps[0].GetStart(), peps[1].GetStart())
                self.assertEqual( peps[0].GetStop(),  peps[1].GetStop())
                self.assertEqual( peps[0].Strand(),   peps[1].Strand())
                self.assertEqual( peps[0].name,       peps[1].name)

        self.assertEqual( chrom1Peps, chrom.numORFsWithPeptides())

        os.remove( gffOut )

    def testChromsomeGBInput(self):
        'Reading and mapping ORFs to a genbank chromosome'

        gbkFile  = '../PGP.Regression.Test/NC_004088.gbk'
        sixFrame = '../PGP.Regression.Test/NC_004088.6frame.RS.fasta'
        chromReader = PGPeptide.GenbankGenomeReader(gbkFile,sixFrame)
        genome = chromReader.makeGenomeWithProteinORFs()
        self.assertEqual( 1, genome.numChromosomes())
        acc = 'NC_004088'
        self.assertEqual( acc, genome.chromosomes.keys()[0] )
        self.assertEqual( 4086, len( genome.chromosomes[acc].simpleOrfs ))
        self.assertEqual( 0, genome.chromosomes[acc].numORFsWithPeptides())
        
    def testProteinLocation(self):
        'testProteinLocation: reading and mapping of protein locations from gbk onto 6frame'
        gbkFile = "NC_004837.gbk"
        sixFrame = "NC_004837.6frame.trie"
        Accession = "NC_004837"
        chromReader = PGPeptide.GenbankGenomeReader(gbkFile, sixFrame)
        genome = chromReader.makeGenomeWithProteinORFs()
        self.assertEqual( 1, genome.numChromosomes())
        #now we do the hard stuff.  Let's check exact locations of proteins
        self.SetUpProteinLocs()
        SimpleORFs = genome.chromosomes[Accession].simpleOrfs
        for (Name, Seq) in self.Proteins.items():
            ORFName = self.ProteinORFLocation[Name]
            ORF = SimpleORFs[ORFName]# a PGPeptide.OpenReadingFrame object
            LocatedProtein = ORF.annotatedProtein
            self.assertEqual(LocatedProtein.GetStart(), self.ProteinStarts[Name])
            self.assertEqual(LocatedProtein.GetStop(), self.ProteinStops[Name])
            self.assertEqual(LocatedProtein.GetORFName(), self.ProteinORFLocation[Name])
            
            


    def SetUpProteinLocs(self):
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
        self.ProteinORFLocation[Name1] = "NC_004837.Protein87"
        Name2 = "gi|39980762|ref|NP_951040.1| hypothetical protein YPKp06 [Yersinia pestis KIM]"
        Seq2 = "MKFHFCDLNHSYKNQEGKIRSRKTAPGNIRKKQKGDNVSKTKSGRHRLSKTDKRLLAALVVAGYEERTARDLIQKHVYTLTQADLRHLVSEISNGVGQSQAYDAIYQARRIRLARKYLSGKKPEGVEPREGQEREDLP"
        self.Proteins[Name2] = Seq2
        self.ProteinStarts[Name2] = 6006
        self.ProteinStops[Name2] = 6422
        self.ProteinORFLocation[Name2] = "NC_004837.Protein189"
        Name3 = "gi|39980763|ref|NP_951041.1| putative transcriptional regulator [Yersinia pestis KIM]"
        Seq3 = "MRTLDEVIASRSPESQTRIKEMADEMILEVGLQMMREELQLSQKQVAEAMGISQPAVTKLEQRGNDLKLATLKRYVEAMGGKLSLDVELPTGRRVAFHV"
        self.Proteins[Name3] = Seq3
        self.ProteinStarts[Name3] = 7790
        self.ProteinStops[Name3] = 8089
        self.ProteinORFLocation[Name3] = "NC_004837.Protein352"


if __name__ == "__main__":
    unittest.main()
