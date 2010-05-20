#!/usr/bin/env python

'''
Created on Oct 8, 2009

@author: eli

Some unit tests for the database setup code.
'''
import unittest
import filecmp
import os

import bioseq

class Test(unittest.TestCase):

    IN = 'NC_004837.fa'
    OUT = 'test.fa'
    SIX = 'NC_004837.6frame.fa'
    TRIE = 'NC_004837.6frame.trie'

    def testTrieReader(self):
        "The trie sequence is correctly parsed into acc, seq and desc."
        reader = bioseq.SequenceIO(Test.TRIE)
        i=0
        for seq in reader:
            i+=1
            if i == 1:
                self.assertEqual('Protein0.Chr:NC_004837.Frame1.StartNuc1.Strand+',seq.acc)
                self.assertEqual('CNERCNSDPHPTPEIRSRG',seq.seq)
            elif i == 2:
                self.assertEqual('Protein1.Chr:NC_004837.Frame2.StartNuc2.Strand+',seq.acc)
                self.assertEqual('VTNGAIVIHTQRLKSDPGGNLLS',seq.seq)
            elif i == 596:
                self.assertEqual('Protein595.Chr:NC_004837.Frame3.StartNuc74.Strand-',seq.acc)
                self.assertEqual('IRRADYPLDLISGVGCGSLLHRSL',seq.seq)

    def testTrieIndexReader(self):
        "The trie is indexed and finds peptides."
        trieIndex = bioseq.TrieIndexSeqs(Test.TRIE)
        trieIndex.index()
        acc,index = trieIndex.accessionWherePeptideFound('TPEIRSR')
        self.assertEqual('Protein0.Chr:NC_004837.Frame1.StartNuc1.Strand+',acc)
        self.assertEqual('CNERCNSDPHPTPEIRSRG',trieIndex.seqs[index])

        acc,index = trieIndex.accessionWherePeptideFound('VTNGAIV')
        self.assertEqual('Protein1.Chr:NC_004837.Frame2.StartNuc2.Strand+',acc)
        self.assertEqual('VTNGAIVIHTQRLKSDPGGNLLS',trieIndex.seqs[index])

        acc,index = trieIndex.accessionWherePeptideFound('RADYPLD')
        self.assertEqual('Protein595.Chr:NC_004837.Frame3.StartNuc74.Strand-',acc)
        self.assertEqual('IRRADYPLDLISGVGCGSLLHRSL',trieIndex.seqs[index])

    def testFastaReader(self):
        "Fasta sequence correctly parsed into acc, seq and desc."
        reader = bioseq.SequenceIO(Test.IN)
        for seq in reader:
            self.assertEqual('gi|31795327|ref|NC_004837.1|',seq.acc)
            self.assertEqual('TGTAACGAACGGTG',seq.seq[0:14])
            self.assertEqual('TACCCCGACCCCTG',seq.seq[-14:])
            self.assertEqual('Yersinia pestis KIM plasmid pPCP1, complete sequence',seq.desc)

    def testFastaMultiReader(self):
        "All of the sequences in the multiFasta file are read."
        reader = bioseq.SequenceIO(Test.SIX)
        i = 0
        for seq in reader:
            i += 1

        self.assertEqual( 596, i )

    def testFastaOut(self):
        "No changes on the round trip through the fasta reader and writer."
        reader = bioseq.SequenceIO(Test.IN)
        output = bioseq.SequenceIO(Test.OUT,'w')
        output.set('linesize',70)
        for seq in reader:
            output.write(seq)

        output.close()
        self.assert_(filecmp.cmp(Test.OUT,Test.IN))
        os.remove(Test.OUT)

if __name__ == "__main__":
    unittest.main()
