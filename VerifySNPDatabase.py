"""
Verify the integrity of the SNP database:
For each of the SNPs, verify that ONE of the alleles corresponds to the genomic sequence.
Run this after running ParseSNPDatabase.
"""
import os
import sys
import struct

def VerifySNPDatabase(ChromosomeNumber):
    GenomeFile = open("e:\\chromosome\\chr%s.trie"%ChromosomeNumber, "rb")
    # A big read:
    GenomicSequence = GenomeFile.read()
    SNPFile = open("SNP\\%s.snp"%ChromosomeNumber, "rb")
    RecordNumber = 0
    while (1):
        Data = SNPFile.read(4)
        if not Data:
            break # eof
        try:
            GenomePos = struct.unpack("<i", Data)[0]
        except:
            print RecordNumber, Data, len(Data)
            raise
        PolyType = ord(SNPFile.read(1))
        FormA = SNPFile.read(1)
        FormB = SNPFile.read(1)
        GenomicNuc = GenomicSequence[GenomePos].upper()
        if GenomicNuc != FormA and GenomicNuc != FormB:
            print "** Error in record %d pos %d: Reported snp alleles are %s and %s but genomic sequence is %s"%(\
                RecordNumber, GenomePos, FormA, FormB, GenomicNuc)
        RecordNumber += 1
        if RecordNumber % 1000 == 0:
            print "Record %s..."%RecordNumber
for X in range(1, 25):
    print "Verify SNP database for chrom %s..."%X
    VerifySNPDatabase(X)
        