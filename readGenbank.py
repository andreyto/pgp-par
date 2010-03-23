#!/usr/bin/env python

import sys, os
from Bio import SeqIO

import PGPeptide
import bioseq


gbFile   = sys.argv[1]
orfFasta = sys.argv[2]
pepGFF   = sys.argv[3]

def checkTimeAndMem(oldtime,cnt):
    t = os.times()[4]
    fd = open("/proc/self/status")
    mem = ''
    for line in fd.xreadlines():
        if line.startswith('VmSize:'):
            mem = line.strip()
            break
    print "Count %d took %d seconds, %s " % (cnt,t-oldtime,mem)
    fd.close()
    return t

endToCDS = {}

count = 0
timecheck = checkTimeAndMem( os.times()[4], count )

chromReader = PGPeptide.GenbankChromosomeReader(gbFile,orfFasta) 
genome = chromReader.locateOrfs()

timecheck = checkTimeAndMem(timecheck,count)

for chrom in genome.chromosomes.values():
    print "Chromosome %s" % chrom.accession
    print "Number of loci with simple ORF end hit %d" % len(chrom.simpleOrfs)
    print "Number of loci as other(complex) ORFs %d" % len(chrom.otherOrfs)
    print "Complex orf loci:"

    for orf in chrom.otherOrfs.values():
        cds = chrom.endToCDS[ orf.location.GetThreePrime() ]
        print cds.qualifiers['locus_tag'][0]

timecheck = checkTimeAndMem(timecheck,count)


gffReader = PGPeptide.GFFPeptide( pepGFF )
gffReader.generateORFs( orfFasta, genome )

for chromName,chrom in genome.chromosomes.items():
    simpleOrfNames = set(chrom.simpleOrfs.keys())
    otherOrfNames  = set(chrom.otherOrfs.keys())

    print "Num orfs from peptides %d, num simpleOrfs %d num otherOrfs %d " % (
        chrom.numORFsWithPeptides(), len(simpleOrfNames), len(otherOrfNames))

    print otherOrfNames

timecheck = checkTimeAndMem(timecheck,count)
