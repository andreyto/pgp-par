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

timecheck = checkTimeAndMem( os.times()[4], 0 )

chromReader = PGPeptide.GenbankGenomeReader(gbFile,orfFasta) 
genome = chromReader.makeGenomeWithProteinORFs()

timecheck = checkTimeAndMem(timecheck,1)

for chrom in genome.chromosomes.values():
    print "Chromosome %s" % chrom.accession
    print "Number of loci with simple ORF end hit %d" % len(chrom.simpleOrfs)
    print "Number of loci as complex ORFs %d" % len(chrom.complexOrfs)
    print "Complex orf loci:"

    for orf in chrom.complexOrfs.values():
        print orf.CDS.qualifiers['locus_tag'][0]

timecheck = checkTimeAndMem(timecheck,2)


gffReader = PGPeptide.GFFPeptide( pepGFF )
gffReader.generateORFs( orfFasta, genome )

for chromName,chrom in genome.chromosomes.items():
    pepOrfNames  = chrom.pepOnlyOrfs.keys()

    print "Chrom %s num orfs from peptides %d, num simpleOrfs %d num pepOnlyOrfs %d " % (
        chromName,chrom.numORFsWithPeptides(), len(chrom.simpleOrfs), len(pepOrfNames))

    print pepOrfNames

timecheck = checkTimeAndMem(timecheck,3)
