#!/usr/bin/env python

import sys, os
from Bio import SeqIO

import PGPeptide
import bioseq


gbFile   = sys.argv[1]
orfFasta = sys.argv[2]

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
chromDict = chromReader.locateOrfs()

timecheck = checkTimeAndMem(timecheck,count)


for chrom in chromDict.values():
    print "Chromosome %s" % chrom.accession
    print "Number of loci without simple ORF end hit %d" % len(chrom.otherOrfs)
    print "Number of loci with simple ORF end hit %d" % len(chrom.simpleOrfs)
    print "Complex orf loci:"
    for orf in chrom.otherOrfs:
        cds = chrom.endToCDS[ orf.location.GetThreePrime() ]
        print cds.qualifiers['locus_tag'][0]

timecheck = checkTimeAndMem(timecheck,count)
