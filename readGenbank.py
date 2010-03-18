#!/usr/bin/env python

import sys, os
from Bio import SeqIO

import PGPeptide
import bioseq


gbFile   = sys.argv[1]
orfFasta = sys.argv[2]

def checkTimeAndMem(oldtime,cnt):
    t = os.times()[4]
    fd = open("/proc/%d/status"%os.getpid())
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

nucs = ''
acc  = ''
for gb_rec in SeqIO.parse(open(gbFile), 'genbank'):
    nucs = gb_rec.seq
    acc  = gb_rec.name
    count += 1
    if count % 1000 == 0:
        timecheck = checkTimeAndMem(timecheck,count)

    for feat in gb_rec.features:
        if feat.type == 'CDS':
            if feat.strand == 1:
                endToCDS[ feat.location.end.position ] = feat
            else:
                # biopython 1.53 seems to use 0, or space based coords
                # so start is 1 less then what's in the genbank file
                endToCDS[ feat.location.start.position + 1 ] = feat

timecheck = checkTimeAndMem(timecheck,count)

hitByOrf = set()
orfReader = bioseq.SequenceIO(orfFasta)
for orfSeq in orfReader:
    tmpOrf    = PGPeptide.OpenReadingFrame(orfSeq.acc, orfSeq.seq)
    if tmpOrf.name.startswith('XXX'): # Skip the decoy ORFs
        continue
    if acc != tmpOrf.chromosome:
      raise ValueError("Wrong chromosome have %s and %s"%(acc,tmpOrf.chromosome))
    
    orfThreePrime = tmpOrf.location.GetThreePrime()

    if endToCDS.has_key( orfThreePrime ):
        f = endToCDS[ orfThreePrime ]
        hitByOrf.add( f.qualifiers['locus_tag'][0] )
        print f
        print orfSeq.seq
        print "Orf begin %d end %d" % (tmpOrf.location.start, tmpOrf.location.stop )
        print f.extract( nucs )
    else:
        print "No CDS end found at %d" % tmpOrf.location.stop

timecheck = checkTimeAndMem(timecheck,count)

allLocus = set([c.qualifiers['locus_tag'][0] for c in endToCDS.values()])
notHit = allLocus - hitByOrf
print "len(all) %d len(hits) %d" % (len(allLocus), len(hitByOrf))
print "Number of loci without an ORF end hit %d" % len(notHit)
for locus in notHit:
    print locus

timecheck = checkTimeAndMem(timecheck,count)
