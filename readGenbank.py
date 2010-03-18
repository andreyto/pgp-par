#!/usr/bin/env python

import sys, os
from bisect import bisect_left
from Bio import SeqIO

import PGPeptide
import bioseq


gbFile   = sys.argv[1]
orfFasta = sys.argv[2]

def binary_search(a, x, lo=0, hi=None):   # can't use a to specify default for hi
    hi = hi if hi is not None else len(a) # hi defaults to len(a)   
    pos = bisect_left(a,x,lo,hi)          # find insertion position
    return (pos if pos != hi and a[pos] == x else -1) # don't walk off the end

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

#cds  = []
#ends = []
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
                #ends.append( feat.location.end.position )
                endToCDS[ feat.location.end.position ] = feat
            else:
                #ends.append( feat.location.start.position + 1 )
                # biopython uses 0 based coords
                endToCDS[ feat.location.start.position + 1 ] = feat
            #cds.append( feat )

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
 #   orfPos = binary_search( ends, orfThreePrime )

    if tmpOrf.name == 'Protein115787':
        print "Locus y3406 Orf begin %d end %d 3' %d" % (tmpOrf.location.start, tmpOrf.location.stop, orfThreePrime )
 #       pos = bisect_left(ends, orfThreePrime)
 #       print "OrfPos %d, end-1 %d end %d end+1 %d" % (pos,ends[pos-1],ends[pos],ends[pos+1])

        
#    if orfPos != -1:
    if endToCDS.has_key( orfThreePrime ):
#        f = cds[orfPos]
        f = endToCDS[ orfThreePrime ]
        hitByOrf.add( f.qualifiers['locus_tag'][0] )
        print f
        print orfSeq.seq
        print "Orf begin %d end %d" % (tmpOrf.location.start, tmpOrf.location.stop )
        print f.extract( nucs )
    else:
        print "No CDS end found at %d" % tmpOrf.location.stop

timecheck = checkTimeAndMem(timecheck,count)

#allLocus = set([c.qualifiers['locus_tag'][0] for c in cds])
allLocus = set([c.qualifiers['locus_tag'][0] for c in endToCDS.values()])
notHit = allLocus - hitByOrf
print "Number of loci without an ORF end hit %d" % len(notHit)
for locus in notHit:
    print locus

timecheck = checkTimeAndMem(timecheck,count)
