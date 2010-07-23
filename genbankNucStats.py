#!/usr/bin/env python

import bioseq
from optparse import OptionParser
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature,FeatureLocation
from NucleotideStatistics import GetGC
from setTblStartsFromGFF import StartsGFF

class SequenceStats(bioseq.FlatFileIO):
    def __init__(self,gbk,startGff):
        bioseq.FlatFileIO.__init__(self,gbk)
        self.gff = StartsGFF(startGff)

    def genStats(self, dna, context, locus):
        gc = GetGC(dna)
        print "GC is %f for %s %s" % (gc, context, locus)
#        codonVec = CodonUsageFractions(dna)

    def processGBK(self):
        novelCount = 0
        for gb_rec in SeqIO.parse(self.io, 'genbank'):
            # First check for novel genes on this sequence
            for (orfName, start, end) in self.gff.novelIter():
                acc = orfName[ 0:orfName.index('.') ]
                if acc == gb_rec.name:
                    print "Novel %d %d" % (start,end)
                    if start < end:
                        f = SeqFeature(FeatureLocation( start, end), strand=1)
                    else:
                        f = SeqFeature(FeatureLocation( end, start), strand=-1)
                    dna = f.extract(gb_rec.seq)
                    self.genStats(dna,'Novel',"Novel%d"%novelCount)
                    novelCount+=1

            # Now do the regular genes and extensions
            for feat in gb_rec.features:
                if feat.type == 'CDS':
                    dna = feat.extract(gb_rec.seq)
                    locus = feat.qualifiers['locus_tag'][0]
                    self.genStats(dna, 'Original',locus)
                    if self.gff.starts.has_key( locus ):
                        gffRec = self.gff.starts[locus]
                        print "%s gene %s start %d %d" % ( gffRec.strand,
                                    feat.location, gffRec.start,gffRec.end)
                        if gffRec.strand == '+':
                            f = SeqFeature(FeatureLocation( gffRec.start, feat.location.start),
                                           strand=1)
                        else:
                            f = SeqFeature(FeatureLocation( feat.location.end, gffRec.end),
                                           strand=-1)
                        dna = f.extract(gb_rec.seq)
                        self.genStats(dna,'Extension',locus)


    def Main(self):
        self.gff.readGFF()
        self.processGBK()

def ParseCommandLine():
    Desc = 'Reads in a genbank file and starts GFF and dumps stats on the genes.'
    opts = OptionParser(description=Desc)
    opts.add_option('-s','--gff',
                   help='Input file path of GFF file with new start codons.')
    opts.add_option('-k','--gbk',
                   help='Input file path of genbank file with gene annotations.')
    opts.add_option('-o','--output',
               help='Output file path of altered tbl file to write, can be .gz')
    (options,args) = opts.parse_args()

    if options.gff and options.gbk: # and options.output
        pass
    else:
        opts.error('gff, gbk must all be set.')

    return options

if __name__ == '__main__':
    options = ParseCommandLine()
    driver = SequenceStats(options.gbk, options.gff)
    driver.Main()
