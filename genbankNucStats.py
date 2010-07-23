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

    def processGBK(self):
        for gb_rec in SeqIO.parse(self.io, 'genbank'):
            for feat in gb_rec.features:
                if feat.type == 'CDS':
                    dna = feat.extract(gb_rec.seq)
                    locus = feat.qualifiers['locus_tag'][0]
                    self.genStats(dna, 'Full',locus)
                    if self.gff.starts.has_key( locus ):
                        gffRec = self.gff.starts[locus]
                        if gffRec.strand == '+':
                            f = SeqFeature(FeatureLocation( gffRec.start, feat.location.start),
                                           strand=1)
                        else:
                            f = SeqFeature(FeatureLocation( feat.location.end, gffRec.end),
                                           strand=-1)
                        dna = f.extract(gb_rec.seq)
                        self.genStats(dna,'Extension',locus)


#                    codonVec = CodonUsageFractions(dna)

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
