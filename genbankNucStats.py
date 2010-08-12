#!/usr/bin/env python

import bioseq
from optparse import OptionParser
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
#from NucleotideStatistics import GetGC, CodonUsageFractions, CodonCount #make these functions locally scoped
import NucleotideStatistics
from setTblStartsFromGFF import StartsGFF

# Need R_HOME and LD_LIBRARY_PATH. Set I used:
# export R_HOME=/usr/local/packages/R-2.10.1/lib64/R
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$R_HOME/lib
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr

class SequenceStats(bioseq.FlatFileIO):
    def __init__(self,gbk,startGff):
        bioseq.FlatFileIO.__init__(self,gbk)
        self.gff = StartsGFF(startGff)
        self.gcContent = {}
        self.geneLen   = {}
        self.CodonMatrix = {}

    def genStats(self, dna, context, locus):
        for (vecs,func, kludge) in ((self.gcContent, NucleotideStatistics.GetGC, "a"), (self.geneLen,len, "b"), (self.CodonMatrix, NucleotideStatistics.CodonUsageFractions, "c")):
            data = func(dna) # apply the function from the list of tuples
            
            if not vecs.has_key(context):
                vecs[ context ] = []
                if kludge == "c":
                    print "I am making a context %s "%context
                    print vecs.keys()
            vecs[ context ].append( data )
#        print "stats for %s %s %d" % (context,locus,len(vecs[context]))


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
                    if self.gff.starts.has_key( locus ):
                        self.genStats(dna, 'Original',locus)
                        gffRec = self.gff.starts[locus]
                        print "%s gene %s start %d %d" % ( gffRec.strand,
                                    feat.location, gffRec.start,gffRec.end)
                        if gffRec.strand == '+':
                            f = SeqFeature(FeatureLocation( gffRec.start, feat.location.start),
                                           strand=1)
                        else:
                            f = SeqFeature(FeatureLocation( feat.location.end, gffRec.end),
                                           strand=-1)
                        dna = f.extract(gb_rec.seq)# the extracted dna sequence
                        self.genStats(dna,'Extension',locus)
                    else:
                        self.genStats(dna, 'Unchanged',locus)


    def writeDataVecs(self,vecs,outName):
        for (key,value) in vecs.items():
            out = open("%s.%s.vec"%(key,outName),'w')
            for v in value:
                out.write("%d\n"%v)
            out.close()
            
    def writeDataMatrix(self,matrix, outName):
        print matrix.keys()
        for (key, value) in matrix.items():
            print "Now using key %s"%key
            out = open("%s.%s.matrix"%(key,outName), 'w')
            #this value is a list of lists.  We need to get each list tab separated onto a line
            for List in value:
                ListAsAString =  "\t".join(map(str, List))#tab separate this list, cast numbers to string
                out.write("%s\n"%ListAsAString)
            out.close()

    def plotVecs(self,data,dataName):
        r = robjects.r
        rdevoff = r['dev.off']
        rVecs = {}
        self.writeDataVecs( data, dataName)
        for (key,value) in data.items():
            rVecs[key] = robjects.IntVector(value)

#        r.png("geneTypes_GC.png", width=1024, height=768)
        for key in data:
            r.png("%s_%s.png"%(key,dataName), width=1024, height=768)
            r.hist(rVecs[key], freq=False, breaks=10,main=key,xlab=dataName)
            lens = ['Extension','Original']
            orig = ['Unchanged','Novel']
            if key in lens:
                lens.remove(key)
                alt = lens[0]
            else:
                orig.remove(key)
                alt = orig[0]

            r.lines(r.density(rVecs[alt]), col='Blue')
            rdevoff()
#        colors = ['grey','blue','green','red']
        #r.require("plotrix",lib_loc="~eventer/R")
#        plotrix = importr("plotrix", lib_loc='/home/eventer/R')
#        plotrix.multhist(r.list(rVecs.values()), freq=False, col=colors,
#                   xlab="GC %",ylab="Density",main="Gene Type GC")
#        r.legend(0,0.08,['Unchanged','Novel','Original','Extension'],
#                text_col=colors)

#        rdevoff()

    def Main(self):
        self.gff.readGFF()
        self.processGBK()
        self.plotVecs(self.geneLen,'Length')
        self.plotVecs(self.gcContent,'GC')
        print self.CodonMatrix.keys()
        self.writeDataMatrix(self.CodonMatrix, 'CodonMatrix')

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
