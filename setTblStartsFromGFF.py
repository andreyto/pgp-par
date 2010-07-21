#!/usr/bin/env python

from optparse import OptionParser

import re, sys, bioseq, GFFIO

class AlterTblStarts:
    def __init__(self, startGFF, tblInput, outputPath):
        self.gff    = startGFF
        self.tbl    = tblInput
        self.output = bioseq.FlatFileIO(outputPath,'w')
        self.starts = {}
        self.newGenes = []
        self.curRec = ''
        self.curBeg = None
        self.curEnd = None
        self.curType = ''

    def readGFF(self):
        lt = 'locus_tag'
        for gffRec in GFFIO.File(self.gff):
            if gffRec.attributes.has_key( lt ):
                self.starts[gffRec.attributes[ lt ] ] = gffRec
            else:
                self.newGenes.append( gffRec )

    def writeRec(self):
        if self.curRec:
            self.output.write( "%d\t%d\t%s\n%s" % (self.curBeg,self.curEnd,
                                               self.curType, self.curRec) )
    def processTbl(self):
        tbl = bioseq.FlatFileIO( self.tbl )
        # RE to match the start of a feature record
        recRE  = re.compile('^(\d+)\t(\d+)\t(\w+)')
        # RE to match the qualifiers associated in the record
        qualRE = re.compile('^\t\t\t(\w+)\t(.*)')
        # Accession of current sequence in tbl file
        curSeq = None
        modifyCDS = False
        for line in tbl.io:
            if line[0] == '>':
                accStart = line.index('|') + 1
                curSeq = line[accStart:-2]
                self.writeRec()
                self.output.write( line )
                continue

            match = recRE.match( line )
            if match:
                # Match feature start
                self.writeRec()
                self.curRec = ''
                (b,e,self.curType) = match.groups()
                if modifyCDS and self.curType == 'CDS':
                    # Keep the same begin & end as for the previous gene
                    modifyCDS = False
                else:
                    (self.curBeg,self.curEnd) =(int(b),int(e))
                continue

            match = qualRE.match( line )
            if match:
                (qual,val) = match.groups()
                if qual == 'locus_tag' and self.starts.has_key(val):
                    # We need to modify this records start
                    gff = self.starts[val]
                    seqId = gff.seqid
                    if seqId != curSeq:
                        print "Wrong seq for %s, %s, %s" % (val,seqId,curSeq)
                        self.curRec += line
                        continue
                    if gff.strand == '+':
                        self.curBeg = gff.start
                    else:
                        self.curEnd = gff.end
                    modifyCDS = True
                else:
                    pass
                self.curRec += line
            elif line == "\t\t\tpseudo\n":
                # special case, pseudo lines don't match our qualRE
                self.curRec += line
            else:
                print "Bad line: %s" % line
        # write out the last record
        self.writeRec()

    def Main(self):
        self.readGFF()
        self.processTbl()

def ParseCommandLine():
    Desc = 'Alter the starts of genes in a tbl file using a GFF file of starts.'
    opts = OptionParser(description=Desc)
    opts.add_option('-g','--gff',dest='startGFF',
                   help='Input file path of GFF file with new start codons.')
    opts.add_option('-t','--tbl',dest='inTbl',
                   help='Input file path of tbl file with gene annotations.')
    opts.add_option('-o','--output',
               help='Output file path of altered tbl file to write, can be .gz')
    (options,args) = opts.parse_args()

    if not options.startGFF and not options.inTbl and not options.output:
        opts.error('gff, tbl, and output must all be set.')

    return options

if __name__ == '__main__':
    options = ParseCommandLine()
    driver = AlterTblStarts( options.startGFF, options.inTbl, options.output)
    driver.Main()
