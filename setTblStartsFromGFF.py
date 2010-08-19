#!/usr/bin/env python

from optparse import OptionParser

import re, sys, bioseq, GFFIO

class StartsGFF(object):
    def __init__(self,startGFF):
        self.gff    = startGFF
        self.starts = {}
        self.newGenes = {}

    def checkStart(self,gffRec):
        if gffRec.end - gffRec.start != 2:
            raise ValueError("Start codon wrong size: %s" % gffRec)

    def readGFF(self):
        lt = 'locus_tag'
        for gffRec in GFFIO.File(self.gff):
            if gffRec.attributes.has_key( lt ):
                self.checkStart( gffRec )
                if self.starts.has_key( gffRec.attributes[ lt ] ):
                    print "Warning multiple starts for %s, skipping" % gffRec.attributes[lt]
                    self.starts[gffRec.attributes[ lt ] ] = None
                    continue
                self.starts[gffRec.attributes[ lt ] ] = gffRec

            elif gffRec.attributes['Name'] == 'Observed2Stop':
                self.newGenes[ gffRec.attributes['Parent']] = (gffRec,[])
            else:
                self.checkStart( gffRec )
                self.newGenes[ gffRec.attributes['Parent']][1].append(gffRec)

    def novelIter(self):
        for (parent,value) in self.newGenes.items():
            observedRec = value[0]
            if len(value[1]) == 0:
                continue
            elif len(value[1]) > 1:
                print "Warning multiple starts for %s, skipping" % parent
                continue
            startRec = value[1][0]
            start = startRec.start
            end   = observedRec.end
            if observedRec.strand == '-':
                start = startRec.end
                end   = observedRec.start
            yield( (parent, start, end))

class AlterTblStarts(object):
    def __init__(self, startGFF, tblInput, outputPath):
        self.startGFF = StartsGFF(startGFF)
        self.tbl    = tblInput
        self.output = bioseq.FlatFileIO(outputPath,'w')
        self.excludes = []
        self.curRec = ''
        self.curBeg = None
        self.curEnd = None
        self.curType = ''
        self.locusPrefix = 'locus_'
        self.locusOffset = 5000

    def starts(self):
        return self.startGFF.starts

    def readExcludeFile(self, excludeFile):
        for line in open(excludeFile):
            l = line.strip()
            # Tbl file has it's protein_id's wrapped in ref| |
            # so make sure the exclude list has them as such
            if l[0:4] == 'ref|':
                self.excludes.append( l )
            else:
                self.excludes.append("ref|%s|"%l)

    def writeRec(self):
        if self.curRec:
            self.output.write( "%d\t%d\t%s\n%s" % (self.curBeg, self.curEnd,
                                               self.curType, self.curRec) )
    def processTbl(self):
        tbl = bioseq.FlatFileIO( self.tbl )
        # RE to match the start of a feature record
        recRE  = re.compile('^(\d+)\t(\d+)\t(\w+)')
        # RE to match the qualifiers associated in the record
        qualRE = re.compile('^\t\t\t(\w+)\t?(.*)')
        # Accession of current sequence in tbl file
        curSeq = None
        modifyCDS = False
        origBeg = None
        # Iterate through the tbl file a line at at time, reading it's features
        # and modifing the begin coordinates of both gene and CDS records
        for line in tbl.io:
            # Match start of a contig sequence and get it's ID
            if line[0] == '>':
                accStart = line.index('|') + 1
                curSeq = line[accStart:-2]
                self.writeRec()
                self.output.write( line )
                continue

            # Match beginning of a record
            match = recRE.match( line )
            if match:
                # Write the previous record
                self.writeRec()
                self.curRec = ''
                # begin, end and keep track of the current record type
                (b,e,self.curType) = match.groups()
                if modifyCDS and self.curType == 'CDS':
                    # Keep the same begin & end as for the preceeding gene record
                    modifyCDS = False
                else:
                    # Use the begin and end from the tbl Record
                    (self.curBeg,self.curEnd) =(int(b),int(e))
                    origBeg = self.curBeg
                continue

            match = qualRE.match( line )
            # Match the qualifiers in the tbl records
            if match:
                (qual,val) = match.groups()
                if qual == 'locus_tag' and self.starts().has_key(val):
                    # We need to modify this records start
                    gff = self.starts()[val]
                    if not gff: # Get none if multiple starts set, so skip
                        continue
                    seqId = gff.seqid
                    # Sometimes GFF sequence don't have the .# so check for that
                    tSeq = curSeq
                    if -1 == seqId.find('.'):
                        print "Warning, unversioned sequence from GFF %s."%seqId
                        tSeq = curSeq[0:curSeq.index('.')]

                    if seqId != tSeq:
                        print "Wrong seq for %s, %s, %s." % (val,seqId,tSeq)
                        self.curRec += line
                        continue
                    # Save the original begin, in case we end up excluding this protein
                    origBeg = self.curBeg
                    if gff.strand == '+':
                        self.curBeg = gff.start
                    else:
                        self.curBeg = gff.end
                    modifyCDS = True
                elif qual == 'protein_id':
                    for exclude in self.excludes:
                        if -1 != val.find( exclude ):
                            # This protein matches our exclude list, so don't change it
                            self.curBeg = origBeg
                            break
                else:
                    pass
                self.curRec += line
            else:
                print "Bad line: %s" % line
        # write out the last record
        self.writeRec()

    def writeNovel(self):
        for (orfName,start,end) in self.startGFF.novelIter():
            self.output.write( "%d\t%d\t%s\n" % (start,end,'gene'))
            self.output.write("\t\t\tlocus_tag\t%s%d\n" % (self.locusPrefix,self.locusOffset))
            self.locusOffset += 1
            self.output.write( "%d\t%d\t%s\n" % (start,end,'CDS'))
            self.output.write("\t\t\tproduct\t%s\n" % orfName)

    def Main(self):
        self.startGFF.readGFF()
        self.processTbl()
        self.writeNovel()

def ParseCommandLine():
    Desc = 'Alter the starts of genes in a tbl file using a GFF file of starts.'
    opts = OptionParser(description=Desc)
    opts.add_option('-g','--gff',dest='startGFF',
                   help='Input file path of GFF file with new start codons.')
    opts.add_option('-t','--tbl',dest='inTbl',
                   help='Input file path of tbl file with gene annotations.')
    opts.add_option('-o','--output',
               help='Output file path of altered tbl file to write, can be .gz')
    opts.add_option('-l','--locusPrefix',
               help='Prefix for locus_tag written out for novel genes.')
    opts.add_option('-s','--locusStart', type='int',
               help='Start number to use for locus counting.')
    opts.add_option('-x','--excludeFile',
               help='File of protein_ids, 1 per line, to exclude from modification.')
    (options,args) = opts.parse_args()

    if options.startGFF and options.inTbl and options.output and \
       options.locusPrefix and options.locusStart:
        pass
    else:
        opts.error('gff, tbl, output, locusPrefix and locusStart must all be set.')

    return options

if __name__ == '__main__':
    options = ParseCommandLine()
    driver = AlterTblStarts( options.startGFF, options.inTbl, options.output)
    driver.locusPrefix = options.locusPrefix
    driver.locusOffset = options.locusStart
    if options.excludeFile:
        driver.readExcludeFile( options.excludeFile )
    driver.Main()
