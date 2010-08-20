#!/usr/bin/env python

from optparse import OptionParser

import os, re, sys, bioseq, GFFIO

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
        self.outputPath = outputPath
        self.output = bioseq.FlatFileIO(outputPath,'w')
        self.excludes = []
        self.cdsRec = ''
        self.otherRec = ''
        self.notFirst = False
        self.curBeg = None
        self.curEnd = None
        self.curType = ''
        self.prevType = ''
        self.locusPrefix = 'locus_'
        self.locusOffset = 5000
        self.excludeCurrent = False

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

    def writeRec(self,last=False):
        if self.curType == 'CDS' and not last:
            return # Don't write out CDS features, unless it's the last
        elif self.prevType == 'CDS' or last:
            if not self.excludeCurrent:
                self.output.write( "%d\t%d\tgene\n%s" % ( self.curBeg, self.curEnd, self.otherRec) )
                self.output.write( "%d\t%d\tCDS\n%s" % ( self.curBeg, self.curEnd, self.cdsRec) )
            # Reset everything
            self.excludeCurrent = False
            self.otherRec = ''
            self.cdsRec = ''
        elif self.notFirst and self.otherRec:
            # output other records like RNA's that still have gene's
            self.output.write( "%d\t%d\t%s\n%s" % (
                self.curBeg, self.curEnd, self.prevType, self.otherRec) )
            self.otherRec = ''
        self.notFirst = True

    def saveLine(self,line):
        if self.curType == 'CDS':
            self.cdsRec += line
        else:
            self.otherRec += line

    def processTbl(self):
        tbl = bioseq.FlatFileIO( self.tbl )
        # RE to match the start of a feature record
        recRE  = re.compile('^(\d+)\t(\d+)\t(\w+)')
        # RE to match the qualifiers associated in the record
        qualRE = re.compile('^\t\t\t(\w+)\t?(.*)')
        # Accession of current sequence in tbl file
        curSeq = None
        # Flag set when reading a gene to indicate changing the coords
        modifyCDS = False
        # Iterate through the tbl file a line at at time, reading it's features
        # and modifing the begin coordinates of both gene and CDS records
        # Flow is, read 1st record and store. Read next record and output previous
        # unless new record is a CDS then save it until the record after when
        # both the gene and CDS can be output (or rejected). Also, when we get
        # a gene record check for coordinate modification
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
                self.prevType = self.curType
                # begin, end and keep track of the current record type
                (b,e,self.curType) = match.groups()
                # Write the previous record
                self.writeRec()
                if modifyCDS:
                    # Use the begin and end from the gff, and reset modify
                    modifyCDS = False
                else:
                    # Use the begin and end from the tbl Record
                    (self.curBeg, self.curEnd) = (int(b),int(e))
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
#                        print "Warning, unversioned sequence from GFF %s."%seqId
                        tSeq = curSeq[0:curSeq.index('.')]

                    if seqId != tSeq:
                        print "Wrong seq for %s, %s, %s." % (val,seqId,tSeq)
                        self.saveLine( line )
                        continue
                    if gff.strand == '+':
                        self.curBeg = gff.start
                    else:
                        self.curBeg = gff.end
                    modifyCDS = True
                elif qual == 'protein_id':
                    for exclude in self.excludes:
                        if -1 != val.find( exclude ):
                            # This protein matches our exclude list, don't output it
                            self.excludeCurrent = True
                            break
                else:
                    pass
                self.saveLine( line )
            else:
                print "Unknown line: %s" % line
                self.saveLine( line )
        # write out the last record
        self.writeRec(True)

    def writeNovel(self):
        for (orfName,start,end) in self.startGFF.novelIter():
            self.output.write( "%d\t%d\t%s\n" % (start,end,'gene'))
            self.output.write("\t\t\tlocus_tag\t%s%d\n" % (self.locusPrefix,self.locusOffset))
            self.locusOffset += 1
            self.output.write( "%d\t%d\t%s\n" % (start,end,'CDS'))
            self.output.write("\t\t\tproduct\thypothetical protein\n" )

    def Main(self):
        self.startGFF.readGFF()
        self.processTbl()
        self.writeNovel()
        # Record the differences between the original and our alternate
        self.output.close()
        os.system("diff -u %s %s > %s.diff" % (self.tbl, self.outputPath, self.outputPath))

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
