#!/usr/bin/env python

'''
Created on Feb 18, 2010

@author: eli

Some unit tests for the gff code
'''
import unittest
import os
import StringIO

import GFFIO

class Test(unittest.TestCase):

    IN = 'test.gff'

    def testGFFInput(self):
        "Reading of a GFF file."
        recCount = 0
        for gff in GFFIO.File( open(self.IN) ):
            recCount += 1
            if recCount == 1:
                self.assertEqual( 'ctg123', gff.seqid )
                self.assertEqual( 'gene', gff.type )
                self.assertEqual( 1000, gff.start )
                self.assertEqual( 9000, gff.end )
                self.assertEqual( '+', gff.strand )
                self.assertEqual( 'gene00001', gff.attributes['ID'] )
                self.assertEqual( 'EDEN', gff.attributes['Name'] )
            elif recCount == 7:
                self.assertEqual( 'ctg123', gff.seqid )
                self.assertEqual( 'exon', gff.type )
                self.assertEqual( 1050, gff.start )
                self.assertEqual( 1500, gff.end )
                self.assertEqual( '+', gff.strand )
                self.assertEqual( 'exon00002', gff.attributes['ID'] )
                self.assertEqual( 'mRNA00001,mRNA00002', gff.attributes['Parent'] )
            elif recCount == 17:
                self.assertEqual( 'ctg123', gff.seqid )
                self.assertEqual( 'CDS', gff.type )
                self.assertEqual( 7000, gff.start )
                self.assertEqual( 7600, gff.end )
                self.assertEqual( '+', gff.strand )
                self.assertEqual( 0, gff.phase )
                self.assertEqual( 'cds00002', gff.attributes['ID'] )
                self.assertEqual( 'mRNA00002', gff.attributes['Parent'] )
                self.assertEqual( 'edenprotein.2', gff.attributes['Name'] )

        self.assertEqual(17, recCount )

    def testGFFOutput(self):
        "Writing a gff record to a string file."
        rec = GFFIO.Record()
        shandle = StringIO.StringIO()
        gffout = GFFIO.File(shandle)
        gffout.write(rec)
        self.assertEqual("\t".join(list('.'*9)), shandle.getvalue().strip())

    def testGFFRoundTrip(self):
        "GFF roundtrip from input to output" 
        OUT = '.test.gff'
        gffout = GFFIO.File( open(OUT,'w'))
        gffs = []
        for gff in GFFIO.File( open(self.IN) ):
            gffs.append( gff )
            gffout.write( gff )

        gffout.close()

        i = 0
        for gff in GFFIO.File( open(OUT) ):
            self.assertEqual(gffs[i].seqid, gff.seqid)
            self.assertEqual(gffs[i].source, gff.source)
            self.assertEqual(gffs[i].type, gff.type)
            self.assertEqual(gffs[i].start, gff.start)
            self.assertEqual(gffs[i].end, gff.end)
            self.assertEqual(gffs[i].score, gff.score)
            self.assertEqual(gffs[i].strand, gff.strand)
            self.assertEqual(gffs[i].phase, gff.phase)
            self.assertEqual(gffs[i].attributes, gff.attributes)
            i += 1

        os.remove( OUT )

if __name__ == "__main__":
    unittest.main()
