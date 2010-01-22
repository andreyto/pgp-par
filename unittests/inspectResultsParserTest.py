#!/usr/bin/env python

'''
Created on Oct 8, 2009

@author: eli

Some unit tests for the InspectResults Parser class.
'''
import unittest

import InspectResults

class Test(unittest.TestCase):

    IN = 'inspect.txt.bz2'
    ANNO = [
    'R.GGAAPA.P',
    'G.GGAAAP.G',
    'H.GGALAG.I',
    'A.GGAIAG.A',
    'Q.GGALGA.G',
    'Q.GGAIGA.N',
    'L.GGIAAG.A',
    'G.GGLAAG.E',
    'S.GGLAGA.R']
    PVAL = [
    0.99972,
    0.99990,
    0.99999,
    0.99999,
    0.99999,
    0.99999,
    0.99999,
    0.99999,
    0.99999]


    def testResultParserPrint(self):
        parser = InspectResults.Parser( self.IN )
        for result in parser:
            continue
#            print result
#            print "%s\n" % result

    def testResultParser2Cols(self):
        "Test that the parser is reading the ANNO and PVAL columns correctly."
        parser = InspectResults.Parser( self.IN )
        i = 0
        self.assertEqual( 9, len(self.ANNO) )
        self.assertEqual( 9, len(self.PVAL) )

        for result in parser:
            self.assertEqual( self.ANNO[i], result.Annotation )
            self.assertEqual( self.PVAL[i], result.PValue )

            i += 1

        self.assertEqual( 9, i)

if __name__ == "__main__":
    unittest.main()
