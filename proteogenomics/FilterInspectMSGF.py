#!/usr/bin/env python

from optparse import OptionParser

import sys, InspectResults, bioseq

class FilterInspect:
    def __init__(self, pCutOff, justBest=None, inputPath=None, outputPath=None):
        self.pCutOff    = pCutOff
        self.justBest   = justBest
        self.inputPath  = inputPath
        self.output     = bioseq.FlatFileIO(outputPath,'w')

    def filter(self):
        parser = InspectResults.Parser(self.inputPath)
        best = None
        bestPvalue = 1.1
        lastScan = -1
        for result in parser:
            # LFDR is the 21st column, also the MSGF score
            pval = result.LFDR
            scan = result.ScanNumber
            # we have a new set of scan results
            if scan != lastScan:
                if self.pCutOff > bestPvalue:
                    self.output.write( str(best) )
                best = result
                bestPvalue = pval

            if pval <= self.pCutOff:
                if self.justBest:
                    # Want just the best result per scan
                    if pval < bestPvalue:
                        best = result
                        bestPvalue = pval
                else:
                    # Want all results below the cutoff
                    self.output.write( str(result) )
            lastScan = scan
        # After for loop print out best result for last scan
        if self.pCutOff > bestPvalue:
            self.output.write( str(best) )

    def Main(self):
        self.filter()

def ParseCommandLine():
    Desc = 'Filter a MSGFInspectOutput file given a maximum P-value cutoff.'
    opts = OptionParser(description=Desc)
    opts.add_option('-a','--all',dest='justBest',default=True, action='store_false',
                   help='Display all values instead of just the best.')
    opts.add_option('-p','--pvalue',dest='pcutoff',type='float',
                   help='Maximum P-Value cutoff for filtering.' )
    opts.add_option('-r','--input','-i',dest='inPath',
                   help='Input file path of InspectResults.')
    opts.add_option('-w','--output','-o',dest='outPath',
                   help='Output file path of filtered results')
    (options,args) = opts.parse_args()

    if not options.inPath and not options.pcutoff and not options.outPath:
        opts.error('Input, Output, and P-Value must all be set.')

    return options

if __name__ == '__main__':
    options = ParseCommandLine()
    driver = FilterInspect( options.pcutoff, options.justBest,
                            options.inPath, options.outPath)
    driver.Main()
