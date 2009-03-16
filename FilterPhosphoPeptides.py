UsageInfo = """CountPhosphoPeptides.py - 
  This simple script takes an input file and prints to output only the
  spectra with a phosphorylation

Required options:
 -r [FileName] - The name of the results file to parse.  If a directory is
    specified, then all .txt files will be read.
"""

import os
import sys
import getopt
import ResultsParser
from Utils import *
Initialize()

class CountClass(ResultsParser.ResultsParser):
    def __init__(self):
        self.ResultsDir = None
        self.AllPeptides = {}
        self.ModifiedPeptides = {}
        self.ModifiedStrippedPeptides = {}
        self.UnmodifiedPeptides = {}        
        ResultsParser.ResultsParser.__init__(self)

    def Main(self):
        self.ProcessResultsFiles(self.ResultsDir, self.ParseFileCallback)

        

    def ParseFileCallback(self, FilePath):
        print "Parse %s..."%FilePath
        File = open(FilePath, "rb")
        (Path, Stub) = os.path.split(FilePath)
        (ShortName, Ext) = os.path.splitext(Stub)
        OutFileName = "%s.PhosOnly.txt"%ShortName
        OutPath = os.path.join(Path, OutFileName)
        OutHandle = open(OutPath, "wb")
        LineNumber = 0
        for FileLine in File.xreadlines():
            LineNumber += 1
            if LineNumber % 100000 == 0:
                print "%s %s..."%(Stub, LineNumber)
            if FileLine[0] == "#":
                OutHandle.write(FileLine)
                continue
            if not FileLine.strip():
                continue                        
            Bits = FileLine.split("\t")
            Annotation = Bits[self.Columns.Annotation]
            if not Annotation.find("phos") == -1:
                OutHandle.write(FileLine)

            
        


    def ParseCommandLine(self, Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-r":
                self.ResultsDir = Value
            else:
                print "** Warning: Option %s not understood"%Option
        if not self.ResultsDir:
            print UsageInfo
            sys.exit(-1)




if __name__ =="__main__":
    Abacus = CountClass()
    Abacus.ParseCommandLine(sys.argv[1:])
    Abacus.Main()