UsageInfo = """CountPhosphoPeptides.py - 
  This simple script takes an input file and counts the number of phospho
  peptides.  Since this is actually a rather complex issue, depending on the
  placement of the phosphate group, I report a number of things and then
  dump the list for your pleasure.

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
        print "Total number of peptides: %d"%(len(self.AllPeptides))
        self.PrintModdedOverlap()

    def PrintModdedOverlap(self):
        SpectraCount =0
        for Annotation in self.ModifiedPeptides.keys():
            Peptide = GetPeptideFromModdedName(Annotation)
            Aminos = Peptide.Aminos
            SpectraCount += self.ModifiedPeptides[Annotation]
            #print "%s\t%s"%(Aminos, Annotation)
        print "%d total spectra were identified with a PTM"%SpectraCount
        print "there were %d phosphopeptides"%len(self.ModifiedPeptides)
        print "There were %d phosphopeptides, when just looking at amino acid sequence"%len(self.ModifiedStrippedPeptides)
        

    def ParseFileCallback(self, FileName):
        print "Parse %s..."%FileName
        File = open(FileName, "rb")
        Stub = os.path.split(FileName)[1]
        LineNumber = 0
        for FileLine in File.xreadlines():
            LineNumber += 1
            if LineNumber % 100000 == 0:
                print "%s %s..."%(Stub, LineNumber)
            if FileLine[0] == "#":
                continue
            if not FileLine.strip():
                continue                        
            Bits = FileLine.split("\t")
            Annotation = Bits[self.Columns.Annotation]
            self.AddPeptide(Annotation)

    def AddPeptide(self, Annotation):
        Peptide = GetPeptideFromModdedName(Annotation)
        Aminos = Peptide.Aminos

        ##1. to the all
        if not self.AllPeptides.has_key(Aminos):
            self.AllPeptides[Aminos] = 0
        self.AllPeptides[Aminos] += 1
        ## 2. to the modified or not list
        if not len(Peptide.Modifications) == 0: # has mod
            if Annotation.find("phos") == -1:
                # I hate modified peptides, which are not phosphorylated
                # the purpose of this script is to focus on phosphopeptides!!!!
                return
            if not self.ModifiedPeptides.has_key(Annotation):
                self.ModifiedPeptides[Annotation] = 0
            self.ModifiedPeptides[Annotation] += 1
            if not self.ModifiedStrippedPeptides.has_key(Aminos):
                self.ModifiedStrippedPeptides[Aminos] = 0
            self.ModifiedStrippedPeptides[Aminos] += 1
        else: #no mod
            if not self.UnmodifiedPeptides.has_key(Aminos):
                self.UnmodifiedPeptides[Aminos] =0
            self.UnmodifiedPeptides[Aminos] += 1
            
        


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