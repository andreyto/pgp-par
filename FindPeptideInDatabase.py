"""FindPeptideInDatabase.py
This is a quicky script to search for our novel peptides in TAIR8

"""
UsageInfo = """FindPeptideInDatabase.py
This is a quicky script to search for our novel peptides in TAIR8

Required Options
 -r [Directory] Directory where the reference results files are kept.
 -d [Trie file] Database

"""

import sys
import os
import getopt
import traceback
import ResultsParser
import SelectProteins
from Utils import *
Initialize()


class CompileClass(ResultsParser.ResultsParser):
    def __init__(self):
        self.ReferenceResults = None
        self.OutputPath = None # the training set
        self.ReferenceDatabasePath = [] #possibly multiple
        ResultsParser.ResultsParser.__init__(self)
        #hopefully all programs will annotate scans with the same peptide, but it's no guarentee
        self.AllReferenceAnnotations = {} # {Aminos} => spectrum count
        self.GlobalCount = 0
        self.GlobalSpliceCount = 0
        self.ProteinCoverageHash = {}  # {ProteinID} => [0,0,0,....] len of the sequence
        self.UniqueProteinSequences = {}
        self.SharedProteinSequences = {}
        self.ProteinSpectra = {}

        
    def Main(self):
        self.ProteinPicker = SelectProteins.ProteinSelector()
        self.ProteinPicker.LoadMultipleDB(self.ReferenceDatabasePath)
        self.ProcessResultsFiles(self.ReferenceResults, self.ParseReferenceFileCallback)
        self.FindInDatabase()



    def FindInDatabase(self):
        FoundCount = 0
        FoundAminos = []
        for Aminos in self.AllReferenceAnnotations.keys():
            Locations = self.ProteinPicker.FindPeptideLocations(Aminos)
            if len(Locations) > 0:
                FoundCount += 1
                print Aminos
                FoundAminos.append(Aminos)
        print "Out of %d novel peptides, I found %s"%(len(self.AllReferenceAnnotations), FoundCount)

    def ParseReferenceFileCallback(self, FilePath):
        """SInce these are novel results, we will be dealing with Natalie's format
        """
        Handle = open(FilePath, "rb")
        for Line in Handle.xreadlines():
            #print "Parsing line %s"%Line
            if Line[0] in ["#", "\n"]:
                continue # comment line or blank
            if not Line.strip():
                continue
            Bits = list(Line.split("\t"))
            try:
                Aminos = Bits[4]
            except:
                traceback.print_exc()
                continue # SNAFU
            #print Spectrum
            self.AllReferenceAnnotations[Aminos] = 1  # keep track of spectrum count, if I want to use it
                
        Handle.close()


    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:d:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-r":
                # -r results file(s)
                if not os.path.exists(Value):
                    print "** Error: couldn't find results file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.ReferenceResults = Value
            if Option == "-d":
                if not os.path.exists(Value):
                    print "** Error: couldn't find database file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.ReferenceDatabasePath.append( Value)
            if Option == "-t":
                self.TrackName = Value
        if not OptionsSeen.has_key("-r") or not OptionsSeen.has_key("-d"):
            print UsageInfo
            sys.exit(1)


if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(psyco not found - running in non-optimized mode)"
    Gather = CompileClass()
    Gather.ParseCommandLine(sys.argv[1:])
    Gather.Main()