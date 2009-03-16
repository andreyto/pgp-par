UsageInfo = """MapNTerminalPeptides.py
Given a set of annotation for some spectra, and a protein database
we are going to determine how far from the initial methionine
the n-terminal most peptide is, and then give that out in a 
histogram.

Required Options
 -r [FileName] File or directory of results.
 -d [Trie file] Proteome Database used for search
 -w [FileName]  Output file
 
Additional Options
 -v   Flag to print a more verbose trace of the program progression
     through various methods and whatnot
 -p [float] Pvalue cutoff (Default 0.05)
 -t Tryptic flag.  When set, this excluded fully tryptic peptides from the result.
     Default includes any peptide

"""

import sys
import os
import getopt
import traceback
import ResultsParser
import SelectProteins
import BasicStats
import ProteinStatistics
from Utils import *
Initialize()


class FinderClass(ResultsParser.ResultsParser):
    def __init__(self):
        self.ReferenceResults = None
        self.OutputPath = "RenameYourOutput.txt"
        self.DatabasePaths = [] #possibly multiple
        self.AllPeptides = {} # AminoSequence -> sPeptide Object
        self.ProteinsAndFirstObservedAA = {} #proteinID-> list of start position of closest peptide
        self.SpectrumCount = 0
        self.PValueLimit = 0.05 #pvalues are 0.00 (good) to 1.0 (bad) 
        self.MinPeptidesPerProtein = 2
        self.Verbose = 0
        self.LimitTryptic = 0
        self.NameToID = {}
        ResultsParser.ResultsParser.__init__(self)
        
    def Main(self):
        
        self.ProteinPicker = SelectProteins.ProteinSelector()
        self.ProteinPicker.LoadMultipleDB(self.DatabasePaths)
        self.ProcessResultsFiles(self.ReferenceResults, self.ParseInspectCallback)
        print "I found %s peptides from %s spectra"%(len(self.AllPeptides), self.SpectrumCount)
        self.PutPeptidesOnProteins()
        #print self.ProteinsAndFirstObservedAA.values()
        self.PrintSpecificDetails()
        Histogram = BasicStats.CreateHistogram(self.ProteinsAndFirstObservedAA.values(), 0, 5)
        Handle = open(self.OutputPath, "wb")
        BasicStats.PrettyPrintHistogram(Histogram, Handle) 
        Handle.close()
        
    def PrintSpecificDetails(self):
        """quicky"""
        HPlot = ProteinStatistics.HyrdopathyPlot()
        for (ProteinName, PeptideStart) in self.ProteinsAndFirstObservedAA.items():
            MaxLen = 50
            if PeptideStart == 0:
                ID = self.NameToID[ProteinName]
                BeginningSequence = self.ProteinPicker.ProteinSequences[ID][:5]
                String = "%s\t%s\t%s"%(ProteinName, "Mcapture", BeginningSequence)
                print String
            if PeptideStart > 15 and PeptideStart < MaxLen:
                ID = self.NameToID[ProteinName]
                PrefixSequence = self.ProteinPicker.ProteinSequences[ID][:PeptideStart]
                AfterCut = self.ProteinPicker.ProteinSequences[ID][PeptideStart:PeptideStart+2]
                #if PrefixSequence[-1] == "M":
                #    print "%s\t%s\t%s\t%s"%(ProteinName, PeptideStart, PrefixSequence, AfterCut)
                #continue
                #now some trickery
                Plot = HPlot.MakePlot(PrefixSequence)
                LongerSignal = HPlot.IsConsistentSignal(Plot) #default to 10
                ShorterSignal = HPlot.IsConsistentSignal(Plot, 0.5, 8)
                #if ShorterSignal and not FoundSignal:
                #    print "SHORTER ! %s\t%s\t%s\t%s\t"%(ProteinName, PeptideStart, PrefixSequence, AfterCut)
                String = "%s\t%s\t%s\t%s\t%s\t"%(ProteinName, PeptideStart, PrefixSequence, AfterCut, ShorterSignal)
                PlotString =""
                #freaking python, making me cast all of this crap
                # need to front pad these with zeros
                PadLen = MaxLen - len(Plot)
                for i in range(PadLen):
                    PlotString += "\t0"
                
                for Item in Plot:
                    PlotString += "\t%s"%Item
                #print "%s%s"%(String, PlotString) #plotString has an explicit \t at start
                


    def PutPeptidesOnProteins(self):
        """Get the location for the start of the peptide, and then just keep track of the smallest one
        per protein
        """
        PeptideObjectInProtein = {} #protein->list of starts
        PeptideStartInProtein = {}
        PeptidePerProteinCounter = {}
        for Aminos in self.AllPeptides.keys():
            DBLocations = self.ProteinPicker.FindPeptideLocations(Aminos)
            for (ProteinID, PeptideStartAA) in DBLocations:
                ProteinName = self.ProteinPicker.ProteinNames[ProteinID]
                self.NameToID[ProteinName] = ProteinID
                if not PeptidePerProteinCounter.has_key(ProteinName):
                    PeptidePerProteinCounter[ProteinName] = 1
                    PeptideObjectInProtein[ProteinName] = self.AllPeptides[Aminos]
                    PeptideStartInProtein[ProteinName] = PeptideStartAA
                else:
                    PeptidePerProteinCounter[ProteinName] += 1
                    if PeptideStartAA < PeptideStartInProtein[ProteinName]:
                        PeptideObjectInProtein[ProteinName] = self.AllPeptides[Aminos]
                        PeptideStartInProtein[ProteinName] = PeptideStartAA
                        
        #now that we're done going through all our peptides, we limit it to proteins passing our min peptides
        for (ProteinName, Count) in PeptidePerProteinCounter.items():
            if Count < self.MinPeptidesPerProtein:
                continue
            #now we see if this passes our tryptic concerns
            PeptideObject = PeptideObjectInProtein[ProteinName]
            if self.LimitTryptic and (PeptideObject.Prefix in ["R", "K"]):
                continue
            #okay, apparently we care about this protein.  Let's print stuff out
            self.ProteinsAndFirstObservedAA[ProteinName] = PeptideStartInProtein[ProteinName]

        
    def ParseInspectCallback(self, FilePath):
        """Called by the ResultsParser.ProcessResultsFiles
        Here I parse out Inspect Results to get peptide annotations, 
        Putting them in a hash for safe keeping
        """
        Handle = open(FilePath, "rb")
        for Line in Handle.xreadlines():
            #print "Parsing line %s"%Line
            if Line[0] == "#":
                continue # comment line
            if not Line.strip():
                continue
            Bits = list(Line.split("\t"))
            try:
                Annotation = Bits[self.Columns.Annotation]
                Peptide = GetPeptideFromModdedName(Annotation)
                Aminos = Peptide.Aminos
                PValue = float(Bits[self.Columns.PValue])
            except:
                traceback.print_exc()
                continue # SNAFU
            #here's some hacking that needs to be fixed.  I currently want to filter stuff by lfdr, but
            #that may not always be present
            if len(Bits) < self.Columns.LFDR: #meaning that I don't have the column in question
                if PValue > self.PValueLimit: #substitute pvalue for LFDR
                    continue
            else:
                LFDR = float(Bits[self.Columns.LFDR])
                if LFDR > self.PValueLimit:
                    continue
            if not self.AllPeptides.has_key(Peptide.Aminos):
                self.AllPeptides[Aminos] = Peptide
            self.SpectrumCount += 1
        Handle.close()


    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:d:w:vtp:")
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
                self.DatabasePaths.append( Value)
            if Option == "-w":
                self.OutputPath = Value
            if Option == "-p":
                #undocumented pvalue filter
                self.PValueLimit = float(Value)
            if Option == "-v":
                self.Verbose = 1
            if Option == "-t":
                self.LimitTryptic =1
        if not OptionsSeen.has_key("-r")  or not OptionsSeen.has_key("-w") or not OptionsSeen.has_key("-d"):
            print UsageInfo
            sys.exit(1)




if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(psyco not found - running in non-optimized mode)"
    Gumshoe = FinderClass()
    Gumshoe.ParseCommandLine(sys.argv[1:])
    Gumshoe.Main()