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
 -s [FileName] Signal peptide predictions from signalp

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
        self.PredictedSignalPeptideFile = None
        self.OutputPath = "RenameYourOutput.txt"
        self.DatabasePaths = [] #possibly multiple
        self.AllPeptides = {} # AminoSequence -> sPeptide Object
        self.ProteinsAndFirstObservedAA = {} #proteinID-> list of start position of closest peptide
        self.FirstPeptideInProtein = {}
        self.SpectrumCount = 0
        self.PValueLimit = 0.05 #pvalues are 0.00 (good) to 1.0 (bad) 
        self.MinPeptidesPerProtein = 2
        self.Verbose = 0
        self.NameToID = {}
        ResultsParser.ResultsParser.__init__(self)
        
    def Main(self):
        
        
        self.ProteinPicker = SelectProteins.ProteinSelector()
        self.ProteinPicker.LoadMultipleDB(self.DatabasePaths)
        self.ProcessResultsFiles(self.ReferenceResults, self.ParseInspectCallback)
        print "I found %s peptides from %s spectra"%(len(self.AllPeptides), self.SpectrumCount)
        self.PutPeptidesOnProteins()
        self.SignalPeptideHunt()
        #self.NMEHunt()
        #self.MCaptureHunt()
        #if self.PredictedSignalPeptideFile:
        #    self.DebunkSignalP()
        
    def NMEHunt(self):
        """Now we go through all the peptide sequences, and see which look like NME (M.peptide start)
        We print off whether these are or are not the first peptide in their respective protein
        """
        for Aminos in self.AllPeptides.keys():
            PeptideObject = self.AllPeptides[Aminos]
            if not PeptideObject.Prefix == "M":
                continue
            #so it looks promising.  Let's see where it lands on the protein, and it it's 
            #we see if it mapps to the same place as the first reported amino acid in it's protein
            DBLocations = self.ProteinPicker.FindPeptideLocations(Aminos)
            for (ProteinID, PeptideStartAA) in DBLocations:
                ProteinName = self.ProteinPicker.ProteinNames[ProteinID]
                ThisProteinStart = self.ProteinsAndFirstObservedAA.get(ProteinName, None) #none as the default
                if not ThisProteinStart:
                    #this occurs where we did not record the first peptide start for a protein.  This is 
                    #perfectly reasonable, because it might have too few peptides
                    continue
                #so now we know what the first start is.  If we're the same, then we're the start
                ## if we're not, then were a look alike NME
                if PeptideStartAA > ThisProteinStart:
                    #look alike
                    print "%s\t%s\tInterior NME like"%(ProteinName, PeptideObject.GetFullModdedName())
                elif PeptideStartAA == ThisProteinStart:
                    #potential NME for real
                    if PeptideStartAA == 1:
                        print "%s\t%s\tTrue NME "%(ProteinName, PeptideObject.GetFullModdedName())
                    else:
                        print "%s\t%s\tPotential NME correction\t%s"%(ProteinName, PeptideObject.GetFullModdedName(), PeptideStartAA)
                else:
                    print "Snafu mapping %s\t%s\t"%(ProteinName, PeptideObject.GetFullModdedName())

    def MCaptureHunt(self):
        for Aminos in self.AllPeptides.keys():
            PeptideObject = self.AllPeptides[Aminos]
            if not Aminos[0] == "M":
                continue
            DBLocations = self.ProteinPicker.FindPeptideLocations(Aminos)
            for (ProteinID, PeptideStartAA) in DBLocations:
                ProteinName = self.ProteinPicker.ProteinNames[ProteinID]
                ThisProteinStart = self.ProteinsAndFirstObservedAA.get(ProteinName, None) #none as the default
                if not ThisProteinStart:
                    #this occurs where we did not record the first peptide start for a protein.  This is 
                    #perfectly reasonable, because it might have too few peptides
                    continue
                #so now we know what the first start is.  If we're the same, then we're the start
                ## if we're not, then were a look alike NME
                if PeptideStartAA > 0:
                    #see if it is a tryptic cut
                    if not PeptideObject.Prefix in ["R", "K"]:
                        print "%s\t%s\t Interior M-capture like"%(ProteinName, PeptideObject.GetFullModdedName())
                else :
                    print "%s\t%s\tTue M-capture "%(ProteinName, PeptideObject.GetFullModdedName())

    def DebunkSignalP(self):
        """This function takes as input a file with predictions from signal p (list of proteins)
        and then looks for proteins that we find evidence for in the first 25 amino acids of peptides
        from our proteomics data
        """
        SignalPPredictions = self.ParseSignalPredictionFile()
        #now we go through all our protiens and see if their first residue is within the first 25 base pairs
        ## thus disproving the call, we should really keep track of the position they said, but that's more work
        ## perhaps for this afternoon
        for (Name, Probability, PredictedStart) in SignalPPredictions:
            ThisProteinStart = self.ProteinsAndFirstObservedAA.get(Name, None)
            if not ThisProteinStart:
                continue #we've go nothing to comare with
            if ThisProteinStart < PredictedStart:
                PredictedMotif = self.GetMotifArea(Name, PredictedStart)
                print "%s\t%s\t%s\t0\t%s"%(Name, ThisProteinStart, PredictedStart, PredictedMotif)
            elif ThisProteinStart == PredictedStart:
                print "%s\t%s\t%s\t1"%(Name, ThisProteinStart,PredictedStart)


    def GetMotifArea(self, ProteinName, PlusOneResidue):
        """This function is to get the sequence surrounding the cleavage site of a protein
        we want -3, -2, -1, +1, +2  that's five letters
        """
        ID = self.NameToID[ProteinName]
        Sequence = self.ProteinPicker.ProteinSequences[ID]
        MotifStart = PlusOneResidue -3
        MotifStop = PlusOneResidue +2
        MotifArea = Sequence[MotifStart:MotifStop]
        return MotifArea
        
    def ParseSignalPredictionFile(self):
        """This function parses the signal p prediction file, which is a list of gi numbers.
        I just have to map those the the proteins in our database, and then return the IDs.
        gi|22125314|ref|NP_668737.1|    0.554    30
        gi|22127555|ref|NP_670978.1|    0.831    21
        gi|22127683|ref|NP_671106.1|    0.792    20
        gi|22125472|ref|NP_668895.1|    0.866    28
        ## warning this number is 1 based, not zero based like we often use.  that's the reason for
        subtracting zero
        """
        Names = []
        Handle = open(self.PredictedSignalPeptideFile, "rb")
        for Line in Handle.xreadlines():
            Line = Line.strip()
            (NameStub, Probability, Residue) = Line.split("\t") #and some other crap that should also match.
            FirstObservableResidue  = int(Residue) -1
            #print GiNum
            #now go through our database and find matches to the name
            for (ID, Name) in self.ProteinPicker.ProteinNames.items():
                Location = Name.find(NameStub)
                if not Location == -1: #it's found
                    #print "I found %s in %s"%(GiNum, Name)
                    Names.append((Name, Probability, FirstObservableResidue)) #swap in Name (not NameStub) for dictionary referencing later
                    continue
        return Names

    
    def SignalPeptideHunt(self):
        """look at the first observed peptide and see if we have evidence of 
        signal peptide cleavage"""
        HPlot = ProteinStatistics.HyrdopathyPlot()
        for ProteinName in self.ProteinsAndFirstObservedAA.keys():
            PeptideStart = self.ProteinsAndFirstObservedAA[ProteinName]
            FirstPeptideObject = self.FirstPeptideInProtein[ProteinName]
            ## we first filter by trypticness.  A signal peptide can't be tryptic,
            ## or at least we don't have confidence in the invivo cleavage if it appears tryptic
            if FirstPeptideObject.Prefix in ["R", "K"]:
                continue
            #now we expect signal peptides to fall within a certain length range
            MaxLen = 50
            MinLen = 15
            if PeptideStart > MinLen and PeptideStart < MaxLen:
                ID = self.NameToID[ProteinName]
                PrefixSequence = self.ProteinPicker.ProteinSequences[ID][:PeptideStart]
                AfterCut = self.ProteinPicker.ProteinSequences[ID][PeptideStart:PeptideStart+2]
                #now some trickery
                Plot = HPlot.MakePlot(PrefixSequence)
                LongerSignal = HPlot.IsConsistentSignal(Plot) #default to 10
                ShorterSignal = HPlot.IsConsistentSignal(Plot, 0.5, 8)
                AxBMotif = ProteinStatistics.HasSignalPeptidaseMotif(PrefixSequence)
                #I decided to go with the shorter hydrophobic patch length, but I kept the other var in there just incase
                String = "%s\t%s\t%s\t%s\t%s\t%s\t"%(ProteinName, PeptideStart, PrefixSequence, AfterCut, ShorterSignal, AxBMotif)
                PlotString =""
                # need to front pad these with zeros
                PadLen = MaxLen - len(Plot)
                for i in range(PadLen):
                    PlotString += "\t0"
                
                for Item in Plot:
                    PlotString += "\t%s"%Item
                print "%s%s"%(String, PlotString) #plotString has an explicit \t at start
                


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
            #okay, apparently we care about this protein.  Let's print stuff out
            self.ProteinsAndFirstObservedAA[ProteinName] = PeptideStartInProtein[ProteinName]
            self.FirstPeptideInProtein[ProteinName] = PeptideObjectInProtein[ProteinName]

        
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
        (Options, Args) = getopt.getopt(Arguments, "r:d:w:vp:s:")
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
            if Option == "-s":
                if not os.path.exists(Value):
                    print "** Error: couldn't find results file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.PredictedSignalPeptideFile = Value
                    
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
    