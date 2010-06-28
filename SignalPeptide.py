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
import BasicStats
import ProteinStatistics



class FinderClass():
    def __init__(self, OutputPath):
        self.PredictedSignalPeptideFile = None
        self.OutputPath = "RenameYourOutput.txt"
        (Path, Ext) = os.path.splitext(OutputPath)
        self.InfoPath = "%s.%s"%(Path, "signalpeptide.info")
        


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
        
    
    def EvaluateSignalPeptide(self, ORF):
        """Parameters: an openreadingframe object
        Return: true/false as for our call on whether it has a signal peptide or not
        look at the first observed peptide and see if we have evidence of 
        signal peptide cleavage
        """
        if ORF.numPeptides() == 0:
            return "NO EVIDENCE"

        HPlot = ProteinStatistics.HyrdopathyPlot()
        FirstObservedPeptide = ORF.GetFivePrimePeptide(Unique=1)
        # level 1 - is it non-tryptic on the N-terminus
        if FirstObservedPeptide.IsTrypticNterm():
            return "FAIL: N-term tryptic"
        #Level 2 - the length filter, which is somewhat stupid, because a 
        ##protein could be mispredicted, but it's what we have to go with 
        ##and we expect signal peptides to fall within a certain length range
        PeptideOffsetIntoProtein = None ### FIX THIS
        MaxLen = 50
        MinLen = 15
        if PeptideOffsetIntoProtein < MinLen:
            return "FAIL: too short"
        if PeptideOffsetIntoProtein > MaxLen:
            return "FAIL: too long"
        #meet's level 2
        PrefixSequence = ORF.GetProteinSequence(0, PeptideOffsetIntoProtein) # zero is the start.
        AfterCut = ORF.GetProteinSequence(PeptideOffsetIntoProtein, PeptideOffsetIntoProtein + 2)
        ProteinName = ORF.GetProteinName()
        #now some trickery
        Plot = HPlot.MakePlot(PrefixSequence)
        LongerSignal = HPlot.IsConsistentSignal(Plot) #default to 10
        ShorterSignal = HPlot.IsConsistentSignal(Plot, 0.5, 8)
        AxBMotif = ProteinStatistics.HasSignalPeptidaseMotif(PrefixSequence)
        #I decided to go with the shorter hydrophobic patch length, but I kept the other var in there just incase
        String = "%s\t%s\t%s\t%s\t%s\t%s\t"%(ProteinName, PeptideOffsetIntoProtein, PrefixSequence, AfterCut, ShorterSignal, AxBMotif)
        PlotString =""
        # need to front pad these with zeros
        PadLen = MaxLen - len(Plot)
        for i in range(PadLen):
            PlotString += "\t0"
        
        for Item in Plot:
            PlotString += "\t%s"%Item
        print "%s%s"%(String, PlotString) #plotString has an explicit \t at start
        
        return ShorterSignal #it has a hyrophobic patch or not
                
    