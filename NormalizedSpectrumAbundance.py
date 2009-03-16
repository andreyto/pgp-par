"""NormalizedSpectrumAbundance.py
Given a dataset of Inspect results, we calculate the normalized spectrum
abundance factor for every protein in the dataset.  This is a form of
relative quantitation.  The original reference for this metric is
Zybailov B, Mosley A L , Sardiu M E, Coleman M K, Florens L, Washburn M P.
Statistical Analysis of Membrance Proteome Expression Changes in
Saccharomyces cerevisiae. J Prot Research 2006. 5:2339-2347.  The basic
equation is

NSAF = (SpectrumCount_k/Length_k)/SUM_i: (SpectrumCount_i / Length_i)

Perhaps I will also do something fun like replacing SpectrumCount with PeptideCount.

"""

UsageInfo="""NormalizedSpectrumAbundance.py
Calculate all the NSAF for each protein in the proteome.

Required Options:
 -r [FileName] File or directory of filtered Inspect results
 -d [TrieFile] Database of protein sequences, use multiple times if needed
 -w [FileName] Output from program


"""

import os
import sys
import getopt
import SelectProteins
import ResultsParser
from Utils import *
Initialize()

class AbacusClass(ResultsParser.ResultsParser):
    def __init__ (self):
        self.DBPath = []
        self.OutputPath = None
        self.InspectInput = None
        self.Denominator = None # the denominator for the equation
        self.SpectrumCounts = {} # proteinID -> SpectrumCount  ## prehaps extend to hash of peptide and count
        ResultsParser.ResultsParser.__init__(self)
    def Main(self):
        self.ProteinPicker = SelectProteins.ProteinSelector()
        self.ProteinPicker.LoadMultipleDB(self.DBPath)
        self.ProcessResultsFiles(self.InspectInput, self.ParseResults)
        self.CalculateDenominator()
        self.CalculateNSAF()

    def CalculateDenominator(self):
        """Run through the calculation to sum over all proteins"""
        Sum = 0.0
        for ID in self.ProteinPicker.ProteinSequences.keys():
            Length = float(len(self.ProteinPicker.ProteinSequences[ID])) #cast to float for division
            SpectrumCount = self.SpectrumCounts.get(ID, 0)
            Sum += (SpectrumCount / Length)
            
        self.Denominator = Sum

    def CalculateNSAF(self):
        Handle = open(self.OutputPath, "wb")
        for ID in self.ProteinPicker.ProteinSequences.keys():
            Length = float(len(self.ProteinPicker.ProteinSequences[ID]))
            SpectrumCount = self.SpectrumCounts.get(ID, 0)
            Name = self.ProteinPicker.ProteinNames[ID]
            NSAF = (SpectrumCount/Length) / self.Denominator
            Handle.write("%s\t%s\n"%(Name, NSAF))
        Handle.close()
                         

    def ParseResults(self, FileName):
        Handle = open(FileName, "rb")
        for Line in Handle.xreadlines():
            if Line[0] == "#":
                continue
            Bits = Line.strip().split("\t")
            Annotation = Bits[self.Columns.Annotation]
            Peptide = GetPeptideFromModdedName(Annotation)
            Aminos = Peptide.Aminos
            ProteinList = self.ProteinPicker.FindPeptideLocations(Aminos)
            for (ProteinID, Residue) in ProteinList:
                if not self.SpectrumCounts.has_key(ProteinID):
                    self.SpectrumCounts[ProteinID] = 0
                self.SpectrumCounts[ProteinID] += 1
                #print "Added %s to %s:%s"%(Aminos, ProteinID, self.ProteinPicker.ProteinNames[ProteinID])
                
    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "d:w:r:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-d":
                # -r results file(s)
                if not os.path.exists(Value):
                    print "** Error: couldn't find database file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.DBPath.append(Value)
            if Option == "-w":
                # -r results file(s)
                self.OutputPath = Value
            if Option == "-r":
                self.InspectInput = Value
        if not OptionsSeen.has_key("-d") or not OptionsSeen.has_key("-w") or not OptionsSeen.has_key("-r"):
            print UsageInfo
            sys.exit(1)


if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(psyco not found - running in non-optimized mode)"
    DoStuff = AbacusClass()
    DoStuff.ParseCommandLine(sys.argv[1:])
    DoStuff.Main()        
