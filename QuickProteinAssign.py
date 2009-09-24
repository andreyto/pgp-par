UsageInfo = """QuickProteinAssign.py
Given a set of annotation for some spectra, we map peptides
onto proteins, and then output the proteins, with peptide counts.
This is call the Quick protein assignment, because it makes
NO attempt at parsimony or any other protein inference.  It
merely takes the annotated proteins as given. As such it does
not require a database.  It relies on the Inspect identifications.

Required Options
 -r [FileName] File or directory of results.
 -w [FileName]  Output file
 
Required Options
 -v   Flag to print a more verbose trace of the program progression
     through various methods and whatnot
 -p [float] Pvalue cutoff (Default 0.05)

"""

import sys
import os
import getopt
import traceback
import ResultsParser
from Utils import *
Initialize()


class FinderClass(ResultsParser.ResultsParser):
    def __init__(self):
        self.ReferenceResults = None
        self.OutputPath = "RenameYourOutput.txt"
        self.AllPeptides = {} # AminoSequence -> spectrumCount, nulled out in MapAllPeptides
        self.ProteinsWithTheirPeptides = {} # proteinName -> list of peptide sequences
        self.SpectrumCount = 0
        self.PValueLimit = 0.05 #pvalues are 0.00 (good) to 1.0 (bad) 
        self.Verbose = 0
        ResultsParser.ResultsParser.__init__(self)
        
    def Main(self):
        self.ProcessResultsFiles(self.ReferenceResults, self.ParseInspectCallback)
        print "I found %s peptides from %s spectra in %s proteins"%(len(self.AllPeptides), self.SpectrumCount, len(self.ProteinsWithTheirPeptides.keys()))
        ##now map all the peptides
        self.PrettyPrint()
        
    def PrettyPrint(self):
        """This function is supposed to print out proteins, and their associated peptides
        """
        Handle = open (self.OutputPath, "wb")
        for (Protein, PeptideList) in self.ProteinsWithTheirPeptides.items():
            String = "%s\t"%Protein
            String += "\t".join(PeptideList)
            String += "\t%s\n"%len(PeptideList)
            #print String
            Handle.write(String)
        Handle.close()
        
        
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
                AminosWithMod = Peptide.GetModdedName()
                PValue = float(Bits[self.Columns.PValue])
                Protein = Bits[self.Columns.ProteinName]
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
            # whew ! we made it through all kinds of checks, let's store the data
            if not self.AllPeptides.has_key(AminosWithMod):
                self.AllPeptides[AminosWithMod] = 0
            self.AllPeptides[AminosWithMod] += 1
            self.SpectrumCount += 1
            if not self.ProteinsWithTheirPeptides.has_key(Protein):
                self.ProteinsWithTheirPeptides[Protein] = []
            self.ProteinsWithTheirPeptides[Protein].append(AminosWithMod)
        Handle.close()


    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:w:vp:")
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
            if Option == "-w":
                self.OutputPath = Value
            if Option == "-p":
                #undocumented pvalue filter
                self.PValueLimit = float(Value)
            if Option == "-v":
                self.Verbose = 1
        if not OptionsSeen.has_key("-r")  or not OptionsSeen.has_key("-w"):
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