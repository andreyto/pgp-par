"""GatherSpectraParis.py
This  program takes the output of Inspect (Pvalued and protein assigned),
and gets the list of peptides that need to be reported according to the
Paris guidelines.
1. Any annotation with a PTM
2. Any annotation for a protein which is a one hit wonder

With this list, it then calls LabelItAll, and creates all the images.

"""

UsageInfo = """GatherSpectraParis.py
This  program takes the output of Inspect (Pvalued and protein assigned),
and gets the list of peptides that need to be reported according to the
Paris guidelines.

Required Options
 -r [FileName] - The Inspect output file or directory
 -w [FileName] - The output of this program, a subset of the input
 -D [Directory] - Directory to store all the labeled output
 -m [Directory] - Directory for the spectra
"""

import os
import getopt
import sys
import string
import ResultsParser
import LabelItAll
from Utils import *
Initialize()

class ParisClass(ResultsParser.ResultsParser):
    def __init__(self):
        self.LabelSpewageDir = None
        self.SpectraDir = None
        self.InspectFile = None
        self.OutputFile = None
        self.DatabasePath = None
        self.ProteinPeptide = {} # Protein => (Pep1, Pep2, ..)
        self.BestScores = {} #annotation > score
        self.OneHitWonders = {}
        self.OutHandle = None
        self.Header = None
        self.SavedID = None
        ResultsParser.ResultsParser.__init__(self)

    def ParseForParis(self, FilePath):
        """In each file I am going to grep out
        filename, scan number, annotation
        """
        
        ##in the filename are the scan number and mzxml filename
        Handle= open(FilePath, "rb")
        for Line in Handle.xreadlines():
            #print Line
            Line = Line.strip()
            if Line[0] == "#":
                self.Header = Line
                continue
            #parse crap
            Bits = Line.split("\t")
            SpectrumFile = Bits[self.Columns.SpectrumFile]
            ScanNumber = Bits[self.Columns.ScanNumber]
            ID = (SpectrumFile, ScanNumber)
            Annotation = Bits[self.Columns.Annotation]
            Protein = Bits[self.Columns.ProteinName]
            Score = float(Bits[self.Columns.MQScore])
            # add everybody to the best scores hash
            if not self.BestScores.has_key(Annotation):
                self.BestScores[Annotation] = Score
            else:
                CurrScore = self.BestScores[Annotation]
                if Score > CurrScore:
                    self.BestScores[Annotation] = Score
            # get annotations fixed onto a protein
            self.AddPeptideToProtein(Protein, Annotation, Score)
        Handle.close()

    def CheckForMods(self, Annotation):
        "check for non caps"
        Found = 0
        if Annotation[1] == ".":
            Annotation = Annotation[2:-2]
        for Index in range(len(Annotation)):
            Letter = Annotation[Index]
            if not Letter in string.uppercase:
                Found = 1
                break
        return Found


    def AddPeptideToProtein(self, Protein, Annotation, Score):
        """Note that here we are only using aminos, because
        the paris guidelines say that SAM+16PAYNE and SAMPAYNE
        count as a single peptide for the purpose of counting
        peptides in a protein
        """
        Peptide = GetPeptideFromModdedName(Annotation)
        Aminos = Peptide.Aminos
        if not self.ProteinPeptide.has_key(Protein):
            self.ProteinPeptide[Protein] = []
        AllPeptides = self.ProteinPeptide[Protein]
        Found = 0
        for Peptide in AllPeptides:
            if Aminos == Peptide:
                Found = 1
                break
        if not Found:
            self.ProteinPeptide[Protein].append(Aminos)
                
            
    def WriteResults(self, FilePath):
        """ We are going to reparse the output file
        and print to the output file results that
        we feel like keeping.
        """
        Handle = open(FilePath, "rb")
        for Line in Handle.xreadlines():
            if not Line.strip():
                continue
            Line = Line.strip() #assigning. 
            if Line[0] == "#":
                continue
            WriteLine = 0 
            Bits = Line.split("\t")
            SpectrumFullPath = Bits[self.Columns.SpectrumFile]
            Scan = Bits[self.Columns.ScanNumber]
            ProteinName = Bits[self.Columns.ProteinName]
            Annotation = Bits[self.Columns.Annotation]
            MQScore = float(Bits[self.Columns.MQScore])
            ## 1. check for mod
            if self.CheckForMods(Annotation):
                ## is this the best score for this peptide.
                if MQScore == self.BestScores[Annotation]:
                    WriteLine = 1
            ## 2. check for one hit wonder
            PeptideObject = GetPeptideFromModdedName(Annotation)
            Aminos = PeptideObject.Aminos
            if self.OneHitWonders.has_key((ProteinName, Aminos)):
                if MQScore == self.BestScores[Annotation]:
                    WriteLine = 1
            if WriteLine:
                String = "\t".join(Bits)
                self.OutHandle.write("%s\n"%String)
        Handle.close()
            
    def FindOneHitWonders(self):
        for (Protein, PeptideList) in self.ProteinPeptide.items():
            if len(PeptideList) == 1:
                Peptide = PeptideList[0]
                self.OneHitWonders[(Protein, Peptide)] = 1
        
    def Main(self):
        self.ProcessResultsFiles(self.InspectFile, self.ParseForParis)
        self.FindOneHitWonders()
        self.OutHandle = open(self.OutputFile, "wb")
        self.OutHandle.write("%s\n"%self.Header)
        self.ProcessResultsFiles(self.InspectFile, self.WriteResults)
        self.OutHandle.close()
        self.GenerateImages(self.OutputFile)

    def GenerateImages(self, FileName):
        """This function wraps the LabelItAll script, 
        """
        Args = " -r %s -w %s -m %s"%(FileName,self.LabelSpewageDir, self.SpectraDir)
        ArgsList = Args.split()
        #print "Parsing Results for %s, scan %s, charge %s"%(FileName, Scan, Charge)
        #print "Args: %s"%Args
        Dymo = LabelItAll.DetectiveClass()
        Dymo.ParseCommandLine(ArgsList)
        Dymo.Main()        
        
        

    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:w:D:m:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-r":
                # -r results file(s)
                if not os.path.exists(Value):
                    print "** Error: couldn't find results file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.InspectFile = Value
            elif Option == "-w":
                self.OutputFile = Value
            elif Option == "-D":
                self.LabelSpewageDir = Value
            elif Option == "-m":
                self.SpectraDir = Value
            else:
                print "Unknown option %s"%Option
                print UsageInfo
                sys.exit(1)
        if not OptionsSeen.has_key("-r") or not OptionsSeen.has_key("-w") or not OptionsSeen.has_key("-D") or not OptionsSeen.has_key("-m") :
            print "**********************\nRequired options missing\n***********************"
            print UsageInfo
            sys.exit(1)
                
if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(psyco not found - running in non-optimized mode)"
    Fix = ParisClass()
    Fix.ParseCommandLine(sys.argv[1:])
    Fix.Main()        