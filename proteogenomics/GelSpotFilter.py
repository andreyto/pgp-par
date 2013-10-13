import os
import sys
import getopt
import traceback
from Utils import *
import ResultsParser
import SelectProteins

UsageInfo = """GelSpotFilter.py
Given results files that correspond to a single gel spot, we filter results (already run through pvalue)
and require peptides correspond to proteins with multiple peptides in this single spot.

-r [Directory] Directory containing PValued Results
-w [Directory] Directory to write out Results
-d [TrieFile] Protein sequence database
-p [float] PValue cutoff. Default 0.05
"""


class BinMaker(ResultsParser.ResultsParser):

    def __init__(self):
        self.PvalueCutoff = 0.05
        self.FileNames = []
        self.InputDir = None
        self.OutputDir = None
        self.DatabasePaths = []
        ResultsParser.ResultsParser.__init__(self)

    def Main(self):
        self.ProteinPicker = SelectProteins.ProteinSelector()
        self.ProteinPicker.LoadMultipleDB(self.DatabasePaths)
        self.ProcessResultsFiles(self.InputDir, self.FirstPass)

    def FirstPass(self, FilePath):
        """Parameters: path to an inspect results file (representing a single gel spot)
        Return: none
        Description: map peptides on to proteins and then add them to the list
        """
        ProteinPeptides = {} #protein name=> peptidelist (you can do a len to get the count
        Handle = open(FilePath, "rb")
        for Line in Handle.xreadlines():
            if Line[0] == "#":
                continue # comment line
            
            if not Line.strip():
                continue
            Bits = list(Line.split("\t"))
            try:
                Annotation = Bits[self.Columns.Annotation]
                Peptide = GetPeptideFromModdedName(Annotation)
                Aminos = Peptide.Aminos
            except:
                traceback.print_exc()
                continue # SNAFU
            DBLocations = self.ProteinPicker.FindPeptideLocations(Aminos)
            for (ID, Residue) in DBLocations:
                if not ProteinPeptides.has_key(ID):
                    ProteinPeptides[ID] = []
                #now test if the amino sequence is in the list of results
                if not Aminos in ProteinPeptides[ID]:
                    ProteinPeptides[ID].append(Aminos)
                    
        #now we're done running through a file.  Let's pass this off to the second pass
        CountDict= {}
        for (ID, List) in ProteinPeptides.items():
            CountDict[ID] = len(List)
        self.SecondPass(FilePath, CountDict)
        
    def SecondPass(self, FilePath, ProteinPeptideCount):
        """Parameters: path to an inspect results file , count dictionary
        Return: none
        Description: filter result lines to those with at least two peptides/protein/spot
        """
        Handle = open(FilePath, "rb")
        (Path, FileName) = os.path.split(FilePath)
        OutPath = os.path.join(self.OutputDir, FileName)
        OutHandle = open(OutPath, "wb")
        for Line in Handle.xreadlines():
            if Line[0] == "#":
                OutHandle.write(Line)                
                continue # comment line
            
            if not Line.strip():
                continue
            Bits = list(Line.split("\t"))
            try:
                Annotation = Bits[self.Columns.Annotation]
                Peptide = GetPeptideFromModdedName(Annotation)
                Aminos = Peptide.Aminos
            except:
                traceback.print_exc()
                continue # SNAFU
            #here's some hacking that needs to be fixed.  I currently want to filter stuff by lfdr, but
            #that may not always be present
            if len(Bits) < self.Columns.LFDR: #meaning that I don't have the column in question
                if PValue > self.PvalueCutoff: #substitute pvalue for LFDR
                    continue
            else:
                LFDR = float(Bits[self.Columns.LFDR])
                if LFDR > self.PvalueCutoff:
                    continue
            #### end hack
            DBLocations = self.ProteinPicker.FindPeptideLocations(Aminos)
            Keep = 0
            for (ID, Residue) in DBLocations:
                if ProteinPeptideCount[ID] >= 2:
                    Keep = 1
                    break
            if Keep:
                OutHandle.write(Line)
        Handle.close()
        OutHandle.close()

    def ParseCommandLine(self,Arguments):
        (Options,Args) = getopt.getopt(Arguments, "r:w:d:p:")
        OptionsSeen = {}
        for(Option, Value) in Options:
            OptionsSeen[Option] = 1

            if Option == "-r":
                if not os.path.isdir(Value):
                    print "**Error: %s is not a valid directory"%Value
                    print UsageInfo
                    sys.exit(1)
                self.InputDir = Value
            elif Option == "-w":
                if not os.path.isdir(Value):
                    os.mkdir(Value)
                self.OutputDir = Value
            elif Option == "-d":
                self.DatabasePaths.append(Value)
            elif Option == "-p":
                self.PvalueCutoff = float(Value)
            else:
                print "ERROR: Argument %s is not known. Program Exit"%Option
                print UsageInfo
                sys.exit(1)

        if not OptionsSeen.has_key("-r") or not OptionsSeen.has_key("-w") or not OptionsSeen.has_key("-d"):
            print "**Error: required options -r, -w, -d"
            print UsageInfo
            sys.exit(1)
                    


if __name__ == "__main__":
    Dummy = BinMaker()
    Dummy.ParseCommandLine(sys.argv[1:])
    Dummy.Main()
