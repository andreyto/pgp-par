UsageInfo = """NotInDatabase.py
This is a quicky script to search for novel peptides

Required Options
 -r [Directory] Directory where the results files are kept.
 -d [Trie file] Database (can use multiple times)
 
Additional Options
 -w [FileName] File for writing out all peptides (best scoring example)
 -n [FileName] File for all novel peptides
 -m [FileName] File for all mundane (not novel) peptides

"""

import sys
import os
import getopt
import traceback
import InspectResults
import SelectProteins
from Utils import *
Initialize()


class CompileClass():
    def __init__(self):
        self.ReferenceResults = None
        self.AllPeptideOutPath = None
        self.NovelPeptideOutPath = None
        self.MundanePeptideOutPath = None
        self.ReferenceDatabasePath = [] #possibly multiple
        self.PeptideScores = {} # {Aminos} => pvalue
        self.PValueCutoff = None
        self.LFDRCutoff = None
        self.PeptideORFs = {} #aminos => ORFs
        

        
    def Main(self):
        self.ProteinPicker = SelectProteins.ProteinSelector()
        self.ProteinPicker.LoadMultipleDB(self.ReferenceDatabasePath)
        self.ParseInspect(self.ReferenceResults)
        self.FindInDatabase()

    def FindInDatabase(self):
        FoundCount = 0
        TotalPeptides =0
        if self.AllPeptideOutPath:
            AllOutHandle = open(self.AllPeptideOutPath, "wb")
        if self.NovelPeptideOutPath:
            NovelOutHandle = open(self.NovelPeptideOutPath, "wb")
        if self.MundanePeptideOutPath:
            MundaneOutHandle = open(self.MundanePeptideOutPath, "wb")
        print "Searching %d peptides"%len(self.PeptideScores)
        for Aminos in self.PeptideScores.keys():
            TotalPeptides += 1
            if TotalPeptides%1000 == 0:
                print "finished with %d peptides"%TotalPeptides
            Locations = self.ProteinPicker.FindPeptideLocations(Aminos)
            if self.AllPeptideOutPath:
                AllOutHandle.write("%s\n"%Aminos)
            if len(Locations) == 0:
                FoundCount += 1
                if self.NovelPeptideOutPath:
                    NovelOutHandle.write("%s\t%s\n"%(Aminos, self.PeptideORFs[Aminos]))
            else:
                if self.MundanePeptideOutPath:
                    MundaneOutHandle.write("%s\n"%Aminos)
                #print Aminos
        print "Out of %d total peptides, I found %s"%(TotalPeptides, FoundCount)
        if self.AllPeptideOutPath:
            AllOutHandle.close()
        if self.NovelPeptideOutPath:
            NovelOutHandle.close()
        if self.MundanePeptideOutPath:
            MundaneOutHandle.close()

    def ParseInspect(self, FilePath):
        inspectParser = InspectResults.Parser( FilePath )
        for result in inspectParser:
            try:
                Annotation = result.Annotation
                Peptide = GetPeptideFromModdedName(Annotation)
                Aminos = Peptide.Aminos
                PValue = result.PValue
                ORF = result.ProteinName
            except:
                traceback.print_exc()
                continue # SNAFU
            if not self.PeptideScores.has_key(Peptide.Aminos):
                self.PeptideScores[Peptide.Aminos] = PValue
                self.PeptideORFs[Peptide.Aminos] = ORF
            else:
                if PValue < self.PeptideScores[Peptide.Aminos]:
                    self.PeptideScores[Peptide.Aminos] = PValue


    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:d:w:n:m:p:l:")
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
            if Option == "-w":
                self.AllPeptideOutPath = Value
            if Option == "-n":
                self.NovelPeptideOutPath = Value
            if Option == "-m":
                self.MundanePeptideOutPath = Value
            #two secret undocumented options to limit stuff by pvalue, in case I need that option
            if Option == "-p":
                self.PValueCutoff = float(Value)
            if Option == "-l":
                self.LFDRCutoff = float(Value)
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