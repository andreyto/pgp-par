"""GetGrossGenomeCoverage.py
This script takes in a set of annotation for some spectra.  And we comput the
genome coverage, meaning of the XXX residues in the input database, we get some
annotation for XXX of them.  For this we count a peptide every place that it occurs
in the proteome.

Maybe I will print it out gene by gene too.
"""
UsageInfo = """GetGrossGenomeCoverage.py
Takes annotations from a searches of a dataset and
computes the genome coverage.

Required Options
 -r [Directory] Directory where the reference results files are kept.
 -d [Trie file] Database


Additional Options

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
        self.SetUpCoverageArrays()
        self.SetCoverage()
        self.WriteResults()



    def WriteResults(self):
        """On the gene level, and the whole genome 
        """
        TotalAminos = 0
        TotalCoveredAminos = 0
        Histogram = {}
        for I in range (11):
            Histogram[I] = 0
        for ID in self.ProteinCoverageHash.keys():
            CoverageArray = self.ProteinCoverageHash[ID]
            Length = len(CoverageArray)
            Covered = 0.0
            for Index in range(Length):
                if CoverageArray[Index] > 0:
                    Covered += 1
            CoverageThisProtein = Covered/ Length
            TotalAminos += Length
            TotalCoveredAminos += Covered
            Percent = int (CoverageThisProtein * 100)
            Bin = Percent / 10 # 0-10, 11-20, ...
            Histogram[Bin] += 1
            if Covered > 0:
                print "Coverage for %s is %s"%(self.ProteinPicker.ProteinNames[ID][:12], CoverageThisProtein)
                print "\t%d unique peptides and %d shared peptides"%(self.UniqueProteinSequences[ID], self.SharedProteinSequences[ID])
        print "Residues covered / Total  %s / %s "%(TotalCoveredAminos, TotalAminos)
        for I in range(11):
            print ".%d\t%d"%(I, Histogram[I])
        self.PlotUniqueAndShared()
        self.PlotSpectrumCount()

    def PlotSpectrumCount(self):
        """ just print out the number of spectra for each protein, we'll sort and do stuff later"""
        print "\nSpectrumCount By Protein"
        for (ID, Count) in self.ProteinSpectra.items():
            UniquePeptideCount = self.UniqueProteinSequences[ID]
            SharedPeptideCount = self.SharedProteinSequences[ID]
            if UniquePeptideCount == 0:
                #no unique peptides, let's not print this one out.
                continue
            if SharedPeptideCount + UniquePeptideCount < 2:
                continue
            print "%s\t%s"%(self.ProteinPicker.ProteinNames[ID], Count)

    def PlotUniqueAndShared(self):
        """Make some histogram plots of the unique and shared peptides
        """
        UniqueHistogram = {}
        SharedHistogram = {}
        SharedUZero = {}
        SharedUOne = {}
        SharedUTwoFour = {}
        SharedUFivePlus = {}
        for ID in self.ProteinCoverageHash.keys():
            Shared = self.SharedProteinSequences[ID]
            Unique = self.UniqueProteinSequences[ID]
            if Shared + Unique == 0:
                ## let's just not include these in our display
                continue
            ## simple histograms
            if not UniqueHistogram.has_key(Unique):
                UniqueHistogram[Unique] = 0
            UniqueHistogram[Unique] += 1
            if not SharedHistogram.has_key(Shared):
                SharedHistogram[Shared] = 0
            SharedHistogram[Shared] += 1
            ## complex, ANOVAish
            if Unique == 0:
                if not SharedUZero.has_key(Shared):
                    SharedUZero[Shared] = 0
                SharedUZero[Shared] += 1
            elif Unique == 1:
                if not SharedUOne.has_key(Shared):
                    SharedUOne[Shared] = 0
                SharedUOne[Shared] += 1
            elif Unique < 5:
                if not SharedUTwoFour.has_key(Shared):
                    SharedUTwoFour[Shared] = 0
                SharedUTwoFour[Shared] += 1
            else:
                if not SharedUFivePlus.has_key(Shared):
                    SharedUFivePlus[Shared] = 0
                SharedUFivePlus[Shared] += 1
        # now print it out
        Max = self.GetMaxKey(SharedHistogram, UniqueHistogram)
        print "Key\tUnique\tShared\tShared U=0\tShared U=1\tShared U=2-4\tShared U=5+"
        for Key in range(Max + 1):
            U = UniqueHistogram.get(Key, 0)
            S = SharedHistogram.get(Key, 0)
            S0 = SharedUZero.get(Key, 0)
            S1 = SharedUOne.get(Key, 0)
            S2 = SharedUTwoFour.get(Key, 0)
            S5 = SharedUFivePlus.get(Key, 0)
            print "%s\t%s\t%s\t%s\t%s\t%s\t%s"%(Key, U, S, S0, S1, S2, S5)


    def GetMaxKey(self, Hash1, Hash2):
        Max = 0
        for Key in Hash1.keys():
            if Key > Max:
                Max = Key
        for Key in Hash2.keys():
            if Key > Max:
                Max = Key
        return Max
                

    def SetUpCoverageArrays(self):
        """start the self.ProteinCoverageHash
        """
        for ID in self.ProteinPicker.ProteinNames.keys():
            self.ProteinCoverageHash[ID] = [0]*len(self.ProteinPicker.ProteinSequences[ID])
            self.UniqueProteinSequences[ID] = 0
            self.SharedProteinSequences[ID] = 0
            self.ProteinSpectra[ID] = 0
        

    def SetCoverage(self):
        print "setting coverage for %d peptides"%len(self.AllReferenceAnnotations)
        for Aminos in self.AllReferenceAnnotations.keys():
            #1. Get all the locations of these aminos in the proteinpicker.db
            Locations = self.ProteinPicker.FindPeptideLocations(Aminos)
            Length =  len(Aminos)
            NumSpectra = self.AllReferenceAnnotations[Aminos]
            #print "%s has %d locations"%(Aminos, len(Locations))
            #2. write to the coverage arrays
            for (ID, FirstResidue) in Locations: # (ID, Residue Number)
                self.ProteinSpectra[ID] += NumSpectra
                if len(Locations) > 1:
                    self.SharedProteinSequences[ID] += 1
                else:
                    self.UniqueProteinSequences[ID] += 1
                for Index in range(FirstResidue, FirstResidue + Length):
                    self.ProteinCoverageHash[ID][Index] +=1
                    

    def ParseReferenceFileCallback(self, FilePath):
        """Called by the ResultsParser.ProcessResultsFiles
        Here I parse out lines and get annotations, Putting them in a list
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
                Score = float(Bits[self.Columns.MQScore])
            except:
                traceback.print_exc()
                continue # SNAFU
            #print Spectrum
            if not self.AllReferenceAnnotations.has_key(Aminos):
                self.AllReferenceAnnotations[Aminos] = 0
            self.AllReferenceAnnotations[Aminos] += 1  # keep track of spectrum count, if I want to use it
                
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