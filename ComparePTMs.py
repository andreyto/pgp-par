"""
Perform relative quantitation of post-translational modifications between
samples, based (for now, at least!) upon spectrum counts.

At present, this analysis seems to generate mostly noise: The modifications
with the best p-values are all chemical damage (e.g. pyroglutamate formation)
where any variation that DOES exist between samples is probably not
biologically relevant.  We need a sample with a known-valid difference
in modification states to be certain of our techniques, but...it appears that
spectrum counts may not be the way to go for differential modification rates.
Note that spectrum counts are discontinuous, that we have on average only ONE
peptide to measure (as compared with protein quantitation where we can count
many peptides), and that peptides may not be fragmented and identified at the
same rate in both samples.  If we had differential labeling (via ICAT or
different isotopes), we could compare signal strengths, and that might perform
better (note again, we require a known sample to test the methods)

"""
import os
import sys
import traceback
import getopt
import struct
import ChiSquare
from Utils import *
from TrainPTMFeatures import FormatBits
Initialize()

class PTMQuantifier:
    def __init__(self):
        self.InputFileName = None
        self.ClusterRootDir = None
        self.OutputFileName = "PTMQuantitation.txt"
    def LoadCoverage(self, Coverage, ModdedCoverage, FilePath):
        File = open(FilePath, "rb")
        BlockSize = struct.calcsize("<II")
        while 1:
            Str = File.read(BlockSize)
            if not Str:
                break
            (Count, CountModded) = struct.unpack("<II", Str)
            Coverage.append(Count)
            ModdedCoverage.append(CountModded)
            #self.CoverageA.append(Coverage)
            #self.ModdedCoverageA.append(Coverage)
        File.close()
    def Main(self):
        """
        Parse coverage levels for each sample.
        Parse the table of modded peptides.
        For each modded peptide, read the cluster-members file.  Keep track of
        how many cluster-members come from control (spectrum file name contains NameA)
        and treatment (NameB) consditions.  Determine whether the modificaiton rate
        is significantly different.
        """
        ###########################################
        # Load all coverage:
        print "Load coverage..."
        self.CoverageA = []
        self.ModdedCoverageA = []
        self.LoadCoverage(self.CoverageA, self.ModdedCoverageA, self.CoveragePathA)
        self.CoverageB = []
        self.ModdedCoverageB = []
        self.LoadCoverage(self.CoverageB, self.ModdedCoverageB, self.CoveragePathB)
        ###########################################
        # Iterate over peptides, and compare modification rates:
        self.OutputFile = open(self.OutputFileName, "wb")
        HeaderString = "DBPos\tModMass\tPeptideExample\tSpeciesCount\tModdedA\tModlessA\tRateA\tModdedB\tModlessB\tRateB\tChi\tPValue\tNote\t\n"
        self.OutputFile.write(HeaderString)
        File = open(self.InputFileName, "rb")
        CurrentSite = None
        CurrentSpeciesList = []
        LineNumber = 0
        for FileLine in File.xreadlines():
            LineNumber += 1
            if LineNumber % 1000 == 0:
                print "Line %s..."%LineNumber
            if FileLine[0] == "#":
                continue
            Bits = FileLine.split("\t")
            if len(Bits) < 2:
                continue
            try:
                DBPos = int(Bits[FormatBits.DBPos])
                ModMass = int(Bits[FormatBits.ModificationMass])
            except:
                traceback.print_exc()
                print Bits
                continue
            Site = (DBPos, ModMass)
            if Site != CurrentSite:
                self.ComparePTMRate(CurrentSpeciesList)
                CurrentSpeciesList = []
                CurrentSite = Site
            CurrentSpeciesList.append(Bits)
        File.close()
        self.ComparePTMRate(CurrentSpeciesList)
    def ComparePTMRate(self, CurrentSpeciesList):
        if not CurrentSpeciesList:
            return
        Bits = CurrentSpeciesList[0]
        DBPos = int(Bits[FormatBits.DBPos])
        ModMass = int(Bits[FormatBits.ModificationMass])
        # Get the UNMODIFIED coverage:
        ModlessCoverageA = self.CoverageA[DBPos]
        ModlessCoverageB = self.CoverageB[DBPos]
        # Get the MODIFIED coverage:
        ModdedCoverageA = 0
        ModdedCoverageB = 0        
        for Bits in CurrentSpeciesList:
            Annotation = Bits[FormatBits.Peptide]
            Charge = int(Bits[FormatBits.Charge])
            Peptide = GetPeptideFromModdedName(Annotation)
            # Find the cluster-members file:
            MemberFileName = "%s.%s.txt"%(Annotation, Charge)
            Path = os.path.join(self.ClusterMembersAdjusted, Peptide.Aminos[0], MemberFileName)
            if not os.path.exists(Path):
                Path = os.path.join(self.ClusterMembersDir, Peptide.Aminos[0], MemberFileName)
                if not os.path.exists(Path):
                    print "** Error: Couldn't find cluster members at %s"%Path
                    return
            ##############################################
            # Parse the cluster-members file:
            File = open(Path, "rb")
            for FileLine in File.xreadlines():
                if FileLine.find(self.NameA) != -1:
                    ModdedCoverageA += 1
                if FileLine.find(self.NameB) != -1:
                    ModdedCoverageB += 1
            File.close()
        ##############################################
        # Chi-square test:
        if ModdedCoverageA + ModlessCoverageA + ModdedCoverageB + ModlessCoverageB == 0:
            Chi = ""
            PValue = ""
            Notes = "C"
        if ModdedCoverageA + ModlessCoverageA == 0:
            Chi = ""
            PValue = ""
            Notes = "B"
        elif ModdedCoverageB + ModlessCoverageB == 0:
            Chi = ""
            PValue = ""
            Notes = "A"
        else:
            Chi = ChiSquare.ComputeChiSquare(ModdedCoverageA, ModlessCoverageA, ModdedCoverageB, ModlessCoverageB)
            PValue = ChiSquare.GetPValue(Chi)
            Notes = ""
        # Verbose output:
        Str = "%s\t%s\t%s\t%s\t"%(DBPos, ModMass, Annotation, len(CurrentSpeciesList))
        RateA = 100 * ModdedCoverageA / float(max(1, ModdedCoverageA + ModlessCoverageA))
        Str += "%s\t%s\t%s\t"%(ModdedCoverageA, ModlessCoverageA, RateA)
        RateB = 100 * ModdedCoverageB / float(max(1, ModdedCoverageB + ModlessCoverageB))
        Str += "%s\t%s\t%s\t"%(ModdedCoverageB, ModlessCoverageB, RateB)
        Str += "%s\t%s\t%s\t"%(Chi, PValue, Notes)
        self.OutputFile.write(Str + "\n")
        print Str
    def ParseCommandLine(self, Arguments):
        (Options, Args) = getopt.getopt(Arguments, "C:c:S:s:r:d:w:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-r":
                # -r results file
                self.InputFileName = Value
                if not self.ClusterRootDir:
                    self.ClusterRootDir = os.path.split(Value)[0]
            elif Option == "-C":
                self.CoveragePathA = Value
            elif Option == "-c":
                self.CoveragePathB = Value
            elif Option == "-S":
                self.NameA = Value
            elif Option == "-s":
                self.NameB = Value
            elif Option == "-d":
                self.ClusterRootDir = Value
            elif Option == "-w":
                self.OutputFileName = Value
                
        self.ClusterMembersAdjusted = os.path.join(self.ClusterRootDir, "ClusterMembersAdjusted")
        self.ClusterMembersDir = os.path.join(self.ClusterRootDir, "ClusterMembers")
if __name__ == "__main__":
    Quan = PTMQuantifier()
    Quan.ParseCommandLine(sys.argv[1:])
    Quan.Main()