"""
Simple summary script: Count the PRMs found in data-sets, by type.
"""
import os
import sys
import random
random.seed(1)
import types
import traceback
from TrainPTMFeatures import FormatBits
from Utils import *
Initialize()

class PTMCounter:
    def __init__(self):
        self.MassPeptideCounts = {}
        self.MassAAPeptideCounts = {}
        self.MassSiteCounts = {}
        self.MassAASiteCounts = {}
        self.SitePValueCutoff = 0.1
        self.ReconcileFlag = 1
        self.FDRCutoff = 0.05
    def FindPValueCutoff(self, FileName):
        File = open(FileName, "rb")
        SortedList = []
        for FileLine in File.xreadlines():
            if FileLine[0] == "#":
                continue
            Bits = list(FileLine.strip().split("\t"))
            PValue = float(Bits[FormatBits.SitePValue])
            try:
                KnownPValue = float(Bits[FormatBits.KnownPTMSitePValue])
            except:
                KnownPValue = 999
            PValue = min(PValue, KnownPValue)
            if Bits[FormatBits.ProteinName][:3] == "XXX":
                ValidFlag = 0
            else:
                ValidFlag = 1
            SortedList.append((PValue, random.random(), ValidFlag))
        Cutoff = 0.0
        SortedList.sort()
        TPCount = 0
        FPCount = 0
        FPScalingFactor = 1.0
        self.FDRByPValue = {}
        Name = os.path.split(FileName)[1]
##        if Name == "ST1LensSitesFullK.txt":
##            FPScalingFactor = 1.0 / 99.0
        for (PValue, Dummy, ValidFlag) in SortedList:
            if ValidFlag:
                TPCount += 1
            else:
                FPCount += 1
            FDR = min(1.0, (FPCount * FPScalingFactor) / TPCount)
            if FDR <= self.FDRCutoff:
                Cutoff = PValue
            self.FDRByPValue[PValue] = FDR
        print "FDR %s%% cutoff is %s"%(100 * self.FDRCutoff, Cutoff)
        self.SitePValueCutoff = Cutoff
    def ParseFile(self, FileName, OutputFileName):
        """
        Parse a PTM table.  Ignore known-unmodified rows.  Add the rest of the
        rows to our cumulative dictionaries.
        """
        File = open(FileName, "rb")
        OutFile = open(OutputFileName, "wb")
        PrevSite = None
        LineNumber = 0
        SiteCount = 0
        PeptideCount = 0
        OKSiteCount = 0
        OKPeptideCount = 0
        for FileLine in File.xreadlines():
            OKFlag = 1
            if FileLine[0] == "#":
                OutFile.write(FileLine)
                continue
            LineNumber += 1
            Bits = list(FileLine.strip().split("\t"))
            if len(Bits) <= FormatBits.SitePValue:
                OutFile.write(FileLine)
                continue
            if len(Bits) < 51 or not Bits[50]:
                continue # this site is considered to be unmodified!
            PValue = float(Bits[FormatBits.SitePValue])
            ProteinName = Bits[FormatBits.ProteinName]
            if ProteinName[:3] == "XXX":
                OKFlag = 0
            try:
                KnownPValue = float(Bits[FormatBits.KnownPTMSitePValue])
            except:
                KnownPValue = 999
            PValue = min(PValue, KnownPValue)
            if PValue > self.SitePValueCutoff:
                OKFlag = 0
            if self.ReconcileFlag:
                if len(Bits) < 51 or not Bits[50].strip():
                    continue # unmodified peptide!
                ModifiedAA = Bits[50]
                Mass = int(Bits[51])
            else:
                Mass = int(Bits[FormatBits.ModificationMass])
                ModifiedAA = Bits[FormatBits.ModifiedAA]
            Site = (Bits[FormatBits.ModificationMass], Bits[FormatBits.ModifiedResidueNumber])
            Key = (Mass, ModifiedAA)
            if Site != PrevSite:
                # It's a new site; update the site counts
                SiteCount += 1
                if OKFlag:
                    self.MassSiteCounts[Mass] = self.MassSiteCounts.get(Mass, 0) + 1
                    self.MassAASiteCounts[Key] = self.MassAASiteCounts.get(Key, 0) + 1
                    OKSiteCount += 1
            PeptideCount += 1
            if OKFlag:
                OKPeptideCount += 1
                self.MassPeptideCounts[Mass] = self.MassPeptideCounts.get(Mass, 0) + 1
                self.MassAAPeptideCounts[Key] = self.MassAAPeptideCounts.get(Key, 0) + 1
            PrevSite = Site
            while len(Bits) < 54:
                Bits.append("")
            Bits[52] = str(PValue)
            Bits[53] = str(self.FDRByPValue[PValue])
            OutFile.write(string.join(Bits, "\t") + "\n")
        File.close()
        print "%s lines"%LineNumber
        print "%s peptides, %s sites"%(PeptideCount, SiteCount)
        print "%s valid peptides, %s valid sites"%(OKPeptideCount, OKSiteCount)
        print "Cutoff was %s"%self.SitePValueCutoff
    def ReportPTMDict(self, Title, DictA, DictB):
        print
        print Title
        print "SiteCount\tPeptideCount\tMass\tAA\t"
        SortedList = []
        for (Key, Count) in DictA.items():
            SortedList.append((Count, Key))
        SortedList.sort()
        SortedList.reverse()
        for (Count, Key) in SortedList:
            CountB = DictB[Key]
            if type(Key) in (types.ListType, types.TupleType):
                print "%s\t%s\t%s\t%s\t"%(Count, CountB, Key[0], Key[1])
            else:
                print "%s\t%s\t%s\t"%(Count, CountB, Key)
    def ReportAdductCount(self, Dict):
        KnownAdductCount = 0
        KnownGoodCount = 0
        UnknownCount = 0
        SortedUnknown = []
        for (Key, Count) in Dict.items():
            (Mass, AA) = Key
            KnownAdduct = 0
            KnownGood = 0
            if Mass == 57:
                KnownAdduct = 1
            elif Mass == -17 and AA in ("N", "Q"):
                KnownAdduct = 1
            elif Mass in (16, 32) and AA in ("M", "W"):
                KnownAdduct = 1
            elif Mass == -48 and AA == "M":
                KnownAdduct = 1
            elif Mass == 58 and AA == "M":
                KnownAdduct = 1 # it's both, so take credit for the +42
            elif Mass == 161 and AA == "Q":
                KnownGood = 1
            elif Mass == -43 and AA == "C":
                KnownGood = 1
            elif Mass == -57 and AA == "C":
                KnownAdduct = 1
            elif Mass == -18 and AA in ("D", "E", "S", "T"):
                KnownAdduct = 1
            elif Mass == 22 and AA in ("D", "E", "S", "T"):
                KnownAdduct = 1
            elif Mass == 38 and AA in ("D", "E", "S", "T"):
                KnownAdduct = 1
            elif Mass in (80, 160) and AA in ("S", "T", "Y"):
                KnownGood = 1
            elif Mass == 42:
                KnownGood = 1
            elif Mass == 72 and AA == "K":
                KnownGood = 1
            elif Mass == 55 and AA == "R":
                KnownGood = 1
            elif Mass == 58 and AA == "K":
                KnownGood = 1
            elif Mass == 43:
                KnownAdduct = 1
            elif Mass == 14 and AA in ("H", "K"):
                KnownGood = 1
            elif Mass == 12 and AA == "W":
                KnownAdduct = 1
            elif Mass == 28 and AA in ("H", "K"):
                KnownGood = 1
            elif Mass == 28 and AA in ("S", "T"):
                KnownAdduct = 1
            elif Mass == 40 and AA == "Q":
                KnownAdduct = 1
            elif Mass == 4 and AA == "W":
                KnownGood = 1
            if KnownAdduct == 1:
                KnownAdductCount += Count
            elif KnownGood == 1:
                KnownGoodCount += Count
            else:
                UnknownCount += Count
                SortedUnknown.append((Count, Mass, AA))
        print "%s unknown, %s known adduct, %s known"%(UnknownCount, KnownAdductCount, KnownGoodCount)
        SortedUnknown.sort()
        SortedUnknown.reverse()
        for (Count, Mass, AA) in SortedUnknown[:25]:
            print "Unknown:", Mass, AA, Count

    def Main(self):
        ParseFileName = sys.argv[1]
        OutputFileName = sys.argv[2]
        self.FindPValueCutoff(ParseFileName)
        self.ParseFile(ParseFileName, OutputFileName)
        self.ReportPTMDict("Sites/peptides by mass", self.MassSiteCounts, self.MassPeptideCounts)
        self.ReportPTMDict("Sites/peptides by mass and aa", self.MassAASiteCounts, self.MassAAPeptideCounts)
        print "Adduct summary for SITES:"
        self.ReportAdductCount(self.MassAASiteCounts)
        print "Adduct summary for PEPTIDES:"
        self.ReportAdductCount(self.MassAAPeptideCounts)
        ####self.ReportPTMDict("Peptides by mass", self.MassPeptideCounts)
        ####self.ReportPTMDict("Peptides by mass and aa", self.MassAAPeptideCounts)

if __name__ == "__main__":
    Counter = PTMCounter()
    Counter.Main()