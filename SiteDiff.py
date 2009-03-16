"""
Compare sites between two sets of annotations.
"""
import os
import sys
import psyco
import string
import traceback
from Utils import *
Global.FixedMods = {"C":57.0518} 
Initialize()

MinSpectraForSite = 5

class SiteDiffClass:
    def __init__(self):
        self.AnnotationsA = []
        self.AnnotationsB = []
    def ReadDB(self, FileName):
        self.ProteinNames = []
        self.Sequences = []
        File = open(FileName, "r")
        for FileLine in File.xreadlines():
            FileLine = FileLine.strip()
            if FileLine[0] == ">":
                self.ProteinNames.append(FileLine[1:])
                self.Sequences.append("")
            else:
                self.Sequences[-1] += FileLine
        File.close()
    def ReadAnnotations(self, FileName):
        File = open(FileName, "r")
        Annotations = []
        # Annotations[Protein#][Residue][Mass] = count
        for X in range(len(self.ProteinNames)):
            List = []
            for Y in range(len(self.Sequences[X])):
                List.append({})
            Annotations.append(List)
        LineNumber = 0
        for FileLine in File.xreadlines():
            LineNumber += 1
            if LineNumber % 1000 == 0:
                print LineNumber
            Bits = FileLine.split("\t")
            try:
                Annot = Bits[2]
            except:
                continue # dummy line
            if Annot[2] == ".":
                Annot = Annot[2:-2]
            Peptide = GetPeptideFromModdedName(Annot)
            FoundFlag = 0
            for SequenceIndex in range(len(self.Sequences)):
                Pos = self.Sequences[SequenceIndex].find(Peptide.Aminos)
                if Pos==-1:
                    continue
                FoundFlag = 1
                for X in range(len(Peptide.Aminos)):
                    Annotations[SequenceIndex][Pos + X][None] = Annotations[SequenceIndex][Pos + X].get(None, 0) + 1
                    if Peptide.Modifications.has_key(X):
                        Mass = int(round(Peptide.Modifications[X][0].Mass))
                        Annotations[SequenceIndex][Pos + X][Mass] = Annotations[SequenceIndex][Pos + X].get(Mass, 0) + 1
                        #print "%s %s %s"%(SequenceIndex, Pos+X, Mass)
                    else:
                        Annotations[SequenceIndex][Pos + X][0] = Annotations[SequenceIndex][Pos + X].get(0, 0) + 1
            if not FoundFlag:
                #print "wtf?", Annot
                pass
        File.close()
        return Annotations
    def SummarizeAnnotations(self, Annotations):
        print "Protein\tResidue\tAA\tMass\tCountA\tModA\tRateA\tCountB\tModB\tRateB\tMinC\tMinRate\tCall"
        for ProteinIndex in range(len(self.ProteinNames)):
            ProteinName = self.ProteinNames[ProteinIndex]
            Sequence = self.Sequences[ProteinIndex]
            #print "Find sites on %s (seqlen %s)"%(ProteinName, len(Sequence))
            for Pos in range(len(Sequence)):
                CoverageA = Annotations[ProteinIndex][Pos]
                CountA = CoverageA.get(None, 0)
                #print CountA, CountB, CoverageA, CoverageB
                for (Key, ModCountA) in CoverageA.items():
                    if Key and ModCountA >= MinSpectraForSite:
                        Str = "%s\t%s\t%s\t%s\t"%(ProteinName, Pos+1, Sequence[Pos], Key)
                        RateA = ModCountA / float(max(CountA, 1))
                        Str += "%s\t%s\t%.1f\t"%(CountA, ModCountA, 100*RateA)
                        print Str
    def CompareAnnotations(self, AnnotFileA, AnnotFileB):
        self.AnnotationsA = self.ReadAnnotations(AnnotFileA)
        self.AnnotationsB = self.ReadAnnotations(AnnotFileB)
        # Output bits:
        # Protein, Residue, ModMass, CountA, ModCount, Rate, CountB, ModCount, Rate, Call
        print "Protein\tResidue\tAA\tMass\tCountA\tModA\tRateA\tCountB\tModB\tRateB\tMinC\tMinRate\tCall"
        for ProteinIndex in range(len(self.ProteinNames)):
            ProteinName = self.ProteinNames[ProteinIndex]
            Sequence = self.Sequences[ProteinIndex]
            #print "Find sites on %s (seqlen %s)"%(ProteinName, len(Sequence))
            for Pos in range(len(Sequence)):
                CoverageA = self.AnnotationsA[ProteinIndex][Pos]
                CountA = CoverageA.get(None, 0)
                CoverageB = self.AnnotationsB[ProteinIndex][Pos]
                CountB = CoverageB.get(None, 0)
                #print CountA, CountB, CoverageA, CoverageB
                for (Key, ModCountA) in CoverageA.items():
                    if Key and ModCountA >= MinSpectraForSite:
                        Str = "%s\t%s\t%s\t%s\t"%(ProteinName, Pos+1, Sequence[Pos], Key)
                        RateA = ModCountA / float(max(CountA, 1))
                        ModCountB = CoverageB.get(Key, 0)
                        RateB = ModCountB / float(max(CountB, 1))
                        if RateA > RateB*2:
                            Call = "A"
                        elif RateB > RateA*2:
                            Call = "B"
                        else:
                            Call = "Both"
                        Str += "%s\t%s\t%.1f\t%s\t%s\t%.1f\t%s\t%s\t%.1f\t"%(CountA, ModCountA, 100*RateA, CountB, ModCountB, 100*RateB, Call, min(ModCountA, ModCountB), min(100*RateA, 100*RateB))
                        print Str
                for (Key, ModCountB) in CoverageB.items():
                    if Key and ModCountB >= MinSpectraForSite:
                        Str = "%s\t%s\t%s\t%s\t"%(ProteinName, Pos+1, Sequence[Pos], Key)
                        RateB = ModCountB / float(max(CountB, 1))
                        ModCountA = CoverageA.get(Key, 0)
                        if ModCountA >= MinSpectraForSite:
                            continue # we already printed a line for this site :)
                        RateA = ModCountA / float(max(CountA, 1))
                        if RateA > RateB*2:
                            Call = "A"
                        elif RateB > RateA*2:
                            Call = "B"
                        else:
                            Call = "Both"
                        Str += "%s\t%s\t%.1f\t%s\t%s\t%.1f\t%s\t%s\t%.1f\t"%(CountA, ModCountA, 100*RateA, CountB, ModCountB, 100*RateB, Call, min(ModCountA, ModCountB), min(100*RateA, 100*RateB))
                        print Str

psyco.full()                        
Differ = SiteDiffClass()
Differ.ReadDB("Database\\DrosProtAbundant.fasta")
#Differ.ReadDB("Database\\Lens.fasta")
Annotations = Differ.ReadAnnotations("HXAnnotations.txt")
Differ.SummarizeAnnotations(Annotations)
#Differ.CompareAnnotations("JYAnnotations.txt", "LDAnnotations.txt")