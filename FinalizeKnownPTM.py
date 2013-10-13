"""
Compare p-values of raw and known-ptm annotations, produce a "final" annotation
"""
import os
import sys
import math
from Utils import *
from TrainPTMFeatures import FormatBits
Initialize()

class ReconcilerClass:
    def Main(self):
        File = open(sys.argv[1], "rb")
        OutputFile = open(sys.argv[2], "wb")
        for FileLine in File.xreadlines():
            if FileLine[0] == "#":
                OutputFile.write(FileLine)
                continue
            Bits = FileLine.strip().split("\t")
            PValueOriginal = float(Bits[FormatBits.SitePValue])
            try:
                PValueKnown = float(Bits[FormatBits.KnownPTMSitePValue])
            except:
                PValueKnown = 9999
            KnownFlag = 0
            ScoreDiff = math.log(PValueKnown) - math.log(PValueOriginal)
            if ScoreDiff < 0.5:
                # Use the KNOWN modification type:
                Peptide = GetPeptideFromModdedName(Bits[FormatBits.KnownPTMAnnotation])
                Keys = Peptide.Modifications.keys()
                KnownFlag = 1
                if Keys:
                    ModIndex = Peptide.Modifications.keys()[0]
                    ModAA = Peptide.Aminos[ModIndex]
                    ModMass = Peptide.Modifications[ModIndex][0].Mass
                else:
                    ModAA = ""
                    ModMass = ""
            else:
                ModAA = Bits[FormatBits.ModifiedAA]
                ModMass = Bits[FormatBits.ModificationMass]
            Bits = list(Bits)
            while len(Bits) < 52:
                Bits.append("")
            Bits[49] = str(KnownFlag)
            Bits[50] = ModAA
            Bits[51] = str(ModMass)
            #OutputFile.write("%s\t%s\t"%(ModAA, ModMass))
            Str = string.join(Bits, "\t")
            OutputFile.write(Str + "\n")

if __name__ == "__main__":
    Bob = ReconcilerClass()
    Bob.Main()