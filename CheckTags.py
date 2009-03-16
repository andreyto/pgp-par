"""
Read tagger output, and determine whether one or more tags are correct.
"""
import os
import sys
from Utils import *
Initialize()

FLANKING_MASS_EPSILON = 4.0

class TagChecker:
    def __init__(self):
        pass
    def Main(self, TagFileName, TruePeptide, ScanNumber):
        ValidTagCount = 0
        FirstValidTag = None
        FilteredAminos = TruePeptide.Aminos.replace("Q", "K").replace("I", "L")
        File = open(TagFileName, "rb")
        for FileLine in File.xreadlines():
            Bits = FileLine.split("\t")
            try:
                TagScan = int(Bits[1])
                TagIndex = int(Bits[3])
                PrefixMass = float(Bits[4])
                SuffixMass = float(Bits[6])
            except:
                continue
            if ScanNumber!=None and TagScan != ScanNumber:
                continue
            Tag = GetPeptideFromModdedName(Bits[5])
            FilteredTagAminos = Tag.Aminos.replace("Q", "K").replace("I", "L")
            for AminoIndex in range(len(TruePeptide.Aminos)):
                if FilteredAminos[AminoIndex:AminoIndex + 3] != FilteredTagAminos:
                    continue
                PrefixDiff = PrefixMass - TruePeptide.Masses[AminoIndex]
                if abs(PrefixDiff) > FLANKING_MASS_EPSILON:
                    continue
                ValidTagCount += 1
                if FirstValidTag==None:
                    FirstValidTag = TagIndex
                ThisScanFirstCorrectTag = TagIndex
                print "Tag #%s '%s' correct (prefix %s delta %s)"%(TagIndex, Bits[5], PrefixMass, PrefixDiff)
            
if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(Warning: psyco not found; running without optimization)"
    TagFileName = sys.argv[1]
    TruePeptide = GetPeptideFromModdedName(sys.argv[2])
    if len(sys.argv)>3:
        ScanNumber = int(sys.argv[3])
    else:
        ScanNumber = None
    Bob = TagChecker()
    Bob.Main(TagFileName, TruePeptide, ScanNumber)