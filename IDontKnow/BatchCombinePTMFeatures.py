"""
Call CombinePTMFeatures several times, once for each database-block.
"""
import sys
import os
import math
import getopt
import string
import CombinePTMFeatures
from TrainPTMFeatures import FormatBits

class CombineBatcher:
    def __init__(self):
        self.DatabaseBlockSize = 100000
    def ParseCommandLine(self, Arguments):
        self.MergeFlag = 0
        (Options, Args) = getopt.getopt(Arguments, "d:r:w:mb:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-r":
                self.PTMFeatureDirectory = Value
            elif Option == "-d":
                self.DBPath = Value
            elif Option == "-w":
                self.OutputDir = Value
            elif Option == "-m":
                self.MergeFlag = 1
            elif Option == "-b":
                self.DatabaseBlockSize = int(Value)
            else:
                print "* Error: Unrecognized option %s"%Option
    def MergeBlocks(self):
        """
        Parse all the block files.  Then merge them into one large file, updating the
        sites-per-species and spectra-for-species columns.
        MergeBlocks is called AFTER all the blocks, in a separate script run.
        """
        print "MERGING all blocks..."
        self.ModTypeSpectrumCount = {}
        self.ModTypeSiteCount = {}
        # Iterate over files once to get site counts:
        for FileName in os.listdir(self.OutputDir):
            if FileName.find("Block") != 0:
                continue
            Path = os.path.join(self.OutputDir, FileName)
            File = open(Path, "rb")
            for FileLine in File.xreadlines():
                if FileLine[0] == "#":
                    continue
                Bits = FileLine.strip().split("\t")
                try:
                    SpectrumCount = int(Bits[FormatBits.SpectrumCount])
                except:
                    print Bits
                    traceback.print_exc()
                    continue
                ModAA = Bits[FormatBits.ModifiedAA]
                ModMass = int(Bits[FormatBits.ModificationMass])
                ModTypeKey = (ModAA, ModMass)
                self.ModTypeSiteCount[ModTypeKey] = self.ModTypeSiteCount.get(ModTypeKey, 0) + 1
                self.ModTypeSpectrumCount[ModTypeKey] = self.ModTypeSpectrumCount.get(ModTypeKey, 0) + SpectrumCount
            File.close()
        # Iterate over files again to write them out:
        OutputFile = open(os.path.join(self.OutputDir, "PTMFeatures.txt"), "wb")
        FirstFileFlag = 1
        for FileName in os.listdir(self.OutputDir):
            if FileName.find("Block") != 0:
                continue
            Path = os.path.join(self.OutputDir, FileName)
            File = open(Path, "rb")
            for FileLine in File.xreadlines():
                if FileLine[0] == "#":
                    if FirstFileFlag:
                        OutputFile.write(FileLine)
                    continue
                Bits = list(FileLine.strip().split("\t"))
                ModAA = Bits[FormatBits.ModifiedAA]
                ModMass = int(Bits[FormatBits.ModificationMass])
                ModTypeKey = (ModAA, ModMass)
                Bits[FormatBits.SpectraWithThisModType] = str(self.ModTypeSpectrumCount[ModTypeKey])
                Bits[FormatBits.SitesWithThisModType] = str(self.ModTypeSiteCount[ModTypeKey])
                Bits[FormatBits.LogSpectraThisModType] = str(math.log(self.ModTypeSpectrumCount[ModTypeKey]))
                Bits[FormatBits.LogSitesThisModType] = str(math.log(self.ModTypeSiteCount[ModTypeKey]))
                String = string.join(Bits, "\t")
                OutputFile.write(String + "\n")
            File.close()
            FirstFileFlag = 0
    def Run(self):
        """
        Main entry point.  For each database block, call CombinePTMFeatures.
        Note: ComputePTMFeatures handles one block of *spectra* at once, and CombinePTMFeatures
        handles one block of the *databas* at once; these steps allow us to keep
        memory usage reasonably low!
        """
        if self.MergeFlag:
            self.MergeBlocks()
            return
        DBSize = os.stat(self.DBPath).st_size
        StartPos = 0
        BlockIndex = 0
        while StartPos < DBSize:
            Merger = CombinePTMFeatures.PTMFeatureMerger()
            BlockPath = os.path.join(self.OutputDir, "Block%s.txt"%BlockIndex)
            ArgumentList = ["-d", self.DBPath, "-w", BlockPath, "-r", self.PTMFeatureDirectory]
            ArgumentList.append("-s")
            ArgumentList.append(str(StartPos))
            ArgumentList.append("-e")
            ArgumentList.append(str(StartPos + self.DatabaseBlockSize))
            if BlockIndex == 0:
                ArgumentList.append("-x")
            print ArgumentList
            Merger.ParseCommandLine(ArgumentList)
            Merger.MergeResults()
            # Iterate to next:
            BlockIndex += 1
            StartPos += self.DatabaseBlockSize

if __name__ == "__main__":
    Batcher = CombineBatcher()
    Batcher.ParseCommandLine(sys.argv[1:])
    Batcher.Run()