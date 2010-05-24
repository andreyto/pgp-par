#!/usr/bin/env python
"""
Count the number of level-2 scans with 10+ peaks from an mgf or mzxml file.
Iterate over a directory tree and produce sums for each file, and totals by directory.
"""
import os
import sys
import getopt

MEG = 1024 * 1024

# Keep track of the files for which we've written the scan count out.  When we come to
# the end, if there are any known scan counts that we haven't yet written out,
# then we write them out.
FilesWrittenOut = {}

class CountScanBits:
    FileName = 0
    Stub = 1
    Count = 2
    MaxScan = 3
    FileSize = 4

class ScanCounter:
    def __init__(self):
        self.KnownCounts = {} # filename -> number
        self.KnownMaxScans = {} # filename -> number
        self.FileSizes = {}
        self.CountFileName = "ScanCount.txt"  # default
        self.CountDirectory = None
    def SaveKnownCounts(self, CountFileName):
        Keys = self.KnownCounts.keys()
        Keys.sort()
        File = open(CountFileName, "wb")
        for FileName in Keys:
            Stub = os.path.splitext(FileName)[0]
            Str = "%s\t%s\t"%(FileName, Stub)
            Str += "%s\t%s\t"%(self.KnownCounts[FileName], self.KnownMaxScans[FileName])
            Str += "%s\t"%self.FileSizes[FileName]
            File.write(Str + "\n")
        File.close()
    def LoadKnownCounts(self, CountFileName):
        if not os.path.exists(CountFileName):
            return
        File = open(CountFileName, "rb")
        for FileLine in File.xreadlines():
            Bits = FileLine.strip().split("\t")
            if FileLine[0] == "#" or len(Bits) < 4:
                continue
            Count = int(Bits[CountScanBits.Count])
            MaxScan = int(Bits[CountScanBits.MaxScan])
            FileSize = int(Bits[CountScanBits.FileSize])
            FileName = Bits[CountScanBits.FileName]
            self.KnownCounts[FileName] = Count
            self.KnownMaxScans[FileName] = MaxScan
            self.FileSizes[FileName] = FileSize
        File.close()
    def ParseCommandLine(self, Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:w:R")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-r":
                self.CountDirectory = Value
            elif Option == "-w":
                self.CountFileName = Value
    def Main(self):
        # Load counts:
        self.LoadKnownCounts(self.CountFileName)
        # Count our directory:
        self.CountScansInDirectory(self.CountDirectory)
        # Save counts:
        self.SaveKnownCounts(self.CountFileName)
    def CountScansInDirectory(self, InputFilePath):
        pathLen = len(InputFilePath) + 1
        for root,dirs,files in os.walk( InputFilePath ):
            if root == InputFilePath:
                subDir = ''
            else:
                subDir = root[pathLen:]
            for fileStr in files:
                # Handle a single file:
                FileName = os.path.join(subDir,fileStr)
                FilePath = os.path.join(root,fileStr)
                if self.KnownCounts.has_key(FileName):
#                    print "(Skip %s - already counted it)"%FileName
                    continue
                print "(Count %s...)"%FileName
                (Stub, Extension) = os.path.splitext(fileStr)
                Extension = Extension.lower()
                if Extension == ".mgf":
                    Result = self.CountScansMGF(FilePath)
                elif Extension == ".mzxml":
                    Result = self.CountScansMZXML(FilePath)
                elif Extension == ".ms2":
                    Result = self.CountScansMS2(FilePath)
                else:
                    Result = None
                if not Result:
                    continue
                self.KnownCounts[FileName] = Result[0]
                self.KnownMaxScans[FileName] = Result[1]
                self.FileSizes[FileName] = os.stat(FilePath).st_size
    def CountScansMZXML(self, FilePath):
        FileScanCount = 0
        File = open(FilePath, "r")
        Text = ""
        MZXMLScanNumbers = {}
        while 1:
            Block = File.read(MEG)
            if not Block:
                break
            Text += Block
            Pos = -1
            while 1:
                NextPos = Text.find("<scan num=\"", Pos + 1)
                if NextPos == -1 or NextPos > len(Text) - 100:
                    Text = Text[Pos + 1:]
                    break
                MSLevelPos = Text.find("msLevel=\"", Pos + 1)
                PeakCountPos = Text.find("peaksCount=\"", Pos + 1)
                QuotePos = Text.find('"', NextPos + 11)
                ScanNumber = int(Text[NextPos + 11:QuotePos])
                QuotePos = Text.find('"', PeakCountPos + 12)
                PeakCount = int(Text[PeakCountPos + 12:QuotePos])
                # Check the number of peaks AND the MSLevel:
                if PeakCount > 10 and Text[MSLevelPos + 9] in ("2", "3"):
                    #ScanCount += 1
                    FileScanCount += 1
                    MZXMLScanNumbers[ScanNumber] = 1
                Pos = NextPos
        File.close()
        ScanNumbers = MZXMLScanNumbers.keys()
        ScanNumbers.sort()
        FileScanCount = len(ScanNumbers)
#        MinScanNumber = ScanNumbers[0]
        MaxScanNumber = ScanNumbers[-1]
        return (FileScanCount, MaxScanNumber)
        #return (MinScanNumber, MaxScanNumber, FileScanCount)
    def CountScansMGF(self, FilePath):
        "Count the number of SCAN= lines in an .mgf file."
        FileScanCount = 0
        File = open(FilePath, "rb")
        Text = ""
        ScanNumbers = {}
        ScanNumber = 0
        while 1:
            Block = File.read(MEG)
            if not Block:
                break
            Text += Block
            Pos = -1
            while 1:
                NextPos = Text.find("BEGIN IONS", Pos + 1)
                if NextPos == -1 or NextPos > len(Text) - 100:
                    Text = Text[Pos + 1:]
                    break
                #ScanNumber = int(Text[NextPos + 5:NextPos + 100].split()[0])
                ScanNumbers[ScanNumber] = 1
                ScanNumber += 1
                Pos = NextPos
            Text = Text[-100:]
        File.close()
        ScanNumbers = ScanNumbers.keys()
        #ScanNumbers.sort()
        FileScanCount = len(ScanNumbers)
        if ScanNumbers:
            MaxScanNumber = max(ScanNumbers)
        else:
            MaxScanNumber = 0
        return (FileScanCount, MaxScanNumber)

    def CountScansMS2(self, FilePath):
        """New version to count scans in an MS2 file.
        """
        file = open(FilePath)
        scanCount = 0
        maxScan = 0
        for line in file:
            if line[0] == 'S':
                scanCount += 1
                (s,scanNum,num2,other) = line.split("\t")
                if scanNum > max:
                    maxScan = scanNum
        return (scanCount, maxScan)
            
def CountScansMS2(FilePath):
    """Non working version of MS2 count scans.
    """
    ScanNumbers = {}
    File = open(FilePath, "rb")
    for FileLine in File.xreadlines():
        if FileLine[0] == ":":
            Bits = FileLine[1:].split(".")
            ScanNumber = int(Bits[0])
            ScanNumbers[ScanNumber] = 1
    File.close()
    ScanNumbers = ScanNumbers.keys()
    ScanNumbers.sort()
    FileScanCount = len(ScanNumbers)
    if len(ScanNumbers):
        MinScanNumber = ScanNumbers[0]
        MaxScanNumber = ScanNumbers[-1]
    else:
        MinScanNumber = 0
        MaxScanNumber = 0
    return (MinScanNumber, MaxScanNumber, FileScanCount)

UsageInfo = """
CountScans.py - report how many spectra are contained in an .mzXML or .mgf file
Arguments:
 -r [Directory]: Count scans from here
 -w [FileName]: Save scan counts here.
"""

def main(args):
    Counter = ScanCounter()
    Counter.ParseCommandLine(args)
    if not Counter.CountDirectory:
        raise Exception( UsageInfo )
    Counter.Main()

if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(no psyco)"

    main(sys.argv[1:])
