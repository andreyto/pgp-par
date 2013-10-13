"""
When we run a search on the grid, we might search several .mzxml files
in a single run.  It's convenient to have the results split into
one file per .mzxml file.
- Iterate over results-files from NEWEST to OLDEST
- If you've already seen results for a file, IGNORE any further results for it!
"""
import os
import sys

def Main():
    ScansSeen = {}
    #BlockResultsDir = "E:\ms\PeptideAtlas\BlockResultsX"
    #ResultsDir = "E:\ms\PeptideAtlas\ResultsX"
    BlockResultsDir = "E:\ms\PeptideAtlas\BlockResultsIPI"
    ResultsDir = "E:\ms\PeptideAtlas\ResultsIPI"
    SortedFiles = []
    for ResultsFileName in os.listdir(BlockResultsDir):
        BlockResultsPath = os.path.join(BlockResultsDir, ResultsFileName)
        if os.path.isdir(BlockResultsPath):
            continue
        ModTime = os.stat(BlockResultsPath).st_mtime
        SortedFiles.append((ModTime, ResultsFileName))
    SortedFiles.sort()
    SortedFiles.reverse()
    CurrentScanKey = None
    SkippingFlag = 0
    for (ModTime, ResultsFileName) in SortedFiles:
        print ModTime, ResultsFileName
        BlockResultsPath = os.path.join(BlockResultsDir, ResultsFileName)
        if os.path.isdir(BlockResultsPath):
            continue
        File = open(BlockResultsPath, "rb")
        CurrentXMLFileName = None
        CurrentOutputFile = None
        for FileLine in File.xreadlines():
            if not FileLine or FileLine[0] == "#":
                continue
            Bits = FileLine.split("\t")
            try:
                ScanNumber = int(Bits[1])
            except:
                continue
            XMLFileName = Bits[0].replace("/","\\").split("\\")[-1]
            # Check for duplicate results:
            ScanKey = (XMLFileName, ScanNumber)
            if ScanKey != CurrentScanKey:
                # We've encountered a new spectrum!  We'll REJECT it, if we already saw
                # search hits for this spectrum:
                if ScansSeen.has_key(ScanKey):
                    SkippingFlag = 1
                else:
                    SkippingFlag = 0
                ScansSeen[ScanKey] = 1
                CurrentScanKey = ScanKey
            if SkippingFlag:
                continue
            if XMLFileName != CurrentXMLFileName:
                if CurrentOutputFile:
                    CurrentOutputFile.close()
                (Stub, Extension) = os.path.splitext(XMLFileName)
                OutputPath = os.path.join(ResultsDir, "%s.txt"%Stub)
                print "%s\t%s"%(BlockResultsPath, OutputPath)
                CurrentOutputFile = open(OutputPath, "a+")
                CurrentXMLFileName = XMLFileName
            CurrentOutputFile.write(FileLine)
        File.close()

if __name__ == "__main__":
    import psyco
    psyco.full()
    Main()