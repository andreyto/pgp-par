"""
Count the number of spectra searched.
"""
import os
import sys

class Counter:
    def __init__(self):
        # XMLFileName -> dictionary of scans
        self.AllSpectraDict = {}
        self.MaxHitScanNumber = {}
        self.LineCount = 0
    def CountSpectraFromFile(self, Path):
        try:
            File = open(Path, "rb")
        except:
            print "* Can't open '%s', not counting hits."%Path
            return
        PrevSpectrum = None
        for FileLine in File.xreadlines():
            self.LineCount += 1
            Bits = FileLine.split("\t")
            try:
                ScanNumber = int(Bits[1])
            except:
                continue
            MZXMLFileName = Bits[0].replace("/","\\").split("\\")[-1]
            if not self.AllSpectraDict.has_key(MZXMLFileName):
                self.AllSpectraDict[MZXMLFileName] = {}
            Spectrum = (Bits[0], Bits[1])
            # If we see a spectrum that we've seen before, and it's not from the immediately
            # preceding rows, that's a BAD thing and indicates overlap in the results.
            if Spectrum != PrevSpectrum and self.AllSpectraDict[MZXMLFileName].has_key(ScanNumber):
                print "** Warning: Apparent duplication!", Path, MZXMLFileName, ScanNumber, self.AllSpectraDict[MZXMLFileName][ScanNumber]
            self.AllSpectraDict[MZXMLFileName][ScanNumber] = Path
            self.MaxHitScanNumber[MZXMLFileName] = max(ScanNumber, self.MaxHitScanNumber.get(MZXMLFileName, 0))
            PrevSpectrum = Spectrum
        File.close()
    def Main(self):
        OutputFile = open("IPISearchScanCounts.txt", "wb")
        #Dirs = [r"E:\ms\PeptideAtlas\ResultsX", "e:\\ms\\Briggs\\ResultsX"]
        Dirs = [r"E:\ms\PeptideAtlas\ResultsIPI", "e:\\ms\\Briggs\\ResultsIPI"]
        for Dir in Dirs:
            for FileName in os.listdir(Dir):
                Path = os.path.join(Dir, FileName)
                if not os.path.isdir(Path): 
                    print "Count %s..."%Path
                    self.CountSpectraFromFile(Path)
        # Summarize:
        print "LineCount:", self.LineCount
        print "XML files:", len(self.AllSpectraDict.keys())
        AllScanCount = 0
        Keys = self.AllSpectraDict.keys()
        Keys.sort()
        for XMLFileName in Keys:
            ScanCount = len(self.AllSpectraDict[XMLFileName].keys())
            AllScanCount += ScanCount
            Str = "File\t%s\t%s\t%s\t"%(XMLFileName, ScanCount, self.MaxHitScanNumber[XMLFileName])
            OutputFile.write(Str+"\n")
            print Str
Bob = Counter()
Bob.Main()