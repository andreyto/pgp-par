"""
Simple filter: Drop second-and-later hits with poor p-value.
"""
import os
import sys
import traceback
import getopt
import ResultsParser

class ResultsFilter(ResultsParser.ResultsParser):
    def __init__(self):
        self.OutputPath = "FilteredResults.txt"
        self.PValueCutoff = 0.1
        ResultsParser.ResultsParser.__init__(self)
    def ParseCommandLine(self, Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:w:p:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-r":
                self.ParsePath = Value
            elif Option == "-w":
                self.OutputPath = Value
            elif Option == "-r":
                self.PValueCutoff = float(Value)
    def Main(self):
        self.OutputFile = open(self.OutputPath, "wb")
        self.FirstFileFlag = 1
        self.ProcessResultsFiles(self.ParsePath, self.FilterOneFile)
    def FilterOneFile(self, FilePath):
        File = open(FilePath, "rb")
        OldSpectrum = None
        for FileLine in File.xreadlines():
            if FileLine[0] == "#":
                if self.FirstFileFlag:
                    self.OutputFile.write(FileLine)
                continue
            Bits = FileLine.split("\t")
            Spectrum = (Bits[0], Bits[1])
            PValue = float(self.Columns.PValue)
            if (Spectrum != OldSpectrum) or (PValue <= self.PValueCutoff):
                self.OutputFile.write(FileLine)
            OldSpectrum = Spectrum
        self.FirstFileFlag = 0
if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "psyco not found - running without optimization"
    Filter = ResultsFilter()
    Filter.ParseCommandLine(sys.argv[1:])
    Filter.Main()
    