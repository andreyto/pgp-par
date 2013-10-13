"""RelativeQuantitation.py
Given the NSAF of several datasets (see NormalizedSpectrumAbundance.py),
compute the mean and stdev of each dataset.  Then report genes with very
different abundance.  Following this formula

Define expression category
0. not expressed        =            0  <= Expression <= mean - 3stdev
1. low expression       = mean - 3stdev <= Expression <= mean - stdev
2. average expression   = mean - stdev  <= Expression <= mean + stdev
3. high expression      = mean + stdev  <= Expression <= mean + 3stdev
4. very high expression = mean + 3stdev <= Expression <= inf.

report all genes for which the cagegory(Sample0)- category(Sample1) >=2


"""

UsageInfo="""RelativeQuantitation.py
Given the NSAF of several datasets (see NormalizedSpectrumAbundance.py),
compute the mean and stdev of each dataset.  Then report genes with very
different abundance.

Required Options:
 -r [FileName] File of NSAF, col 1-n are NSAF values.  col 0 = gene
 -w [FileName] Output from program


"""

import os
import sys
import getopt
import SelectProteins
import ResultsParser
from Utils import *
Initialize()

class AbacusClass(ResultsParser.ResultsParser):
    def __init__ (self):
        self.DBPath = []
        self.OutputPath = None
        self.InspectInput = None
        self.NSAFLists = {} # key= column header. value = list of NSAF
        self.NSAFStats = {} #key = column header. value = mean, median, stdev
        self.Headers = None # list of column headers
        ResultsParser.ResultsParser.__init__(self)
    def Main(self):
        self.ProcessResultsFiles(self.InspectInput, self.ParseResults)
        self.GetNSAFStats()
        self.ProcessResultsFiles(self.InspectInput, self.CategorizeResults)

    def CategorizeResults(self, FileName):
        """Get each gene line, and the category of the NSAF values.  Print to file
        those that are 2 categories away
        """
        Handle = open(Filename, "rb")
        for Line in Handle.xreadlines():
            if Line[0] == "#":
                continue
            Bits = Line.strip().split("\t")
            for Index in range(1, len(Bits)):
                Value = float(Bits[Index])
                NSAFKey = self.Headers[Index]
                Category = self.GetCategory(Value, NSAFKey)

    def GetCategory(self, Value, NSAFKey):
        """see how many stdev out it is
        """
        
            
        

    def ParseResults(self, FileName):
        Handle = open(FileName, "rb")
        for Line in Handle.xreadlines():
            if Line[0] == "#":
                self.Headers = Line.strip().split("\t")
                for Header in self.Headers:
                    self.NSAFLists[Header] = []
                continue
            Bits = Line.strip().split("\t")
            #cludge for gene name, col 0
            for Index in range(1, len(Bits)):
                Value = float(Bits[Index])
                NSAFKey = self.Headers[Index]
                self.NSAFLists[NSAFKey].append(Value)
        Handle.close()

    def GetNSAFStats(self):
        """Go through the list and get mean, median, stdev
        """
        NumHeaders = len(self.Headers)
        for Index in range(1, NumHeaders):
            HeaderKey = self.Headers[Index]
            (Average, Median) = self.AverageList(self.NSAFLists[HeaderKey])
            Stdev = self.StdevList(self.NSAFLists[HeaderKey], Average)
            self.NSAFStats[HeaderKey] = (Average, Median, Stdev)

    def AverageList(self, List1, List2 = []):
        """get the average, median, min, max of a list or two"""
        List1.extend(List2)
        Min = 1000000
        Max = -1
        Sum = 0.0
        Average = 0.0
        for Measurement in List1:
            Sum += Measurement
            if Measurement > Max:
                Max = Measurement
            if Measurement < Min:
                Min = Measurement
        Average = Sum / len(List1)
        MedianIndex = len(List1) / 2
        Median = List1[MedianIndex]
        return (Average, Median)

    def StdevList(self, List, average):
        """Given the average (mean, median, random number), get the stdev"""
        SumSquaredDev = 0.0
        for Item in List:
            Deviation = Item - Average
            SumSquaredDev += (Deviation * Deviation)
        Variance = SumSquaredDev / Len
        Stdev = math.sqrt(Variance)
        return Stdev


                
    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "d:w:r:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-d":
                # -r results file(s)
                if not os.path.exists(Value):
                    print "** Error: couldn't find database file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.DBPath.append(Value)
            if Option == "-w":
                # -r results file(s)
                self.OutputPath = Value
            if Option == "-r":
                self.InspectInput = Value
        if not OptionsSeen.has_key("-d") or not OptionsSeen.has_key("-w") or not OptionsSeen.has_key("-r"):
            print UsageInfo
            sys.exit(1)


if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(psyco not found - running in non-optimized mode)"
    DoStuff = AbacusClass()
    DoStuff.ParseCommandLine(sys.argv[1:])
    DoStuff.Main()        
