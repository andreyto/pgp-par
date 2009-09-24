UsageInfo = """ MakeHistogram.py
Takes an input file(s) and makes a histogram on the
input attribute, between true and false annotations


Required Options
 -r [FileName] Filename or directory of annotations
 -w [FileName] single output filename
 -a [AttributeName] An attribute from Inspect's output format
 -f [Bin Size] The size of a bin.  enter a float
 -b [Bin Size] Integer
 -c [Flag] Cumulative Histogram
"""


import os
import getopt
import sys
import ResultsParser

class AssassinClass(ResultsParser.ResultsParser):
    def __init__(self):
        self.InputFile = None
        self.OutputFile = None
        self.Attribute = None
        self.AttributeColumn = None
        self.BinSize = None
        self.FloatBin = 0
        self.TrueDistribution = [] # attribute values for peps hitting the true db
        self.FalseDistribution = []
        self.CumulativeFlag =0
        ResultsParser.ResultsParser.__init__(self)
        
    def Main(self):
        self.GetAttributeColumn()
        self.ProcessResultsFiles(self.InputFile, self.ParseAnnotations)
        if self.CumulativeFlag:
            self.SetUpCumulativeHistogram()
        else:
            self.RegularHistogram()

    def RegularHistogram(self):
        self.TrueDistribution.sort()
        self.FalseDistribution.sort()
        Max = int(self.TrueDistribution[-1])
        Handle = open(self.OutputFile, "wb")
        Handle.write("%s\tTrueHits\tFalseHits\n"%self.Attribute)
        print "trying to use range with %s as max and %s as bin"%(Max, self.BinSize)
        for Bin in range(0,Max, self.BinSize):
            # find all the answers with this score 
            TrueCount =0
            FalseCount =0
            for Value in self.TrueDistribution:
                if Value == Bin:
                    TrueCount +=1
            for Value in self.FalseDistribution:
                if Value == Bin:
                    FalseCount +=1
            #now print out this bin
            # I'm printing out the true-false, because that's the number of truepositives within the true db hits
            Handle.write("%s\t%s\t%s\n"%(Bin, TrueCount-FalseCount, FalseCount))
        Handle.close()
        

    def SetUpCumulativeHistogram(self):
        """now sort on score, and then get the values for each bin"""
        ## get them in a top down
        self.FalseDistribution.sort()
        self.TrueDistribution.sort()
        ## now tally
        TrueCount = 0
        FalseCount = 0
        Max = self.TrueDistribution[-1]
        Handle = open(self.OutputFile, "wb")
        Handle.write("%s\tTrueHits\tFalseHits\n"%self.Attribute)
        print "trying to use range with %s as max and %s as bin"%(Max, self.BinSize)
        ##Range does not do well with floats, so let's get our own list, and we'll pop it out
        if self.FloatBin:
            BinSet = self.GetRange(0,Max, self.BinSize)
        else:
            BinSet = range(0,Max, self.BinSize)
        for Bin in BinSet:
            # find all the answers with this score or greater.
            if len(self.TrueDistribution) > 0:
                PeekTrue = self.TrueDistribution[-1]
                while PeekTrue > Bin:
                    TrueCount += 1
                    self.TrueDistribution.pop()
                    # now peek at the next one for the while loop
                    if len(self.TrueDistribution) > 0:
                        PeekTrue = self.TrueDistribution[-1]
                    else:
                        break

            if len(self.FalseDistribution) > 0:
                PeekFalse = self.FalseDistribution[-1]
                while PeekFalse > Bin:
                    FalseCount += 1
                    self.FalseDistribution.pop()
                    if len(self.FalseDistribution) > 0:
                        PeekFalse = self.FalseDistribution[-1]
                    else:
                        break
            #now print out this bin
            # I'm printing out the true-false, because that's the number of truepositives within the true db hits
            Handle.write("%s\t%s\t%s\n"%(Bin, TrueCount-FalseCount, FalseCount))
        Handle.close()


    def GetRange(self, Min, Max, Step):
        """This is a substitute for the python range, which can't handle float steps
            You hand me the min and max.  MIN < MAX !!!!!!!!!!!!!!! if you want a reverse list
            then supply a negative step and I'll reverse it for you
        """
        List = []
        Iter = Min
        ReverseIt =0
        if Step < 0:
            Step *= -1
            ReverseIt = 1
        while(Iter < Max):
            List.append(Iter)
            Iter += Step
        if ReverseIt:
            List.reverse()
        return List
        
        

    def GetAttributeColumn(self):
        """parse user input and try to get the actual column"""
        if self.Attribute == "FScore":
            self.AttributeColumn = self.Columns.FScore
        elif self.Attribute == "Length":
            self.AttributeColumn = self.Columns.Length
        else:
            print "sam you need to set the rest of the attribute columns"
            sys.exit(1)
        
    def ParseAnnotations(self, FileName):
        Handle = open(FileName, "rb")
        True = 0
        False = 0
        for Line in Handle.xreadlines():
            if not Line.strip():
                continue
            if Line[0] == "#":
                self.Header = Line
                continue
            Bits = Line.strip().split("\t")
            Value = float(Bits[self.AttributeColumn])
            Protein = Bits[self.Columns.ProteinName]
            if Protein[:3] == "XXX":
                self.FalseDistribution.append(Value)
                False +=1
            else:
                self.TrueDistribution.append(Value)
                True += 1
        Handle.close()
        
                
    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:w:a:b:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-r":
                # -r results file(s)
                if not os.path.exists(Value):
                    print "** Error: couldn't find results file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.InputFile = Value
            elif Option == "-w":
                self.OutputFile = Value
            elif Option == "-a":
                self.Attribute = Value
            elif Option == "-f":
                self.BinSize = float(Value)
                self.FloatBin = 1
            elif Option == "-b":
                self.BinSize = int(Value)
            elif Option == "-c":
                self.CumulativeFlag = 1
        if not OptionsSeen.has_key("-r") or not OptionsSeen.has_key("-w") or not OptionsSeen.has_key("-a") or not OptionsSeen.has_key("-b"):
            print UsageInfo
            sys.exit(1)
            

if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(psyco not found - running in non-optimized mode)"
    Guido = AssassinClass()
    Guido.ParseCommandLine(sys.argv[1:])
    Guido.Main()                