UsageInfo="""Claudio.Quantitation.py
Given some quantitation data from Claudio, we will create a distribution for
each protein (median of peptides).  he informs me that the quantitation
is always plotted in ratio to light.  e.g. peptide1 ratio 0.88 means
light peptide1 1: heavy peptide1 0.88

Required Options:
 -r [FileName] File or directory of Input
 -w [FileName] Output from program for good data (n>=4)
 -b [FileName] output for bad data n<4
 -p [FileName] Peptide values output
 -m [int] Minimum number of observations for being good data
 -a [FileName] Absolute quantitation file for proteins
"""

"""For the moment
repeat measurements of the sampe peptide are averaged
heavy and light versions of the peptide are counted separately
N observations refers to the number of peptides (heavy and light separate) observed for a protein.

"""


import os
import sys
import getopt
import string
import math
import ResultsParser

class AbacusClass(ResultsParser.ResultsParser):
    def __init__ (self):
        self.OutputPath = None
        self.InspectInput = None
        self.AbsoluteQuantitaionFile = None
        self.AbsoluteQuantitaionDict = {} # protein->quant
        self.ProteinDict = {} # protein -> {} # second key is peptide, value is dummy
        self.PeptideToProtein = {} # peptide ->Protein
        self.LightPeptideQuantitation = {} # peptide -> [quantvalue1, quantvalue2, ..]
        self.HeavyPeptideQuantitation = {}
        self.ProteinCol = 5
        self.QuantCol = 6
        self.PeptideCol = 4
        self.MinimumObservations = 4
        ResultsParser.ResultsParser.__init__(self)
        self.TestProtein = "YAL041W"
    def Main(self):
        self.ProcessResultsFiles(self.InspectInput, self.ParseResults)
        self.ParseAbsoluteQuantitationFile()
        self.PrintPeptideFile()
        #self.GetProteinStats()


    def PrintPeptideFile(self):
        """Here we are going to print out stuff that looks like
        protein | peptide | nHeavy | nlight | min | max | average 
        """
        Handle = open(self.PeptideOutputFile, "wb")
        Handle.write("Protein\tPeptide\tnHeavy\tnLight\tMin\tMax\tAverage\tMedian\n")
        for Peptide in self.PeptideToProtein.keys(): # a list of all peptides
            HeavyQuantMeasurements = self.HeavyPeptideQuantitation.get(Peptide, [])
            LightQuantMeasurements = self.LightPeptideQuantitation.get(Peptide, [])
            nHeavy = len(HeavyQuantMeasurements)
            nLight = len(LightQuantMeasurements)
            (Average, Median, Min, Max) = self.AverageList(HeavyQuantMeasurements, LightQuantMeasurements)
            Handle.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(self.PeptideToProtein[Peptide], Peptide, nHeavy, nLight, Min, Max, Average, Median))
        Handle.close()

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
        return (Average, Median, Min, Max)
            

    def ParseAbsoluteQuantitationFile(self):
        """This is a file with the absolute protein molecules/cell
        defined in some experimental way, like an Erin O'Shea paper
        """
        Handle = open(self.AbsoluteQuantitaionFile, "rb")
        for Line in Handle.xreadlines():
            if Line[0] == "#":
                continue # header
            Bits = Line.strip().split("\t")
            ORF = Bits[0]
            MolsPerCellString = Bits[2]
            try:
                MolsPerCell = float(MolsPerCellString)
                self.AbsoluteQuantitaionDict[ORF] = MolsPerCell
            except:
                pass
        Handle.close()

    def GetProteinStats(self):
        Handle = open(self.OutputPath, "wb")
        Handle.write("Protein\tObservations\tAverage\tStdev\tMedian\tMolsPerCell\tMinimum\tMaximum\n")
        BadHandle = open(self.BadOutputPath, "wb")
        BadHandle.write("Protein\tObservations\tAverage\tStdev\tMedian\tMolsPerCell\tMinimum\tMaximum\n")
        for Protein in self.ProteinDict.keys():
            ListOfQuantitations= []
            #if Protein == self.TestProtein:
            #    print "GetProteinStats, %s"%self.TestProtein
            #    print "Peptides %s"%self.ProteinDict[Protein]
            for Peptide in self.ProteinDict[Protein].keys():
                ListOfQuantitations.append(self.PeptideQuant[Peptide])
                #if Protein == self.TestProtein:
                #    print "GetProteinStats. Peptide quantitation"
                #    print Peptide, self.PeptideQuant[Peptide]
            (Average, StandardDeviation, Median, Min, Max) = self.GetAvgStdevMedian(ListOfQuantitations)
            MolsPerCell = self.AbsoluteQuantitaionDict.get(Protein, "N/A")
            if len(ListOfQuantitations) < self.MinimumObservations:
                BadHandle.write("%s\t%s\t%s\t%s\t%s\t%s"%(Protein, len(ListOfQuantitations), Average, StandardDeviation, Median, MolsPerCell))
                BadHandle.write("\t%s\t%s\n"%( Min, Max))
            else:
                Handle.write("%s\t%s\t%s\t%s\t%s\t%s"%(Protein, len(ListOfQuantitations), Average, StandardDeviation, Median, MolsPerCell))
                Handle.write("\t%s\t%s\n"%( Min, Max))
        BadHandle.close()            
        Handle.close()
            
    def GetAvgStdevMedian(self, List):
        "for all the measurements avg, median and stdev"
        List.sort()
        Len = len(List)
        Min = 100000
        Max = 0
        Sum = 0.0 # float, always a float if you want to divide
        for Item in List:
            if Item > Max:
                Max = Item
            if Item < Min:
                Min = Item
            Sum += Item
        Average = Sum / Len
        SumSquaredDev = 0.0
        for Item in List:
            Deviation = Item - Average
            SumSquaredDev += (Deviation * Deviation)
        Variance = SumSquaredDev / Len
        Stdev = math.sqrt(Variance)
        # now median
        Median = -999999999999999999
        if Len %2 == 0: #even
            # average the two middle ones
            MedianIndex = Len /2
            Median = (List[MedianIndex] + List[MedianIndex - 1]) / 2.0
        else:
            MedianIndex = Len / 2 # intentional int division
            Median = List[MedianIndex]
        return (Average, Stdev, Median, Min, Max)

    def ParseResults(self, FileName):
        """Gets the annotation, puts the quantitation value in the peptide list (heavy or light)
        then keeps track of the peptide<->protein connection
        """
        Handle = open(FileName, "rb")
        for Line in Handle.xreadlines():
            if Line[0] == "#":
                continue
            Bits = Line.strip().split("\t")
            QuantValue = float(Bits[self.QuantCol])  # cast at the beginning
            Annotation = Bits[self.PeptideCol][2:-2]
            Annotation = Annotation.replace("C[160.03]", "C")  # get rid of that stupid bracket
            CleanedAminos = self.CleanAnnotation(Annotation)
            if not CleanedAminos:
                #this annotation is light!, that's why we did not clean it.
                Aminos = Annotation
                if not self.LightPeptideQuantitation.has_key(Aminos):
                    self.LightPeptideQuantitation[Aminos] = []
                self.LightPeptideQuantitation[Aminos].append(QuantValue)
            else:
                Aminos = CleanedAminos
                if not self.HeavyPeptideQuantitation.has_key(Aminos):
                    self.HeavyPeptideQuantitation[Aminos] = []
                self.HeavyPeptideQuantitation[Aminos].append(QuantValue)
            Protein = Bits[self.ProteinCol] # Claudio ensures me that this is a single entry
            ## map peptide to protein
            self.PeptideToProtein[Aminos] = Protein
            ## map protein to peptide
            if not self.ProteinDict.has_key(Protein):
                self.ProteinDict[Protein] = {}
            if not self.ProteinDict[Protein].has_key(Aminos):
                self.ProteinDict[Protein][Aminos] = 1
                    
    def CleanAnnotation(self, Annotation):
        "get rid of any brackets and numbers.  only works when prefix and suffix havae been stripped!!!!!!"
        if Annotation.find("[") == -1:
            #this annotation did not need any cleaning
            return None
        ToReturn = ""
        Index = 0
        while Index < len(Annotation):
            if Annotation[Index] in string.letters:
                ToReturn += Annotation[Index]
            Index += 1
        return ToReturn
        
    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "w:r:b:m:a:p:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-w":
                # -r results file(s)
                self.OutputPath = Value
            elif Option == "-r":
                self.InspectInput = Value
            elif Option == "-b":
                self.BadOutputPath = Value
            elif Option == "-m":
                self.MinimumObservations = int (Value)
            elif Option == "-a":
                self.AbsoluteQuantitaionFile = Value
            elif Option == "-p":
                self.PeptideOutputFile = Value
            else:
                print "Option %s is not understood"%Option
                print UsageInfo
                sys.exit(1)
        if not OptionsSeen.has_key("-w") or not OptionsSeen.has_key("-r"):
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
