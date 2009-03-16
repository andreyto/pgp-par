"""
Takes Post PValue results (1 hit per scan, decoy hits reported) and determines
instantaneous FDR at different pvalue cutoffs.  Writes out each forward DB hit with added column of localFDR based on pvalue.
"""
import os
import sys
import getopt
import ResultsParser


UsageInfo = """ProduceFDRBins.py
Takes Post PValue results (1 hit per scan, decoy hits reported) and determines
instantaneous FDR at different pvalue cutoffs.  Writes out each forward DB hit with added column of localFDR based on pvalue.


-r [Directory] Directory containing PValued Results
-w [Directory] Directory to write out Results
-p [FileName] File to write tabular summary

-b [Number] Size of bin in pvalue Units (default 0.01)
-v Verbose
"""


class BinMaker:

    def __init__(self):
        self.BinWidth = 0.01
        self.VerboseFlag = 0
        self.PValueCutoff = 0.01
        self.Header = ""

        self.IScoreCounts_1_2 = {} # ScoreCounts[XX.X] = (#TP,#FP)
        self.IScoreCounts_3 = {} # ScoreCounts[XX.X] = (#TP,#FP)
        self.LFDR_1_2 = {}
        self.LFDR_3 = {}
        self.FileNames = []

    def Main(self):

        FullFileNames = []
        for FileName in os.listdir(self.InputDir):
            self.FileNames.append(FileName)
            FullFileNames.append(os.path.join(self.InputDir,FileName))

        self.FileNames.sort()
        FullFileNames.sort()

        for i in range(0,len(FullFileNames)):
            FileName = FullFileNames[i]

            print "(%s/%s) Reading Values in %s"%(i,len(FullFileNames)-1,FileName)
                
            if os.path.isdir(FileName):
                continue
            

            File = open(FileName,'r')

            for Line in File:
                Line = Line.strip()
                if Line == "" or Line[0] == "#":
                    if Line[0] == "#":
                        self.Header = Line
                        #print self.Header
                    continue 

                Bits = Line.split("\t")

                Charge = int(Bits[ResultsParser.Columns.Charge])
                #Score = float(Bits[ResultsParser.Columns.FScore])
                Score = float(Bits[ResultsParser.Columns.PValue])
                Index = int(Score/(self.BinWidth))

                if Charge <= 2:
                    (TP,FP) = self.IScoreCounts_1_2.setdefault(Index,(0,0))
                else:
                    (TP,FP) = self.IScoreCounts_3.setdefault(Index,(0,0))
                    
                Protein = Bits[ResultsParser.Columns.ProteinName]
                if self.VerboseFlag:
                    print Bits
                    print Protein
                    print self.IScoreCounts_1_2[Index]
                    print self.IScoreCounts_3[Index]
                if Protein[0:3] == "XXX":
                    if Charge <= 2:
                        self.IScoreCounts_1_2[Index] = (TP,FP+1)
                    else:
                        self.IScoreCounts_3[Index] = (TP,FP+1)
                else:
                    if Charge <= 2:
                        self.IScoreCounts_1_2[Index] = (TP+1,FP)
                    else:
                        self.IScoreCounts_3[Index] = (TP+1,FP)

                if self.VerboseFlag:
                    print self.IScoreCounts_1_2[Index]
                    print self.IScoreCounts_3[Index]
                    raw_input()
                
            File.close()

        if self.VerboseFlag:
            self.PrintResults()

        self.ComputeAndWriteTable()
        self.WriteResults()
        
    def WriteResults(self):

        Scans = 0
        GoodScans = 0
        for i in range(0,len(self.FileNames)):
            FileName = self.FileNames[i]
    
            print "(%s/%s) Writing to %s"%(i,len(self.FileNames)-1,FileName)
            InFile = open(os.path.join(self.InputDir,FileName),'r')
            OutFile = open(os.path.join(self.OutputDir,FileName),'w')
            OutFile.write(self.Header)
            for Line in InFile:
                Line = Line.strip()
                if Line == "" or Line[0] == "#":
                    continue
                Bits = Line.split("\t")
                #if len(Bits) != 25:
                #    print Bits
                #    print len(Bits)
                #    raw_input()
                PValue = float(Bits[ResultsParser.Columns.PValue])
                Index = int(PValue/self.BinWidth)
                Charge = int(Bits[ResultsParser.Columns.Charge])
                Protein = Bits[ResultsParser.Columns.ProteinName]
                if Protein[0:3] == "XXX":
                    continue

                if Charge <= 2:
                    Scans += 1
                    OutFile.write(Line + "\t")
                    
                    EntryCount = len(Bits)
                    while(EntryCount < ResultsParser.Columns.LFDR):
                        OutFile.write("\t")
                        EntryCount += 1
                    OutFile.write(str(self.LFDR_1_2[Index]) + "\n")
                    
                    if self.LFDR_1_2[Index] < self.PValueCutoff:
                        GoodScans += 1
                else:
                    Scans += 1
                    OutFile.write(Line + "\t")
                    EntryCount = len(Bits)
                    while(EntryCount < ResultsParser.Columns.LFDR):
                        OutFile.write("\t")
                        EntryCount += 1
                    OutFile.write(str(self.LFDR_3[Index]) + "\n")
                    
                    if self.LFDR_3[Index] < self.PValueCutoff:
                        GoodScans += 1
                        
            OutFile.close()
            InFile.close()
        print "%s of %s (%s) forward hits retained"%(GoodScans,Scans,float(GoodScans)/float(Scans))
                
    def ComputeAndWriteTable(self):

        File = open(self.OutputFileName,'w')

        Keys = self.IScoreCounts_1_2.keys()
        Keys.sort()
        Keys = Keys[::-1]

        File.write("ChargeState 1 or 2\n")
        File.write("StartBin\tEndBin\tI_FDR\tC_FDR\tTP\tFP\tTotalScans\tTotalFP\n")
        TotalScans = 0
        TotalFP = 0
        #PrevKey = -1
        for k in Keys:
            (TP,FP) = self.IScoreCounts_1_2[k]

            TotalScans += TP + FP
            TotalFP += FP

            self.LFDR_1_2[k] = float(FP)/float(TP+FP)
            File.write(str(k*self.BinWidth) + "\t"+ str((k+1)*self.BinWidth) + "\t" + str(float(FP)/float(TP+FP)) + "\t" +  str(float(TotalFP)/float(TotalScans)) + "\t" + str(TP) + "\t" + str(FP) + "\t" + str(TotalScans) + "\t" + str(TotalFP) + "\n")


        Keys = self.IScoreCounts_3.keys()
        Keys.sort()
        Keys = Keys[::-1]
        File.write("ChargeState 3\n")
        File.write("StartBin\tEndBin\tI_FDR\tC_FDR\tTP\tFP\tTotalScans\tTotalFP\n")
        TotalScans = 0
        TotalFP = 0
        #PrevKey = -1
        for k in Keys:
            (TP,FP) = self.IScoreCounts_3[k]

            TotalScans += TP + FP
            TotalFP += FP
            self.LFDR_3[k] = float(FP)/float(TP+FP)
            
            File.write(str(k*self.BinWidth) + "\t"+ str((k+1)*self.BinWidth) + "\t" + str(float(FP)/float(TP+FP)) + "\t" +  str(float(TotalFP)/float(TotalScans)) + "\t" + str(TP) + "\t" + str(FP) + "\t" + str(TotalScans) + "\t" + str(TotalFP) + "\n")

        File.close()

    def PrintResults(self):
        Keys = self.IScoreCounts.keys()
        Keys.sort()

        print("StartBin EndBin I_FDR C_FDR TP FP")
        TotalScans = 0
        TotalFP = 0
        for k in Keys:
            (TP,FP) = self.IScoreCounts[k]

            TotalScans += TP + FP
            TotalFP += FP
            print(str(k*self.BinWidth) + " "+ str((k+1)*self.BinWidth) + " " + str(float(FP)/float(TP+FP)) + " " +  str(float(TotalFP)/float(TotalScans)) + " " + str(TP) + " " + str(FP))


    def ParseCommandLine(self,Arguments):
        (Options,Args) = getopt.getopt(Arguments, "r:w:b:p:v")
        OptionsSeen = {}
        for(Option, Value) in Options:
            OptionsSeen[Option] = 1

            if Option == "-r":
                if not os.path.isdir(Value):
                    print "**Error: %s is not a valid directory"%Value
                    print UsageInfo
                    sys.exit(1)
                self.InputDir = Value
            elif Option == "-p":
                
                self.OutputFileName = Value
            elif Option == "-w":
                if not os.path.isdir(Value):
                    os.mkdir(Value)
                self.OutputDir = Value
            elif Option == "-b":
                self.BinWidth = float(Value)

            elif Option == "-v":
                self.VerboseFlag = 1

        if not OptionsSeen.has_key("-r") or not OptionsSeen.has_key("-w") or not OptionsSeen.has_key("-p"):
            print "**Error: Not enough arguments given"
            print UsageInfo
            sys.exit(1)
                    


if __name__ == "__main__":
    Dummy = BinMaker()
    Dummy.ParseCommandLine(sys.argv[1:])
    Dummy.Main()
