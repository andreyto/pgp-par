#!/usr/bin/env python

"""
Takes Post PValue results (1 hit per scan, decoy hits reported) and determines
instantaneous FDR at different pvalue cutoffs.  Writes out each forward DB hit with added column of localFDR based on pvalue.
"""
import os
import sys
import getopt
import InspectResults


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

        inspectParser = InspectResults.Parser( self.InputDir )

        for result in inspectParser:
            Charge = result.Charge
            Score = result.PValue
            Index = int(Score/(self.BinWidth))

            if Charge <= 2:
                (TP,FP) = self.IScoreCounts_1_2.setdefault(Index,(0,0))
            else:
                (TP,FP) = self.IScoreCounts_3.setdefault(Index,(0,0))
                    
            Protein = result.ProteinName
            if self.VerboseFlag:
                print result
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
                
        self.Header = inspectParser.header

        if self.VerboseFlag:
            self.PrintResults()

        self.ComputeAndWriteTable()
        self.WriteResults()
        
    def WriteResults(self):

        Scans = 0
        GoodScans = 0
        inspectParser = InspectResults.Parser( self.InputDir, inputMirrorTo=self.OutputDir)
        for result in inspectParser:
            if not Scans:
                inspectParser.mirrorOutHandle.write(self.Header)

            PValue = result.PValue
            Index = int(PValue/self.BinWidth)
            Charge = result.Charge
            Protein = result.ProteinName
            if Protein[0:3] == "XXX":
                continue

            if Charge <= 2:
                result.LFDR = self.LFDR_1_2[Index]
            else:
                result.LFDR = self.LFDR_3[Index]
                    
            Scans += 1
            if result.LFDR < self.PValueCutoff:
                GoodScans += 1
            inspectParser.mirrorOutHandle.write( str(result) )
                        
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
