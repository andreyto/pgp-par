"""
Test charge-correction.  Run this *after* training parent mass correction.
"""
import os
import sys
import random
import getopt
import traceback
import time
import ResultsParser
import PyInspect
from Utils import *
Initialize()

OracleDir = r"e:\ms\TrainingCorpus"
QuickParseLineCount = 100

class ChargeCorrector(ResultsParser.ResultsParser):
    def __init__(self):
        ResultsParser.ResultsParser.__init__(self)
        self.QuickParseFlag = 0
        self.SVMFlag = 1
        self.ComputeFeaturesFlag = 0
        #self.OracleFiles = ("OMICS2002.pv", "OMICS2004.pv", "HEK.pv")
    def ComputeAndReportStats(self):
        """
        Compute cutoffs such that "if spectrum score is above X, then spectrum has a 75% of being charge n"
        """
        self.RightCallScores = ([], [], [], [])
        self.WrongCallScores = ([], [], [], [])
        self.OutputFile = open("CCStats.txt", "wb")
        for OracleFileName in os.listdir(OracleDir):
            OraclePath = os.path.join(OracleDir, OracleFileName)
            self.ComputeStatsOnFile(OraclePath)
        for Charge in (1, 2, 3):
            SortedScores = []
            for Score in self.RightCallScores[Charge]:
                SortedScores.append((Score, random.random(), 1))
            for Score in self.WrongCallScores[Charge]:
                SortedScores.append((Score, random.random(), 0))
            SortedScores.sort()
            SortedScores.reverse()
            CumulativeFalse = 0
            CumulativeTrue = 0
            for Tuple in SortedScores:
                if Tuple[-1]:
                    CumulativeTrue += 1
                else:
                    CumulativeFalse += 1
                OddsTrue = CumulativeTrue / float(CumulativeTrue + CumulativeFalse)
                Str = "%s\t%s\t%s\t%s\t%s\t"%(Charge, Tuple[0], CumulativeTrue, CumulativeFalse, OddsTrue)
                Str += "%s\t"%(CumulativeTrue / float(len(self.RightCallScores[Charge])))
                self.OutputFile.write(Str + "\n")
            self.OutputFile.write("\n\n")
    def ComputeStatsOnFile(self, OracleFileName):
        File = open(OracleFileName, "rb")
        for FileLine in File.xreadlines():
            Bits = FileLine.split("\t")
            if FileLine[0] == "#" or not FileLine.strip():
                continue
            try:
                PValue = float(Bits[self.Columns.PValue])
            except:
                traceback.print_exc()
                print Bits
                continue
            Peptide = GetPeptideFromModdedName(Bits[self.Columns.Annotation])
            Charge = int(Bits[self.Columns.Charge])
            FilePos = int(Bits[self.Columns.FileOffset])
            FileName = Bits[self.Columns.SpectrumFile]
            FileName = self.FixPath(FileName)
            Spectrum = PyInspect.Spectrum(FileName, FilePos)
            Tuple = Spectrum.CorrectCharge()
            # Set values:
            for TestCharge in (1, 2, 3):
                if TestCharge == Charge:
                    ValidFlag = 1
                else:
                    ValidFlag = 0
                Score = Tuple[TestCharge - 1]
                if ValidFlag:
                    self.RightCallScores[TestCharge].append(Score)
                else:
                    self.WrongCallScores[TestCharge].append(Score)
    def ComputeAndReportAccuracy(self):
        self.CorrectCount = 0
        self.IncorrectCount = 0
        self.OutputFile = open("CCResults.txt", "wb")
        for OracleFileName in os.listdir(OracleDir):
            OraclePath = os.path.join(OracleDir, OracleFileName)
            if os.path.isdir(OraclePath):
                continue
            self.ComputeAccuracyOnFile(OraclePath)
        print "In all: %s correct, %s incorrect"%(self.CorrectCount, self.IncorrectCount)
    def ComputeAccuracyOnFile(self, OracleFileName):
        File = open(OracleFileName, "rb")
        for FileLine in File.xreadlines():
            Bits = FileLine.split("\t")
            if FileLine[0] == "#" or not FileLine.strip():
                continue
            try:
                PValue = float(Bits[self.Columns.PValue])
            except:
                traceback.print_exc()
                print Bits
                continue
            Peptide = GetPeptideFromModdedName(Bits[self.Columns.Annotation])
            Charge = int(Bits[self.Columns.Charge])
            FilePos = int(Bits[self.Columns.FileOffset])
            FileName = Bits[self.Columns.SpectrumFile]
            FileName = self.FixPath(FileName)
            Spectrum = PyInspect.Spectrum(FileName, FilePos)
            CorrectedCharge = Spectrum.CorrectCharge()
            if CorrectedCharge == Charge:
                ValidFlag = 1
                self.CorrectCount += 1
            else:
                ValidFlag = 0
                self.IncorrectCount += 1
            Str = "%s\t%s\t%s\t%s\t%s\t"%(FileName, FilePos, Charge, CorrectedCharge, ValidFlag)
            self.OutputFile.write(Str + "\n")
    def ComputeAllCCFeatures(self):
        self.OutputFiles = [None]
        for Charge in (1, 2):
            FileName = "CCFeatures%d.txt"%Charge
            self.OutputFiles.append(open(FileName, "wb"))
        for OracleFileName in os.listdir(OracleDir):
            OraclePath = os.path.join(OracleDir, OracleFileName)
            if os.path.isdir(OraclePath):
                continue
            self.ComputeCCFeatures(OraclePath)
        # Close our files and tidy up:
        for File in self.OutputFiles:
            if File:
                File.close()
        self.OutputFiles = []
    def ComputeCCFeatures(self, OracleFileName):
        File = open(OracleFileName, "rb")
        LineNumber = 0
        for FileLine in File.xreadlines():
            LineNumber += 1
            if self.QuickParseFlag and LineNumber > QuickParseLineCount:
                break
            Bits = FileLine.split("\t")
            if FileLine[0] == "#" or not FileLine.strip():
                continue
            try:
                PValue = float(Bits[self.Columns.PValue])
            except:
                traceback.print_exc()
                print Bits
                continue
            Peptide = GetPeptideFromModdedName(Bits[self.Columns.Annotation])
            Charge = int(Bits[self.Columns.Charge])
            FilePos = int(Bits[self.Columns.FileOffset])
            FileName = Bits[self.Columns.SpectrumFile]
            FileName = self.FixPath(FileName)
            Spectrum = PyInspect.Spectrum(FileName, FilePos)
            TupleCh1 = Spectrum.GetCCFeatures(1)
            TupleCh2 = Spectrum.GetCCFeatures(0)
            # Write out the charge-1 tuple:
            if (Charge == 1):
                ValidFlag = 1
            else:
                ValidFlag = 0
            Str = "%s:%s (%d)\t"%(FileName, FilePos, Charge)
            Str += "%s\t"%ValidFlag
            for Value in TupleCh1:
                Str += "%s\t"%Value
            self.OutputFiles[1].write(Str + "\n")
            # Write out the charge-2 tuple:
            if (Charge == 2):
                ValidFlag = 1
            else:
                ValidFlag = 0
            Str = "%s:%s (%d)\t"%(FileName, FilePos, Charge)
            Str += "%s\t"%ValidFlag
            for Value in TupleCh2:
                Str += "%s\t"%Value
            self.OutputFiles[2].write(Str + "\n")
    def CCVerbose(self):
        FileNames = [r"C:\ms\OMICS04\0700-1010mz_01_dta.mgf"]
        ByteOffsets = [2235446]
        for Index in range(len(FileNames)):
            FileName = FileNames[Index]
            FilePos = ByteOffsets[Index]
            Spectrum = PyInspect.Spectrum(FileName, FilePos)
            Tuple = Spectrum.CorrectCharge()
            print Tuple
    def FixPath(self, Path):
        if sys.platform == "win32":
            return Path
        else:
            return Path.replace("\\", "/").split("/")[-1]
        #return "c" + Path[1:]
        
    def DetermineScoreCutoffs(self):
        """
        Determine a score confidence interval:
        First model:
        - For what score can we be 95% confident that a spectrum is *not* singly-charged?
        - For what score can we be 95% confident that a spectrum is *not* multiply-charged?
        Second model:
        - For what score can we be 95% confident that a spectrum is *not* doubly-charged?
        - For what score can we be 95% confident that a spectrum is *not* triply-charged?
        """
        # entries of the form (score, actual charge):
        Model1SortedResults = []
        Model2SortedResults = []
        # Iterate over the spectra.  Compute the scores from model 1 and model 2.  Accumulate
        # entries in the sorted lists:
        for OracleFileName in os.listdir(OracleDir):
            OraclePath = os.path.join(OracleDir, OracleFileName)
            if os.path.isdir(OraclePath):
                continue
            print ">>Parse %s..."%OraclePath
            File = open(OraclePath, "rb")
            LineNumber = 0
            for FileLine in File.xreadlines():
                LineNumber += 1
                if LineNumber % 1000 == 0:
                    print "%s..."%LineNumber
                Bits = FileLine.split("\t")
                if FileLine[0] == "#" or not FileLine.strip():
                    continue
                try:
                    PValue = float(Bits[self.Columns.PValue])
                except:
                    traceback.print_exc()
                    print Bits
                    continue
                Peptide = GetPeptideFromModdedName(Bits[self.Columns.Annotation])
                Charge = int(Bits[self.Columns.Charge])
                FilePos = int(Bits[self.Columns.FileOffset])
                FileName = Bits[self.Columns.SpectrumFile]
                FileName = self.FixPath(FileName)
                Spectrum = PyInspect.Spectrum(FileName, FilePos)
                (ModelScore1, ModelScore2) = Spectrum.CorrectCharge(1)
                Model1SortedResults.append((ModelScore1, Charge))
                Model2SortedResults.append((ModelScore2, Charge))
        Model1SortedResults.sort()
        Model2SortedResults.sort()
        print ">> Now determine cutoffs:"
        # Now that the sorted lists are populated, let's determine the score cutoffs to use:
        Cutoffs1 = self.DetermineCCScoreCutoffs(Model1SortedResults, 1, 2)
        #print ">>> Charge correction cutoffs, model 1:", Cutoffs1
        Cutoffs2 = self.DetermineCCScoreCutoffs(Model2SortedResults, 2, 3)
        #print ">>> Charge correction cutoffs, model 2:", Cutoffs2
    def DetermineCCScoreCutoffs(self, SortedList, ChargeA, ChargeB):
        # Plot a histogram:
        HistogramA = {}
        HistogramB = {}
        MinBin = 9999
        MaxBin = -9999
        for (Score, Charge) in SortedList:
            Bin = int(round(Score * 10))
            MinBin = min(MinBin, Bin)
            MaxBin = max(MaxBin, Bin)
            if Charge == ChargeA:
                HistogramA[Bin] = HistogramA.get(Bin, 0) + 1
            else:
                HistogramB[Bin] = HistogramB.get(Bin, 0) + 1
        PlotFile = open("CCScoreHistogram.%s.txt"%(ChargeA), "wb")
        for Bin in range(MinBin, MaxBin + 1):
            Str = "%s\t%s\t%s\t"%(Bin, HistogramA.get(Bin, 0), HistogramB.get(Bin, 0))
            PlotFile.write(Str + "\n")
        # Determine thresholds:
        SortedListA = []
        SortedListB = []
        for (Score, Charge) in SortedList:
            if Charge == ChargeA:
                SortedListA.append(Score)
            else:
                SortedListB.append(Score)
        SortedListA.sort()
        SortedListB.sort()
        print ">>DetermineCCScoreCutoffs(%s, %s):"%(ChargeA, ChargeB)
        Percentile5 = SortedListA[int(round(len(SortedListA) * 0.05))]
        Percentile95 = SortedListA[int(round(len(SortedListA) * 0.95))]
        Percentile1 = SortedListA[int(round(len(SortedListA) * 0.01))]
        Percentile99 = SortedListA[int(round(len(SortedListA) * 0.99))]
        print ">>ChargeA: 5th percentile is %s, 95th percentile is %s"%(Percentile5, Percentile95)
        print ">>ChargeA: 1th percentile is %s, 99th percentile is %s"%(Percentile1, Percentile99)
        Percentile5 = SortedListB[int(round(len(SortedListB) * 0.05))]
        Percentile95 = SortedListB[int(round(len(SortedListB) * 0.95))]
        Percentile1 = SortedListB[int(round(len(SortedListB) * 0.01))]
        Percentile99 = SortedListB[int(round(len(SortedListB) * 0.99))]
        print ">>ChargeB: 5th percentile is %s, 95th percentile is %s"%(Percentile5, Percentile95)
        print ">>ChargeB: 1th percentile is %s, 99th percentile is %s"%(Percentile1, Percentile99)
    def ParseCommandLine(self, Arguments):
        (Options, Args) = getopt.getopt(Arguments, "flQ")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-f":
                self.ComputeFeaturesFlag = 1
            elif Option == "-l":
                self.SVMFlag = 0
            elif Option == "-Q":
                self.QuickParseFlag = 1
    def Main(self):
        print "-- Compute features..."
        if self.ComputeFeaturesFlag:
            self.ComputeAllCCFeatures()
        print "-- Train models..."
        if SVMFlag:
            ModelName = "SVM"
        else:
            ModelName = "LDA"
        for Charge in (1, 2):
            Command = "TrainMachineLearner.py -r CCFeatures%s.txt -m %s -w CC%s%s"%(Charge, ModelName, ModelName, Charge)
            os.system(Command)
        print "-- Determine score cutoffs..."
        self.DetermineScoreCutoffs()
        StartComputeAccuracy = time.clock()
        print "-- Compute accuracy..."
        self.ComputeAndReportAccuracy()
        EndTime = time.clock()
        print "Elapsed time to compute accuracy: %s"%(EndTime - StartComputeAccuracy)
        print "Overall elapsed time: %s"%(EndTime - StartAll)
        
if __name__ == "__main__":
    SVMFlag = 1
    StartAll = time.clock()
    CC = ChargeCorrector()
    CC.ParseCommandLine(sys.argv[1:])
    CC.Main()
##    CC.ComputeAndReportStats()