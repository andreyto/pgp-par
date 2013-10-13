"""
Compute features for scoring of peptide matches.
"""

import os
import sys
import getopt
import traceback
import struct
import random
import ResultsParser
import PyInspect
from Utils import *
Initialize()

class MQScoreTrainer(ResultsParser.ResultsParser, ResultsParser.SpectrumOracleMixin):
    def __init__(self):
        self.CorpusDir = r"e:\ms\TrainingCorpusUnfiltered"
        self.OneMatchPerSpectrumFlag = 1
        self.MinimumPeptideLength = 6
        self.MinOldFScore = -5
        self.QuickParseFlag = 0
        self.QuickParseLineLimit = 10
        ResultsParser.ResultsParser.__init__(self)
        ResultsParser.SpectrumOracleMixin.__init__(self)
    def ParseCommandLine(self, Arguments):
        (Options, Args) = getopt.getopt(Arguments, "c:Q")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-c":
                self.CorpusDir = Value
            elif Option == "-Q":
                self.QuickParseFlag = 1
    def FixSpectrumPath(self, Path):
        if Path[1] == ":":
            return Path
        FileName = Path.replace("/", "\\").split("\\")[-1]
        return os.path.join(r"e:\ms\lens\spectra", FileName)
    def Main(self):
        self.ComputeFeaturesForFiles()
    def ComputeFeaturesForFiles(self):
        Header = "#File\tValidFlag\t"
        Header += "LegacyMQScore\tLegacyFScore\t"
        Header += "Length\tTotalCutScore\tMeanCutScore\tTotalCutScore1\tMeanCutScore1\tTotalCutScore2\tMeanCutScore2\tMedianCutScore2\tMedianCutScore1\tMedianCutScore\t"
        Header += "FractionY\tFractionB\tFractionY1\tFractionB1\tBYStrong\tYstrong\tBstrong\tBYStrongWeighted\tYStrongWeighted\tBStrongWeighted\tIntensityBY\tIntensityY\tIntensityB\t"
        Header += "IntensityBYSeries\tIntensityYSeries\tIntensityBSeries\tConvolve0\tConvolve+\tNTT\t"
        Header += "LogLength\tLogLengthB\tLogLengthC\t"
        self.FeatureFile2 = open("MQScoreFeatures2.txt", "wb")
        self.FeatureFile2.write(Header + "\n")
        self.FeatureFile3 = open("MQScoreFeatures3.txt", "wb")
        self.FeatureFile3.write(Header + "\n")
        self.ComputeFeaturesForDir(self.CorpusDir)
        self.FeatureFile2.close()
        self.FeatureFile3.close()        
    def ComputeFeaturesForDir(self, CorpusDir):
        FileNames = []
        for FileName in os.listdir(CorpusDir):
            if FileName[:5] == "H293b":
                continue #%%% TEMP
            FileNames.append(FileName)            
        random.shuffle(FileNames)
        if self.QuickParseFlag:
            FileNames = FileNames[:10]
        for FileName in FileNames:
            FilePath = os.path.join(CorpusDir, FileName)
            if os.path.isdir(FilePath):
                self.ComputeFeaturesForDir(FilePath)
            else:
                self.ComputeFeaturesForFile(FilePath)
    def ComputeFeaturesForFile(self, FilePath):
        File = open(FilePath, "rb")
        print FilePath
        OldSpectrum = None
        LineNumber = 0
        for FileLine in File.xreadlines():
            LineNumber += 1
            if self.QuickParseFlag and LineNumber > self.QuickParseLineLimit:
                break
            if LineNumber % 1000 == 0:
                print "Line %s..."%LineNumber
            if FileLine[0] == "#":
                continue
            Bits = FileLine.split("\t")
            if len(Bits) < 2:
                continue
            Spectrum = (Bits[0], Bits[1])
            if self.OneMatchPerSpectrumFlag:
                if Spectrum == OldSpectrum:
                    continue
            OldSpectrum = Spectrum
            OldMQScore = float(Bits[self.Columns.MQScore])
            OldFScore = 0.3 * OldMQScore + 1.5 * float(Bits[self.Columns.DeltaScore]) / 2.194
            # skip crappy matches:
            if OldFScore < self.MinOldFScore:
                continue
            Charge = int(Bits[self.Columns.Charge])
            Annotation = Bits[self.Columns.Annotation]
            Peptide = GetPeptideFromModdedName(Annotation)
            if len(Peptide.Aminos) < self.MinimumPeptideLength:
                continue
            ProteinName = Bits[self.Columns.ProteinName]
            SpectrumPath = self.FixSpectrumPath(Bits[self.Columns.SpectrumFile])
            SpectrumFilePos = int(Bits[self.Columns.FileOffset])
            Spec = PyInspect.Spectrum(SpectrumPath, SpectrumFilePos)
            Features = Spec.GetMatchFeatures(Annotation)
            ValidFlag = 1
            if ProteinName[:3] == "XXX":
                ValidFlag = 0
            Str = "%s:%s %s\t%s\t"%(SpectrumPath, SpectrumFilePos, Annotation, ValidFlag)
            Str += "%s\t%s\t"%(OldMQScore, OldFScore)
            for Feature in Features:
                Str += "%s\t"%Feature
            if Charge < 3:
                self.FeatureFile2.write(Str + "\n")
            else:
                self.FeatureFile3.write(Str + "\n")
        
if __name__ == "__main__":
    try:
        Trainer = MQScoreTrainer()
        Trainer.ParseCommandLine(sys.argv[1:])
        Trainer.Main()
    except:
        traceback.print_exc()
        print ">>> ENTER <<<"
        sys.stdin.readline()