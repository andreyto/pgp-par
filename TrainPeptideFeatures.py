"""
TrainPeptideFeatures:
Given the features produced by ComputePeptideFeatures, let's make a list of
peptides which we're confident in.
"""
import os
import sys
import struct
import traceback
import getopt
import MSSpectrum
import PyInspect
import random
import shutil
import Learning
import BasicStats
import ResultsParser
import SpectralSimilarity
from ComputePeptideFeatures import FormatBits
from Utils import *
Initialize()

try:
    from numpy import *
    import numpy.linalg
    FloatType = float
    MatrixMultiply = dot
except:
    print "** Warning: Unable to import Numpy.  Logit training not available"

random.seed(1)

UsageInfo = """
TrainPeptideFeatures: Train a model to differentiate between CORRECT and
INCORRECT peptides identified in a database search.
-m [TYPE] Model type
-u [FILENAME] Feature filename
-v [FILENAME] Output feature path

Optional:
-R [FILENAME]: Filename to report ROC curve
-f [FLAG] Feature selection mode
"""

class PeptideModeler:
    def __init__(self):
        self.ModelType = None
        self.TrainingSetDBRatio = 1.0
        self.TestingSetDBRatio = 1.0
        self.ReadModelFilePath = None
        self.WriteModelFilePath = "Peptide.mdl" #default
        self.InputFeaturePath = "PeptideFeatures.txt" # default
        self.OutputFeaturePath = "ScoredPeptides.txt" # default
        self.ModelTestFilePath = None
        self.FeatureSelectionFlag = 0
        self.ReportROCPath = None
        self.SelectedFeatures = None
        self.HeaderLines = []
        self.FilterFactor = None
        self.MinFScore = 0.5
    def ParseFeatureFile(self, FilePath, FeatureSet2, FeatureSet3, DBRatio):
        """
        Initialize the FeatureSet objects, by parsing features from the specified FilePath.
        """
        File = open(FilePath, "rb")
        # Parse the top header lines specially:
        HeaderLine = File.readline()
        self.HeaderLines.append(HeaderLine)
        Bits = HeaderLine.strip().split("\t")
        print "Grab feature names:", Bits, FormatBits.FirstFeature
        for BitIndex in range(len(Bits)):
            FeatureNumber = BitIndex - FormatBits.FirstFeature
            if FeatureNumber in self.FeaturesAll:
                self.FeatureNames[FeatureNumber] = Bits[BitIndex]
            #if BitIndex >= FormatBits.FirstFeature:
            #    self.FeatureNames[BitIndex - FormatBits.FirstFeature] = Bits[BitIndex]
        # Iterate over non-header lines:
        LineNumber = 0
        for FileLine in File.xreadlines():
            LineNumber += 1
            if LineNumber % 5000 == 0:
                print "Line %s: %s+%s vectors so far"%(LineNumber, len(FeatureSet2.AllVectors), len(FeatureSet3.AllVectors))
            if FileLine[0] == "#":
                self.HeaderLines.append(FileLine)
                continue # skip comment line
            if not FileLine.strip():
                continue # skip blank line
            Bits = FileLine.replace("\r","").replace("\n","").split("\t")
            # If there are TOO MANY bits, then discard the extras:
            Bits = Bits[:FormatBits.LastFeature + 1]
            try:
                TrueFlag = int(Bits[FormatBits.TrueProteinFlag])
                Charge = int(Bits[FormatBits.Charge])
            except:
                continue # skip; not a valid instance line
            FScore = float(Bits[FormatBits.BestFScore])
            if self.MinFScore != None and FScore < self.MinFScore:
                continue
            # If our filter-factor is set, then SKIP some percentage of the input lines.
            if self.FilterFactor != None:
                if random.random() > self.FilterFactor:
                    continue
            Vector = Learning.FeatureVector()
            if Charge > 2:
                FeatureSet = FeatureSet3
            else:
                FeatureSet = FeatureSet2
            try:
                for FeatureIndex in range(FormatBits.FirstFeature, FormatBits.LastFeature + 1):
                    FeatureNumber = FeatureIndex - FormatBits.FirstFeature
                    if FeatureNumber not in self.FeaturesAll:
                        continue # skip features we're not allowed to use
                    if FeatureIndex < len(Bits) and Bits[FeatureIndex].strip() and Bits[FeatureIndex] != "None":
                        Vector.Features.append(float(Bits[FeatureIndex]))
                    else:
                        Vector.Features.append(0)
                Vector.FileBits = Bits
                Vector.TrueFlag = TrueFlag
                if TrueFlag:
                    FeatureSet.TrueVectors.append(Vector)
                else:
                    FeatureSet.FalseVectors.append(Vector)
                FeatureSet.AllVectors.append(Vector)
            except:
                traceback.print_exc()
                print "** Error on line %s column %s of feature file"%(LineNumber, FeatureIndex)
                print Bits
        File.close()
        # Initialize counts:
        FeatureSet2.SetCounts()
        FeatureSet2.GetPriorProbabilityFalse(DBRatio)
        FeatureSet3.SetCounts()
        FeatureSet3.GetPriorProbabilityFalse(DBRatio)
        print "CHARGE 1,2: Read in %s true and %s false features"%(FeatureSet2.TrueCount, FeatureSet2.FalseCount)
        print "CHARGE  3+: Read in %s true and %s false features"%(FeatureSet3.TrueCount, FeatureSet3.FalseCount)
    def ParseTrainingAndTestingSets(self):
        # Load in features for a collection of TRUE and FALSE instances.
        File = open(self.InputFeaturePath, "rb")
        self.FeatureNames = {}
        FeatureCount = FormatBits.LastFeature - FormatBits.FirstFeature + 1
        self.FeaturesAll = range(FeatureCount)
        if self.SelectedFeatures != None:
            self.FeaturesAll = list(self.SelectedFeatures)
        #self.FeaturesAll = [0, 1, 2, 3, 4, 5] #%%%
        print "Permitted features:", self.FeaturesAll
        # Parse the features from the TRAINING and TESTING files.  We generate
        # training sets for the FACULTATIVE (F) and for CONSTITUTIVE (C) sites.
        self.TrainingSet2 = Learning.FeatureSetClass()
        self.TrainingSet2.Type = "All"
        self.TrainingSet3 = Learning.FeatureSetClass()
        self.TrainingSet3.Type = "All"
        self.ParseFeatureFile(self.InputFeaturePath, self.TrainingSet2, self.TrainingSet3, self.TrainingSetDBRatio)
        if self.ModelTestFilePath:
            self.TestingSet2 = FeatureSetClass()
            self.TestingSet3 = FeatureSetClass()
            self.ParseFeatureFile(self.ModelTestFilePath, self.TestingSet2, self.TestingSet3, self.TestingSetDBRatio)        
    def Train(self):
        """
        Our training data-set is in self.InputFeaturePath.
        Let's train a model to predict which entries come from the true database.
        """
        self.ParseTrainingAndTestingSets()
        # Call TrainMachineLearner:
        (Stub, Extension) = os.path.splitext(self.WriteModelFilePath)
        WriteModelFilePath2 = "%s.2%s"%(Stub, Extension)
        WriteModelFilePath3 = "%s.3%s"%(Stub, Extension)
        self.TrainModel(0, WriteModelFilePath2)
        self.TrainModel(1, WriteModelFilePath3)
        if self.ModelType == "feature":
            return
        # Load our model objects:
        self.Model2 = Learning.LoadGeneralModel(WriteModelFilePath2)
        self.Model3 = Learning.LoadGeneralModel(WriteModelFilePath3)
        # Iterate over our vectors, and write the model score for each one.
        OutputFile = open(self.OutputFeaturePath, "wb")
        #self.OutputScoredVectors(OutputFile, self.TrainingSet2, self.Model2)
        #self.OutputScoredVectors(OutputFile, self.TrainingSet3, self.Model3)
        self.OutputScoredVectors(self.InputFeaturePath, OutputFile)
    def OutputScoredVectors(self, FeaturePath, OutputFile):
        File = open(FeaturePath, "rb")
        for FileLine in File.xreadlines():
            if FileLine[0] == "#":
                continue
            Bits = FileLine.replace("\r","").replace("\n","").split("\t")
            # If there are TOO MANY bits, then discard the extras:
            Bits = Bits[:FormatBits.LastFeature + 1]
            try:
                TrueFlag = int(Bits[FormatBits.TrueProteinFlag])
                Charge = int(Bits[FormatBits.Charge])
            except:
                continue # skip; not a valid instance line
            Features = []
            try:
                for FeatureIndex in range(FormatBits.FirstFeature, FormatBits.LastFeature + 1):
                    FeatureNumber = FeatureIndex - FormatBits.FirstFeature
                    if FeatureNumber not in self.FeaturesAll:
                        continue # skip features we're not allowed to use
                    if FeatureIndex < len(Bits) and Bits[FeatureIndex].strip() and Bits[FeatureIndex] != "None":
                        Features.append(float(Bits[FeatureIndex]))
                    else:
                        Features.append(0)
            except:
                traceback.print_exc()
                print "** Error on line %s column %s of feature file"%(LineNumber, FeatureIndex)
                print Bits
                continue
            if Charge > 2:
                Score = self.Model3.ScoreInstance(Features)
                PValue = self.Model3.GetPValue(Score)
            else:
                Score = self.Model2.ScoreInstance(Features)
                PValue = self.Model2.GetPValue(Score)
            String = FileLine.strip() + "\t%s\t%s\t"%(Score, PValue)
            OutputFile.write(String + "\n")
    def xxOutputScoredVectors(self, OutputFile, TrainingSet, Model):
        for Vector in TrainingSet.AllVectors:
            String = string.join(Vector.FileBits[:FormatBits.FirstFeature], "\t")
            String += "\t"
            for Feature in Vector.Features:
                String += "%s\t"%Feature
            #String += string.join(Vector.Features, "\t")
            ModelScore = Model.ScoreInstance(Vector.Features)
            PValue = Model.GetPValue(ModelScore)
            String += "%s\t%s\t"%(ModelScore, PValue)
            OutputFile.write(String + "\n")
    def TrainModel(self, Charge3Flag, WriteModelPath):
        print "Train model  (%s, %s)"%(Charge3Flag, WriteModelPath)
        FeatureFileName = "PeptideFeatures.%s.txt"%Charge3Flag
        TempFeatureFile = open(FeatureFileName, "wb")
        TempFeatureFile.write("#Index\tTrueFlag\t")        
        Keys = self.FeatureNames.keys()
        Keys.sort()
        for Key in Keys:
            TempFeatureFile.write(self.FeatureNames[Key] + "\t")
        TempFeatureFile.write("\n")
        if Charge3Flag:
            self.TrainingSet3.SaveTabDelimited(TempFeatureFile)
        else:
            self.TrainingSet2.SaveTabDelimited(TempFeatureFile)
        TempFeatureFile.close()
        Command = "TrainMachineLearner.py -r %s -m %s -w %s"%(FeatureFileName, self.ModelType, WriteModelPath)
        if self.ReportROCPath:
            (Stub, Extension) = os.path.splitext(self.ReportROCPath)
            Command += " -R %s.%s%s"%(Stub, Charge3Flag, Extension)
        print Command
        os.system(Command)
    def ParseCommandLine(self):
        (Options, Args) = getopt.getopt(sys.argv[1:], "m:u:v:r:w:f:e:R:D:a:F:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-m":
                self.ModelType = Value
            elif Option == "-D":
                self.TrainingSetDBRatio = float(Value)
            elif Option == "-F":
                self.FilterFactor = float(Value)
            elif Option == "-r":
                if not os.path.exists(Value):
                    print "** Error: Model file '%s' not found for reading.\n"%Value
                    return 0
                self.ReadModelFilePath = Value
            elif Option == "-w":
                self.WriteModelFilePath = Value
            elif Option == "-u":
                if not os.path.exists(Value):
                    print "** Error: Feature file '%s' not found for reading.\n"%Value
                    return 0
                self.InputFeaturePath = Value
            elif Option == "-v":
                self.OutputFeaturePath = Value
            elif Option == "-e":
                self.ModelTestFilePath = Value
            elif Option == "-f":
                self.FeatureSelectionFlag = int(Value)
            elif Option == "-R":
                self.ReportROCPath = Value
            elif Option == "-a":
                self.SelectedFeatures = self.ParseCommaDelimitedInts(Value)
            else:
                print "* Error: Unrecognized option %s"%Option
                return 0
        return 1 # success
    def ParseCommaDelimitedInts(self, String):
        ValueList = []
        for Bit in String.split(","):
            FeatureIndex = int(Bit)
            ValueList.append(FeatureIndex)
        return ValueList        
    def GetModelObject(self, Features):
        if self.ModelType == "lda":
            return Learning.LDAModel(Features)
        elif self.ModelType == "svm":
            return Learning.SVMModel(Features)
        elif self.ModelType == "logit":
            return Learning.LogitModel(Features)
        else:
            print "** Model type NOT KNOWN!", self.ModelType
            return
    def ReportAccuracy(self, SortedList, ROCCurvePlotPath = None):
        """
        The list should have entries of the form (ModelScore, TrueFlag)
        We'll sort them from high model scores to low, and report how many
        TRUE positives we have for a given FALSE DISCOVERY RATE.
        """
        SortedList.sort()
        SortedList.reverse()
        AllTrueCount = 0
        for Tuple in SortedList:
            AllTrueCount += Tuple[-1]
        AllFalseCount = len(SortedList) - AllTrueCount
        print "SortedList has %s entries, %s true"%(len(SortedList), AllTrueCount)
        # Iterate through the list from best to worst.  Report the number of hits
        # before false positive rate rises above 1%, and before it rises above 5%.
        # ALSO: Compute the area under the ROC curve!
        TrueCount = 0
        FalseCount = 0
        Cutoffs = (0.01, 0.03, 0.05, 0.07, 0.1)
        HitFlags = [0] * len(Cutoffs)
        Thresholds = [0] * len(Cutoffs)
        BestCounts = [0] * len(Cutoffs)
        BestCountsGenerous = [0] * len(Cutoffs)
        PrevStuff = None
        TopCount = 0
        TopCountFalse = 0
        if ROCCurvePlotPath:
            ROCCurvePlotFile = open(ROCCurvePlotPath, "wb")
        # Count true and false instances:
        OverallTrueCount = 0
        OverallFalseCount = 0
        for Tuple in SortedList:
            if Tuple[-1]:
                OverallTrueCount += 1
            else:
                OverallFalseCount += 1
        ROCTPForFP = {}
        ROCTPForFPCount = {}
        # Find the cutoff that gives a particular DISCOVERY RATE:
        for Index in range(len(SortedList)):
            Tuple = SortedList[Index]
            if Tuple[-1]:
                TrueCount += 1
            else:
                FalseCount += 1
            if (TrueCount + FalseCount) <= 200:
                TopCount = (TrueCount + FalseCount)
                TopCountFalse = FalseCount
            OverallTPRate = TrueCount / float(max(1, OverallTrueCount))
            OverallFPRate = FalseCount / float(max(1, OverallFalseCount))
            Bin = int(round(OverallFPRate * 100))
            ROCTPForFP[Bin] = ROCTPForFP.get(Bin, 0) + OverallTPRate
            ROCTPForFPCount[Bin] = ROCTPForFPCount.get(Bin, 0) + 1
            if ROCCurvePlotPath:
                ROCCurvePlotFile.write("%s\t%s\t%s\t%s\t%s\t\n"%(Index, TrueCount, FalseCount, OverallFPRate, OverallTPRate))
            #print Index, Tuple[0], TrueCount, FalseCount, OverallTrueCount, OverallFalseCount, OverallTPRate, OverallFPRate
            if Tuple[0] == PrevStuff:
                if TopCount == (TrueCount + FalseCount - 1):
                    TopCount = (TrueCount + FalseCount)
                    TopCountFalse = FalseCount
                continue
            PrevStuff = Tuple[0]
            FDRate = FalseCount / float(max(1, TrueCount))
            FDRate = min(1.0, FDRate)            
            for CutIndex in range(len(Cutoffs)):
                if FDRate > Cutoffs[CutIndex]:
                    HitFlags[CutIndex] = 1
                if not HitFlags[CutIndex]:
                    BestCounts[CutIndex] = max(BestCounts[CutIndex], TrueCount)
                    Thresholds[CutIndex] = Tuple[0]
                if FDRate <= Cutoffs[CutIndex]:
                    BestCountsGenerous[CutIndex] = max(BestCountsGenerous[CutIndex], TrueCount)
        # Compute the area under the ROC curve.
        for Bin in range(0, 100):
            if ROCTPForFP.has_key(Bin):
                ROCTPForFP[Bin] /= float(ROCTPForFPCount[Bin])
        ROCArea = 0
        for Bin in range(0, 100):
            if ROCTPForFP.has_key(Bin):
                ROCArea += 0.01 * ROCTPForFP[Bin]
                #print "%s: %s"%(Bin, ROCTPForFP[Bin])
            else:
                # Interpolate between points:
                PrevX = 0 # default
                PrevY = 0 # default
                for PrevBin in range(Bin - 1, -1, -1):
                    if ROCTPForFP.has_key(PrevBin):
                        PrevX = PrevBin
                        PrevY = ROCTPForFP[PrevBin]
                        break
                NextX = 100
                NextY = 1
                for NextBin in range(Bin + 1, 101):
                    if ROCTPForFP.has_key(NextBin):
                        NextX = NextBin
                        NextY = ROCTPForFP[NextBin]
                        break
                InterpolatedValue = PrevY + (Bin - PrevX) * float(NextY - PrevY) / (NextX - PrevX)
                ROCArea += 0.01 * InterpolatedValue
##                print "%s: %s (interpolate (%s, %s) to (%s, %s)"%(Bin, InterpolatedValue,
##                                                                  PrevX, PrevY, NextX, NextY)
        for CutIndex in range(len(Cutoffs)):
            Sensitivity = 100 * BestCounts[CutIndex] / float(max(1, AllTrueCount))
            print "  At %.1f%% FDRate (cutoff %.5f), got %s PTMs (sensitivity %.2f%%)"%(Cutoffs[CutIndex] * 100, Thresholds[CutIndex],
                BestCounts[CutIndex], Sensitivity)
            print "  ->True sensitivity: %.4f%%"%(100 * BestCounts[CutIndex] / float(max(1, AllTrueCount - AllFalseCount)))
##        print "Generous-count version:"
##        for CutIndex in range(len(Cutoffs)):
##            print "  At %.1f%% FDRate, got %s PTMs (sensitivity %.2f%%)"%(Cutoffs[CutIndex] * 100, BestCounts[CutIndex],
##                100 * BestCounts[CutIndex] / float(max(1, AllTrueCount)))
        print "False positive rate amoung top %s sites: %s"%(TopCount, 100*TopCountFalse/float(max(1, TopCount)))
        print "Overall, %s true and %s false features."%(TrueCount, FalseCount)
        print "ROC curve area: %.5f"%ROCArea
        # The 'score' we return is a tuple giving the best accuracy at several cutoffs:
        return (BestCounts[2], BestCounts[0], BestCounts[4], BestCounts[3], BestCounts[2])
    def xxxWriteScoredFeatureSet(self, FeatureSet):
        # Write out the features with their model-scores:
        if not self.OutputFeaturePath:
            return
        File = open(self.OutputFeaturePath, "wb")
        for FileLine in self.HeaderLines:
            File.write(FileLine)
        # Iterate over all vectors, write them all out:
        for Vector in FeatureSet.AllVectors:
            Bits = Vector.FileBits
            while len(Bits) <= FormatBits.ModelPValue:
                Bits.append("")
            Bits[FormatBits.ModelScore] = str(Vector.Score)
            Bits[FormatBits.ModelPValue] = str(self.Model.GetPValue(Vector.Score))
            Str = string.join(Bits, "\t")
            File.write(Str + "\n")
        File.close()


if __name__ == "__main__":
    Modeler = PeptideModeler()
    Modeler.ParseCommandLine()
    if Modeler.ModelType:
        Modeler.Train()
    else:
        print UsageInfo