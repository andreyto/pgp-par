"""
General-purpose machine learner driver.
"""
import os
import sys
import random
import traceback
import getopt
import Learning
import BasicStats

class FileColumns:
    Comment = 0
    TrueFlag = 1
    FirstFeature = 2

UsageInfo = """
TrainMachineLearner:
 -r [FileName]: Train on these features
 -m [Type]: Model type.  Options are 'feature', 'lda', 'logit', 'svm'
 -w [FileName]: Save model to file

Optional:
 -f [Flag]: Feature selection type: 1 for additive, 2 for subtractive
 -a [FeatureIndices]: Use this list of features.  For additive selection, begin with
    this list of features.
 -x [FeatureIndices]: For feature indices, forbid this list of features.
 -v: Verbose training
"""

class MachineTrainer:
    def __init__(self):
        self.FeatureNames = []
        self.ForbiddenFeatureList = []
        self.ModelType = None
        self.InputFileName = None
        self.Columns = FileColumns
        self.FeatureSelectionFlag = 0
        self.ReadModelFilePath = None
        self.WriteModelFilePath = None
        self.ModelTestFilePath = None
        self.SpecifyFeatureList = None
        self.TrainVerboseFlag = 0
        self.ReportROCPath = None
        self.MaxFeatureSetSize = 100000
    def ParseCommaDelimitedInts(self, String):
        ValueList = []
        for Bit in String.split(","):
            FeatureIndex = int(Bit)
            ValueList.append(FeatureIndex)
        return ValueList        
    def ParseCommandLine(self, Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:m:f:w:a:x:vR:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-m":
                self.ModelType = Value.lower()
            elif Option == "-r":
                self.InputFileName = Value
            elif Option == "-f":
                self.FeatureSelectionFlag = int(Value)
            elif Option == "-w":
                self.WriteModelFilePath = Value
            elif Option == "-a":
                self.SpecifyFeatureList = self.ParseCommaDelimitedInts(Value)
            elif Option == "-x":
                self.ForbiddenFeatureList = self.ParseCommaDelimitedInts(Value)
            elif Option == "-v":
                self.TrainVerboseFlag = 1
            elif Option == "-R":
                self.ReportROCPath = Value
                
            else:
                print "???", Option, Value
                raise ValueError
        if not self.InputFileName and not self.ModelType:
            print UsageInfo
            sys.exit(-1)
    def ParseFeaturesFromFile(self, FileName):
        """
        Parse one Learning.FeatureSet from the specified file.
        Format is tab-delimited, with format:
        Comment\tTrueFlag\tFeature0\tFeature1\t...
        """
        FeatureSet = Learning.FeatureSetClass()
        File = open(FileName, "rb")
        LineNumber = 0
        HeaderLineSeen = 0
        self.FeatureNames = []
        for FileLine in File.xreadlines():
            LineNumber += 1
            Bits = FileLine.strip().split("\t")
            if len(Bits) < 2:
                continue
            if FileLine[0] == "#":
                if not HeaderLineSeen:
                    # Special: Assume the first line has column headings for the features.
                    self.FeatureNames = []
                    if self.SpecifyFeatureList and not self.FeatureSelectionFlag:
                        for FeatureIndex in self.SpecifyFeatureList:
                            BitIndex = self.Columns.FirstFeature + FeatureIndex
                            if BitIndex < len(Bits):
                                self.FeatureNames.append(Bits[BitIndex])
                            else:
                                self.FeatureNames.append("FeatureIndex%s"%FeatureIndex)
                        
                    else:
                        for BitIndex in range(self.Columns.FirstFeature, len(Bits)):
                            self.FeatureNames.append(Bits[BitIndex])
                    HeaderLineSeen = 1
                continue
            try:
                Vector = Learning.FeatureVector()
                Vector.FileBits = Bits
                Vector.TrueFlag = int(Bits[self.Columns.TrueFlag])
                if self.SpecifyFeatureList and not self.FeatureSelectionFlag:
                    for FeatureIndex in self.SpecifyFeatureList:
                        BitIndex = self.Columns.FirstFeature + FeatureIndex
                        Vector.Features.append(float(Bits[BitIndex]))
                        
                else:
                    for BitIndex in range(self.Columns.FirstFeature, len(Bits)):
                        Vector.Features.append(float(Bits[BitIndex]))
            except:
                traceback.print_exc()
                print "ERROR on line %s"%LineNumber
                print Bits
                continue
            # Skip truncated lines (e.g. if the feature-file is still being
            # written out)
            if HeaderLineSeen and len(Vector.Features) < len(self.FeatureNames):
                continue
            if Vector.TrueFlag:
                FeatureSet.TrueVectors.append(Vector)
            else:
                FeatureSet.FalseVectors.append(Vector)
            FeatureSet.AllVectors.append(Vector)
            if len(FeatureSet.AllVectors) >= self.MaxFeatureSetSize:
                break
        FeatureSet.SetCounts()
        self.Features = range(FeatureSet.Size)
        while len(self.FeatureNames) < FeatureSet.Size:
            self.FeatureNames.append("Feature%s"%len(self.FeatureNames))
        print "Feature count:", FeatureSet.Size        
        return FeatureSet
    def Main(self):
        self.TrainingSet = self.ParseFeaturesFromFile(self.InputFileName)
        if self.ModelType == "feature":
            self.TrainOneFeature(self.TrainingSet)
            return
        self.TestingSet = self.TrainingSet # default
        self.Model = self.GetModelObject(self.ModelType, self.Features)
        if self.FeatureSelectionFlag == 0:
            # Load a pre-trained model, if we received a path:
            if self.ReadModelFilePath:
                self.Model.LoadModel(self.ReadModelFilePath)
            else:
                self.Model.Train(self.TrainingSet, self.TrainVerboseFlag)
                self.Model.Test(self.TrainingSet)
                self.Model.ReportAccuracy(self.TrainingSet, self.ReportROCPath)
                if self.WriteModelFilePath:
                    self.Model.SaveModel(self.WriteModelFilePath)
                    if self.ModelType == "lda":
                        self.Model.SaveBinaryModel(self.WriteModelFilePath + ".model")
                    if self.ModelType == "svm":
                        Stub = os.path.splitext(self.WriteModelFilePath)[0]
                        self.Model.SaveTextModel(Stub)
        #######################################################################
        # Accumulation feature selection: Start with no features, and iteratively
        # add the BEST marginal feature.
        if self.FeatureSelectionFlag == 1:
            self.TrainingResultsFile = open("FeatureSelectionAdd.txt", "wb")
            self.TrainingResultsFile.write("Feature selection for %s, accumulating features:\n"%self.ModelType)
            self.PerformFeatureSelectionAdditive(self.Features)
            self.TrainingResultsFile.close()
        #######################################################################
        # Subtractive feature selection: Iteratively REMOVE features
        if self.FeatureSelectionFlag == 2:
            self.TrainingResultsFile = open("FeatureSelectionRemove.txt", "wb")
            self.TrainingResultsFile.write("Feature selection for %s, pruning features:\n"%self.ModelType)
            self.PerformFeatureSelectionSubtractive(self.Features)
            self.TrainingResultsFile.close()
    def PerformFeatureSelectionAdditive(self, Features):
        FeatureList = []
        OrderedFeatureList = []
        if self.SpecifyFeatureList:
            FeatureList = self.SpecifyFeatureList[:]
        while (1):
            BestNewFeature = None
            BestNewScore = None
            for NewFeature in Features:
                # consider adding NewFeature to FeatureList:
                if NewFeature in FeatureList:
                    continue
                if NewFeature in self.ForbiddenFeatureList:
                    continue
                TempFeatureList = FeatureList[:]
                TempFeatureList.append(NewFeature)
                TempFeatureList.sort()
                print
                print ">>>>>>>>>>>>>>>>>>>>"
                print "GetModel(%s)"%TempFeatureList
                self.Model = self.GetModelObject(self.ModelType, TempFeatureList)
                # Train on these features, and report accuracy on the testing set:
                try:
                    self.Model.Train(self.TrainingSet)
                    if self.ModelTestFilePath:
                        self.Model.Test(self.TestingSet)
                        Score = self.Model.ReportAccuracy(self.TestingSet)
                    else:
                        self.Model.Test(self.TrainingSet)
                        Score = self.Model.ReportAccuracy(self.TrainingSet)
                except:
                    traceback.print_exc()
                    print "** Erorr: Unable to train model on feature list:", TempFeatureList
                    Score = -9999
                #print TempFeatureList
                #sys.stdin.readline()
                self.TrainingResultsFile.write("  %s -> %s\n"%(TempFeatureList, Score))
                self.TrainingResultsFile.flush()
                if Score > BestNewScore:
                    BestNewScore = Score
                    BestNewFeature = NewFeature
            if BestNewScore == None:
                break # no features left to add!
            FeatureList.append(BestNewFeature)
            OrderedFeatureList.append(BestNewFeature)
            FeatureList.sort()
            self.TrainingResultsFile.write("\nBest new feature: %s\n"%BestNewFeature)
            self.TrainingResultsFile.write("Feature list is now: %s\n"%FeatureList)
            self.TrainingResultsFile.write("Ordered list: %s\n"%OrderedFeatureList)
            self.TrainingResultsFile.write("Accuracy: %s\n\n"%str(BestNewScore))
            self.TrainingResultsFile.flush()
    def PerformFeatureSelectionSubtractive(self, Features):
        FeatureList = []
        for Feature in Features:
            if Feature not in self.ForbiddenFeatureList:
                FeatureList.append(Feature)
        self.Model = self.GetModelObject(self.ModelType, FeatureList)
        try:
            self.Model.Train(self.TrainingSet)
            if self.ModelTestFilePath:
                self.Model.Test(self.TestingSet)
                Score = self.Model.ReportAccuracy(self.TestingSet)
            else:
                self.Model.Test(self.TrainingSet)
                Score = self.Model.ReportAccuracy(self.TrainingSet)
        except:
            print "** Erorr: Unable to train model on feature list:", FeatureList
            Score = -9999
        self.TrainingResultsFile.write("Initial accuracy on all features: %s\n"%str(Score))
        while (1):
            if len(FeatureList) <= 1:
                break # we've pruned our features down to the stump
            BestNewFeature = None
            BestNewScore = None
            for NewFeature in FeatureList[:]:
                # Consider removing NewFeature from FeatureList:
                TempFeatureList = FeatureList[:]
                TempFeatureList.remove(NewFeature)
                TempFeatureList.sort()
                # Train on these features, and report accuracy on the testing set:
                try:
                    self.Model.Train(self.TrainingSet)
                    if self.ModelTestFilePath:
                        self.Model.Test(self.TestingSet)
                        Score = self.Model.ReportAccuracy(self.TestingSet)
                    else:
                        self.Model.Test(self.TrainingSet)
                        Score = self.Model.ReportAccuracy(self.TrainingSet)
                except:
                    traceback.print_exc()
                    print "** Erorr: Unable to train model on feature list:", TempFeatureList
                    Score = -9999
                if Score > BestNewScore:
                    BestNewScore = Score
                    BestNewFeature = NewFeature
            if BestNewScore == None:
                break # no features left to remove!
            FeatureList.remove(BestNewFeature)
            FeatureList.sort()
            self.TrainingResultsFile.write("\nRemoved feature: %s\n"%BestNewFeature)
            self.TrainingResultsFile.write("Feature list is now: %s\n"%FeatureList)
            self.TrainingResultsFile.write("Accuracy: %s\n"%str(BestNewScore))
            self.TrainingResultsFile.flush()
    def GetModelObject(self, ModelType, Features):
        if ModelType == "lda":
            return Learning.LDAModel(Features)
        elif ModelType == "svm":
            return Learning.SVMModel(Features)
        elif ModelType == "logit":
            return Learning.LogitModel(Features)
        else:
            print "** Model type NOT KNOWN!", self.ModelType
            sys.exit(-1)
    def TrainOneFeature(self, FeatureSet):
        """
        Compute accuracy for a very simple-minded model:
        Rank sites by the value of a SINGLE FEATURE (descending order)
        """
        File = open("SingleFeatureResults.txt", "wb")
        Header = "FeatureIndex\tFeatureName\tMeanTrue\tMeanFalse\tFCS\tROCArea\tROCValue\t\n"
        File.write(Header)
        for FeatureIndex in range(FeatureSet.Size):
            SortedList = []
            MaxValue = -9999
            MinValue = 9999
            TrueValues = []
            FalseValues = []
            for Vector in FeatureSet.TrueVectors:
                if FeatureIndex < len(Vector.Features):
                    Value = Vector.Features[FeatureIndex]
                else:
                    # Padding value:
                    Value = 0
                MaxValue = max(MaxValue, Value)
                MinValue = min(MinValue, Value)
                SortedList.append((Value, random.random(), 1))
                TrueValues.append(Value)
            for Vector in FeatureSet.FalseVectors:
                if FeatureIndex < len(Vector.Features):
                    Value = Vector.Features[FeatureIndex]
                else:
                    # Padding value:
                    Value = 0
                MaxValue = max(MaxValue, Value)
                MinValue = min(MinValue, Value)
                SortedList.append((Value, random.random(), 0))
                FalseValues.append(Value)
            # If the feature doesn't vary at all, then don't bother
            # trying to compute accuracy:
            if (MaxValue == MinValue):
                return
            (AverageTrue, StdDevTrue) = BasicStats.GetMeanStdDev(TrueValues)
            (AverageFalse, StdDevFalse) = BasicStats.GetMeanStdDev(FalseValues)
            # And report the accuracy of this lonely feature:
            print
            print "Feature %s (%s):"%(FeatureIndex, self.FeatureNames[FeatureIndex])
            print "Average on true  vectors: %.6f"%(AverageTrue)
            print "Average on false vectors: %.6f"%(AverageFalse)
            FCS = (AverageTrue - AverageFalse)**2 / (StdDevTrue + StdDevFalse)
            print "FCS: %.6f"%FCS
            ROCArea = self.ReportAccuracy(SortedList)        
            Str = "%s\t%s\t"%(FeatureIndex, self.FeatureNames[FeatureIndex])
            Str += "%s\t%s\t%s\t"%(AverageTrue, AverageFalse, FCS)
            ROC = 2 * abs(0.5 - ROCArea)
            Str += "%s\t%s\t"%(ROCArea, ROC)
            File.write(Str + "\n")
        # All features are complete:
        File.close()
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
            OverallTPRate = TrueCount / float(max(1, AllTrueCount))
            OverallFPRate = FalseCount / float(max(1, AllFalseCount))
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
        for CutIndex in range(len(Cutoffs)):
            Sensitivity = 100 * BestCounts[CutIndex] / float(max(1, AllTrueCount))
            print "  At %.1f%% FDRate (cutoff %.5f), got %s true (sensitivity %.2f%%)"%(Cutoffs[CutIndex] * 100, Thresholds[CutIndex],
                BestCounts[CutIndex], Sensitivity)
            print "  ->True sensitivity: %.4f%%"%(100 * BestCounts[CutIndex] / float(max(1, AllTrueCount - AllFalseCount)))
        print "False positive rate amoung top %s sites: %s"%(TopCount, 100*TopCountFalse/float(max(1, TopCount)))
        print "Overall, %s true and %s false features."%(TrueCount, FalseCount)
        print "ROC curve area: %.5f"%ROCArea
        return ROCArea

if __name__ == "__main__":
    Trainer = MachineTrainer()
    Trainer.ParseCommandLine(sys.argv[1:])
    Trainer.Main()