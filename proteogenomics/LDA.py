"""
Linear discriminant analysis
Assumes the input file is tab-delimited, with category in the first column, and float
values in the remaining columns.
"""

USE_NUMPY = 1
import traceback
try:
    if USE_NUMPY:
        from numpy import *
        import numpy.linalg
        MatrixMulitply = dot
        InvertMatrix = numpy.linalg.inv
        FloatType = float
    else:
        from Numeric import *
        import LinearAlgebra
        InvertMatrix = LinearAlgebra.inverse
        MatrixMulitply = matrixmultiply
        FloatType = Float
except:
    print "\n* Warning: Unable to import numpy.  LDA training not available."
    print "  Please install NumPy (see http://numpy.scipy.org/ for details)"
    print "  Error details are shown here:"
    traceback.print_exc()
    
import math
import os
import sys
import random
import struct
import traceback

ForbiddenFeatures = [2, 3, 4, 5, 13, 31, 32, 33, 34, 43, 47] #[2,3,4]

def PrintHistogram(Histogram, HistoFile):
    Bins = Histogram.keys()
    Bins.sort()
    #Bins.reverse()
    TotalBads = 0
    TotalGoods = 0
    print "\nHistogram results:"
    for Bin in Bins:
        TotalBads += Histogram[Bin][0]
    Bads = TotalBads
    for Bin in range(Bins[0], Bins[-1]):
        if Histogram.has_key(Bin):
            Bads -= Histogram[Bin][0]
        PValue = Bads / float(TotalBads)
        print "%s\t%s\t%s\t%s\t"%(Bin, PValue, Bads, TotalBads)
        if HistoFile and Bin >=- 70 and Bin < 150:
            PValue = min(0.99, max(0.0001, PValue))
            Str = struct.pack("<f", PValue)
            #print "PValue struct:", Str
            HistoFile.write(Str)
            
class LDAClassifier:
    def __init__(self):
        pass
    def GetCovarianceArray(self, VectorList):
        Size = len(VectorList[0])
        VectorCount = float(len(VectorList))
        C = zeros((Size, Size), FloatType)
        for Vector in VectorList:
            for X in range(Size):
                for Y in range(Size):
                    C[X][Y] += Vector[X] * Vector[Y] / VectorCount
        return C
    def LoadVectors(self, FileName, CategoryBit, FeatureList):
        Size = len(FeatureList)
        self.GoodVectors = []
        self.BadVectors = []
        # Iterate over file lines, and read vetors in:
        File = open(FileName, "r")
        for FileLine in File.xreadlines():
            if FileLine[0] == "#":
                continue # comment
            Bits = FileLine.split("\t")
            try:
                Category = int(Bits[CategoryBit])
            except:
                continue
            # Turn -1 vs 1 into 0 vs 1:
            if Category < 0:
                Category = 0
            #for X in range(len(Bits)):
            #    print "%s: %s"%(X, Bits[X])
            Vector = []
            try:
                for Index in FeatureList:
                    if Index >= len(Bits) or not Bits[Index].strip():
                        Vector.append(0)
                    else:
                        Vector.append(float(Bits[Index]))
            except:
                traceback.print_exc()
                print Bits
                continue
            if Category:
                self.GoodVectors.append(Vector)
            else:
                self.BadVectors.append(Vector)
        print "First GoodVector:\n", self.GoodVectors[0]
        print "First BadVector:\n", self.BadVectors[0]
    def ScaleVectors(self):
        """
        Scale all vectors so that 90% of all values lie in the range [-1, 1]
        """
        Values = []
        MinValues = []
        MaxValues = []
        FeatureCount = len(self.GoodVectors[0])
        for X in range(FeatureCount):
            Values.append([])
        for Vector in self.GoodVectors:
            for X in range(FeatureCount):
                Values[X].append(Vector[X])
        for Vector in self.BadVectors:
            for X in range(FeatureCount):
                Values[X].append(Vector[X])
        print "Value count:", len(Values[0])
        for X in range(FeatureCount):
            Values[X].sort()
            ValueCount = len(Values[X])
            MinValues.append(Values[X][int(round(ValueCount * 0.05))])
            MaxValues.append(Values[X][int(round(ValueCount * 0.95))])
        print "Range:"
        for X in range(FeatureCount):
            print "%s: %.4f ... %.4f"%(X, MinValues[X], MaxValues[X])
        for X in range(FeatureCount):
            HalfRange = (MaxValues[X] - MinValues[X]) / 2.0
            if not HalfRange:
                continue
            for Vector in self.BadVectors:            
                Vector[X] = (Vector[X] - MinValues[X]) / HalfRange - 1.0
                #Vector[X] = max(-1.0, min(Vector[X], 1.0))
            for Vector in self.GoodVectors:            
                Vector[X] = (Vector[X] - MinValues[X]) / HalfRange - 1.0
                #Vector[X] = max(-1.0, min(Vector[X], 1.0))
    def PerformLDA(self, FileName, CategoryBit, FeatureList, ScaleVectors = 1, FoldValidation = 0):
        VerboseFlag = 1
        Size = len(FeatureList)
        self.LoadVectors(FileName, CategoryBit, FeatureList)
        if ScaleVectors:
            self.ScaleVectors()
        if FoldValidation:
            random.seed(1)
            random.shuffle(self.GoodVectors)
            random.shuffle(self.BadVectors)
        # n-fold validation:
        self.MasterGoodVectors = self.GoodVectors
        self.MasterBadVectors = self.BadVectors
        WorstAccuracy = 1.0
        for Fold in range(max(1, FoldValidation)):
            # Slice the master lists of good and bad vectors to separate 1/FoldValidation of
            # them into a test set.  FoldValidation can be 0, to do no such splitting.
            self.GoodVectors = []
            self.GoodTestVectors = []
            for X in range(len(self.MasterGoodVectors)):
                if FoldValidation and X % FoldValidation == Fold:
                    self.GoodTestVectors.append(self.MasterGoodVectors[X])
                else:
                    self.GoodVectors.append(self.MasterGoodVectors[X])
            self.BadVectors = []
            self.BadTestVectors = []
            for X in range(len(self.MasterBadVectors)):
                if FoldValidation and X % FoldValidation == Fold:
                    self.BadTestVectors.append(self.MasterBadVectors[X])
                else:
                    self.BadVectors.append(self.MasterBadVectors[X])
            ############################################################
            # Compute the mean vectors:
            GoodCount = float(len(self.GoodVectors))
            BadCount = float(len(self.BadVectors))
            AllCount = GoodCount + BadCount
            self.MeanGood = [0]*Size
            self.MeanBad = [0]*Size
            self.MeanGlobal = [0]*Size
            for Vector in self.GoodVectors:
                for Index in range(Size):
                    self.MeanGood[Index] += Vector[Index] / GoodCount
                    self.MeanGlobal[Index] += Vector[Index] / AllCount
            for Vector in self.BadVectors:
                for Index in range(Size):
                    self.MeanBad[Index] += Vector[Index] / BadCount
                    self.MeanGlobal[Index] += Vector[Index] / AllCount
            print "MeanGood:\n  ", self.MeanGood
            print "MeanBad:\n  ", self.MeanBad
            print "MeanGlobal:\n  ", self.MeanGlobal
            ############################################################
            # Compute the mean-corrected vectors:
            MeanCorrectedGoodVectors = []
            MeanCorrectedBadVectors = []
            for Vector in self.GoodVectors:
                NewVector = []
                for X in range(Size):
                    NewVector.append(Vector[X] - self.MeanGlobal[X])
                MeanCorrectedGoodVectors.append(NewVector)
            for Vector in self.BadVectors:
                NewVector = []
                for X in range(Size):
                    NewVector.append(Vector[X] - self.MeanGlobal[X])
                MeanCorrectedBadVectors.append(NewVector)
            ############################################################
            # Compute covariance matrices:
            CovarArrayGood = self.GetCovarianceArray(MeanCorrectedGoodVectors)
            if VerboseFlag:
                print "CovarArrayGood:", CovarArrayGood
            CovarArrayBad = self.GetCovarianceArray(MeanCorrectedBadVectors)
            if VerboseFlag:
                print "CovarArrayBad:", CovarArrayBad
            # CovarArrayFull is the pooled within-group covariance matrix, it's
            # computed componentwise as weighted sum of CovarArrayGood and CovarArrayBad.
            CovarArrayFull = zeros((Size, Size), FloatType)
            for X in range(Size):
                for Y in range(Size):
                    CovarArrayFull[X][Y] += CovarArrayGood[X][Y] * GoodCount / AllCount
                    CovarArrayFull[X][Y] += CovarArrayBad[X][Y] * BadCount / AllCount
            if VerboseFlag:
                print "CovarArrayFull:", CovarArrayFull
            ############################################################
            # Invert the covariance array:
            try:
                self.CI = InvertMatrix(CovarArrayFull)
            except:
                traceback.print_exc()
                print "Unable to invert covariance matrix!  Invalid feature set."
                return 0
            if VerboseFlag:
                print "CI:", self.CI
            self.GoodMuC = MatrixMulitply(self.CI, self.MeanGood)
            if VerboseFlag:
                print "GoodMuC:", self.GoodMuC
            self.BadMuC = MatrixMulitply(self.CI, self.MeanBad)
            if VerboseFlag:
                print "BadMuC:", self.BadMuC
            self.ConstantGood = -MatrixMulitply(self.MeanGood, self.GoodMuC) / 2.0
            self.ConstantBad = -MatrixMulitply(self.MeanBad, self.BadMuC) / 2.0
            if VerboseFlag:
                print "Constant good %.4f constant bad %.4f"%(self.ConstantGood, self.ConstantBad)
            #######################################################
            if VerboseFlag:
                # Print C initializers:
                for X in range(Size):
                    Str = "double CovInv%s[] = {"%chr(ord("A") + X)
                    for Y in range(Size):
                        Str += "%.3f,"%self.CI[X][Y]
                    Str = Str[:-1] + "};"
                    print Str
                Str = "double MeanVectorTrue[] = {"
                for X in range(Size):
                    Str += "%.3f,"%self.MeanGood[X]
                Str = Str[:-1] + "};"
                print Str
                Str = "double MeanVectorFalse[] = {"
                for X in range(Size):
                    Str += "%.3f,"%self.MeanBad[X]
                Str = Str[:-1] + "};"
                print Str
                print "double SubProdTrue = (float)%.3f;"%self.ConstantGood
                print "double SubProdFalse = (float)%.3f;"%self.ConstantBad
                print "CG and CB:", self.ConstantGood, self.ConstantBad
            #######################################################
            Weights = []
            for X in range(Size):
                Weights.append(self.GoodMuC[X] - self.BadMuC[X])
            Str = "*-*->Weights:"
            for X in range(Size):
                Str += " %s: %.4f"%(X, Weights[X])
            print Str
            #######################################################
            # Compute our accuracy on the testing set:
            if FoldValidation:
                CorrectCount = 0
                IncorrectCount = 0
                for Vector in self.GoodTestVectors:
                    NewVector = []
                    for X in range(Size):
                        NewVector.append(Vector[X] - self.MeanGlobal[X])
                    Reading = self.GetReading(NewVector)
                    if Reading > 0:
                        CorrectCount += 1
                    else:
                        IncorrectCount += 1
                for Vector in self.BadTestVectors:
                    NewVector = []
                    for X in range(Size):
                        NewVector.append(Vector[X] - self.MeanGlobal[X])
                    Reading = self.GetReading(NewVector)
                    if Reading > 0:
                        IncorrectCount += 1
                    else:
                        CorrectCount += 1
                TotalCount = CorrectCount + IncorrectCount
                Accuracy = CorrectCount / float(TotalCount)
                print "Cross-validation accuracy: %d of %d (%.3f%%)"%(CorrectCount, TotalCount, Accuracy*100)
                WorstAccuracy = min(Accuracy, WorstAccuracy)
            else:
                # Compute accuracy on all vectors:
                CorrectCount = 0
                IncorrectCount = 0
                for Vector in MeanCorrectedGoodVectors:
                    Reading = self.GetReading(NewVector)
                    if Reading > 0:
                        CorrectCount += 1
                    else:
                        IncorrectCount += 1
                for Vector in MeanCorrectedBadVectors:
                    Reading = self.GetReading(NewVector)
                    if Reading > 0:
                        IncorrectCount += 1
                    else:
                        CorrectCount += 1
                TotalCount = CorrectCount + IncorrectCount
                Accuracy = CorrectCount / float(TotalCount)
                print "Accuracy: %d of %d (%.3f%%)"%(CorrectCount, TotalCount, Accuracy*100)
                WorstAccuracy = min(Accuracy, WorstAccuracy)
        print "Min. cross-validation accuracy: %.3f%%"%(WorstAccuracy*100)
        return WorstAccuracy
    def GetReading(self, Vector):
        CIProduct = MatrixMulitply(self.CI, Vector)
        ReadingGood = MatrixMulitply(self.MeanGood, CIProduct) + self.ConstantGood
        ReadingBad = MatrixMulitply(self.MeanBad, CIProduct) + self.ConstantBad
        print
        print "Vector:", Vector
        print "CIProduct:", CIProduct
        print "ReadingGood %s ReadingBad %s Net %s"%(ReadingGood, ReadingBad, ReadingGood - ReadingBad)
        return (ReadingGood - ReadingBad)
    def ReportROCCurve(self):
        SortedList = []
        MeanCorrectedGoodVectors = []
        MeanCorrectedBadVectors = []
        PositiveCount = len(self.GoodVectors)
        NegativeCount = len(self.BadVectors)
        Size = len(self.GoodVectors[0])
        for Vector in self.GoodVectors:
            NewVector = []
            for X in range(Size):
                NewVector.append(Vector[X] - self.MeanGlobal[X])
            CIProduct = MatrixMulitply(self.CI, NewVector)
            ReadingGood = MatrixMulitply(self.MeanGood, CIProduct) + self.ConstantGood
            ReadingBad = MatrixMulitply(self.MeanBad, CIProduct) + self.ConstantBad
            SortedList.append((ReadingGood - ReadingBad, 1))
        for Vector in self.BadVectors:
            NewVector = []
            for X in range(Size):
                NewVector.append(Vector[X] - self.MeanGlobal[X])
            CIProduct = MatrixMulitply(self.CI, NewVector)
            ReadingGood = MatrixMulitply(self.MeanGood, CIProduct) + self.ConstantGood
            ReadingBad = MatrixMulitply(self.MeanBad, CIProduct) + self.ConstantBad
            SortedList.append((ReadingGood - ReadingBad, 0))
        SortedList.sort()
        SortedList.reverse()
        TPCount = 0
        FPCount = 0
        Area = 0
        ROCCurveFile = open("ROCCurve.txt", "wb")
        for (Reading, TrueFlag) in SortedList:
            #print Reading, TrueFlag
            if (TrueFlag):
                TPCount += 1
            else:
                FPCount += 1
                Area += (TPCount / float(PositiveCount))
            TPRate = TPCount / float(PositiveCount)
            FPRate = FPCount / float(NegativeCount)
            ROCCurveFile.write("%s\t%s\t\n"%(FPRate, TPRate))
            if FPRate < 0.05:
                HappyTPRate = TPRate
        Area /= float(FPCount)
        print "ROC curve area:", Area
        print "TP rate for FP < 0.05: %s"%HappyTPRate
        return Area
    def ProducePValueCurve(self):
        Histogram = {}
        HistogramShort = {}
        HistogramMedium = {}
        HistogramLong = {}
        MediumCutoff = 9
        LongCutoff = 13
        for X in range(len(self.BadVectors)):
            Vector = self.BadVectors[X]
            MQScore = self.GetReading(Vector)
            Bin = int(round(MQScore * 10))
            if not Histogram.has_key(Bin):
                Histogram[Bin] = [0,0]
            if not HistogramShort.has_key(Bin):
                HistogramShort[Bin] = [0,0]
            if not HistogramMedium.has_key(Bin):
                HistogramMedium[Bin] = [0,0]
            if not HistogramLong.has_key(Bin):
                HistogramLong[Bin] = [0,0]
            #Len = self.BadVectorPepLengths[X]
            Histogram[Bin][0] += 1
        PrintHistogram(Histogram, None)

def FeatureSelectMain():
    BestAccuracy = 0
    FeatureList = [6, 49]
    FeatureCount = 55
    while len(FeatureList)<12:
        BestAccuracy = 0
        BestList = FeatureList
        for FeatureA in range(4, FeatureCount):
            if FeatureA in FeatureList:
                continue
            if FeatureA in ForbiddenFeatures:
                continue
            AugmentedFeatureList = FeatureList[:]
            AugmentedFeatureList.append(FeatureA)
            LDA = LDAClassifier()
            LDA.PerformLDA("TrainingSet.Table.txt", 0, AugmentedFeatureList, 1, 0)
            Accuracy = LDA.ReportROCCurve()
            if Accuracy > BestAccuracy:
                BestAccuracy = Accuracy
                BestList = AugmentedFeatureList
            print "Feature set %s has accuracy %.4f%%"%(AugmentedFeatureList, Accuracy*100)
            print "So far...best accuracy %.4f%%, feature set %s"%(BestAccuracy*100, BestList)
        FeatureList = BestList
        print "Best accuracy %s, feature list %s"%(BestAccuracy, FeatureList)
            

def Main():
    LDA = LDAClassifier()
    LDA.PerformLDA("LDATrainingSet.txt", 0, [1, 5], 0, 0)
    LDA.ReportROCCurve()
    #LDA.ProducePValueCurve()
        
if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(psyco not loaded)"
    Main()
