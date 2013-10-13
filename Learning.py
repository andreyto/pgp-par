"""
The LearnerClass is an abstract machine learner.  It can be trained,
saved, and loaded.
"""
import os
import sys
import struct
import random
import math
import traceback
import cPickle
import LDA
import RunPySVM
try:
    import PySVM
except:
    print "(Warning: PySVM not present!)"

try:
    from numpy import *
    import numpy.linalg
    FloatType = float
    MatrixMultiply = dot
    InvertMatrix = numpy.linalg.inv
except:
    print "\n* Warning: Unable to import NumPy.  Logit training not available"
    print "  Please install NumPy (see http://numpy.scipy.org/ for details)"
    print "  Error details are shown here:"
    traceback.print_exc()
    
random.seed(1)

MaxSVMFeatureCount = 500

if sys.platform == "win32":
    PATH_SVMSCALE = r"C:\libsvm\windows\svmscale.exe"
    PATH_SVMTRAIN = r"C:\libsvm\windows\svmtrain.exe"
else:
    PATH_SVMSCALE = os.path.join(os.environ["HOME"], "libsvm", "svm-scale")
    PATH_SVMTRAIN = os.path.join(os.environ["HOME"], "libsvm", "svm-train")

SQRT2PI = math.sqrt(2 * math.pi)
SQRT2 = math.sqrt(2)
Cof = [76.18009172947146, -86.50532032941677,
    24.01409824083091, -1.231739572450155, 
    0.1208650973866179e-2, -0.5395239384952e-5]

def Gamma(Z):
    X = Z
    Y = Z
    Temp = X + 5.5
    Temp -= (X + 0.5) * math.log(Temp)
    Ser = 1.000000000190015
    for J in range(6):
        Y += 1
        Ser += Cof[J] / Y
    Z = -Temp + math.log(2.5066282746310005 * Ser / X)
    return math.exp(Z)




class MixtureModelClass:
    def __init__(self, BinMultiplier = 10.0):
        self.BinMultiplier = BinMultiplier
    def Model(self, Values, Histogram = None):
        if Values:
            print "Model scores.  Range is %s...%s"%(min(Values), max(Values))
        else:
            if not Histogram.keys():
                # There's nothing to model!
                self.MinBin = 0
                self.MaxBin = 0
                self.OddsTrue = {}
                return
            print "Model scores.  Range is %s...%s"%(min(Histogram.keys()), max(Histogram.keys()))
        self.MaxCycleCount = 300
        self.VerboseFlag = 0     
        if Histogram:
            self.ScoreHistogram = Histogram
        else:
            self.ScoreHistogram = {}
            for Value in Values:
                Bin = int(round(Value * self.BinMultiplier))
                self.ScoreHistogram[Bin] = self.ScoreHistogram.get(Bin, 0) + 1
        Keys = self.ScoreHistogram.keys()
        self.MinBin = min(Keys)
        self.MaxBin = max(Keys) + 1
        self.InitializeModel()
        try:
            self.ModelDistribution()
        except:
            print "* Warning: Unable to compute p-values via mixture model"
            print "* Error trace follows:"
            traceback.print_exc()
            print "self.VarianceFalse:", self.VarianceFalse
            print "self.VarianceTrue:", self.VarianceTrue
            print "MeanFalse:", self.MeanFalse
            print "GammaOffset:", self.GammaOffset
            print "ThetaFalse:", self.ThetaFalse
            print "KFalse:", self.KFalse
    def GetOddsTrue(self, X):
##        if self.CumulativeFlag:
##            # Set our odds using the cumulative probability p(score >= X) instead of
##            # the odds that the score is in this bin.
##            ErfArg = (X - self.MeanTrue) / (self.StdDevTrue * SQRT2)
##            NormalCDF = 0.5 + 0.5 * PyInspect.erf(ErfArg)
##            GX = max(0.01, X + self.GammaOffset)
##            GammaCDF = PyInspect.GammaIncomplete(self.KFalse, GX / self.ThetaFalse) #/ Gamma(self.KFalse)
##            TrueNormal = 1.0 - NormalCDF
##            FalseGamma = 1.0 - GammaCDF
##        else:
        Pow = - ((X - self.MeanTrue)**2) / (2 * self.VarianceTrue)
        TrueNormal = math.exp(Pow) / (self.StdDevTrue * SQRT2PI)
        GX = max(0.01, X + self.GammaOffset)
        FalseGamma = math.pow(GX, self.KFalse - 1) * math.exp(-GX / self.ThetaFalse) / self.GammaDemonFalse
        # Special patch-up code:
        # Toward the edges of the mixture model, odd behavior may occur where one curve falls off
        # slower than the other.  We force very low scores to get a bad odds-true, and very
        # high scores to get a good odds-true.
        if X < self.MeanTrue - self.VarianceTrue:
            FalseGamma = max(FalseGamma, 0.001)
        if X > self.MeanTrue + self.VarianceTrue:
            TrueNormal = max(TrueNormal, 0.001)
        OddsTrue = (TrueNormal * self.PriorProbabilityTrue) / (TrueNormal * self.PriorProbabilityTrue + FalseGamma * (1 - self.PriorProbabilityTrue))
        return OddsTrue
    def InitializeModel(self):
        # Initialize mixture model:
        MinValue = self.MinBin / self.BinMultiplier
        MaxValue = self.MaxBin / self.BinMultiplier
        self.MeanFalse = MinValue + (MaxValue - MinValue) * 0.25
        self.MeanTrue = MaxValue - (MaxValue - MinValue) * 0.25
        self.VarianceFalse = (MaxValue - MinValue) * 0.1
        self.VarianceTrue = (MaxValue - MinValue) * 0.1
        self.PriorProbabilityTrue = 0.1
        if MinValue < 0:
            self.GammaOffset = -MinValue
        elif MinValue > 0.1:
            self.GammaOffset = -MinValue
        else:
            self.GammaOffset = 0
        self.OddsTrue = {}
    def ModelDistribution(self):
        self.ThetaFalse = self.VarianceFalse / (self.MeanFalse + self.GammaOffset)
        self.StdDevTrue = math.sqrt(self.VarianceTrue)
        self.KFalse = (self.MeanFalse + self.GammaOffset) / self.ThetaFalse
        self.GammaDemonFalse = math.pow(self.ThetaFalse, self.KFalse) * Gamma(self.KFalse)
        for Cycle in range(self.MaxCycleCount):
            self.Cycle = Cycle
            self.EstimateOddsTrue()
            self.ComputeDistributionParameters()
    def EstimateOddsTrue(self):
        """
        One half of the E/M cycle: Estimate the probability true for each bin.
        """
        # For each bin, compute the probability that it's true:
        BestOddsTrue = 0
        for Bin in range(self.MinBin, self.MaxBin):
            X = Bin / self.BinMultiplier
            self.OddsTrue[Bin] = self.GetOddsTrue(X)
            # Somewhat hacky: If the left tail of the normal distribution falls off more slowly
            # than that of the gamma distribution, the value of OddsTrue is often rather
            # high for these (very bad!) bins.  We fix that.
            if Bin < 0:
                self.OddsTrue[Bin] = min(self.OddsTrue[Bin], 1 / float(-Bin))
            # Somewhat hacky: If the right tail of the normal distribution falls off too quickly,
            # then the odds true will decay:
            BestOddsTrue = max(BestOddsTrue, self.OddsTrue[Bin])
            if X >= self.MeanTrue:
                self.OddsTrue[Bin] = max(BestOddsTrue, self.OddsTrue[Bin])
            #print "%s: %s"%(X, self.OddsTrue[Bin])
    def ComputeDistributionParameters(self):
        """
        One half of the E/M cycle: Optimize the distribution parameters.
        """
        # Compute the new mean and variance for the true and the false distributions:
        self.CountTrue = 0
        self.MeanTrue = 0
        self.CountFalse = 0
        self.MeanFalse = 0
        for Bin in range(self.MinBin, self.MaxBin):
            X = Bin / self.BinMultiplier
            Count = self.ScoreHistogram.get(Bin, 0)
            self.MeanTrue += X * self.OddsTrue[Bin] * Count
            self.CountTrue += self.OddsTrue[Bin] * Count
            self.MeanFalse += X * (1.0 - self.OddsTrue[Bin]) * Count
            self.CountFalse += (1.0 - self.OddsTrue[Bin]) * Count
        if self.CountTrue <= 0 or self.CountFalse <= 0:
            print "** Error: Unable to fit mixture model.  Appears to be %s true and %s false matches."%(self.CountTrue, self.CountFalse)
            return 0
        self.MeanTrue /= self.CountTrue
        self.MeanFalse /= self.CountFalse
        self.PriorProbabilityTrue = self.CountTrue / (self.CountTrue + self.CountFalse)
        # Adjust GammaOffset, if the false distribution's mean is getting close to 0:
        if self.MeanFalse + self.GammaOffset < 0.1:
            print "False distribution mean is small; BUMP gamma offset up"
            self.GammaOffset += 0.5
        ##################################
        # Compute the new variation for the true and the false distributions:
        self.VarianceTrue = 0
        self.VarianceFalse = 0
        for Bin in range(self.MinBin, self.MaxBin):
            X = Bin / self.BinMultiplier
            Count = self.ScoreHistogram.get(Bin, 0)
            try:
                self.VarianceTrue += (X - self.MeanTrue)**2 * Count * self.OddsTrue[Bin]
                self.VarianceFalse += (X - self.MeanFalse)**2 * Count * (1.0 - self.OddsTrue[Bin])
            except:
                print X
                print self.MeanTrue
                print self.MeanFalse
                print self.OddsTrue[Bin]
                raise
        self.VarianceTrue /= self.CountTrue
        self.StdDevTrue = math.sqrt(self.VarianceTrue)
        self.VarianceFalse /= self.CountFalse
        #print "  True mean %.4f var %.4f"%(self.MeanTrue, self.VarianceTrue)
        #print " False mean %.4f var %.4f"%(self.MeanFalse, self.VarianceFalse)
        self.ThetaFalse = self.VarianceFalse / (self.MeanFalse + self.GammaOffset)
        self.KFalse = (self.MeanFalse + self.GammaOffset) / self.ThetaFalse
        self.GammaDemonFalse = math.pow(self.ThetaFalse, self.KFalse) * Gamma(self.KFalse)
        if self.VerboseFlag:
            print "-----------------------"
            print "Cycle %s report:"%self.Cycle
            print "Theta %.4f K %.4f GammaDenominator %.8f GammaOffset %.2f"%(self.ThetaFalse, self.KFalse, self.GammaDemonFalse, self.GammaOffset)
            print "True: Count %s mean %s variance %s"%(self.CountTrue, self.MeanTrue, self.VarianceTrue)
            print "False: Count %s mean %s variance %s"%(self.CountFalse, self.MeanFalse, self.VarianceFalse)
            print "Prior probability true: %s"%self.PriorProbabilityTrue
    def PlotDistribution(self, FileName):
        File = open(FileName, "wb")
        Header = "Bin\tValue\tHistogram\tOddsTrue\tTrueNormal\tFalseGamma\tMixture\t"
        File.write(Header + "\n")
        for Bin in range(self.MinBin, self.MaxBin):
            Str = "%s\t%s\t"%(Bin, Bin / self.BinMultiplier)
            Str += "%s\t%s\t"%(self.ScoreHistogram.get(Bin, 0), self.OddsTrue[Bin])
            X = Bin / self.BinMultiplier
            # Plot gamma and normal curves:
            Pow = - ((X - self.MeanTrue)**2) / (2 * self.VarianceTrue)
            TrueNormal = math.exp(Pow) / (self.StdDevTrue * SQRT2PI)
            GX = max(0.01, X + self.GammaOffset)
            FalseGamma = math.pow(GX, self.KFalse - 1) * math.exp(-GX / self.ThetaFalse) / self.GammaDemonFalse
            Str += "%s\t%s\t"%(TrueNormal, FalseGamma)
##            # Plot gamma and noraml CDF:
##            ErfArg = (X - self.MeanTrue) / (self.StdDevTrue * SQRT2)
##            NormalCDF = 0.5 + 0.5 * PyInspect.erf(ErfArg)
##            GX = max(0.01, X + self.GammaOffset)
##            GammaCDF = PyInspect.GammaIncomplete(self.KFalse, GX / self.ThetaFalse) #/ Gamma(self.KFalse)
##            Str += "%s\t%s\t"%(NormalCDF, GammaCDF)
            MergedMixture = TrueNormal * self.PriorProbabilityTrue
            MergedMixture += FalseGamma * (1.0 - self.PriorProbabilityTrue)
            Str += "%s\t"%(MergedMixture)
            File.write(Str + "\n")
    def PickleSelf(self, File):
        cPickle.dump(self.BinMultiplier, File)
        cPickle.dump(self.PriorProbabilityTrue, File)
        cPickle.dump(self.MeanTrue, File)
        cPickle.dump(self.VarianceTrue, File)
        cPickle.dump(self.StdDevTrue, File)
        cPickle.dump(self.GammaOffset, File)
        cPickle.dump(self.KFalse, File)
        cPickle.dump(self.ThetaFalse, File)
        cPickle.dump(self.GammaDemonFalse, File)

def UnpickleMixtureModel(File):
    Model = MixtureModelClass()
    Model.BinMultiplier = cPickle.load(File)
    Model.PriorProbabilityTrue = cPickle.load(File)
    Model.MeanTrue = cPickle.load(File)
    Model.VarianceTrue = cPickle.load(File)
    Model.StdDevTrue = cPickle.load(File)
    Model.GammaOffset = cPickle.load(File)
    Model.KFalse = cPickle.load(File)
    Model.ThetaFalse = cPickle.load(File)
    Model.GammaDemonFalse = cPickle.load(File)
    return Model

class FeatureVector:
    def __init__(self):
        self.FileBits = []
        self.Features = []
        self.ScaledFeatures = None
        self.TrueFlag = 0
        self.Score = 0 # as assigned by an owning model
        
class FeatureSetClass:
    """
    A feature-set is a list of TRUE tuples, and a list of FALSE tuples.  Normally there
    is a FeatureSetClass for testing, and one for training.
    """
    def __init__(self):
        self.TrueVectors = []
        self.FalseVectors = []
        self.AllVectors = []
        self.TrueCount = 0
        self.FalseCount = 0
        self.Count = 0
        self.PriorProbabilityFalse = 0.5
    def SetCounts(self):
        self.Count = len(self.AllVectors)
        self.TrueCount = len(self.TrueVectors)
        self.FalseCount = len(self.FalseVectors)
        if len(self.AllVectors):
            self.Size = len(self.AllVectors[0].Features)
    def FindFeatureRanges(self):
        """
        Simple scaling function: Find min and max values to push features into [-1, 1]
        """
        Values = []
        Vector = self.AllVectors[0]
        Size = len(Vector.Features)
        print "SIZE:", Size
        for X in range(Size):
            Values.append([])
        for Vector in self.AllVectors:
            for X in range(Size):
                Values[X].append(Vector.Features[X])
        self.MinValues = []
        self.MaxValues = []
        for X in range(Size):
            Values[X].sort()
            ValIndex = int(round(len(Values[X]) * 0.025))
            MinValue = Values[X][ValIndex]
            self.MinValues.append(MinValue)
            ValIndex = int(round(len(Values[X]) * 0.975))
            MaxValue = Values[X][ValIndex]
            self.MaxValues.append(MaxValue)
##        print "Range:"
##        for X in range(Size):
##            print "%s: %.4f-%.4f"%(X, self.MinValues[X], self.MaxValues[X])
        pass
    def ScaleFeatures(self):
        """
        Simple scaling function: Pushes 90% of feature values into the range [-1, 1].
        Assumes that self.MinValues and self.MaxValues are set!
        """
        self.Size = len(self.AllVectors[0].Features)
        for Vector in self.AllVectors:
            Vector.ScaledFeatures = [0]*self.Size
        for X in range(self.Size):
            HalfRange = (self.MaxValues[X] - self.MinValues[X]) / 2.0
            if not HalfRange:
                continue
            #print "Feature %s: Range %s...%s"%(X, self.MinValues[X], self.MaxValues[X])
            for Vector in self.AllVectors:
                Vector.ScaledFeatures[X] = (Vector.Features[X] - self.MinValues[X]) / HalfRange - 1.0
    def __str__(self):
        return "<%sT %sF>"%(len(self.TrueVectors), len(self.FalseVectors))
    def GetPriorProbabilityFalse(self, DBTrueToFalseRatio):
        """
        Compute the prior probability that an arbitrary peptide is false.
        """
        # In a 1:1 database, there's 1 bogus peptide in a valid protein
        # for each bogus peptide in an invalid protein; in that case, DBTrueToFalseRatio is 1.0
        # In a 1:99 database, the ratio is 1/99.
        FalseWithinTrue = self.FalseCount * DBTrueToFalseRatio
        if FalseWithinTrue >= self.TrueCount:
            # Uh-oh; there are FEWER peptides from valid proteins than we would expect
            # to see by chance!  Let's (arbitrarily) cap the prior probability false
            # at 99%.
            print "Warning: FalseWithinTrue = %s >= %s!"%(FalseWithinTrue, self.TrueCount)
            self.PriorProbabilityFalse = 0.99
            return
        VectorCount = len(self.AllVectors)
        self.PriorProbabilityFalse = (VectorCount - (self.TrueCount - FalseWithinTrue)) / float(VectorCount)
        print "==>>PriorProbabilityFalse: %s"%(self.PriorProbabilityFalse)
    def SaveTabDelimited(self, File):
        if type(File) == type(""):
            File = open(File, "wb")
            CloseFlag = 1
        else:
            CloseFlag = 0
        #File = open(FileName, "wb")
        for VectorIndex in range(len(self.AllVectors)):
            Vector = self.AllVectors[VectorIndex]
            String = "%s\t%s\t"%(VectorIndex, Vector.TrueFlag)
            for Value in Vector.Features:
                String += "%s\t"%Value
            File.write(String + "\n")
        if CloseFlag:
            File.close()

def LoadGeneralModel(FileName):
    File = open(FileName, "rb")
    ModelType = cPickle.load(File)
    File.close()
    if ModelType == "LDA":
        Model = LDAModel()
        Model.LoadModel(FileName)
    elif ModelType == "SVM":
        Model = SVMModel()
        Model.LoadModel(FileName)
    elif ModelType == "LOGIT":
        Model = LogitModel()
        Model.LoadModel(FileName)
    else:
        print "** Error: Unable to load model type '%s'"%ModelType
        return None
    return Model

class LearnerClass:
    def __init__(self, FeatureList = None):
        # The entries in FeatureList are indices into the
        # available features of our training and testing sets.
        self.FeatureList = FeatureList
        # OddsTrue[Bin] = probability that an instance with a score
        # in this bin or HIGHER is correct.
        self.OddsTrue = {}
        self.PValue = {}
        # Bin = int(round(Score * self.BinScalingFactor))
        self.BinScalingFactor = 10
        self.MixtureModel = None
    def SaveModel(self, FileName):
        raise ValueError, "Abstract method - override in subclass!"
    def LoadModel(self, FileName):
        raise ValueError, "Abstract method - override in subclass!"
    def Train(self, FeatureSet):
        raise ValueError, "Abstract method - override in subclass!"
    def Test(self, FeatureSet):
        raise ValueError, "Abstract method - override in subclass!"
    def ReportROC(self, FeatureSet, OutputFileName = None):
        SortedList = []
        for Vector in FeatureSet.AllVectors:
            SortedList.append((Vector.Score, random.random(), Vector))
        SortedList.sort()
        SortedList.reverse()
        OverallTrueCount = 0
        OverallFalseCount = 0
        # If there are many many vectors, then we'll end up with an unwieldy curve that's
        # awkward to plot.  So, consider PHASING things:
        Slice = (len(SortedList) / 30000) + 1 # 2 or larger
        OldSortedList = SortedList
        SortedList = []
        print "SLICE roc-curve list: Take every %sth entry"%Slice
        for X in range(len(OldSortedList)):
            if X % Slice == 0:
                SortedList.append(OldSortedList[X])
                Vector = OldSortedList[X][-1]
                if Vector.TrueFlag:
                    OverallTrueCount += 1
                else:
                    OverallFalseCount += 1
        OldSortedList = None
        TrueCount = 0
        FalseCount = 0
        if OutputFileName:
            ROCCurvePlotFile = open(OutputFileName, "wb")
        RowCount = 0
        ROCTPForFP = {}
        ROCTPForFPCount = {}
        for (Score, Dummy, Vector) in SortedList:
            RowCount += 1
            if Vector.TrueFlag:
                TrueCount += 1
            else:
                FalseCount += 1
            OverallTPRate = TrueCount / float(max(1, OverallTrueCount))
            OverallFPRate = FalseCount / float(max(1, OverallFalseCount))
            Bin = int(round(OverallFPRate * 100))
            ROCTPForFP[Bin] = ROCTPForFP.get(Bin, 0) + OverallTPRate
            ROCTPForFPCount[Bin] = ROCTPForFPCount.get(Bin, 0) + 1
            if OutputFileName:
                ROCCurvePlotFile.write("%s\t%s\t%s\t%s\t%s\t\n"%(RowCount, TrueCount, FalseCount, OverallFPRate, OverallTPRate))
        if OutputFileName:
            ROCCurvePlotFile.close()
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
        print "ROC curve area:", ROCArea
    def ReportAccuracy(self, FeatureSet, ROCFilePath = None):
        """
        Called after Test(FeatureSet), to measure how well we did at separating
        the true and false vectors by score.  Compute OddsTrue, as well.
        """
        SortedList = []
        for Vector in FeatureSet.AllVectors:
            SortedList.append((Vector.Score, Vector))
        # sort from HIGHEST to LOWEST score:
        SortedList.sort()
        SortedList.reverse()
        #self.ComputeOddsTrue(SortedList)
        self.ComputePValues(SortedList)
        #self.ComputePValues(SortedList, FeatureSet.PriorProbabilityFalse)
        Rates = [0.05, 0.01, 0.1, 0.5, 0.005]
        CountsByRate = [0, 0, 0, 0, 0, 0]
        CumulativeTrue = 0
        CumulativeFalse = 0
        Cumulative = 0
        PrevScore = None
        for (Score, Vector) in SortedList:
            if Vector.TrueFlag:
                CumulativeTrue += 1
            else:
                CumulativeFalse += 1
            Cumulative += 1
            FractionFalse = CumulativeFalse / float(Cumulative)
            if Score != PrevScore:
                for RateIndex in range(len(Rates)):
                    if FractionFalse < Rates[RateIndex]:
                        CountsByRate[RateIndex] = CumulativeTrue
                #print Score, FractionFalse, CumulativeTrue, CumulativeFalse
            PrevScore = Score
        print "Counts by FDRate: %d@5%% %d@1%% %d@10%% %d@50%%"%(CountsByRate[0], CountsByRate[1], CountsByRate[2], CountsByRate[3])
        SensitivityByRate = []
        for Count in CountsByRate:
            SensitivityByRate.append(100 * Count / float(max(1, FeatureSet.TrueCount)))
        print "FeatureSet.TrueCount:", FeatureSet.TrueCount
        print "  Sensitivity: %.2f@5%% %.2f@1%% %.2f@10%% %.2f@50%%"%(SensitivityByRate[0], SensitivityByRate[1], SensitivityByRate[2], SensitivityByRate[3])
        if ROCFilePath:
            self.ReportROC(FeatureSet, ROCFilePath)
        return CountsByRate
    def ComputeOddsTrue(self, SortedList):
        "DEPRECATED; use ComputePValue instead"
        if len(SortedList) < 200:
            BlockSize = len(SortedList) / 4
        else:
            BlockSize = 100
        WindowTrueSum = 0
        WindowFalseSum = 0
        for Entry in SortedList[:BlockSize - 1]:
            if Entry[1].TrueFlag:
                WindowTrueSum += 1
            else:
                WindowFalseSum += 1
        for Index in range(len(SortedList)):
            # Add one entry to the window:
            if Index + BlockSize < len(SortedList):
                Entry = SortedList[Index + BlockSize]
                if Entry[1].TrueFlag:
                    WindowTrueSum += 1
                else:
                    WindowFalseSum += 1
            # Compute the probability-true for this window:
            Vector = SortedList[Index][1]
            OddsTrue = WindowTrueSum / float(WindowTrueSum + WindowFalseSum)
            Bin = int(round(Vector.Score * self.BinScalingFactor))
            self.OddsTrue[Bin] = OddsTrue
            # Remove leftmost entry from the window:
            if Index >= BlockSize:
                Entry = SortedList[Index - BlockSize]
                if Entry[1].TrueFlag:
                    WindowTrueSum -= 1
                else:
                    WindowFalseSum -= 1
    def GetPValue(self, Score):
        if self.MixtureModel:
            return 1.0 - self.MixtureModel.GetOddsTrue(Score)
        Bin = int(round(Score * self.BinScalingFactor))
        Keys = self.PValue.keys()
        MinKey = min(Keys)
        MaxKey = max(Keys)
        if Bin < MinKey:
            return self.PValue[MinKey]
        if Bin > MaxKey:
            return self.PValue[MaxKey]
        return self.PValue[Bin]
    def ComputePValuesMixtureModel(self, SortedList):
        """
        Our feature-set has an empirical distribution of scores.  We'll approximate
        this distribution as a mixture of two distributions: gamma (false) and
        normal (true).  Then we'll derive p-value (probability false) for each
        score-bin.
        """
        Scores = []
        for (Score, Vector) in SortedList:
            Scores.append(Score)
        #for Bin in range(Model.MinBin, Model.MaxBin):
        #    self.PValue[Bin] = 1.0 - Model.OddsTrue[Bin]
        self.MixtureModel = MixtureModelClass(self.GetMMBinMultiplier())
        self.MixtureModel.Model(Scores)
        for (Score, Vector) in SortedList:
            Vector.PValue = 1.0 - self.MixtureModel.GetOddsTrue(Score)        
        self.MixtureModel.PlotDistribution("PValues.txt")
    def ComputePValues(self, SortedList):
        self.ComputePValuesMixtureModel(SortedList)
    def ComputePValuesEmpirical(self, SortedList, PriorProbabilityFalse):
        """
        DEPRECATED - called only if mixture model fails.
        
        Input SortedList is a list of the form (ModelScore, FeatureVector), sorted from
        highest to lowest ModelScore.
        
        self.PValue[Bin] is the probability that a peptide P is FALSE, given a score of
        Bin or better.  Formally, it equals P(not P | S(P)>=Bin).  By Bayes' Theorem, this
        is equal to:
        P(S(P)>=Bin | not P) * P(not P) / P(S(P)>=Bin) 
        """
        # Passed in: PriorProbabilityFalse == P(not P)
        # And set up dictionaries PValueTotals / PValueCounts to
        #   compute ProbRatio == P(S(P)>=Bin | not P) / P(S(P)>=Bin)
        #PriorProbabilityFalse = 0
        TrueInstanceCount = 0
        HighScoringCount = 0
        HighScoringFalseInstanceCount = 0
        TotalInstances = len(SortedList)
        PValueTotals = {}
        PValueCounts = {}
        TotalFalseInstances = 0
        for (Score, Vector) in SortedList:
            if not Vector.TrueFlag:
                TotalFalseInstances += 1
        HistoGood = {}
        HistoBad = {}
        TempOutputFile = open("TempPTMPValue.txt", "wb")
        TempOutputFile.write("Bin\tHighFalse\tTotalInstances\tTotalFalseInstances\tHighScoringCount\t\n")
        for (Score, Vector) in SortedList:
            HighScoringCount += 1
            if Vector.TrueFlag:
                TrueInstanceCount += 1
            else:
                HighScoringFalseInstanceCount += 1
            ProbRatio = (HighScoringFalseInstanceCount * TotalInstances) / float(TotalFalseInstances * HighScoringCount)
            #ProbRatio = FalseInstanceCount / float(TrueInstanceCount + FalseInstanceCount)
            Bin = int(round(Score * self.BinScalingFactor))
            PValueTotals[Bin] = PValueTotals.get(Bin, 0) + ProbRatio
            PValueCounts[Bin] = PValueCounts.get(Bin, 0) + 1
            TempOutputFile.write("%s\t%s\t%s\t%s\t%s\t\n"%(Bin, HighScoringFalseInstanceCount, TotalInstances, TotalFalseInstances, HighScoringCount))
            if Vector.TrueFlag:
                HistoGood[Bin] = HistoGood.get(Bin, 0) + 1
            else:
                HistoBad[Bin] = HistoBad.get(Bin, 0) + 1
        #PriorProbabilityFalse = FalseInstanceCount / float(FalseInstanceCount + TrueInstanceCount)
        Keys = PValueTotals.keys()
        Keys.sort()
        for Bin in Keys:
            AverageProbRatio = PValueTotals[Bin] / float(PValueCounts[Bin])
            self.PValue[Bin] = AverageProbRatio * PriorProbabilityFalse
            self.PValue[Bin] = max(self.PValue[Bin], 0.0001)
            print "%s: %s"%(Bin, self.PValue[Bin])
        TempOutputFile.write("\n\n\n")
        for Bin in Keys:
            TempOutputFile.write("%s\t%s\t%s\t\n"%(Bin, HistoGood.get(Bin, 0), HistoBad.get(Bin, 0)))
        ############################################################
        # Interpolate p-values, for missing bins:
        Keys = self.PValue.keys()
        Keys.sort()
        MinKey = min(Keys)
        MaxKey = max(Keys)
        for Bin in range(MinKey, MaxKey):
            if self.PValue.has_key(Bin):
                PrevBin = Bin
                PrevPValue = self.PValue[Bin]
                continue
            # Find the next bin:
            for NextBin in range(Bin + 1, MaxKey + 1):
                if self.PValue.has_key(NextBin):
                    NextPValue = self.PValue[NextBin]
                    break
            # Interpolate from (PrevBin, PrevPValue) to (NextBin, NextPValue):
            Slope = (NextPValue - PrevPValue) / float(NextBin - PrevBin)
            Intermediate = PrevPValue + Slope * (Bin - PrevBin)
            self.PValue[Bin] = Intermediate
    def GetMMBinMultiplier(self):
        return 10.0 #default


class LDAModel(LearnerClass):
    def __init__(self, FeatureList = None):
        if FeatureList:
            self.Size = len(FeatureList)
        LearnerClass.__init__(self, FeatureList)
    def GetCovarianceArray(self, VectorList):
        VectorCount = float(len(VectorList))
        C = zeros((self.Size, self.Size), FloatType)
        for Vector in VectorList:
            for X in range(self.Size):
                for Y in range(self.Size):
                    C[X][Y] += Vector[X] * Vector[Y] / VectorCount
        return C
    def LoadModel(self, FileName):
        File = open(FileName, "rb")
        cPickle.load(File) # model type
        self.FeatureList = cPickle.load(File)
        self.PValue = cPickle.load(File)
        self.MinValues = cPickle.load(File)
        self.MaxValues = cPickle.load(File)
        self.CI = cPickle.load(File)
        self.MeanGood = cPickle.load(File)
        self.ConstantGood = cPickle.load(File)
        self.MeanBad = cPickle.load(File)
        self.ConstantBad = cPickle.load(File)
        self.Size = len(self.FeatureList)
        self.MixtureModel = UnpickleMixtureModel(File)
        # Verbose stuff:
        print "\n>>>PyLoadLDAModel(%s)"%FileName
        print "Features: %s"%self.Size
        print "MinValues: %.4f...%.4f"%(self.MinValues[0], self.MinValues[-1])
        print "MaxValues: %.4f...%.4f"%(self.MaxValues[0], self.MaxValues[-1])
        print "MeanGood: %.4f...%.4f"%(self.MeanGood[0], self.MeanGood[-1])
        print "MeanBad: %.4f...%.4f"%(self.MeanBad[0], self.MeanBad[-1])
        if self.Size > 1:
            print "CI: %.4f, %.4f...%.4f, %.4f"%(self.CI[0][0], self.CI[0][1],
                self.CI[self.Size - 1][self.Size - 2],
                self.CI[self.Size - 1][self.Size - 1]
                                         )
        print "ConstTrue %.4f, ConstFalse %.4f"%(self.ConstantGood, self.ConstantBad)
        File.close()
    def SaveBinaryModel(self, FileName):
        """
        Write out a binary representation of this model.  
        """
        File = open(FileName, "wb")
        File.write(struct.pack("<i", self.Size))
        for FeatureIndex in range(self.Size):
            File.write(struct.pack("<d", self.MinValues[FeatureIndex]))
        for FeatureIndex in range(self.Size):
            File.write(struct.pack("<d", self.MaxValues[FeatureIndex]))
        for FeatureIndex in range(self.Size):
            File.write(struct.pack("<d", self.MeanGood[FeatureIndex]))
        for FeatureIndex in range(self.Size):
            File.write(struct.pack("<d", self.MeanBad[FeatureIndex]))
        File.write(struct.pack("<d", self.ConstantGood))
        File.write(struct.pack("<d", self.ConstantBad))
        for Row in range(self.Size):
            for Column in range(self.Size):
                File.write(struct.pack("<d", self.CI[Row][Column]))
        File.close()
    def SaveModel(self, FileName):
        File = open(FileName, "wb")
        cPickle.dump("LDA", File)
        cPickle.dump(self.FeatureList, File)
        cPickle.dump(self.PValue, File)
        cPickle.dump(self.MinValues, File)
        cPickle.dump(self.MaxValues, File)
        cPickle.dump(self.CI, File)
        cPickle.dump(self.MeanGood, File)
        cPickle.dump(self.ConstantGood, File)
        cPickle.dump(self.MeanBad, File)
        cPickle.dump(self.ConstantBad, File)
        try:
            self.MixtureModel.PickleSelf(File)
        except:
            cPickle.dump(None, File)
        File.close()
    def ScoreInstance(self, RawFeatures):
        Features = []
        for FeatureIndex in self.FeatureList:
            Features.append(RawFeatures[FeatureIndex])
        for FeatureIndex in range(self.Size):
            X = self.FeatureList[FeatureIndex]
            HalfRange = (self.MaxValues[X] - self.MinValues[X]) / 2.0
            if not HalfRange:
                continue
            Features[FeatureIndex] = (Features[FeatureIndex] - self.MinValues[X]) / HalfRange - 1.0
        CIProduct = MatrixMultiply(self.CI, Features)
        ReadingGood = MatrixMultiply(self.MeanGood, CIProduct) + self.ConstantGood
        ReadingBad = MatrixMultiply(self.MeanBad, CIProduct) + self.ConstantBad
        return (ReadingGood - ReadingBad)
    def Test(self, FeatureSet):
        # Scale features, according to our trained scaling
        FeatureSet.MinValues = self.MinValues
        FeatureSet.MaxValues = self.MaxValues
        FeatureSet.ScaleFeatures()
        # Compute scores:
        for Vector in FeatureSet.AllVectors:
            FixedVector = []
            for FeatureIndex in self.FeatureList:
                FixedVector.append(Vector.ScaledFeatures[FeatureIndex])
            CIProduct = MatrixMultiply(self.CI, FixedVector)
            ReadingGood = MatrixMultiply(self.MeanGood, CIProduct) + self.ConstantGood
            ReadingBad = MatrixMultiply(self.MeanBad, CIProduct) + self.ConstantBad
            Vector.Score = (ReadingGood - ReadingBad)
    def Train(self, FeatureSet, VerboseFlag = 0):
        # Get the feature range (training only):
        FeatureSet.FindFeatureRanges()
        self.MinValues = FeatureSet.MinValues
        self.MaxValues = FeatureSet.MaxValues
        # Sanity-checking: If a feature's range is a single point,
        # then it's not useful - AND, it will generate a non-invertible
        # matrix.  So, let's filter out any such features.
        InputFeatureList = self.FeatureList
        self.FeatureList = []
        for FeatureIndex in InputFeatureList:
            if self.MinValues[FeatureIndex] < self.MaxValues[FeatureIndex]:
                self.FeatureList.append(FeatureIndex)
            else:
                print "* Warning: Discarding feature '%s', every entry is %s"%(FeatureIndex, self.MinValues[FeatureIndex])
        self.Size = len(self.FeatureList)
        if self.Size == 0:
            print "<< no features - bailing out >>"
            return
        # Scale features:
        FeatureSet.ScaleFeatures()
        AllVectors = []
        TrueVectors = []
        FalseVectors = []
        for Vector in FeatureSet.AllVectors:
            ScaledVector = []
            for FeatureIndex in self.FeatureList:
                ScaledVector.append(Vector.ScaledFeatures[FeatureIndex])
            AllVectors.append(ScaledVector)
            if Vector.TrueFlag:
                TrueVectors.append(ScaledVector)
            else:
                FalseVectors.append(ScaledVector)
        print "First true vector:", TrueVectors[0]
        print "First false vector:", FalseVectors[0]
##        # Temp: Ensure the vector lists are the same size!
##        VectorCount = min(len(TrueVectors), len(FalseVectors))
##        random.shuffle(FalseVectors)
##        FalseVectors = FalseVectors[:VectorCount]
        ############################################################
        # Compute the mean vectors (training only):
        self.MeanGood = [0] * self.Size
        self.MeanBad = [0] * self.Size
        self.MeanGlobal = [0] * self.Size
        for Vector in TrueVectors:
            for X in range(self.Size):
                self.MeanGlobal[X] += Vector[X] / float(FeatureSet.Count)
                self.MeanGood[X] += Vector[X] / float(FeatureSet.TrueCount)
        for Vector in FalseVectors:
            for X in range(self.Size):
                self.MeanGlobal[X] += Vector[X] / float(FeatureSet.Count)
                self.MeanBad[X] += Vector[X] / float(FeatureSet.FalseCount)
        if VerboseFlag:
            print "MeanGood:\n  ",
            for Value in self.MeanGood:
                print "%.3f"%Value,
            print
            print "MeanBad:\n  ",
            for Value in self.MeanBad:
                print "%.3f"%Value,
            print
            print "MeanGlobal:\n  ",
            for Value in self.MeanGlobal:
                print "%.3f"%Value,
            print
        ############################################################
        # Compute mean-corrected vectors:
        MeanCorrectedGoodVectors = []
        MeanCorrectedBadVectors = []
        for Vector in TrueVectors:
            NewVector = []
            for X in range(self.Size):
                NewVector.append(Vector[X] - self.MeanGlobal[X])
            MeanCorrectedGoodVectors.append(NewVector)
        for Vector in FalseVectors:
            NewVector = []
            for X in range(self.Size):
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
        CovarArrayFull = zeros((self.Size, self.Size), FloatType)
        for X in range(self.Size):
            for Y in range(self.Size):
                CovarArrayFull[X][Y] += CovarArrayGood[X][Y] * FeatureSet.TrueCount / float(FeatureSet.Count)
                CovarArrayFull[X][Y] += CovarArrayBad[X][Y] * FeatureSet.FalseCount / float(FeatureSet.Count)
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
        self.GoodMuC = MatrixMultiply(self.CI, self.MeanGood)
        if VerboseFlag:
            print "GoodMuC:", self.GoodMuC
        self.BadMuC = MatrixMultiply(self.CI, self.MeanBad)
        if VerboseFlag:
            print "BadMuC:", self.BadMuC
        self.ConstantGood = -MatrixMultiply(self.MeanGood, self.GoodMuC) / 2.0
        self.ConstantBad = -MatrixMultiply(self.MeanBad, self.BadMuC) / 2.0
        #if VerboseFlag:
        print "LDA Constant good %.4f constant bad %.4f"%(self.ConstantGood, self.ConstantBad)

class SVMModel(LearnerClass):
    def __init__(self, FeatureList = None):
        self.Scaling = None
        self.SupportVectors = None
        self.PySVMReadyFlag = 0
        LearnerClass.__init__(self, FeatureList)
    def WriteSVMFeaturesToFile(self, FilePath, FeatureSet, ForceEqualCounts = 0):
        #print "TRUE vectors %s, FALSE vectors %s"%(len(FeatureSet.TrueVectors), len(FeatureSet.FalseVectors))
        #print "TRUE count %s, FALSE count %s"%(FeatureSet.TrueCount, FeatureSet.FalseCount)
        File = open(FilePath, "wb")
        if ForceEqualCounts:
            # Shuffle the tuples:
            TrueVectors = FeatureSet.TrueVectors[:]
            random.shuffle(TrueVectors)
            FalseVectors = FeatureSet.FalseVectors[:]
            random.shuffle(FalseVectors)
            # Try writing out equal numbers of true and false tuples:
            MaxIndex = min(FeatureSet.TrueCount, FeatureSet.FalseCount)
            MaxIndex = min(MaxIndex, MaxSVMFeatureCount)
            TrueVectors = TrueVectors[:MaxIndex]
            FalseVectors = FalseVectors[:MaxIndex]
        else:
            TrueVectors = FeatureSet.TrueVectors
            FalseVectors = FeatureSet.FalseVectors
        for Vector in TrueVectors:
            Str = "+1 "
            for FeatureIndex in range(len(self.FeatureList)):
                Str += "%d:%.8f "%(FeatureIndex + 1, Vector.Features[self.FeatureList[FeatureIndex]])
            File.write(Str + "\n")
        for Vector in FalseVectors:
            Str = "-1 "
            for FeatureIndex in range(len(self.FeatureList)):
                Str += "%d:%.8f "%(FeatureIndex + 1, Vector.Features[self.FeatureList[FeatureIndex]])
            File.write(Str + "\n")
        File.close()
    def Train(self, FeatureSet, VerboseFlag = 0):
        print "TRAINSVM()...", FeatureSet
        TempFeaturesFileName = "PTMFeatures.SVM.txt"
        TempScalingFileName = "PTMFeaturesSVMScale.txt"
        TempScaledFeaturesFileName = "PTMFeatures.SVMScaled.txt"
        TempModelFileName = "PTMFeatures.SVMScaled.txt.model"        
        # Write feature vectors, forcing equal true and false instance-counts:
        self.WriteSVMFeaturesToFile(TempFeaturesFileName, FeatureSet, 1)
        ###############################################################
        # SCALE the features, and remember the scaling:
        Command = r"%s -s %s %s > %s"%(PATH_SVMSCALE, TempScalingFileName, TempFeaturesFileName, TempScaledFeaturesFileName)
        print Command
        os.system(Command)
        # Read the scaling limits, for later use:
        File = open(TempScalingFileName, "rb")
        self.Scaling = File.read()
        File.close()
        os.remove(TempScalingFileName)
        print "Train!"
        ###############################################################
        # TRAIN the model.  We won't use cross-validation here, because in the future there'll be
        # a testing-set.
        Command = r"%s %s"%(PATH_SVMTRAIN, TempScaledFeaturesFileName)
        print Command
        os.system(Command)
        File = open(TempModelFileName, "rb")
        self.SupportVectors = File.read()
        File.close()
        ###############################################
        # Clean up temp-files:
        os.remove(TempFeaturesFileName)
        os.remove(TempScaledFeaturesFileName)
        os.remove(TempModelFileName)
    def Test(self, FeatureSet):
        if not self.Scaling or not self.SupportVectors:
            print "Error in SVMModel.Test(): We haven't trained (or loaded) yet!"
            return
        TempFeaturesFileName = "TestFeatures.SVM.txt"
        TempScalingFileName = "PTMFeaturesSVMScale.txt"
        TempScaledFeaturesFileName = "TestFeatures.SVMScaled.txt"
        TempModelFileName = "SVM.model"
        TempOutputFileName = "SVMPrediction.txt"
        TrueFlags = []
        for Tuple in FeatureSet.TrueVectors:
            TrueFlags.append(1)
        for Tuple in FeatureSet.FalseVectors:
            TrueFlags.append(0)        
        ########################################################################
        # WRITE the testing set to file:
        self.WriteSVMFeaturesToFile(TempFeaturesFileName, FeatureSet)
        # Write our scaling-info and our model to files:
        File = open(TempScalingFileName, "wb")
        File.write(self.Scaling)
        File.close()
        File = open(TempModelFileName, "wb")
        File.write(self.SupportVectors)
        File.close()
        # SCALE the testing set:
        Command = r"%s -r %s %s > %s"%(PATH_SVMSCALE, TempScalingFileName, TempFeaturesFileName, TempScaledFeaturesFileName)
        print Command
        os.system(Command)
        os.remove(TempFeaturesFileName)
        os.remove(TempScalingFileName)
        # Ok, now let's run svmpredict on all the instances in the TESTING set:
##        Command = r"%s %s %s %s"%(PATH_SVMPREDICT, TempScaledFeaturesFileName, TempModelFileName, TempOutputFileName)
##        print Command
##        os.system(Command)
        RunPySVM.Predict(TempScaledFeaturesFileName, TempModelFileName, TempOutputFileName)
        # Now read in the results, and assign scores to the vectors of the set:
        File = open(TempOutputFileName, "rb")
        InstanceIndex = 0
        TrueIndex = 0
        FalseIndex = 0
        SortedList = []
        for FileLine in File.xreadlines():
            Score = float(FileLine)
            if TrueFlags[InstanceIndex]:
                Vector = FeatureSet.TrueVectors[TrueIndex]
                TrueIndex += 1
            else:
                Vector = FeatureSet.FalseVectors[FalseIndex]
                FalseIndex += 1
            Vector.Score = float(FileLine)
            InstanceIndex += 1
        File.close()
        ########################################################################
        # Clean up temp-files:
        os.remove(TempScaledFeaturesFileName)
        os.remove(TempModelFileName)
        os.remove(TempOutputFileName)
    def SaveTextModel(self, Stub):
        ModelPath = "%s.model"%Stub
        File = open(ModelPath, "wb")
        File.write(self.SupportVectors)
        File.close()
        ScalingPath = "%s.range"%Stub
        File = open(ScalingPath, "wb")
        File.write(self.Scaling)
        File.close()
    def SaveModel(self, FileName):
        File = open(FileName, "wb")
        cPickle.dump("SVM", File)
        cPickle.dump(self.FeatureList, File)
        cPickle.dump(self.PValue, File)
        cPickle.dump(self.Scaling, File)
        cPickle.dump(self.SupportVectors, File)
        self.MixtureModel.PickleSelf(File)
        File.close()
    def LoadModel(self, FileName):
        File = open(FileName, "rb")
        cPickle.load(File) # model type
        self.FeatureList = cPickle.load(File)
        self.PValue = cPickle.load(File)
        self.Scaling = cPickle.load(File)
        self.SupportVectors = cPickle.load(File)
        self.MixtureModel = UnpickleMixtureModel(File)
        File.close()
    def PreparePySVM(self):
        """
        Prepare PySVM to score some features using our model.
        """
        # Support vectors:
        TempModelFileName = "TempModel.txt"
        File = open(TempModelFileName, "wb")
        File.write(self.SupportVectors)
        File.close()
        PySVM.LoadModel(TempModelFileName)
        os.remove(TempModelFileName)
        # Feature ranges:
        TempScalingFileName = "TempScaling.txt"
        File = open(TempScalingFileName, "wb")
        File.write(self.Scaling)
        File.close()
        PySVM.LoadScaling(TempScalingFileName)
        os.remove(TempScalingFileName)
        # And now, we can score many vectors quickly!
        self.PySVMReadyFlag = 1
    def ScoreInstance(self, Features):
        """
        Compute the score for this instance. 
        """
        if not self.PySVMReadyFlag:
            self.PreparePySVM()
        Vector = []
        for FeatureIndex in self.FeatureList:
            Vector.append(Features[FeatureIndex])
        Score = PySVM.ScaleAndScore(Vector)
        return Score
    
class LogitModel(LearnerClass):
    """
    A maximum-likelihood logistic regression model.  The model's parameters
    are tuned using the Newton-Raphson algorithm.  See p98 in Hastie/Tibshirani,
    The Elements of Statistical Learning.  
    """
    def GetMMBinMultiplier(self):
        return 40.0 #default
    def ComputePValues(self, SortedList):
        pass # The score we output IS a probability!
    def GetFixedTuples(self, FeatureSet):
        # Return fixed-up tuples.  Keep a random selection of
        # true and of false tuples.  
        VectorSize = len(self.FeatureList) + 1 # add one for the CONSTANT input
        AllTuples = [] # entries (TrueFlag, FeatureTuple)
        KeepCount = min(FeatureSet.TrueCount, FeatureSet.FalseCount, 500)
        random.shuffle(FeatureSet.TrueVectors)
        random.shuffle(FeatureSet.FalseVectors)
        for Vector in FeatureSet.TrueVectors[:KeepCount]:
            FixedTuple = [1.0]
            for FeatureIndex in self.FeatureList:
                FixedTuple.append(Vector.Features[FeatureIndex])
            AllTuples.append((1, tuple(FixedTuple)))
        for Vector in FeatureSet.FalseVectors[:KeepCount]:
            FixedTuple = [1.0]
            for FeatureIndex in self.FeatureList:
                FixedTuple.append(Vector.Features[FeatureIndex])
            AllTuples.append((0, tuple(FixedTuple)))        
        return AllTuples
    def Train(self, FeatureSet, VerboseFlag = 0):
        VectorSize = len(self.FeatureList) + 1 # add one for the CONSTANT input
        AllTuples = self.GetFixedTuples(FeatureSet)
        random.shuffle(AllTuples)
        #################################################################
        # Train the model - set the weight-vector self.Beta
        # Initialize vector self.Beta, all zeroes:
        self.Beta = zeros(VectorSize)
        TupleCount = len(AllTuples)
        # Initialize vector Y, indicating which vectors are true:
        Y = zeros(TupleCount, FloatType)
        for I in range(TupleCount):
            if AllTuples[I][0]:
                Y[I] = 1.0
            else:
                Y[I] = 0.0
        # Initialize the input matrix X:
        X = zeros((TupleCount, VectorSize), FloatType)
        for I in range(TupleCount):
            #X[I][0] = 1.0
            for J in range(VectorSize):
                X[I][J] = AllTuples[I][1][J]
        XT = transpose(X)
        PrevLogLikelihood = None
        CycleCount = 0
        while 1:
            # Compute the current log-likelihood:
            LogLikelihood = 0
            for I in range(TupleCount):
                BX = MatrixMultiply(self.Beta, X[I])
                LogLikelihood += Y[I] * BX
                LogLikelihood -= math.log(1 + math.exp(BX))
            if PrevLogLikelihood != None:
                if VerboseFlag:
                    print "Log likelihood: %s (prev %s)"%(LogLikelihood, PrevLogLikelihood)
                Improvement = PrevLogLikelihood - LogLikelihood
                if Improvement < 0.001 and Improvement >= 0:
                    print "Reached optimum: Stop now!"
                    break
            PrevLogLikelihood = LogLikelihood
            # Compute the vector P:
            P = zeros(TupleCount, FloatType)
            for I in range(TupleCount):
                self.BetaSum = 0
                Tuple = AllTuples[I][1]
                for J in range(VectorSize):
                    self.BetaSum += self.Beta[J] * Tuple[J]
                Exp = math.exp(self.BetaSum)
                P[I] = Exp / (1.0 + Exp)
            # Compute the diagonal matrix W:
            W = zeros((TupleCount, TupleCount), FloatType)
            for I in range(TupleCount):
                W[I][I] = P[I] * (1.0 - P[I])
            try:
                WI = numpy.linalg.inv(W)
            except:
                traceback.print_exc()
                print "** Error: Unable to perform logistic regression due to singular matrix."
                print "Feature list was:", self.FeatureList
                return None
            # Compute the "response vector" z:
            z = MatrixMultiply(X, self.Beta)
            Diff = Y - P
            z += MatrixMultiply(WI, Diff)
            # Compute the new self.Beta:
            Product = MatrixMultiply(XT, W)
            Product = MatrixMultiply(Product, X)
            ProdI = numpy.linalg.inv(Product)
            Product = MatrixMultiply(ProdI, XT)
            Product = MatrixMultiply(Product, W)
            NewBeta = MatrixMultiply(Product, z)
            if VerboseFlag:
                print "Old self.Beta:", self.Beta
                print "New self.Beta:", NewBeta
            self.Beta = NewBeta
            CycleCount += 1
            if CycleCount >= 100:
                print "100 cycles performed; stopping now!"
                break
    def ScoreInstance(self, Features):
        FixedFeatures = [1.0,]
        for FeatureIndex in self.FeatureList:
            FixedFeatures.append(Features[FeatureIndex])
        BX = 0
        for I in range(len(FixedFeatures)):
            BX += self.Beta[I] * FixedFeatures[I]
        try:
            Exp = math.exp(BX)
        except:
            print "** exponent unreachable:", BX
            print "Features:", FixedFeatures
            raise ValueError, "Features out-of-range!"
        Score = Exp / (1.0 + Exp)
        return Score
    def Test(self, FeatureSet):
        for Vector in FeatureSet.AllVectors:
            FixedFeatures = [1.0,]
            for FeatureIndex in self.FeatureList:
                FixedFeatures.append(Vector.Features[FeatureIndex])
            BX = 0
            for I in range(len(FixedFeatures)):
                BX += self.Beta[I] * FixedFeatures[I]
            Exp = math.exp(BX)
            Vector.Score = Exp / (1.0 + Exp)
    def SaveModel(self, FileName):
        File = open(FileName, "wb")
        cPickle.dump("LOGIT", File)
        cPickle.dump(self.FeatureList, File)
        cPickle.dump(self.PValue, File)
        cPickle.dump(self.Beta, File)
        File.close()
    def LoadModel(self, FileName):
        File = open(FileName, "rb")
        cPickle.load(File) # model type
        self.FeatureList = cPickle.load(File)
        self.PValue = cPickle.load(File)
        self.Beta = cPickle.load(File)
        File.close()
   
def Test():
    pass

if __name__ == "__main__":
    # Command-line invocation: Test model loading/saving
    Test()