"""
Plan:
Output a large collection of features for each post-translational modification accepted on a
search of a part-bogus database.  All modifications on the bogus proteins are incorrect.
An equivalent number of modifications on the non-bogus proteins are incorrect.  Let's compute
a variety of features for the PTMs observed.

Input:
A collection of annotated spectra, output by SelectSites.py
Output:
A file listing all the observed modification sites, with various features computed.  

Then, we train a model to distinguish between good (correct DB) and bad (incorrect DB)
modifications.  Model types: LDA, logistic regression, SVM, etc.
(Another experiment: Search unmodified spectra against a mutated database,
judge correct precisely those modifications which "undo" the mutations)
"""

UsageInfo = """
TrainPTMFeatures: Train model on PTM features from on a data-set.
Run this AFTER running ComputePTMFeatures.

Arguments:
 -m [model-type]: Train a model and report its accuracy on the specified
    (-u) file
 -u [FeatureFile]: Path to the feature-file written out by ComputePTMFeatures
 -v [FeatureFile]: Write scored features out to the specified file
 -w [FileName]: Write model to the specified file.   (Set either -w OR -r)
 -r [FileName]: Read a model from the specified file (Set either -w OR -r)

Optional:
 -f [Flag]: Perform feature selection (1 for accumulate, 2 for prune)
 -e [TestingFile]: Path to the feature-file that serves as a testing set.
    If not specified, then the same features (-u) will be used for testing.
    For use with -f flag only.
 -R [Path]: Report ROC curve to the specified file
"""
import os
import sys
import struct
import traceback
import getopt
import MSSpectrum
import PyInspect
import random
random.seed(1)
import shutil
import Learning
import BasicStats
import ResultsParser
import SpectralSimilarity
from Utils import *
Initialize()
try:
    from numpy import *
    import numpy.linalg
    FloatType = float
    MatrixMultiply = dot
except:
    print "** Warning: Unable to import Numpy.  Logit training not available"

ValidFeatureIndices = [2,3,5,22,23,24,26]

class FeatureBits:
    SpectrumCount = 0
    ModlessSpectrumCount = 1
    BestMQScore = 2
    BestDeltaScore = 3
    PeptideCount = 4
    ConsensusMQScore = 5
    PeptideLength = 6
    TotalCutScore = 7
    MedianCutScore = 8
    YPresent = 9
    BPresent = 10
    BYIntensity = 11
    NTT = 12
    ModdedFraction = 13
    SpectraThisModType = 15
    SitesThisModType = 16
    Dot = 18
    Shared01 = 19
    Shared11 = 20
    Correlation = 21
    LogSpectrumCount = 22
    LogPeptideLength = 23
    LogSpecThisType = 24
    LogSitesThisType = 25
    DeltaVsBigDB = 26
    
class FormatBits:
    DBPos = 1
    ModificationMass = 2
    ModifiedAA = 3
    ProteinName = 4
    ModifiedResidueNumber = 5
    Peptide = 6
    Charge = 7
    TrueProteinFlag = 8
    SisterAnnotationFlag = 9
    BestSpectrumPath = 10
    BestModlessSpectrumPath = 11
    BestModlessMQScore = 12
    BigDBAnnotation = 13
    BigDBMQScore = 14
    SpectrumCount = 15
    ModlessSpectrumCount = 16
    BestMQScore = 17
    BestDeltaScore = 18
    PeptideCount = 19
    ConsensusMQScore = 20
    NTT = 27
    ModdedFraction = 28
    SpectraWithThisModType = 30
    SitesWithThisModType = 31
    LogSpectrumCount = 37
    LogSpectraThisModType = 39
    LogSitesThisModType = 40
    ConsensusDeltaBigDB = 41
    FirstFeature = 15
    LastFeature = 41
    FeatureCount = LastFeature - FirstFeature + 1
    ModelScore = 42 # score for the PEPTIDE SPECIES
    ModelPValue = 43 # p-value (probability false given this score) for the PEPTIDE SPECIES
    SitePValue = 44 # p-value (probability false given several species) for the SITE 
    KnownPTMName = 45
    KnownPTMAnnotation = 46
    KnownPTMScore = 47
    KnownPTMSitePValue = 48
    
class PTMFeatureTrainer(ResultsParser.ResultsParser):
    def __init__(self):
        self.ResultsFileName = None
        self.DBPath = None
        self.OutputPath = "PTMFeatures.txt"
        self.TempFileDir = "PTMFeatures"
        self.PTMs = {} # keys of the form (DBPos, Mass)
        self.CoverageThreshold = 2 # at least this many spectra to consider a residue 'covered'.
        self.QuickParseFlag = 0 # if true, then parse only the first n lines
        self.SpectrumDir = None
        self.SuperSpectrumDir = None
        self.PoolFlag = 0
        self.ModelType = None
        self.SisterProteins = {} # protein index -> sister protein's index
        self.InputFeaturePath = None
        self.ModelTestFilePath = None
        # Dictionary of unmodified peptides, for computing the coverage level:
        self.UnmodifiedPeptides = {}
        self.FeatureSelectionFlag = None
        self.CachedProteinNames = []
        self.CachedFilePaths = []
        self.CachedFixedFilePaths = []
        self.StartOutputDBPos = 0
        self.HeaderLines = []
        self.ReportROCPath = None
        self.OutputFeaturePath = None
        self.ReadModelFilePath2 = None
        self.ReadModelFilePath3 = None
        self.TrainingSetDBRatio = 1.0
        self.TestingSetDBRatio = 1.0
        self.WriteModelFilePath2 = None
        self.WriteModelFilePath3 = None
        ResultsParser.ResultsParser.__init__(self)
    def TrainFacultative(self):
        """
        Train paired models for CONSTITUTIVE ("always") and FACULTATIVE ("sometimes") PTMs.
        """
        # Train a model on all PTMs, to get initial scores for all PTMs.
        # The initial model uses only CONSTITUTIVE features, and its output
        # is used only to provide an ORACLE for the facultative model:
        print "TRAIN model on all features:"
        self.Model.Train(self.TrainingSetAll)
        print "SCORE all features:"
        self.Model.Test(self.TrainingSetAll)
        ##############################################################
        print "Generate SUB-MODEL of only facultative features:"
        # Sort facultative instances by score:
        SortedList = []
        for Vector in self.TrainingSetAll.AllVectors:
            if not Vector.FileBits[FormatBits.SisterAnnotationFlag]:
                continue
            SortedList.append((Vector.Score, Vector))
        SortedList.sort()
        FacFeatureSet = Learning.FeatureSetClass()
        ChunkSize = min(len(SortedList) / 4, 1000)
        print "Sorted list of %s facultative features, chunk size is %s"%(len(SortedList), ChunkSize)
        for (Score, Vector) in SortedList[:ChunkSize]:
            NewVector = Learning.FeatureVector()
            NewVector.FileBits = Vector.FileBits[:]
            NewVector.Features = Vector.Features[:]
            NewVector.TrueFlag = 0
            FacFeatureSet.AllVectors.append(NewVector)
            FacFeatureSet.FalseVectors.append(NewVector)
        for (Score, Vector) in SortedList[-ChunkSize:]:
            NewVector = Learning.FeatureVector()
            NewVector.FileBits = Vector.FileBits[:]
            NewVector.Features = Vector.Features[:]
            NewVector.TrueFlag = 1
            FacFeatureSet.AllVectors.append(NewVector)
            FacFeatureSet.TrueVectors.append(NewVector)
        FacFeatureSet.SetCounts()
        FacFeatureSet.GetPriorProbabilityFalse(self.TrainingSetDBRatio)
        ##############################################################
        # Write out the FACULTATIVE feature set:
        FacTrainingFile = open("FacultativeTrainingSet.txt", "wb")
        for HeaderLine in self.HeaderLines:
            FacTrainingFile.write(HeaderLine)
        for Vector in FacFeatureSet.AllVectors:
            Bits = Vector.FileBits[:]
            if Vector.TrueFlag:
                Bits[FormatBits.TrueProteinFlag] = "1"
            else:
                Bits[FormatBits.TrueProteinFlag] = "0"
            Str = string.join(Bits, "\t")
            FacTrainingFile.write(Str + "\n")
        FacTrainingFile.close()
        ##############################################################
        # Train the sub-model:
        self.FacModel = self.GetModelObject(self.FeaturesF)
        self.FacModel.Train(FacFeatureSet)
        self.FacModel.Test(FacFeatureSet)
        self.FacModel.ReportAccuracy(FacFeatureSet) # invokes ComputeOddsTrue
        ##############################################################
        # Apply the trained fac-model to *all* facultative features, and
        # train an overall model on all *constitutive* features:
        self.FeatureSetC = Learning.FeatureSetClass()
        self.FeatureSetF = Learning.FeatureSetClass()
        for Vector in self.TrainingSetAll.AllVectors:
            if Vector.FileBits[FormatBits.SisterAnnotationFlag]:
                FeatureSet = self.FeatureSetF
            else:
                FeatureSet = self.FeatureSetC
            FeatureSet.AllVectors.append(Vector)
            if Vector.TrueFlag:
                FeatureSet.TrueVectors.append(Vector)
            else:
                FeatureSet.FalseVectors.append(Vector)
        self.FeatureSetC.SetCounts()
        self.FeatureSetF.SetCounts()
        self.FeatureSetC.GetPriorProbabilityFalse(self.TrainingSetDBRatio)
        self.FeatureSetF.GetPriorProbabilityFalse(self.TrainingSetDBRatio)
        # Score facultative-feature, using facultative-model:
        self.FacModel.Test(self.FeatureSetF)
        # Train constitutive-ONLY model, and score constitutive features:
        self.ConModel = self.GetModelObject(self.FeaturesC)
        self.ConModel.Train(self.FeatureSetC)
        self.ConModel.Test(self.FeatureSetC)
        self.ConModel.ReportAccuracy(self.FeatureSetC) # to invoke ComputeOddsTrue
        ##############################################################
        # Save our models:
        if self.WriteModelFilePath:
            (Stub, Extension) = os.path.splitext(self.WriteModelFilePath)
            ConModelPath = "%s.con"%Stub
            FacModelPath = "%s.fac"%Stub
            self.ConModel.SaveModel(ConModelPath)
            self.FacModel.SaveModel(FacModelPath)
        ##############################################################
        # Write out the scored features:
        OutputFile = open(self.OutputFeaturePath, "wb")
        for Line in self.HeaderLines:
            OutputFile.write(Line)
        for Vector in self.TrainingSetAll.AllVectors:
            if Vector.FileBits[FormatBits.SisterAnnotationFlag]:
                PValue = self.FacModel.GetPValue(Vector.Score)
            else:
                PValue = self.ConModel.GetPValue(Vector.Score)
            while len(Vector.FileBits) <= FormatBits.ModelPValue:
                Vector.FileBits.append("")
            Vector.FileBits[FormatBits.ModelScore] = str(Vector.Score)
            Vector.FileBits[FormatBits.ModelPValue] = str(PValue)
            Str = string.join(Vector.FileBits, "\t")
            OutputFile.write(Str + "\n")
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
    def TrainModel(self):
        """
        Our training data-set is in self.InputFeaturePath.
        Let's train a model to predict which entries come from the true database.
        """
        if not self.InputFeaturePath:
            print "* Please specify an input feature-file."
            print UsageInfo
            sys.exit(-1)
        # Load in features for a collection of TRUE and FALSE instances.
        File = open(self.InputFeaturePath, "rb")
        self.FeatureNames = {}
        FeatureCount = FormatBits.LastFeature - FormatBits.FirstFeature + 1
        # We have one set of features for facultative sites, and one for constitutive.
        # Note that some features (modification rate, correlation with unmodified peptide)
        # are applicable to F but not C.
        #self.FeaturesF = range(FeatureCount)
        # For constitutive modifications: Modification rate, protein coverage,
        # and number of unmodified peptides are all off-limits.  (Those features
        # are "dead giveaways" that we have a non-shuffled protein!)
        #self.FeaturesC = [2, 3, 5, 22, 24, 25, 26]
        self.FeaturesC = ValidFeatureIndices[:]
        #self.FeaturesC = range(FeatureCount)
        self.FeaturesF = self.FeaturesC
        self.FeaturesAll = []
        for FeatureIndex in self.FeaturesF:
            if FeatureIndex in self.FeaturesC:
                self.FeaturesAll.append(FeatureIndex)
        # We can OVERRIDE the list of features here, to forbid the use of some:
        print "Permitted features all:", self.FeaturesAll
        # Parse the features from the TRAINING and TESTING files.  We generate
        # training sets for the FACULTATIVE (F) and for CONSTITUTIVE (C) sites.
        self.TrainingSet2 = Learning.FeatureSetClass()
        self.TrainingSet2.Type = "Charge-2"
        self.TrainingSet3 = Learning.FeatureSetClass()
        self.TrainingSet3.Type = "Charge-3"
        #self.TrainingSetAll = Learning.FeatureSetClass()
        #self.TrainingSetAll.Type = "All"
        self.ParseFeatureFile(self.InputFeaturePath, self.TrainingSet2, self.TrainingSet3,
                              self.TrainingSetDBRatio)
        if self.ModelTestFilePath:
            self.TestingSet2 = FeatureSetClass()
            self.TestingSet3 = FeatureSetClass()
            self.ParseFeatureFile(self.ModelTestFilePath, self.TestingSet2, self.TestingSet3,
                self.TestingSetAll, self.TestingSetDBRatio)
        # SPECIAL values for model, which don't actually cause training:
        if self.ModelType == "feature":
            print "\n\nSINGLE feature:"
            self.TrainOneFeature(self.TrainingSet2)
            self.TrainOneFeature(self.TrainingSet3)
            return
        if self.ModelType == "featurescatter":
            print "\n\nFeature+feature scatter-plots:"
            self.ProduceFeatureScatterPlots(self.TrainingSetAll)
            return
        if self.ModelType == "summary":
            self.PerformFeatureSummary()
            return
        # Instantiate our model:
        self.Model2 = self.GetModelObject(self.FeaturesAll)
        self.Model3 = self.GetModelObject(self.FeaturesAll)
        # Load a pre-trained model, if we received a path:
        if self.ReadModelFilePath2:
            self.Model2.LoadModel(self.ReadModelFilePath2)
            self.Model3.LoadModel(self.ReadModelFilePath3)
        #######################################################################
        # Special value for feature selection (3) means that we train a model on
        # all data, then use it to generate a sub-feature-set for a facultative model!
        if self.FeatureSelectionFlag == 3:
            self.TrainFacultative()
            return
        #######################################################################
        # If we're not doing feature selection: Train on the training set,
        # and then (if we have a testing set) test on the testing set.
        if not self.FeatureSelectionFlag:
            # Train the model (unless we just loaded it in):
            if not self.ReadModelFilePath2:
                self.Model2.Train(self.TrainingSet2)
                self.Model3.Train(self.TrainingSet3)
            # Compute the score of each vector:
            if self.ModelTestFilePath:
                self.Model2.Test(self.TestingSet2)
                self.Model2.ReportAccuracy(self.TestingSet2)
                self.Model3.Test(self.TestingSet3)
                self.Model3.ReportAccuracy(self.TestingSet3)
                self.WriteScoredFeatureSet(self.TestingSet2, self.TestingSet3)
            else:
                self.Model2.Test(self.TrainingSet2)
                self.Model2.ReportAccuracy(self.TrainingSet2)
                shutil.copyfile("Pvalues.txt", "PValues.chg2.txt")
                self.Model3.Test(self.TrainingSet3)
                self.Model3.ReportAccuracy(self.TrainingSet3)
                shutil.copyfile("Pvalues.txt", "PValues.chg3.txt")
                #if self.ReportROCPath:
                #    self.Model.ReportROC(self.TrainingSetAll, self.ReportROCPath)
                self.WriteScoredFeatureSet(self.TrainingSet2, self.TrainingSet3)
            if self.WriteModelFilePath2:
                self.Model2.SaveModel(self.WriteModelFilePath2)
                self.Model3.SaveModel(self.WriteModelFilePath3)
            return
        #######################################################################
        # We're doing feature selection.  We'll need to write out feature files,
        # then call TrainMachineLearner
        print "Feature names:", self.FeatureNames
        print "AllFeatures:", self.FeaturesAll
        self.WriteFeaturesToFile(self.TrainingSet2, "PTMFeatures.2.txt")
        self.WriteFeaturesToFile(self.TrainingSet3, "PTMFeatures.3.txt")
        # *** Additive and subtractive aren't done here, the user can do it!
    def WriteFeaturesToFile(self, TrainingSet, FileName):
        print "Write features to %s..."%FileName
        File = open(FileName, "wb")
        File.write("#Index\tValidFlag\t")
        for Key in self.FeaturesAll:
            File.write("%s\t"%self.FeatureNames[Key])
        File.write("\n")
        TrainingSet.SaveTabDelimited(File)
        File.close()
    def ProduceFeatureScatterPlots(self, FeatureSet):
        """
        Iterate over all pairs of (distinct) features.  For each pair, produce a scatter-plot
        with N true points and N false points.
        """
        OutputFile = open("FeatureScatterPlots.txt", "wb")
        VectorCount = 200
        TrueVectors = FeatureSet.TrueVectors[:]
        random.shuffle(TrueVectors)
        TrueVectors = TrueVectors[:VectorCount]
        FalseVectors = FeatureSet.FalseVectors[:]
        random.shuffle(FalseVectors)
        FalseVectors = FalseVectors[:VectorCount]
        # Write a HEADER:
        HeaderStr = ""
        for FeatureIndex in range(len(self.FeaturesAll)):
            Feature = self.FeaturesAll[FeatureIndex]
            HeaderStr += "T %s\tF %s\t"%(self.FeatureNames[Feature], self.FeatureNames[Feature])
        OutputFile.write(HeaderStr + "\n")
        # Write one row for each pair of vectors:
        for RowIndex in range(len(TrueVectors)):
            Str = ""
            TrueVector = TrueVectors[RowIndex]
            FalseVector = FalseVectors[RowIndex]
            for Feature in self.FeaturesAll:
                Str += "%s\t%s\t"%(TrueVector.Features[Feature], FalseVector.Features[Feature])
            OutputFile.write(Str + "\n")
        return
    def WriteScoredFeatureSet(self, FeatureSet2, FeatureSet3):
        # Write out the features with their model-scores:
        if not self.OutputFeaturePath:
            return
        File = open(self.OutputFeaturePath, "wb")
        for FileLine in self.HeaderLines:
            File.write(FileLine)
        SortedVectors = []
        for Vector in FeatureSet2.AllVectors:
            SortedVectors.append((int(Vector.FileBits[1]), Vector.FileBits[6], int(Vector.FileBits[7]), Vector))
        for Vector in FeatureSet3.AllVectors:
            SortedVectors.append((int(Vector.FileBits[1]), Vector.FileBits[6], int(Vector.FileBits[7]), Vector))
        SortedVectors.sort()
        for Tuple in SortedVectors:
            Vector = Tuple[-1]
            Charge = int(Tuple[2])
            if Charge > 2:
                Model = self.Model3
            else:
                Model = self.Model2
            Bits = Vector.FileBits
            while len(Bits) <= FormatBits.ModelPValue:
                Bits.append("")
            Bits[FormatBits.ModelScore] = str(Vector.Score)
            Bits[FormatBits.ModelPValue] = str(self.Model2.GetPValue(Vector.Score))
            Str = string.join(Bits, "\t")
            File.write(Str + "\n")                
        File.close()
        return
        # Iterate over all vectors, write them all out:
        for Vector in FeatureSet2.AllVectors:
            Bits = Vector.FileBits
            while len(Bits) <= FormatBits.ModelPValue:
                Bits.append("")
            Bits[FormatBits.ModelScore] = str(Vector.Score)
            Bits[FormatBits.ModelPValue] = str(self.Model2.GetPValue(Vector.Score))
            Str = string.join(Bits, "\t")
            File.write(Str + "\n")
        # Iterate over all vectors, write them all out:
        for Vector in FeatureSet3.AllVectors:
            Bits = Vector.FileBits
            while len(Bits) <= FormatBits.ModelPValue:
                Bits.append("")
            Bits[FormatBits.ModelScore] = str(Vector.Score)
            Bits[FormatBits.ModelPValue] = str(self.Model3.GetPValue(Vector.Score))
            Str = string.join(Bits, "\t")
            File.write(Str + "\n")
        File.close()
    def ParseFeatureFile(self, FilePath, FeatureSet2, FeatureSet3, DBRatio):
        """
        Initialize the FeatureSet object, by parsing features from the specified FilePath.
        Facultative features go to FeatureSetF, constitutive features go to FeatureSetC
        """
        File = open(FilePath, "rb")
        # Parse the header line specially:
        HeaderLine = File.readline()
        self.HeaderLines.append(HeaderLine)
        Bits = HeaderLine.strip().split("\t")
        for BitIndex in range(len(Bits)):
            if BitIndex >= FormatBits.FirstFeature:
                self.FeatureNames[BitIndex - FormatBits.FirstFeature] = Bits[BitIndex]
                if BitIndex <= FormatBits.LastFeature:
                    print "Feature %s: %s"%(BitIndex - FormatBits.FirstFeature, Bits[BitIndex])
        # Iterate over non-header lines:
        LineNumber = 0
        for FileLine in File.xreadlines():
            LineNumber += 1
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
            except:
                continue # skip; not a valid instance line
            Charge = int(Bits[FormatBits.Charge])
            SisterAnnotation = Bits[FormatBits.SisterAnnotationFlag]
            Vector = Learning.FeatureVector()
            if Charge > 2:
                FeatureSet = FeatureSet3
            else:
                FeatureSet = FeatureSet2
            try:
                for FeatureBitIndex in range(FormatBits.FirstFeature, FormatBits.LastFeature + 1):
                    FeatureIndex = FeatureBitIndex - FormatBits.FirstFeature
                    #if FeatureIndex not in self.FeaturesAll:
                    #    continue
                    if FeatureBitIndex < len(Bits) and Bits[FeatureBitIndex].strip() and Bits[FeatureBitIndex] != "None":
                        Vector.Features.append(float(Bits[FeatureBitIndex]))
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
        for FeatureSet in (FeatureSet2, FeatureSet3):
            FeatureSet.SetCounts()
            FeatureSet.GetPriorProbabilityFalse(DBRatio)
        print "CHARGE 1,2: Read in %s true and %s false vectors"%(FeatureSet2.TrueCount, FeatureSet2.FalseCount)
        print "CHARGE  3+: Read in %s true and %s false vectors"%(FeatureSet3.TrueCount, FeatureSet3.FalseCount)
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
            print "  At %.1f%% FDRate (cutoff %.5f), got %s PTMs (sensitivity %.2f%%)"%(Cutoffs[CutIndex] * 100, Thresholds[CutIndex],
                BestCounts[CutIndex], Sensitivity)
            print "  ->True sensitivity: %.4f%%"%(100 * BestCounts[CutIndex] / float(max(1, AllTrueCount - AllFalseCount)))
        print "False positive rate amoung top %s sites: %s"%(TopCount, 100*TopCountFalse/float(max(1, TopCount)))
        print "Overall, %s true and %s false features."%(TrueCount, FalseCount)
        print "ROC curve area: %.5f"%ROCArea
        # The 'score' we return is a tuple giving the best accuracy at several cutoffs:
        return (BestCounts[2], BestCounts[0], BestCounts[4], BestCounts[3], BestCounts[2])
    def PerformFeatureSummary(self):
        for FeatureIndex in range(len(self.Features)):
            TrueList = []
            for Tuple in self.TrueTuples:
                TrueList.append(Tuple[FeatureIndex])
            TrueList.sort()
            (TMean, TStdDev) = BasicStats.GetMeanStdDev(TrueList)
            FalseList = []
            for Tuple in self.FalseTuples:
                FalseList.append(Tuple[FeatureIndex])
            FalseList.sort()
            (FMean, FStdDev) = BasicStats.GetMeanStdDev(FalseList)
            print "Feature %s (%s):"%(FeatureIndex, self.FeatureNames[FeatureIndex])
            print "  True: Mean %.4f, stddev %.4f (range %.4f..%.4f)"%(TMean, TStdDev, TrueList[0], TrueList[-1])
            print "  False: Mean %.4f, stddev %.4f (range %.4f..%.4f)"%(FMean, FStdDev, FalseList[0], FalseList[-1])
    def TrainOneFeature(self, TrainingSet):
        """
        Compute accuracy for a very simple-minded model:
        Rank sites by the value of a SINGLE FEATURE (descending order)
        """
        for FeatureIndex in range(FormatBits.FeatureCount):
            SortedList = []
            for Vector in TrainingSet.TrueVectors:
                SortedList.append((Vector.Features[FeatureIndex], random.random(), 1))
            for Vector in TrainingSet.FalseVectors:
                SortedList.append((Vector.Features[FeatureIndex], random.random(), 0))
            # And report the accuracy of this lonely feature:
            print
            print "Feature %s (%s):"%(FeatureIndex, self.FeatureNames[FeatureIndex])
            self.ReportAccuracy(SortedList)
    def ParseCommandLine(self, Arguments):
        (Options, Args) = getopt.getopt(Arguments, "m:u:v:r:w:f:e:R:D:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-m":
                self.ModelType = Value.lower()
            elif Option == "-D":
                self.TrainingSetDBRatio = float(Value)
            elif Option == "-r":
                if not os.path.exists(Value):
                    print "** Error: Model file '%s' not found for reading.\n"%Value
                    return 0
                self.ReadModelFilePath2 = "%s.2"%Value
                self.ReadModelFilePath3 = "%s.3"%Value
            elif Option == "-w":
                #self.WriteModelFilePath = Value
                self.WriteModelFilePath2 = "%s.2"%Value
                self.WriteModelFilePath3 = "%s.3"%Value
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
            else:
                print "* Error: Unrecognized option %s"%Option
                return 0
        return 1 # success

if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(psyco not installed; running unoptimized)"
    Trainer = PTMFeatureTrainer()
    Result = Trainer.ParseCommandLine(sys.argv[1:])
    if not Result:
        sys.exit(-1)
    if Trainer.ModelType:
        Trainer.TrainModel()
        sys.exit()
    print UsageInfo
    sys.exit(-1)
    
    