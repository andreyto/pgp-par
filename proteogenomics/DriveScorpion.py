"""
Code to test the new data-driven Scorpion infrastructure
"""

UsageInfo = """
DriveScorpion.py
This is the code for testing out new scoring models of Scorpion/Inspect.  It takes
a corpus of Inspect results and processes the Bayesian networks on them.  It requires
that the training files be in Inspect/TrainingCorpus
Required options for creating the cut features table:
 -c [Charge] - charge state of the peptides.  Default 2
 -m [SpectraDir] - Directory that holds the mzxml files
 -a [Action] - "table" is the action to make a table of mutual information for
    cut point features.  Other options include "train", "cuttable", "edges".
 -r [Radius] - Intensity radius for peak gathering
 -i [Flag] - Intensity scheme for peak gathering
 -s [Flag] - Sector type
 -e [Score] - Tag edge score multiplier
 -n: Toggle local to global noise model flag 
Additional options:
 -P - phosphorylation flag
"""

import os
import sys
import math
import struct
import random
random.seed(1)
import string
import getopt
import traceback
import PyInspect
import BasicStats
import ResultsParser
from Utils import *
Initialize()

QuickParseMaxLineNumber = 100

class NodeTypes: #Mirrored in IonScoring.h::PRMBayesianNodeType
    PrefixMass = 1
    SuffixMass = 2
    PrefixMassDoublyCharged = 3
    SuffixMassDoublyCharged = 4
    Sector = 5
    FlankAA = 6
    # Simple types: Take on 20 different values, simply reading off the aa
    PrefixAA = 7
    SuffixAA = 8
    PrefixContain = 9
    SuffixContain = 10
    PrefixContainPhos = 11
    SuffixContainPhos = 12

class FragmentTypes:
    Y = 1
    B = 2
    YLoss = 3
    BLoss = 4

class ScorpionDriver(ResultsParser.ResultsParser, ResultsParser.SpectrumOracleMixin):
    def __init__(self):
        self.Charge = 2 # default
        #self.CorpusDir = r"E:\ms\HeadToHead\InspectPV\At-root-3plates-b-2D34-010606-LTQ2-12.txt"
        self.CorpusDir = "TrainingCorpus"
        #self.CorpusDirIncludingFalseMatches = r"e:\ms\TrainingCorpusUnfiltered"
        self.CorpusDirIncludingFalseMatches = "TrainingCorpusUnfiltered"
        # defaults:
        self.IntensityScheme = 2 
        self.IntensityRadius = 0.3 
        self.SectorType = 1 
        self.SpectraDir = None
        self.PhosphorylationFlag = 0
        self.TagEdgeScoreMultiplier = 0
        self.Action = None
        self.NoiseModelFlag = 0
        self.MaxSpectraForPeptideSpecies = 5
        # (Annotation, Charge) -> count
        self.PeptidesSeen = {}
        self.QuickParseFlag = 0
        self.WriteFileName = None
        ResultsParser.ResultsParser.__init__(self)
        ResultsParser.SpectrumOracleMixin.__init__(self)
    def ComputeFeatureTableCut(self):
        self.OutputFile = open("FeatureTableCut.%s.txt"%self.Charge, "w")
        print ">>> Compute feature table for corpus..."
        # Print feature names:
        FeatureNames = PyInspect.GetBNFeatureNames()
        Str = "#"
        for FeatureName in FeatureNames:
            Str += "%s\t"%(FeatureName)
        self.OutputFile.write(Str + "\n")
        self.ProcessCorpus(self.ComputeFeatureTableCallback)
        self.OutputFile.close()
    def ComputeFeatureTable(self, FeatureFileName):
        """
        Called after ConstructPRMNetwork: Computes node values for all the spectra in the training set
        """
        self.OutputFile = open(FeatureFileName, "w")
        print ">>> Compute feature table for corpus..."
        # Print feature names:
        FeatureNames = PyInspect.GetBNFeatureNames()
        Str = ""
        for FeatureName in FeatureNames:
            Str += "%s\t"%(FeatureName)
        self.OutputFile.write(Str + "\n")
        self.ProcessCorpus(self.ComputeFeatureTableCallback)
        self.OutputFile.close()
    def ProcessCorpus(self, CallbackFunction, UnfilteredFlag = 0):
        # Reset peptide count:
        self.PeptidesSeen = {}
        if UnfilteredFlag:
            Dir = self.CorpusDirIncludingFalseMatches
        else:
            Dir = self.CorpusDir
        if not os.path.isdir(Dir):
            self.ProcessCorpusFile(Dir, CallbackFunction)
            return
        FileNames = os.listdir(Dir)
        FileNames.sort()
        for FileName in FileNames:
            FilePath = os.path.join(Dir, FileName)
            if os.path.isdir(FilePath):
                continue
            self.ProcessCorpusFile(FilePath, CallbackFunction)
    def ProcessCorpusFile(self, FilePath, CallbackFunction):
        print FilePath
        File = open(FilePath, "rb")
        LineNumber = 0 
        for FileLine in File.xreadlines():
            LineNumber += 1
            if self.QuickParseFlag and LineNumber >= QuickParseMaxLineNumber:
                break # Speedrun!
            if LineNumber % 100 == 0:
                print "%s line %s..."%(FilePath, LineNumber)
            if FileLine[0] == "#":
                continue
            FileLine = FileLine.strip()
            if not FileLine:
                continue
            Bits = FileLine.split("\t")
            try:
                Charge = int(Bits[self.Columns.Charge])
            except:
                print "* Line %s of file %s is invalid!"%(LineNumber, FilePath)
                traceback.print_exc()
                continue
            if Charge != self.Charge:
                continue
            if self.PhosphorylationFlag:
                #only allow phosphorylation spectra to pass
                Annotation = Bits[self.Columns.Annotation]
                if Annotation.find("phos") == -1:
                    continue
            Key = (Charge, Bits[self.Columns.Annotation])
            OldCount = self.PeptidesSeen.get(Key, 0)
            self.PeptidesSeen[Key] = OldCount + 1
            # Don't process too many occurrences of the same peptide species:
            if OldCount >= self.MaxSpectraForPeptideSpecies:
                continue
            CallbackFunction(Bits)
        File.close()
    def ComputeFeatureTableCallback(self, Bits):
        """
        Called by ComputeFeatureTable via ProcessCorpus: Write one line to the file for each
        PRM / cut, with the model features indicated.
        """
        SpectrumPath = self.FixSpectrumPath(Bits[0])
        Annotation = Bits[self.Columns.Annotation]
        SpectrumFilePos = int(Bits[self.Columns.FileOffset])        
        PySpectrum = PyInspect.Spectrum(SpectrumPath, SpectrumFilePos)
        PySpectrum.SetCharge(self.Charge)
        PySpectrum.PrepareIonScoring(1)
        NodeValueList = PyInspect.ComputeBNValuesForSpectrum(PySpectrum, Annotation)
        # Try stripping the FIRST and LAST values from the list.  Those two values
        # are the least interesting ones (the "cuts" corresponding to the endpoints of
        # the entire peptide).
        # A . V K E A M A P K . A
        #        x x x x
        MinIndex = 2
        MaxIndex = len(NodeValueList) - 1
        for ListIndex in range(MinIndex, MaxIndex):
            ValueList = NodeValueList[ListIndex]
            Str = ""
            for NodeValue in ValueList:
                Str += "%s\t"%NodeValue
            self.OutputFile.write(Str + "\n")
            self.OutputFile.flush()
    def TrainModelCallback(self, Bits):
        try:
            Annotation = Bits[self.Columns.Annotation]
            SpectrumFilePos = int(Bits[self.Columns.FileOffset])
        except:
            print "* Line %s of file %s is invalid!"%(LineNumber, FilePath)
            traceback.print_exc()
            return
        SpectrumPath = self.FixSpectrumPath(Bits[0])
        PySpectrum = PyInspect.Spectrum(SpectrumPath, SpectrumFilePos)
        PySpectrum.SetCharge(self.Charge)
        PySpectrum.PrepareIonScoring(1)
        PyInspect.TrainBNOnSpectrum(PySpectrum, Annotation)
    def TestPRMScoreRanking(self, PySpectrum, Peptide):
        """
        Compute PRM scores for true and false masses.  Check how well the PRMScores
        distinguish true from false PRMs by computing the ROC curve area.  Ideally
        it should be near 1.
        """
        FalseMasses = []
        ParentMass = Peptide.Masses[-1] + 19
        # Generate a list of random (spurious) PRMs to check:
        while len(FalseMasses) < 100:
            FalseMass = random.random() * ParentMass
            FalseMassOK = 1
            for Mass in Peptide.Masses:
                if abs(FalseMass - Mass) <= 0.3:
                    FalseMassOK = 0
                    break
            if FalseMassOK:
                FalseMasses.append(FalseMass)
        # For debugging: Plot spectrum scores:
        #PySpectrum.PlotPRMScores("PRMScores.%s.txt"%Peptide.Aminos, 1)
        # Build a list of the form (Score, ValidFlag):
        SortedScores = []
        #print "Compute TRUE mass PRM scores:"
        for Mass in Peptide.Masses[1:]:
            Score = PySpectrum.GetPRMScore(Mass, 1)
            SortedScores.append((Score, 1, Mass))
        #print "Compute FALSE mass PRM scores:"
        for Mass in FalseMasses:
            Score = PySpectrum.GetPRMScore(Mass, 1)
            SortedScores.append((Score, 0, Mass))
        SortedScores.sort()
        SortedScores.reverse()
        # Compute the ROC curve area:
        OverallPositiveCount = len(Peptide.Masses[1:])
        OverallNegativeCount = len(FalseMasses)
        PositiveCount = 0
        NegativeCount = 0
        ROCArea = 0
        TPRate = 0
        FPRate = 0
        #VerboseROCFileName = "VerboseROC.%s.txt"%Peptide.Aminos
        #VerboseROCFile = open(VerboseROCFileName, "wb")
        for (Score, ValidFlag, Mass) in SortedScores:
            if ValidFlag:
                PositiveCount += 1
                TPRate = PositiveCount / float(OverallPositiveCount)
            else:
                NegativeCount += 1
                FPRate = NegativeCount / float(OverallNegativeCount)
                ROCArea += TPRate / float(OverallNegativeCount)
            Str = "%s\t%s\t%s\t"%(Score, ValidFlag, Mass)
            Str += "%s\t%s\t%s\t%s\t"%(NegativeCount, PositiveCount, FPRate, TPRate)
            #VerboseROCFile.write(Str + "\n")
        #VerboseROCFile.close()
        #print "Wrote results to %s"%VerboseROCFileName
        print "Peptide %s -> ROC area %s"%(Peptide, ROCArea)
        self.PRMROCCount += 1
        self.PRMROCTotal += ROCArea
    def TestCutScoringCallback(self, Bits):
        SpectrumInfo = (Bits[0], Bits[1])
        if SpectrumInfo == self.OldSpectrumInfo:
            #print "(skip repeated annotation for %s:%s)"%(Bits[0], Bits[1])
            return
        try:
            Annotation = Bits[self.Columns.Annotation]
            Peptide = GetPeptideFromModdedName(Annotation)
            Charge = int(Bits[self.Columns.Charge])
            SpectrumFilePos = int(Bits[self.Columns.FileOffset])
        except:
            print "* Line %s of file %s is invalid!"%(LineNumber, FilePath)
            traceback.print_exc()
            return
        self.OldSpectrumInfo = SpectrumInfo
        SpectrumPath = self.FixSpectrumPath(Bits[0])
        PySpectrum = PyInspect.Spectrum(SpectrumPath, SpectrumFilePos)
        PySpectrum.SetCharge(self.Charge)
        PySpectrum.PrepareIonScoring(1)
        CutScores = PySpectrum.GetCutScores(Annotation)
        ValidFlag = 1
        if Bits[self.Columns.ProteinName][:3] == "XXX":
            ValidFlag = 0
        #print "%s:%s %s %s"%(SpectrumPath, SpectrumFilePos, ValidFlag, Annotation)
        #print "CutScores:"
        Str = ""
        for Score in CutScores:
            Str += "%.3f, "%Score
        #print Str
        #CutScores = PySpectrum.GetCutScores(Annotation)
        # Populate self.CutScoreResults:
        CutScoreCount = float(len(CutScores))
        Total = BasicStats.Sum(CutScores)
        Median = BasicStats.GetMedian(CutScores)
        self.CutScoreResults["total"].append((Total / CutScoreCount, ValidFlag))
        self.CutScoreResults["mean"].append((Total / CutScoreCount, ValidFlag))
        self.CutScoreResults["median"].append((Median, ValidFlag))
        Total = BasicStats.Sum(CutScores[1:-1])
        Median = BasicStats.GetMedian(CutScores[1:-1])
        self.CutScoreResults["total1"].append((Total / CutScoreCount, ValidFlag))
        self.CutScoreResults["mean1"].append((Total / CutScoreCount, ValidFlag))
        self.CutScoreResults["median1"].append((Median, ValidFlag))
        Total = BasicStats.Sum(CutScores[2:-2])
        Median = BasicStats.GetMedian(CutScores[2:-2])
        self.CutScoreResults["total2"].append((Total / CutScoreCount, ValidFlag))
        self.CutScoreResults["mean2"].append((Total / CutScoreCount, ValidFlag))
        self.CutScoreResults["median2"].append((Median, ValidFlag))
    def TestScoringCallback(self, Bits):
        try:
            Annotation = Bits[self.Columns.Annotation]
            Peptide = GetPeptideFromModdedName(Annotation)
            Charge = int(Bits[self.Columns.Charge])
            SpectrumFilePos = int(Bits[self.Columns.FileOffset])
        except:
            print "* Line %s of file %s is invalid!"%(LineNumber, FilePath)
            traceback.print_exc()
            return
        SpectrumPath = self.FixSpectrumPath(Bits[0])
        #print "Generate tags for %s\n  %s:%s"%(Annotation, SpectrumPath, SpectrumFilePos)
        PySpectrum = PyInspect.Spectrum(SpectrumPath, SpectrumFilePos)
        ParentMass = Peptide.Masses[-1] + 19
        PySpectrum.SetParentMass(ParentMass, self.Charge)
        PySpectrum.PrepareIonScoring(1)
        #PySpectrum.SetCharge(self.Charge)
        # Test PRM score ranking - rank true and false PRMs by score and report the
        # ROC curve area.
        self.TestPRMScoreRanking(PySpectrum, Peptide)
        #PySpectrum.PlotPRMScores("TestPRM.txt", 1)
        # Test TAG generation:
        TagList = PySpectrum.GenerateTags(Charge, 1, self.TagEdgeScoreMultiplier)
        Index = self.CheckTags(self.TaggingHistogram, TagList, Peptide)
        # Temp - report all cases where the top tag didn't win.
        if Index != 0:
            TagInfoStr = "%s\t%s\t%s\t"%(SpectrumPath, Bits[1], SpectrumFilePos)
            TagInfoStr += "%s\t%s\t%s\t"%(self.Charge, Annotation, Index)
            self.TagResultsFile.write(TagInfoStr + "\n")
        # OLD tag generator:
        ##TagList = PySpectrum.GenerateTags(Charge, 0, self.TagEdgeScoreMultiplier)
        ##self.CheckTags(self.OldTaggingHistogram, TagList, Peptide)
    def TrainAndTestCuts(self):
        """
        Main method for training and testing a CUT scoring model.  We evaluate its
        performance by determing how well cut-scores differentiate between valid and
        invalid matches.
        """
        self.ModelFileName = "TestCut.%s.dat"%self.Charge
        ################################################################
        # Iterate over annotations. Train the model.        
        self.ProcessCorpus(self.TrainModelCallback)
        ################################################################            
        # Generate probability tables, and save model to disk.
        PyInspect.ComputeBNProbabilityTables()
        PyInspect.SaveBNModel(self.ModelFileName)
        PyInspect.DebugPrintBNModel()
        ################################################################
        # Scoring of cut points: Determine the ROC curve if we sort matches
        # by their mean cut score, median cut score, total cut score,
        # mean excluding edges, median excluding edges, etc.  These
        # lists all have entries of the form (Score, ValidFlag)
        self.CutScoreResults = {}
        self.CutScoreResults["mean"] = []
        self.CutScoreResults["median"] = []
        self.CutScoreResults["total"] = []
        self.CutScoreResults["mean1"] = []
        self.CutScoreResults["median1"] = []
        self.CutScoreResults["total1"] = []
        self.CutScoreResults["mean2"] = []
        self.CutScoreResults["median2"] = []
        self.CutScoreResults["total2"] = []
        self.OldSpectrumInfo = None
        self.ProcessCorpus(self.TestCutScoringCallback, 1)
        self.ReportCutTestingResults()
    def ReportCutTestingResults(self):
        self.OutputFileCutTesting = open("CutTestingResults.%s.txt"%(self.Charge), "wb")
        VerboseOutputFile = open("CutTestingResultsVerbose.%s.txt"%(self.Charge), "wb")
        for Key in self.CutScoreResults.keys():
            VerboseOutputFile.write("\n")
            List = self.CutScoreResults[Key]
            List.sort()
            List.reverse()
            ValidCount = 0
            InvalidCount = 0
            for (Score, ValidFlag) in List:
                if ValidFlag:
                    ValidCount += 1
                else:
                    InvalidCount += 1
            ROCArea = 0
            CumulativeTP = 0
            CumulativeFP = 0
            TPRate = 0
            FPRate = 0
            for (Score, ValidFlag) in List:
                if ValidFlag:
                    CumulativeTP += 1
                    TPRate = CumulativeTP / float(ValidCount)
                else:
                    CumulativeFP += 1
                    FPRate = CumulativeFP / float(InvalidCount)
                    ROCArea += TPRate / float(InvalidCount)
                VerboseOutputFile.write("%s\t%s\t%s\t%s\t%s\t\n"%(Key, CumulativeFP, CumulativeTP, FPRate, TPRate))
            Str = "%s\t%s\t"%(Key, ROCArea)
            self.OutputFileCutTesting.write(Str + "\n")
        VerboseOutputFile.close()                                    
    def TrainAndTestPRM(self, TagFlag):
        """
        Call this AFTER the node topology has been set.  We'll train the PRM scorer on a corpus of valid
        annotations.  Then we'll test tagging on the corpus.
        """
        ###########################################
        # Test setting of node values:
        PySpectrum = PyInspect.Spectrum(r"SystemTest\TestSpectrum.dta")
        Peptide = GetPeptideFromModdedName("VKEAMAPK")
        PySpectrum.SetParentMass(Peptide.Masses[-1] + 19)
        PySpectrum.PrepareIonScoring(1)
        for Mass in Peptide.Masses:
            PySpectrum.GetPRMScore(Mass, 1)
        PyInspect.TrainBNOnSpectrum(PySpectrum, "VKEAMAPK")
        ###########################################
        if self.WriteFileName:
            self.ModelFileName = self.WriteFileName
        else:
            if TagFlag:
                self.ModelFileName = "TestBN.Tag.%s.dat"%self.Charge
            else:
                self.ModelFileName = "TestBN.PRM.%s.dat"%self.Charge
        ################################################################
        # Iterate over annotations. Train the model.
        self.ProcessCorpus(self.TrainModelCallback)
        ################################################################            
        # Generate probability tables, and save model to disk.
        PyInspect.ComputeBNProbabilityTables()
        print "SAVE model..."
        PyInspect.SaveBNModel(self.ModelFileName)
        PyInspect.DebugPrintBNModel()
        #print "LOAD model..."
        #PyInspect.LoadBNModel(self.ModelFileName)
        #PyInspect.DebugPrintBNModel()
        # Set up score structures for use by TestScoringCallback:
        self.TaggingHistogram = {}
        self.OldTaggingHistogram = {}
        self.PRMROCTotal = 0
        self.PRMROCCount = 0
        self.TagResultsFile = open("TaggingFailureFile.txt", "wb")
        Header = "#SpectrumFile\tScanNumber\tFilePos\tCharge\tAnnotation\tTrueTagIndex\t"
        self.TagResultsFile.write(Header + "\n")
        ################################################################
        # Iterate over annotations.  Test scoring.
        self.ProcessCorpus(self.TestScoringCallback)
        ################################################################
        # Report results:
        self.FinalOutputFile = open("FinalOutput~%s~%s~%s~%s~%s~%s.txt"%(self.NoiseModelFlag, self.SectorType, self.IntensityScheme, self.IntensityRadius, self.Charge, self.TagEdgeScoreMultiplier), "wb")
        print "\n\n"
        ROC = self.PRMROCTotal / float(max(1, self.PRMROCCount))
        self.FinalOutputFile.write("%s\n"%ROC)
        print "Average ROC area over %s spectra: %s"%(self.PRMROCCount, ROC)
        print "\n\n"
        print "-=- "*4
        print "Tagging results with NEW scorer:"
        self.ReportTaggingResults(self.TaggingHistogram)
    def ReportTaggingResults(self, TaggingHistogram):
        print ">>>PRM scoring results:"
        TotalTaggingTries = TaggingHistogram[None]
        Count = TaggingHistogram.get(0, 0)
        Percent1 = 100 * Count / float(max(1, TotalTaggingTries))
        print " Top tag: %d of %d (%.2f%%)"%(Count, TotalTaggingTries, Percent1)

        Count = 0
        for Index in range(10):
            Count += TaggingHistogram.get(Index, 0)
        Percent10 = 100 * Count / float(max(1, TotalTaggingTries))
        print " Top 10 tags: %d of %d (%.2f%%)"%(Count, TotalTaggingTries, Percent10)

        Count = 0
        for Index in range(25):
            Count += TaggingHistogram.get(Index, 0)
        Percent25 = 100 * Count / float(max(1, TotalTaggingTries))
        print " Top 25 tags: %d of %d (%.2f%%)"%(Count, TotalTaggingTries, Percent25)

        Count = 0
        for Index in range(100):
            Count += TaggingHistogram.get(Index, 0)
        Percent100 = 100 * Count / float(max(1, TotalTaggingTries))
        print " Top 100 tags: %d of %d (%.2f%%)"%(Count, TotalTaggingTries, Percent100)
        
        Score = Percent1 + Percent10 + Percent25 + Percent100
        print "Overall score:", Score
        self.FinalOutputFile.write("%s\t%s\t%s\t%s\t%s\t\n"%(Percent1, Percent10, Percent25, Percent100, Score))
        return Score
    def CheckTags(self, TaggingHistogram, TagList, Peptide, VerboseFlag = 0):
        TaggingHistogram[None] = TaggingHistogram.get(None, 0) + 1
        MassEpsilon = 3.0
        Aminos = Peptide.Aminos.replace("I", "L").replace("Q", "K")
        for TagIndex in range(len(TagList)):
            TagTuple = TagList[TagIndex]
            if VerboseFlag:
                print "Check tag %s for peptide %s: %s"%(TagIndex, Peptide.Aminos, TagTuple)
            (Prefix, Tag, Suffix, Score) = TagTuple[:4]
            Tag = Tag.replace("I", "L").replace("Q", "K")
            for AminoIndex in range(len(Peptide.Masses) - 2):
                if Aminos[AminoIndex:AminoIndex + 3] == Tag and abs(Peptide.Masses[AminoIndex] - Prefix) < MassEpsilon:
                    # Found a valid tag!
                    TaggingHistogram[TagIndex] = TaggingHistogram.get(TagIndex, 0) + 1
                    if VerboseFlag:
                        print ">>> Tag ok!"
                    return TagIndex
        return None
    def TrainTagEdgeScoring(self):
        # Initialize score histograms:
        self.TagHistogramSkew = {}
        self.TagHistogramAbsSkew = {}
        self.TagHistogramAbsTotalSkew = {}
        self.TrueTagHistogramSkew = {}
        self.TrueTagHistogramAbsSkew = {}
        self.TrueTagHistogramAbsTotalSkew = {}
        self.SkewHistogramMultiplier = 50.0
        # Populate score histograms:
        self.ProcessCorpus(self.TrainTagEdgeScoringOnSpectrum)
        # Report statistics for each histogram:
        self.ReportTagEdgeHistogram("TagHistogramTotal.txt", self.TagHistogramSkew, self.TrueTagHistogramSkew)
        self.ReportTagEdgeHistogram("TagHistogramAbs.txt", self.TagHistogramAbsSkew, self.TrueTagHistogramAbsSkew)
        self.ReportTagEdgeHistogram("TagHistogramAbsTotal.txt", self.TagHistogramAbsTotalSkew, self.TrueTagHistogramAbsTotalSkew)
        self.SaveTagEdgeScores("TagSkewScores.dat")
    def SaveTagEdgeScores(self, OutputFileName):
        File = open(OutputFileName, "wb")
        Keys = self.TagHistogramAbsSkew.keys()
        BinCount = 30
        Str = struct.pack("<i", BinCount)
        File.write(Str)
        Scores = self.ComputeEdgeSkewCumulativeHistogram(BinCount, self.TagHistogramAbsTotalSkew, self.TrueTagHistogramAbsTotalSkew)
        for Index in range(BinCount):
            Str = struct.pack("<f", Scores[Index])
            File.write(Str)
        Scores = self.ComputeEdgeSkewCumulativeHistogram(BinCount, self.TagHistogramAbsSkew, self.TrueTagHistogramAbsSkew)
        for Index in range(BinCount):
            Str = struct.pack("<f", Scores[Index])
            File.write(Str)
        File.close()
    def ComputeEdgeSkewCumulativeHistogram(self, BinCount, Histogram, TrueHistogram):
        CumulativeCount = [0] * BinCount
        CumulativeCountTrue = [0] * BinCount
        CountAll = 0
        CountAllTrue = 0
        for (Key, Count) in Histogram.items():
            for Index in range(min(Key, BinCount), BinCount):
                CumulativeCount[Index] += Count
            CountAll += Count
        for (Key, Count) in TrueHistogram.items():
            for Index in range(min(Key, BinCount), BinCount):
                CumulativeCountTrue[Index] += Count
            CountAllTrue += Count
        ScoreList = []
        for Index in range(BinCount):
            OddsTrue = CumulativeCountTrue[Index] / float(CountAllTrue)
            OddsFalse = CumulativeCount[Index] / float(CountAll)
            ScoreList.append(math.log(max(0.001, OddsTrue)) - math.log(max(0.001, OddsFalse)))
            print "%s, %s, %s, %s"%(Index, OddsTrue, OddsFalse, ScoreList[-1])
        return ScoreList
    def TrainTagEdgeScoringOnSpectrum(self, Bits):
        try:
            Annotation = Bits[self.Columns.Annotation]
            Peptide = GetPeptideFromModdedName(Annotation)
            Charge = int(Bits[self.Columns.Charge])
            SpectrumFilePos = int(Bits[self.Columns.FileOffset])
        except:
            print "* Line is invalid!"
            print Bits
            traceback.print_exc()
            return
        SpectrumPath = self.FixSpectrumPath(Bits[0])
        print "Generate tags for %s\n  %s:%s"%(Annotation, SpectrumPath, SpectrumFilePos)
        PySpectrum = PyInspect.Spectrum(SpectrumPath, SpectrumFilePos)
        PySpectrum.SetCharge(self.Charge)
        PySpectrum.PrepareIonScoring(1)
        MassEpsilon = 2.5
        # use zero weight on edge-scoring here, since the edge scores are what we wish to train:
        TagList = PySpectrum.GenerateTags(Charge, 1, 0.0)
        # GenerateTags returns tuples of the form (Prefix, Tag, Suffix, ..., TotalSkew, TotalAbsSkew)
        for TagTuple in TagList:
            (Prefix, Tag, Suffix) = TagTuple[:3]
            (TotalSkew, TotalAbsSkew) = TagTuple[-2:]
            Aminos = Peptide.Aminos.replace("I", "L").replace("Q", "K")
            Tag = Tag.replace("I", "L").replace("Q", "K")
            TagCorrectFlag = 0
            for AminoIndex in range(len(Peptide.Masses) - 2):
                if Aminos[AminoIndex:AminoIndex + 3] == Tag and abs(Peptide.Masses[AminoIndex] - Prefix) < MassEpsilon:
                    # Found a valid tag!
                    TagCorrectFlag = 1
            # Record ABS SKEW in histogram:
            Bin = int(round(TotalAbsSkew / self.SkewHistogramMultiplier))
            if TagCorrectFlag:
                self.TrueTagHistogramAbsSkew[Bin] = self.TrueTagHistogramAbsSkew.get(Bin, 0) + 1
            else:
                self.TagHistogramAbsSkew[Bin] = self.TagHistogramAbsSkew.get(Bin, 0) + 1
            # Return ABS(TOTAL SKEW) in histogram:
            Bin = int(round(abs(TotalSkew) / self.SkewHistogramMultiplier))
            if TagCorrectFlag:
                self.TrueTagHistogramAbsTotalSkew[Bin] = self.TrueTagHistogramAbsTotalSkew.get(Bin, 0) + 1
            else:
                self.TagHistogramAbsTotalSkew[Bin] = self.TagHistogramAbsTotalSkew.get(Bin, 0) + 1
            # Return TOTAL SKEW in histogram:
            Bin = int(round(TotalSkew / self.SkewHistogramMultiplier))
            if TagCorrectFlag:
                self.TrueTagHistogramSkew[Bin] = self.TrueTagHistogramSkew.get(Bin, 0) + 1
            else:
                self.TagHistogramSkew[Bin] = self.TagHistogramSkew.get(Bin, 0) + 1
                
    def ReportTagEdgeHistogram(self, FileName, HistoA, HistoB):
        """
        Helper for TrainTagEdgeScoringOnSpectrum.  How likely is a total skew for
        a true tag versus a false tag?
        """
        File = open(FileName, "wb")
        AllBins = {}
        TotalCountA = 0
        for (Key, Count) in HistoA.items():
            AllBins[Key] = 1
            TotalCountA += Count
        TotalCountB = 0
        for (Key, Count) in HistoB.items():
            AllBins[Key] = 1
            TotalCountB += Count
        Keys = AllBins.keys()
        Keys.sort()
        for Key in Keys:
            Mass = Key * self.SkewHistogramMultiplier * 0.001
            Str = "%s\t%s\t"%(Key, Mass)
            Str += "%s\t%s\t"%(HistoA.get(Key, 0), HistoB.get(Key, 0))
            Str += "%s\t"%Mass
            Str += "%s\t%s\t"%(HistoA.get(Key, 0) / float(TotalCountA), HistoB.get(Key, 0) / float(TotalCountB))
            File.write(Str + "\n")
        File.close()
    def Main(self):
        print ">>> Build PRM scoring network:"
        if self.Action == "prmtable":
            self.ConstructPRMNetwork(0)
            self.ComputeFeatureTable("FeatureTable.%s.txt"%self.Charge)
        elif self.Action == "prmtrain":
            self.ConstructPRMNetwork(0)
            print ">>> Train and test PRM scoring network:"
            self.TrainAndTestPRM(0)
        elif self.Action == "tagtable":    
            self.ConstructPRMNetwork(1)
            self.ComputeFeatureTable("FeatureTable.Tag.%s.txt"%self.Charge)
        elif self.Action == "tagtrain":
            self.ConstructPRMNetwork(1)
            print ">>>Train and test scoring network:"
            self.TrainAndTestPRM(1)
        elif self.Action == "cuttable":
            self.ConstructCutNetwork()
            print ">>> Compute tables to get mutual information, etc. for cut features:"
            self.ComputeFeatureTableCut()
        elif self.Action == "cuttrain":
            self.ConstructCutNetwork()
            print ">>> Train and test cut scoring network:"
            self.TrainAndTestCuts()
        elif self.Action == "edges":
            self.ConstructPRMNetwork(1)
            print ">>> Training for tag edge scoring:"
            self.TrainTagEdgeScoring()
        else:
            print "Action Unknown"
            print UsageInfo
            sys.exit(-1)
    def ConstructCutNetwork(self):
        if self.PhosphorylationFlag:
            self.ConstructPhosphoCutNetwork()                
        elif self.Charge == 2:
            pass
            #self.ConstructCutNetworkCharge2() #for my qstar debugging
        else:
            self.ConstructCutNetworkCharge3()
    def ConstructCutNetworkCharge2(self):
        print "Intensity scheme %s, intensity radius %s"%(self.IntensityScheme, self.IntensityRadius)
        PyInspect.ResetIonScoring(self.IntensityScheme, self.IntensityRadius, 1, 1)
        # Non-peak nodes:
        print "Sector type:", self.SectorType
        IndexSector = PyInspect.AddIonScoringNode("Sector", NodeTypes.Sector, self.SectorType)
        # Flanking amino acids:
        for (Name, TypeFlag) in [("Prefix", NodeTypes.PrefixAA),
                             ("Suffix", NodeTypes.SuffixAA)]:
            PyInspect.AddIonScoringNode("%sA"%Name, TypeFlag, 0)
            PyInspect.AddIonScoringNode("%sC"%Name, TypeFlag, 2)
            PyInspect.AddIonScoringNode("%sD"%Name, TypeFlag, 3)
            PyInspect.AddIonScoringNode("%sE"%Name, TypeFlag, 4)
            PyInspect.AddIonScoringNode("%sF"%Name, TypeFlag, 5)
            PyInspect.AddIonScoringNode("%sG"%Name, TypeFlag, 6)
            PyInspect.AddIonScoringNode("%sH"%Name, TypeFlag, 7)
            PyInspect.AddIonScoringNode("%sI"%Name, TypeFlag, 8)
            PyInspect.AddIonScoringNode("%sK"%Name, TypeFlag, 10)
            PyInspect.AddIonScoringNode("%sL"%Name, TypeFlag, 11)
            PyInspect.AddIonScoringNode("%sM"%Name, TypeFlag, 12)
            PyInspect.AddIonScoringNode("%sN"%Name, TypeFlag, 13)
            PyInspect.AddIonScoringNode("%sP"%Name, TypeFlag, 15)
            PyInspect.AddIonScoringNode("%sQ"%Name, TypeFlag, 16)
            PyInspect.AddIonScoringNode("%sR"%Name, TypeFlag, 17)
            PyInspect.AddIonScoringNode("%sS"%Name, TypeFlag, 18)
            PyInspect.AddIonScoringNode("%sT"%Name, TypeFlag, 19)
            PyInspect.AddIonScoringNode("%sV"%Name, TypeFlag, 21)
            PyInspect.AddIonScoringNode("%sW"%Name, TypeFlag, 22)
            PyInspect.AddIonScoringNode("%sY"%Name, TypeFlag, 24)
        # Prefix/suffix containing residues:
        PyInspect.AddIonScoringNode("PrefixAcidFlag", NodeTypes.PrefixContain, 0)
        PyInspect.AddIonScoringNode("PrefixAcid", NodeTypes.PrefixContain, 1)
        PyInspect.AddIonScoringNode("PrefixBaseFlag", NodeTypes.PrefixContain, 2)
        PyInspect.AddIonScoringNode("PrefixBase", NodeTypes.PrefixContain, 3)
        PyInspect.AddIonScoringNode("SuffixAcidFlag", NodeTypes.SuffixContain, 0)
        PyInspect.AddIonScoringNode("SuffixAcid", NodeTypes.SuffixContain, 1)
        PyInspect.AddIonScoringNode("SuffixBaseFlag", NodeTypes.SuffixContain, 2)
        PyInspect.AddIonScoringNode("SuffixBase", NodeTypes.SuffixContain, 3)        
        # Flank (encoded):
        PyInspect.AddIonScoringNode("FlankB", NodeTypes.FlankAA, 0)
        PyInspect.AddIonScoringNode("FlankY", NodeTypes.FlankAA, 1)
        # And now, the peaks:
        IndexY = PyInspect.AddIonScoringNode("Y", NodeTypes.SuffixMass, 0, 0.0)
        IndexYI = PyInspect.AddIonScoringNode("YI", NodeTypes.SuffixMass, 0, 1.0)
        IndexYII = PyInspect.AddIonScoringNode("YII", NodeTypes.SuffixMass, 0, 2.0)
        IndexYH2O = PyInspect.AddIonScoringNode("Y-H2O", NodeTypes.SuffixMass, 0, -18.0)
        IndexYNH3 = PyInspect.AddIonScoringNode("Y-NH3", NodeTypes.SuffixMass, 0, -17.0)
        #
        IndexB = PyInspect.AddIonScoringNode("B", NodeTypes.PrefixMass, 0, 1.0)
        IndexBI = PyInspect.AddIonScoringNode("BI", NodeTypes.PrefixMass, 0, 2.0)
        IndexBII = PyInspect.AddIonScoringNode("BII", NodeTypes.SuffixMass, 0, 3.0)
        IndexBH2O = PyInspect.AddIonScoringNode("B-H2O", NodeTypes.PrefixMass, 0, -17.0)
        IndexBNH3 = PyInspect.AddIonScoringNode("B-NH3", NodeTypes.PrefixMass, 0, -16.0)
        IndexA = PyInspect.AddIonScoringNode("A", NodeTypes.PrefixMass, 0, -27.0)
        #
        IndexY2 = PyInspect.AddIonScoringNode("Y2", NodeTypes.SuffixMassDoublyCharged, 0, 0.0)
        IndexY2I = PyInspect.AddIonScoringNode("Y2I", NodeTypes.SuffixMassDoublyCharged, 0, 1.0)
        IndexY2H2O = PyInspect.AddIonScoringNode("Y2-H2O", NodeTypes.SuffixMassDoublyCharged, 0, -18.0)
        IndexY2NH3 = PyInspect.AddIonScoringNode("Y2-NH3", NodeTypes.SuffixMassDoublyCharged, 0, -17.0)
        #
        IndexB2 = PyInspect.AddIonScoringNode("B2", NodeTypes.PrefixMassDoublyCharged, 0, 1.0)
        IndexB2I = PyInspect.AddIonScoringNode("B2I", NodeTypes.PrefixMassDoublyCharged, 0, 2.0)
        IndexB2H2O = PyInspect.AddIonScoringNode("B2-H2O", NodeTypes.PrefixMassDoublyCharged, 0, -17.0)
        IndexB2NH3 = PyInspect.AddIonScoringNode("B2-NH3", NodeTypes.PrefixMassDoublyCharged, 0, -16.0)
        #
        PyInspect.FinishIonScoringNetwork()
    def ConstructPhosphoCutNetwork(self):
        #self.SectorType = 3 #have 5 zones
        #self.IntensityScheme = 4 #based on average grass
        ## NOTE: ions are listed here in groups, with Y first, then B then, double charged.
        ## within the group we list by most common first.
        print "Phosphorylation Cut Network, charge 2"
        print "Intensity scheme %s, intensity radius %s"%(self.IntensityScheme, self.IntensityRadius)
        CutFlag = 1
        PyInspect.ResetIonScoring(self.IntensityScheme, self.IntensityRadius, CutFlag, self.NoiseModelFlag)
        # Non-peak nodes:
        print "Sector type:", self.SectorType
        IndexSector = PyInspect.AddIonScoringNode("Sector", NodeTypes.Sector, self.SectorType)
        IndexPrefixContainPhos = PyInspect.AddIonScoringNode("PrefixContainPhos",NodeTypes.PrefixContainPhos,0)
        IndexSuffixContainPhos = PyInspect.AddIonScoringNode("SuffixContainPhos",NodeTypes.SuffixContainPhos,0)
        IndexFlankB = PyInspect.AddIonScoringNode("FlankB", NodeTypes.FlankAA, 0)
        IndexFlankY = PyInspect.AddIonScoringNode("FlankY", NodeTypes.FlankAA, 1)            
        # And now, the peaks:
        IndexY = PyInspect.AddIonScoringNode("Y", NodeTypes.SuffixMass, 0, 0.0, FragmentTypes.Y)
        IndexYI = PyInspect.AddIonScoringNode("YI", NodeTypes.SuffixMass, 0, 1.0, FragmentTypes.Y)
        IndexYII = PyInspect.AddIonScoringNode("YII", NodeTypes.SuffixMass, 0, 2.0, FragmentTypes.Y)
        IndexYH2O = PyInspect.AddIonScoringNode("Y-H2O", NodeTypes.SuffixMass, 0, -18.0, FragmentTypes.YLoss)
        IndexYNH3 = PyInspect.AddIonScoringNode("Y-NH3", NodeTypes.SuffixMass, 0, -17.0, FragmentTypes.YLoss)
        IndexYP = PyInspect.AddIonScoringNode("Y-P", NodeTypes.SuffixMass, 0, -98.0, FragmentTypes.YLoss)
        #
        IndexB = PyInspect.AddIonScoringNode("B", NodeTypes.PrefixMass, 0, 1.0, FragmentTypes.B)
        IndexBI = PyInspect.AddIonScoringNode("BI", NodeTypes.PrefixMass, 0, 2.0, FragmentTypes.B)
        #IndexBII = PyInspect.AddIonScoringNode("BII", NodeTypes.SuffixMass, 0, 3.0, FragmentTypes.B) #don't include charge 2,3
        IndexBP = PyInspect.AddIonScoringNode("B-P", NodeTypes.PrefixMass, 0, -97.0, FragmentTypes.BLoss)
        IndexBH2O = PyInspect.AddIonScoringNode("B-H2O", NodeTypes.PrefixMass, 0, -17.0, FragmentTypes.BLoss)
        IndexBNH3 = PyInspect.AddIonScoringNode("B-NH3", NodeTypes.PrefixMass, 0, -16.0, FragmentTypes.BLoss)
        if self.Charge == 2:
            IndexBPH2O = PyInspect.AddIonScoringNode("B-P-H2O",NodeTypes.PrefixMass, 0, -115.0, FragmentTypes.BLoss)
            IndexBPNH3 = PyInspect.AddIonScoringNode("B-P-NH3",NodeTypes.PrefixMass, 0, -114.0, FragmentTypes.BLoss)
        #
        IndexY2 = PyInspect.AddIonScoringNode("Y2", NodeTypes.SuffixMassDoublyCharged, 0, 0.0, FragmentTypes.Y)
        IndexY2I = PyInspect.AddIonScoringNode("Y2I", NodeTypes.SuffixMassDoublyCharged, 0, 1.0, FragmentTypes.Y)
        IndexY2H2O = PyInspect.AddIonScoringNode("Y2-H2O", NodeTypes.SuffixMassDoublyCharged, 0, -18.0, FragmentTypes.YLoss)
        IndexY2NH3 = PyInspect.AddIonScoringNode("Y2-NH3", NodeTypes.SuffixMassDoublyCharged, 0, -17.0, FragmentTypes.YLoss)
        IndexY2P = PyInspect.AddIonScoringNode("Y2-P", NodeTypes.SuffixMassDoublyCharged, 0, -98.0, FragmentTypes.YLoss)
        
        if not self.Charge == 2:#don't include B2 anything for charge 2
            IndexB2 = PyInspect.AddIonScoringNode("B2", NodeTypes.PrefixMassDoublyCharged, 0, 1.0, FragmentTypes.B)
            IndexB2I = PyInspect.AddIonScoringNode("B2I", NodeTypes.PrefixMassDoublyCharged, 0, 2.0, FragmentTypes.B)
            IndexB2H2O = PyInspect.AddIonScoringNode("B2-H2O", NodeTypes.PrefixMassDoublyCharged, 0, -17.0, FragmentTypes.BLoss)
            IndexB2NH3 = PyInspect.AddIonScoringNode("B2-NH3", NodeTypes.PrefixMassDoublyCharged, 0, -16.0, FragmentTypes.BLoss)
            IndexB2P = PyInspect.AddIonScoringNode("B2-P", NodeTypes.PrefixMassDoublyCharged, 0, -97.0, FragmentTypes.BLoss)
        
        ################### Add parentage ###############################
        PyInspect.SetIonScoringNodeParents(IndexY, [IndexSector,IndexFlankY])
        PyInspect.SetIonScoringNodeParents(IndexYI, [IndexY, IndexSector])
        PyInspect.SetIonScoringNodeParents(IndexYII, [IndexSector, IndexYI])
        PyInspect.SetIonScoringNodeParents(IndexYH2O, [IndexSector, IndexY])
        PyInspect.SetIonScoringNodeParents(IndexYNH3, [IndexSector, IndexYH2O]) 
        PyInspect.SetIonScoringNodeParents(IndexYP, [IndexY, IndexSuffixContainPhos])
        #
        PyInspect.SetIonScoringNodeParents(IndexB, [IndexSector, IndexY, IndexFlankB])
        PyInspect.SetIonScoringNodeParents(IndexBI, [IndexB, IndexSector])
        #PyInspect.SetIonScoringNodeParents(IndexBII, [IndexB, IndexBI]) #don't include
        if self.Charge == 2:
            PyInspect.SetIonScoringNodeParents(IndexBP, [IndexY, IndexPrefixContainPhos]) # yes that's y
        else:
            PyInspect.SetIonScoringNodeParents(IndexBP, [IndexSector, IndexPrefixContainPhos])             
        PyInspect.SetIonScoringNodeParents(IndexBH2O, [IndexSector, IndexB])
        PyInspect.SetIonScoringNodeParents(IndexBNH3, [IndexSector, IndexBH2O])
        if self.Charge == 2:
            PyInspect.SetIonScoringNodeParents(IndexBPH2O, [IndexBP,IndexSector])
            PyInspect.SetIonScoringNodeParents(IndexBPNH3, [IndexBP,IndexSector])
        #
        PyInspect.SetIonScoringNodeParents(IndexY2, [IndexSector, IndexFlankY])
        PyInspect.SetIonScoringNodeParents(IndexY2I, [IndexSector, IndexY2])
        PyInspect.SetIonScoringNodeParents(IndexY2H2O, [IndexSector, IndexY2])
        PyInspect.SetIonScoringNodeParents(IndexY2NH3, [IndexSector, IndexY2])
        if self.Charge == 2:
            PyInspect.SetIonScoringNodeParents(IndexY2P, [IndexSector, IndexY2])
        else:
            PyInspect.SetIonScoringNodeParents(IndexY2P, [IndexSuffixContainPhos, IndexY2])
            
        #
        if not self.Charge == 2:        # don't include B2 for charge 2
            PyInspect.SetIonScoringNodeParents(IndexB2, [IndexSector, IndexY])
            PyInspect.SetIonScoringNodeParents(IndexB2I, [IndexB2, IndexSector])
            PyInspect.SetIonScoringNodeParents(IndexB2H2O, [IndexSector, IndexB2])
            PyInspect.SetIonScoringNodeParents(IndexB2NH3, [IndexSector, IndexB2H2O])
            PyInspect.SetIonScoringNodeParents(IndexB2P, [IndexB2, IndexSector])
        #
        PyInspect.FinishIonScoringNetwork()
        
    def ConstructPRMNetwork(self, TagFlag = 1):
        # Initialize the scoring BN to use intensity threshold scheme 0, and
        # peak-grabbing radius 0.3Da:
        print "Intensity scheme %s, intensity radius %s"%(self.IntensityScheme, self.IntensityRadius)
        PyInspect.ResetIonScoring(self.IntensityScheme, self.IntensityRadius, 0, self.NoiseModelFlag)
        # Add nodes:
        # Non-peak nodes:
        print "Sector type:", self.SectorType
        IndexSector = PyInspect.AddIonScoringNode("Sector", NodeTypes.Sector, self.SectorType)
        if TagFlag:
            IndexFlankB = PyInspect.AddIonScoringNode("FlankB", NodeTypes.FlankAA, 0)
            IndexFlankY = PyInspect.AddIonScoringNode("FlankY", NodeTypes.FlankAA, 1)            
        # Suffix peaks:
        IndexY = PyInspect.AddIonScoringNode("Y", NodeTypes.SuffixMass, 0, 0.0, FragmentTypes.Y)
        IndexYI = PyInspect.AddIonScoringNode("YI", NodeTypes.SuffixMass, 0, 1.0, FragmentTypes.Y)
        IndexYII = PyInspect.AddIonScoringNode("YII", NodeTypes.SuffixMass, 0, 2.0, FragmentTypes.Y)
        IndexYH2O = PyInspect.AddIonScoringNode("Y-H2O", NodeTypes.SuffixMass, 0, -18.0, FragmentTypes.YLoss)
        IndexYNH3 = PyInspect.AddIonScoringNode("Y-NH3", NodeTypes.SuffixMass, 0, -17.0, FragmentTypes.YLoss)
        # Prefix peaks:
        IndexB = PyInspect.AddIonScoringNode("B", NodeTypes.PrefixMass, 0, 1.0, FragmentTypes.B)
        IndexBI = PyInspect.AddIonScoringNode("BI", NodeTypes.PrefixMass, 0, 2.0, FragmentTypes.B)
        #IndexBII = PyInspect.AddIonScoringNode("BII", NodeTypes.SuffixMass, 0, 3.0) # too rare
        IndexBH2O = PyInspect.AddIonScoringNode("B-H2O", NodeTypes.PrefixMass, 0, -17.0, FragmentTypes.BLoss)
        IndexBNH3 = PyInspect.AddIonScoringNode("B-NH3", NodeTypes.PrefixMass, 0, -16.0, FragmentTypes.BLoss)
        IndexA = PyInspect.AddIonScoringNode("A", NodeTypes.PrefixMass, 0, -27.0, FragmentTypes.BLoss)
        # Doubly-charged peaks:
        IndexY2 = PyInspect.AddIonScoringNode("Y2", NodeTypes.SuffixMassDoublyCharged, 0, 0.0, FragmentTypes.Y)
        IndexY2I = PyInspect.AddIonScoringNode("Y2I", NodeTypes.SuffixMassDoublyCharged, 0, 1.0, FragmentTypes.Y)
        IndexY2H2O = PyInspect.AddIonScoringNode("Y2-H2O", NodeTypes.SuffixMassDoublyCharged, 0, -18.0, FragmentTypes.YLoss)
        IndexY2NH3 = PyInspect.AddIonScoringNode("Y2-NH3", NodeTypes.SuffixMassDoublyCharged, 0, -17.0, FragmentTypes.YLoss)
        #
        IndexB2 = PyInspect.AddIonScoringNode("B2", NodeTypes.PrefixMassDoublyCharged, 0, 1.0, FragmentTypes.B)
        IndexB2I = PyInspect.AddIonScoringNode("B2I", NodeTypes.PrefixMassDoublyCharged, 0, 2.0, FragmentTypes.B)
        IndexB2H2O = PyInspect.AddIonScoringNode("B2-H2O", NodeTypes.PrefixMassDoublyCharged, 0, -17.0, FragmentTypes.BLoss)
        IndexB2NH3 = PyInspect.AddIonScoringNode("B2-NH3", NodeTypes.PrefixMassDoublyCharged, 0, -16.0, FragmentTypes.BLoss)
        # Set node parents.  (This needn't be called for nodes without parents)
        if TagFlag:
            PyInspect.SetIonScoringNodeParents(IndexY, [IndexSector,IndexFlankY])
        else:
            PyInspect.SetIonScoringNodeParents(IndexY, [IndexSector,])
        PyInspect.SetIonScoringNodeParents(IndexYI, [IndexY])
        PyInspect.SetIonScoringNodeParents(IndexYII, [IndexY, IndexYI])
        PyInspect.SetIonScoringNodeParents(IndexYH2O, [IndexY,])
        PyInspect.SetIonScoringNodeParents(IndexYNH3, [IndexY, IndexYH2O])
        #
        if TagFlag:
            PyInspect.SetIonScoringNodeParents(IndexB, [IndexSector, IndexY, IndexFlankB])
        else:
            PyInspect.SetIonScoringNodeParents(IndexB, [IndexSector, IndexY])        
        PyInspect.SetIonScoringNodeParents(IndexBI, [IndexB,])
        PyInspect.SetIonScoringNodeParents(IndexA, [IndexB])
        PyInspect.SetIonScoringNodeParents(IndexBH2O, [IndexB])
        PyInspect.SetIonScoringNodeParents(IndexBNH3, [IndexB, IndexBH2O])
        #
        if TagFlag:
            PyInspect.SetIonScoringNodeParents(IndexY2, [IndexSector, IndexB, IndexFlankY])
        else:
            PyInspect.SetIonScoringNodeParents(IndexY2, [IndexSector, IndexB,])
        PyInspect.SetIonScoringNodeParents(IndexY2I, [IndexY2])
        PyInspect.SetIonScoringNodeParents(IndexY2H2O, [IndexSector, IndexY2])
        PyInspect.SetIonScoringNodeParents(IndexY2NH3, [IndexSector, IndexY2H2O])
        #
        if TagFlag:
            PyInspect.SetIonScoringNodeParents(IndexB2, [IndexSector, IndexY, IndexFlankB])
        else:
            PyInspect.SetIonScoringNodeParents(IndexB2, [IndexSector, IndexY,])
        PyInspect.SetIonScoringNodeParents(IndexB2I, [IndexB2])
        PyInspect.SetIonScoringNodeParents(IndexB2H2O, [IndexSector, IndexB2])
        PyInspect.SetIonScoringNodeParents(IndexB2NH3, [IndexSector, IndexB2H2O])
        PyInspect.FinishIonScoringNetwork()
    def ParseCommandLine(self, Arguments):
        (Options, Args) = getopt.getopt(Arguments, "a:c:Qr:s:i:m:Pe:nw:I:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-a":
                self.Action = Value.lower()
            if Option == "-c":
                self.Charge = int(Value)
            if Option == "-m":
                self.SpectraDir = os.path.abspath(Value)
                self.PopulateSpectrumOracle(self.SpectraDir)
            if Option == "-P":
                self.PhosphorylationFlag = 1
            if Option == "-r":
                self.IntensityRadius = float(Value)
            if Option == "-i":
                self.IntensityScheme = int(Value)
            if Option == "-s":
                self.SectorType = int(Value)
            if Option == "-Q":
                self.QuickParseFlag = 1
            if Option == "-e":
                self.TagEdgeScoreMultiplier = float(Value)
            if Option == "-n":
                self.NoiseModelFlag = 1
            if Option == "-w":
                self.WriteFileName = Value
            if Option == "-I":
                if Value == "QTOF":
                    self.InstrumentType = 1 #look in Utils.h for the c #defines
                elif Value == "ESI-ION-TRAP":
                    self.InstrumentType = 0
                elif Value == "FT-HYBRID":
                    self.InstrumentType = 2

if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(psyco not available)"
    try:
        Driver = ScorpionDriver()
        Driver.ParseCommandLine(sys.argv[1:])
        Driver.Main()
    except:
        traceback.print_exc()
        print ">>> Press enter <<<"
        sys.stdin.readline()
    