"""
Invoke parent mass correction code programatically.
"""
import sys
import random
import os
import getopt
import traceback
import Learning
import ResultsParser
import time
import PyInspect
from Utils import *
Initialize()

QuickParseLineCount = 100
ParentMassPPMLimit = 2000

class PMCFeatureGrabber(ResultsParser.ResultsParser):
    """
    Driver for the parent mass correction features of PyInspect.  This class
    is responsible for deriving a training set, and for testing the results
    of the trained model.
    """
    def __init__(self):
        ResultsParser.ResultsParser.__init__(self)
        self.QuickParseFlag = 0
        self.ComputeFeaturesFlag = 0
        self.SVMFlag = 0
        self.TrainFlag = 1
        #self.TrainingCorpusDir = r"e:\ms\TrainingCorpus"
        self.TrainingCorpusDir = "TrainingCorpus"
    def GetModelFileName(self, Charge):
        if self.SVMFlag:
            return "PMCSVM%s"%Charge
        else:
            return "PMCLDA%s"%Charge
    def GetOverallScore(self):
        """
        Score our model, based on its overall success.
        """
        Score = 0
        Score += (self.Success1[1] - self.StartOK1[1]) / float(self.Attempts[1])
        Score += (self.Success1[2] - self.StartOK1[2]) / float(self.Attempts[2])
        Score += (self.Success1[3] - self.StartOK1[3]) / float(self.Attempts[3])
        Score += (self.Success3[1] - self.StartOK3[1]) / float(self.Attempts[1])
        Score += (self.Success3[2] - self.StartOK3[2]) / float(self.Attempts[2])
        Score += (self.Success3[3] - self.StartOK3[3]) / float(self.Attempts[3])
        return Score
    def ComputeAndReportSuccessRate(self):
        """
        Perform parent mass correction.  Report the fraction of the time that we're within 0.1,
        within 0.5, within 0.1 allowing runners-up, within 0.5 allowing runners-up.
        """
        self.Attempts = {}
        self.StartOK1 = {}
        self.StartOK3 = {}
        self.StartOK5 = {}
        self.Success1 = {}
        self.Success3 = {}
        self.Success5 = {}
        self.Success1RU = {}
        self.Success3RU = {}
        self.Success5RU = {}
        PyInspect.ReloadPMC()
        self.PMCOutputFile = open("PMCOutput.txt", "wb")
        Header = "#FileName\tByteOffset\tTruePeptide\tCharge\tTrueParentMass\tFileMass\tCorrectedMass\tOriginalDiff\tFixedDiff\tRunnerUpDiff\t\n"
        self.PMCOutputFile.write(Header)
        self.PMCFailuresFile = open("PMCFailures.txt", "wb")
        Header = "Original\tTrue\tPMCorrected\tCharge\tFile\tPos\tTruePeptide\t"
        self.PMCFailuresFile.write(Header + "\n")
        for OracleFileName in os.listdir(self.TrainingCorpusDir):
            OraclePath = os.path.join(self.TrainingCorpusDir, OracleFileName)
            if os.path.isdir(OraclePath):
                continue
            self.MeasurePMCSuccessRate(OraclePath)
        AllDicts = (("StartOK1", self.StartOK1),
                    ("StartOK3", self.StartOK3),
                    ("StartOK5", self.StartOK5),
                    ("Success1", self.Success1),
                    ("Success3", self.Success3),
                    ("Success5", self.Success5),
                    ("Success1RU", self.Success1RU),
                    ("Success3RU", self.Success3RU),
                    ("Success5RU", self.Success5RU))
        for (Name, Dict) in AllDicts:
            Str = "%s\t"%Name
            for Charge in (1, 2, 3):
                Successes = Dict.get(Charge, 0)
                Attempts = self.Attempts.get(Charge, 0)
                if Attempts:
                    Rate = Successes / float(Attempts)
                    Str += "%.3f\t"%(Rate*100)
                else:
                    Str += "\t"
            print Str
    def MeasurePMCSuccessRate(self, OracleFileName):
        """
        Helper for ComputeAndReportSuccessRate: Measures successes on a single file.
        The "oracle file" is an inspect output file, filtered to a good false discovery
        rate.
        """
        print OracleFileName
        File = open(OracleFileName, "rb")
        LineNumber = 0
        for FileLine in File.xreadlines():
            Bits = FileLine.strip().split("\t")
            LineNumber += 1
            if LineNumber > QuickParseLineCount and self.QuickParseFlag:
                break 
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
            if Charge < 1 or Charge > 3:
                continue
            self.Attempts[Charge] = self.Attempts.get(Charge, 0) + 1
            TrueParentMass = Peptide.Masses[-1] + 19
            FilePos = int(Bits[self.Columns.FileOffset])
            FileName = Bits[self.Columns.SpectrumFile]
            FileName = self.FixPath(FileName)
            Spectrum = PyInspect.Spectrum(FileName, FilePos)
            Spectrum.SetCharge(Charge)
            OriginalMass = Spectrum.GetParentMass()
            # Skip this one if the diff is huge:
            OriginalDiffRaw = (OriginalMass - TrueParentMass)
            OriginalDiff = abs(OriginalMass - TrueParentMass)
            if OriginalDiff > OriginalMass * ParentMassPPMLimit / float(1000000):
                print "(skip %s:%s; difference is %s ppm)"%(Bits[0], FilePos, 1000000 * OriginalDiff / float(OriginalMass))
                continue
            InfoList = Spectrum.CorrectParentMass()
            CorrectMassScore = None
            BestScore = None
            for InfoTuple in InfoList:
                PM = InfoTuple[0]
                Score = InfoTuple[1]
                if (BestScore == None or Score > BestScore):
                    BestPM = PM
                    BestScore = Score
                if abs(PM - TrueParentMass) < 0.1:
                    CorrectMassScore = Score
            RunnerUpScore = None
            RunnerUpPM = 0
            for InfoTuple in InfoList:
                PM = InfoTuple[0]
                Score = InfoTuple[1]                
                if abs(PM - BestPM) < 0.5:
                    continue
                if (RunnerUpScore == None or Score > RunnerUpScore):
                    RunnerUpPM = PM
                    RunnerUpScore = Score
            if OriginalDiff <= 0.1:
                self.StartOK1[Charge] = self.StartOK1.get(Charge, 0) + 1
            if OriginalDiff <= 0.3:
                self.StartOK3[Charge] = self.StartOK3.get(Charge, 0) + 1
            if OriginalDiff <= 0.5:
                self.StartOK5[Charge] = self.StartOK5.get(Charge, 0) + 1
            Diff = abs(BestPM - TrueParentMass)
            DiffRaw = (BestPM - TrueParentMass)
            if Diff <= 0.1:
                self.Success1[Charge] = self.Success1.get(Charge, 0) + 1
            if Diff <= 0.3:
                self.Success3[Charge] = self.Success3.get(Charge, 0) + 1
            if Diff <= 0.5:
                self.Success5[Charge] = self.Success5.get(Charge, 0) + 1
            RUDiff = abs(RunnerUpPM - TrueParentMass)
            RUDiffRaw = (RunnerUpPM - TrueParentMass)
            if Diff <= 0.1 or RUDiff <= 0.1:
                self.Success1RU[Charge] = self.Success1RU.get(Charge, 0) + 1
            if Diff <= 0.3 or RUDiff <= 0.3:
                self.Success3RU[Charge] = self.Success3RU.get(Charge, 0) + 1
            if Diff <= 0.5 or RUDiff < 0.5:
                self.Success5RU[Charge] = self.Success5RU.get(Charge, 0) + 1
            Str = "%s\t%s\t%s\t"%(Bits[self.Columns.SpectrumFile], Bits[self.Columns.FileOffset], Bits[self.Columns.Annotation])
            Str += "%s\t%s\t%s\t%s\t"%(Charge, TrueParentMass, OriginalMass, BestPM)
            Str += "%s\t%s\t%s\t"%(OriginalDiffRaw, DiffRaw, RUDiffRaw)
            Str += "%s\t%s\t"%(BestScore, CorrectMassScore)
            #print Str # verbose!
            self.PMCOutputFile.write(Str + "\n")
            # If the original mass was close and our corrected mass is *not*, then
            # write some info to self.PMCFailuresFile:
            if OriginalDiff <= 0.5 and Diff > 0.5:
                Str = "%.2f\t%.2f\t%.2f\t"%(TrueParentMass, OriginalMass, BestPM)
                Str += "%s\t"%Charge
                Str += "%s\t%s\t%s\t"%(Bits[self.Columns.SpectrumFile],
                    Bits[self.Columns.FileOffset], Bits[self.Columns.Annotation])
                self.PMCFailuresFile.write(Str + "\n")
        File.close()
    def ComputeAllPMCFeatures(self):
        """
        Given some search output with confident annotations, let's write out some
        parent mass correction feature vectors.  Then we'll train a machine learner on
        these feature vectors.
        """
        OutputFile1 = open("PMCFeatures1.txt", "wb")
        OutputFile2 = open("PMCFeatures2.txt", "wb")
        OutputFile3 = open("PMCFeatures3.txt", "wb")
        # Write file-headers:
        Header1 = "#Notes\tValidFlag\tAbsMassOffset\tConv-1\tConv-0.5\tConv0\tConv1\tConv2\t"
        Header2 = "#Notes\tValidFlag\tMassOffset\tConv-1\tConv-0.5\tConv0\tConv1\tConv2\tConv+1.5\tConv+0.5"
        OutputFile1.write(Header1 + "\n")
        OutputFile2.write(Header2 + "\n")
        OutputFile3.write(Header2 + "\n")
        for File in (OutputFile1, OutputFile2, OutputFile3):    
            File.write("#Feature\t\t0\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\t12\t13\t14\t15\t16\t\n")
        self.OutputFiles = [None, OutputFile1, OutputFile2, OutputFile3, OutputFile3]
        for OracleFileName in os.listdir(self.TrainingCorpusDir):
            OraclePath = os.path.join(self.TrainingCorpusDir, OracleFileName)
            if os.path.isdir(OraclePath):
                continue
            self.ComputePMCFeatures(OraclePath)
        for File in (OutputFile1, OutputFile2, OutputFile3):
            File.close()
    def ComputePMCFeatures(self, OracleFileName):
        """
        Helper for ComputeAllPMCFeatures: Compute and write out PMC features for
        each peptide in the given oracle-file.  The output is a tab-delimted file
        suitable for TrainMachineLearner.py.  Important note: Some masses will be correct (within 0.1Da),
        others will be off by more than 0.1da...to achieve a better separation between "right" and "wrong"
        points, we throw out any which are off by 0.1...0.25Da.
        """
        import PyInspect
        print OracleFileName
        File = open(OracleFileName, "rb")
        LineNumber = 0
        for FileLine in File.xreadlines():
            LineNumber += 1
            # Quick parse:
            if LineNumber >= QuickParseLineCount and self.QuickParseFlag:
                break
            if LineNumber % 100 == 0:
                print "  Line %s..."%LineNumber
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
            if Charge < 1 or Charge > len(self.OutputFiles):
                continue
            OutputFile = self.OutputFiles[Charge]
            TrueParentMass = Peptide.Masses[-1] + 19
            FilePos = int(Bits[self.Columns.FileOffset])
            FileName = Bits[self.Columns.SpectrumFile]
            FileName = self.FixPath(FileName)
            Spectrum = PyInspect.Spectrum(FileName, FilePos)
            Spectrum.SetCharge(Charge)
            PMCFeatures = Spectrum.GetPMCFeatures()
            if not PMCFeatures:
                continue
            # GetPMCFeatures returns a list of tuples.  Each tuple consists of
            # a mass, followed by all the features.
            GoodString = None
            BadStringPlus1 = None
            BadStringMinus1 = None
            BadStrings = []
            # Which vectors do we keep?  We keep the correct one, and
            # three randomly-selected others with absolute error >= 0.25
            ####################################################################
            # Keep the correct one, and three randomly-selected others:
            for Tuple in PMCFeatures[1:]:
                Mass = Tuple[0]
                Str = "%s:%s (%.3f T%.3f)\t"%(FileName, FilePos, Mass, TrueParentMass)
                Diff = Mass - TrueParentMass
                if abs(Diff) < 0.1:
                    Str += "1\t"
                    ValidFlag = 1
                elif abs(Diff) < 0.25:
                    # Skip this one - it's wrong but CLOSE.
                    continue
                elif abs(Diff) > 3.5:
                    # Ridiculously far off!
                    continue
                else:
                    Str += "0\t"
                    ValidFlag = 0
                for Feature in Tuple[1:]:
                    Str += "%.6f\t"%Feature
                #print Str
                #self.OutputFile.write(Str + "\n")
                if ValidFlag:
                    GoodString = Str
                else:
                    BadStrings.append(Str)
            random.shuffle(BadStrings)
            if not GoodString:
                print "* Warning: No correct parent mass found for %s:%s"%(FileName, FilePos)
                print "*    True mass %.4f original mass %.4f diff %.4f"%(TrueParentMass, PMCFeatures[0][0], PMCFeatures[0][0] - TrueParentMass)
                continue
            OutputFile.write(GoodString + "\n")
            for BadString in BadStrings[:3]:
                OutputFile.write(BadString + "\n")
    def ComputeAllSelfConvolutionPeaks(self, Charge, TriplyChargedFlag):
        Offsets = []
        X = -20.0
        while X < 5.5:
            Offsets.append(X)
            X += 0.1
        self.SCHistogramRaw = {}
        self.SCHistogramDot = {}
        self.SCHistogramMax = {}
        for Key in Offsets:
            self.SCHistogramRaw[Key] = 0
            self.SCHistogramDot[Key] = 0
            self.SCHistogramMax[Key] = 0
        for OracleFileName in os.listdir(self.TrainingCorpusDir):
            OraclePath = os.path.join(self.TrainingCorpusDir, OracleFileName)
            if os.path.isdir(OraclePath):
                continue            
            self.ComputeSelfConvolutionPeaks(OraclePath, Offsets, Charge, TriplyChargedFlag)
        HistogramFile = open("SelfConvolutionHistogram.%s.%s.txt"%(Charge, TriplyChargedFlag), "wb")
        for Offset in Offsets:
            String = "%s\t%s\t%s\t%s\t"%(Offset, self.SCHistogramRaw[Offset], self.SCHistogramDot[Offset],
                self.SCHistogramMax[Offset])
            HistogramFile.write(String + "\n")
    def ComputeSelfConvolutionPeaks(self, OracleFileName, Offsets, DesiredCharge, TriplyChargedFlag):
        import PyInspect
        SpectrumCount = 1000 # Use few spectra so as to run quickly!
        File = open(OracleFileName, "rb")
        LineNumber = 0
        FileLines = File.readlines()
        random.seed(1)
        random.shuffle(FileLines)
        FileLines = FileLines[:SpectrumCount]
        for FileLine in FileLines:
            LineNumber += 1
            # Ensure this is a valid line:
            Bits = FileLine.split("\t")
            if FileLine[0] == "#" or not FileLine.strip():
                continue
            try:
                PValue = float(Bits[self.Columns.PValue])
            except:
                traceback.print_exc()
                print Bits
                continue
            # Check whether the charge is the state of interest:
            Peptide = GetPeptideFromModdedName(Bits[self.Columns.Annotation])
            Charge = int(Bits[self.Columns.Charge])
            if Charge != DesiredCharge:
                continue
            TrueParentMass = Peptide.Masses[-1] + 19
            FilePos = int(Bits[self.Columns.FileOffset])
            FileName = Bits[self.Columns.SpectrumFile]
            FileName = self.FixPath(FileName)
            Spectrum = PyInspect.Spectrum(FileName, FilePos)
            Spectrum.SetCharge(Charge)
            Spectrum.SetParentMass(TrueParentMass)
            BestConvolution = 0.001
            Convolutions = {}
            # Get the dot product of the spectrum with itself:
            DirectConvolution = Spectrum.BYConvolve(0, -1)
            # Get the self-convolutions:
            for Offset in Offsets:
                Convolution = Spectrum.BYConvolve(Offset, TriplyChargedFlag)
                BestConvolution = max(BestConvolution, Convolution)
                Convolutions[Offset] = Convolution
            # Accumulate self-convolutions in our histograms:
            for Offset in Offsets:
                self.SCHistogramRaw[Offset] += Convolutions[Offset]
                self.SCHistogramDot[Offset] += Convolutions[Offset] / float(max(0.001, DirectConvolution))
                self.SCHistogramMax[Offset] += Convolutions[Offset] / float(max(0.001, BestConvolution))
    def FixPath(self, Path):
        #return "c" + Path[1:]
        return Path
    def CheckTrainingSetSuccess(self, Model, Vectors, Info):
        """
        Helper for ReportTrainingSetSuccessRate.  Score all the vectors, and determine
        whether the high scorer is the correct mass.
        """
        if not Vectors:
            return
        SortedList = []
        for Vector in Vectors:
            Score = Model.ScoreInstance(Vector.Features)
            SortedList.append((Score, Vector))
        SortedList.sort()
        SortedList.reverse()
        TopScoringVector = SortedList[0][1]
        if TopScoringVector.TrueFlag:
            self.ModelSuccessCount += 1
        else:
            self.ModelFailureCount += 1
            self.TrainingSetFailureFile.write("%s\t\n"%(Info))
        #print TopScoringVector.TrueFlag, Info
    def ReportTrainingSetSuccessRate(self):
        """
        Compute our accuracy on the training set:
        - Parse the feature vectors for one spectrum
        - Determine whether the maximum model score is attained on the correct feature
        """
        self.ModelSuccessCount = 0
        self.ModelFailureCount = 0
        PendingVectors = []
        self.TrainingSetFailureFile = open("PMCFailuresTrainingSet.txt", "wb")
        OldSpectrum = None
        for Charge in (1, 2, 3):
            ModelFileName = self.GetModelFileName(Charge)
            if self.SVMFlag:
                Model = Learning.SVMModel()
            else:
                Model = Learning.LDAModel()
            print "LOAD model from %s"%ModelFileName
            Model.LoadModel(ModelFileName)
            File = open("PMCFeatures%s.txt"%Charge, "rb")
            for FileLine in File.xreadlines():
                Bits = FileLine.split("\t")
                if len(Bits) < 10:
                    continue
                try:
                    ValidFlag = int(Bits[1])
                except:
                    continue
                Spectrum = Bits[0].split()[0]
                if Spectrum != OldSpectrum:
                    self.CheckTrainingSetSuccess(Model, PendingVectors, OldSpectrum)
                    PendingVectors = []
                    OldSpectrum = Spectrum
                Vector = Learning.FeatureVector()
                Vector.FileBits = Bits
                Vector.TrueFlag = ValidFlag
                for BitIndex in range(2, len(Bits)):
                    try:
                        Value = float(Bits[BitIndex])
                        Vector.Features.append(Value)
                    except:
                        break
                PendingVectors.append(Vector)
            # Finish the last spectrum:
            self.CheckTrainingSetSuccess(Model, PendingVectors, OldSpectrum)
            File.close()
            print "Charge %s training set results:"%Charge
            TotalCount = self.ModelSuccessCount + self.ModelFailureCount
            SuccessRate = 100 * self.ModelSuccessCount / float(max(1, TotalCount))
            print "%s of %s successful (%.2f%% accuracy)"%(self.ModelSuccessCount, TotalCount, SuccessRate)
    def ParseCommandLine(self, Arguments):
        (Options, Args) = getopt.getopt(Arguments, "fsQc:t")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-f":
                self.ComputeFeaturesFlag = 1
            elif Option == "-s":
                self.SVMFlag = 1
            elif Option == "-Q":
                self.QuickParseFlag = 1
            elif Option == "-c":
                self.TrainingCorpusDir = Value
            elif Option == "-t":
                self.TrainFlag = 0
    def Main(self):
        StartAll = time.clock()
        print "-- Compute features..."
        if self.ComputeFeaturesFlag:
            self.ComputeAllPMCFeatures()
        if self.TrainFlag:
            print "-- Train models..."
            for Charge in (1, 2, 3):
                if self.SVMFlag:
                    ModelType = "SVM"
                else:
                    ModelType = "LDA"
                Command = "TrainMachineLearner.py -r PMCFeatures%s.txt -m %s -w PMC%s%s"%(Charge, ModelType, ModelType, Charge)
                print Command
                os.system(Command)
        print "-- Compute accuracy on the training set..."
        self.ReportTrainingSetSuccessRate()
        print "-- Compute accuracy using production code..."
        GetAccuracyStart = time.clock()
        self.ComputeAndReportSuccessRate()
        Now = time.clock()
        Score = self.GetOverallScore()
        print "Score:", Score
        print "Elasped time: %s for computing PMC for spectra"%(Now - GetAccuracyStart)
        print "Overall running time: %s"%(Now - StartAll)        
            
def MainFindSelfConvolutionPeaks():
    """
    Determine the mass offsets to use when looking for b/y pairs, b2/y pairs, b/y2 pairs.
    """
    Grabber = PMCFeatureGrabber()
    Grabber.ComputeAllSelfConvolutionPeaks(1, 0)
    Grabber.ComputeAllSelfConvolutionPeaks(1, 1)    
    Grabber.ComputeAllSelfConvolutionPeaks(2, 0)
    Grabber.ComputeAllSelfConvolutionPeaks(2, 1)    
    Grabber.ComputeAllSelfConvolutionPeaks(3, 0)
    Grabber.ComputeAllSelfConvolutionPeaks(3, 1)

def MainTestCorrection():
    """
    Perform PMC on spectra with known annotations, and determine how successful we were.
    """
    Grabber = PMCFeatureGrabber()
    Grabber.ComputeAndReportSuccessRate()
    

def MainVerbosePMC():
    import PyInspect
    Stuff = r"""
E:\ms\ISB\sergei_digest_B_full_08.mgf	242966	L.YQEPVLGPVR.G	2	1157.61347	1158.35	1157.25	0.73653	0.36347	0.13653
E:\ms\Briggs\HEK293\H293b-total-try-2nd-digest-c-500ug-2D34-121905-LTQ1\H293b-total-try-2nd-digest-c-500ug-2D34-121905-LTQ1-17.mzXML	70965706	F.NLPSDGSAVDVHINMEQAPIQSEPR.V	2	2704.28164	2705.018	2703.918	0.73636	0.36364	0.83636
E:\ms\OMICS04\0400-2000mz_03_dta.mgf	3923838	K.TKIPAVFK.I	2	903.54835	904.512	903.912	0.96365	0.36365	0.86365
E:\ms\ISB\sergei_digest_B_full_05.mgf	21877	K.VGDANPALQK.V	2	1012.52432	1012.96	1012.16	0.43568	0.36432	0.13568
E:\ms\ISB\sergei_digest_A_full_01.mgf	2880933	G.GFGGAGGFGGPGGFGGSGGFGGPGSLGSPG.G	2	2370.03546	2372.2	2370.4	2.16454	0.36454	0.13546
    """
    for Line in Stuff.split("\n"):
        Bits = Line.split("\t")
        print Bits
        if len(Bits) < 2:
            continue
        TruePeptide = Bits[2]
        FileOffset = int(Bits[1])
        FileName = "C" + Bits[0][1:]
        Peptide = GetPeptideFromModdedName(TruePeptide)
        TruePM = Peptide.Masses[-1] + 19
        print "True PM:", TruePM
        Spectrum = PyInspect.Spectrum(FileName, FileOffset)
        Spectrum.SetCharge(2)
        Tuples = Spectrum.CorrectParentMass()
        for Tuple in Tuples:
            print "Mass %s score %s"%(Tuple[0], Tuple[1])
            Str = ""
            for Feature in Tuple[2:20]:
                Str += "%.4f, "%Feature
            print "  Feature vector:", Str

if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(psyco not loaded - running non-optimized)"
    Grabber = PMCFeatureGrabber()
    Grabber.ParseCommandLine(sys.argv[1:])
    Grabber.Main()
