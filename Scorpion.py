"""
Scorpion Breakpoint Feature Extractor

Given a list of annotations of the form (spectrum file, peptide),
write out - for training purposes - the features for each 'covered'
breakpoint.  A breakpoint is 'covered' if its b *or* y peak lie
between the spectrum's lowest-mass and highest-mass peaks.

Features are gathered in this way: Ions are processed in precedence
order (y, y2, b, b2, y-h, b-h, y-n, b-n, a) For each mass m in the ion
ladder, accumulate all the intensity in the interval (m-delta,
m+delta), and flag those peaks as annotated.  Exception: if a peak is
already annotated by a higher-priority ion type, then you aren't
allowed to "take it away" and re-annotate it (stealing some intensity
from the original), unless you improved the mass skew by at least 0.1

Once all ion types are processed: Then, each ion type for each
breakpoint has an intensity LEVEL assigned (from cutoffs based on the
overall intensity) Then, the feature vectors are written out for each
breakpoint.

The features should hopefully not need to be recalculated much, since they
have only a few parameters:
- Initial window filter (start: 50 daltons, keep top 6)
- Delta for peaks (start: 0.5)
- Intensity level-ing scheme (start: Take bottom third of peaks, and call the
average intensity X; your intensity is based on your total intensity divided
by X; cutoffs are 0.1, 2, 10)

Once the features are computed we can decide on the specifics of the
graph connectivity.  And once the graph is hooked up and the
probability tables saved, we can generate a big stack of features for
training the regression model / SVM / magic 8-ball.
"""

import os
import sys
import math
import getopt
import traceback
import MSSpectrum
from Utils import *
Initialize()

class ScorpionFeatureExtractor:
    Delta = 0.5
    def __init__(self, OutputFileName = None):
        if not OutputFileName:
            OutputFileName = "ScorpionBreakFeatures.txt"
        self.OutputFileName = OutputFileName
        self.OutputFile = None
    def OpenOutputFile(self):
        if self.OutputFile:
            return
        self.OutputFile = open(self.OutputFileName, "w")
        Header = """Peptide\tFile\tBreakIndex\tPrefixAA\tSuffixAA\tNBase\tNAcid\tCBase\tCAcid\tSector\tCharge\tyMass\tyRatio\tySkew\tbMass\tbRatio\tbSkew\ty2Mass\ty2Ratio\ty2Skew\tb2\tb2\tb2\ty-h2o\ty-h2o\ty-h2o\tb-h2o\tb-h2o\tb-h2o\ty-nh3\ty-nh3\ty-nh3\tb-nh3\tb-nh3\tb-nh3\ta\ta\ta\ty2-h2o\ty2-h2o\ty2-h2o\tb2-h2o\tb2-h2o\tb2-h2o\ty2-nh3\ty2-nh3\ty2-nh3\tb2-nh3\tb2-nh3\tb2-nh3\ta-h2o\ta-h2o\ta-h2o\ta-nh3\ta-nh3\ta-nh3\ty-h2o-nh3\ty-h2o-nh3\ty-h2o-nh3\ty-h2o-h2o\ty-h2o-h2o\ty-h2o-h2o\tb-h2o-nh3\tb-h2o-nh3\tb-h2o-nh3\tb-h2o-h2o\tb-h2o-h2o\tb-h2o-h2o"""
        self.OutputFile.write(Header+"\n")
        Str = ""
        for Index in range(len(Header.split("\t"))):
            Str += "%s\t"%Index
        self.OutputFile.write(Str + "\n")
        
    def CorrectParentMass(self, Peptide, Spectrum):
        "Correct the spectrum's parent mass to match our peptide"
        PM = Peptide.GetParentMass()
        BestDiff = 9999
        for Charge in range(1, 5):
            ParentMass = Spectrum.PrecursorMZ * Charge - (Charge - 1)*1.0078
            Diff = abs(ParentMass - PM)
            if Diff < BestDiff:
                BestDiff = Diff
                BestCharge = Charge
                BestMass = ParentMass
        if BestDiff > 3:
            print "\n** WARNING: Parent mass is off by %.2f!\n"%BestDiff
        Spectrum.Charge = BestCharge
    def GetIntensityCutoffs(self, Spectrum):
        "Compute the IntensityCutoffs list based upon ratios with the mean intensity of the weak peaks."
        SortedList = []
        for Peak in Spectrum.Peaks:
            SortedList.append(Peak.Intensity)
        SortedList.sort()
        Total = 0
        GrassCount = len(SortedList) / 3
        for Intensity in SortedList[:GrassCount]:
            Total += Intensity
        AverageGrass = Total / float(GrassCount)
        IntensityCutoffs = []
        IntensityCutoffs.append(AverageGrass * 0.1)
        IntensityCutoffs.append(AverageGrass * 2.0)
        IntensityCutoffs.append(AverageGrass * 10.0)
        self.AverageGrass = AverageGrass
        return IntensityCutoffs
    def Gather(self, FilePath, Peptide, Spectrum):
        """
        Called by GatherFeatures: Handle a single peptide annotation.
        """
        print FilePath #%%%
        self.CorrectParentMass(Peptide, Spectrum)
        PM = Peptide.Masses[-1] + 19
        ############################################################
        # Decide what the first and last valid break points are:
        MinBreakPoint = 999
        MaxBreakPoint = 0
        B = Global.AllIonDict["b"]
        Y = Global.AllIonDict["y"]
        MinPeak = Spectrum.Peaks[0].Mass
        MaxPeak = Spectrum.Peaks[-1].Mass
        for MassIndex in range(1, len(Peptide.Masses) - 1):
            Mass = Peptide.Masses[MassIndex]
            BMass = B.GetPeakMass(Mass, PM)
            YMass = Y.GetPeakMass(Mass, PM)
            #print "Mass index %s: B mass %s, Y mass %s"%(MassIndex, BMass, YMass)
            if (MinPeak <= BMass and BMass <= MaxPeak) or (MinPeak <= YMass and YMass <= MaxPeak):
                # This break point is valid.
                MinBreakPoint = min(MinBreakPoint, MassIndex)
                MaxBreakPoint = max(MaxBreakPoint, MassIndex)
        #print "Use masses %s...%s as break points."%(MinBreakPoint, MaxBreakPoint)
        ############################################################
        IntensityCutoffs = self.GetIntensityCutoffs(Spectrum)
        ############################################################
        # Accumulate the dictionary RawIntensity (iontype, breakindex) -> total intensity found
        # And dictionary Peaks (iontype, breakindex) -> list of all peaks scooped up by that theopeak
        IonTypeNames = ["y","b","y2","b2","y-h2o","b-h2o","y-nh3","b-nh3", "a",
                        "y2-h2o", "b2-h2o", "y2-nh3", "b2-nh3", "a-h2o", "a-nh3", "y-h2o-nh3", "y-h2o-h2o", "b-h2o-nh3", "b-h2o-h2o"]
        RawIntensity = {}
        Peaks = {}
        IonMasses = {} 
        for IonName in IonTypeNames:
            IonType = Global.AllIonDict[IonName]
            for BreakIndex in range(MinBreakPoint, MaxBreakPoint + 1):
                ThisIon = (IonName, BreakIndex)
                Mass = Peptide.Masses[BreakIndex]
                Mass = IonType.GetPeakMass(Mass, PM)
                IonMasses[(IonName, BreakIndex)] = Mass
                MinMass = Mass - self.Delta
                MaxMass = Mass + self.Delta
                for Peak in Spectrum.Peaks:
                    if Peak.Mass > MaxMass:
                        break
                    if Peak.Mass < MinMass:
                        continue
                    if Peak.IonType != None:
                        # Some other ion already tried to claim the peak.  Can we seize it back?
                        OldSkew = abs(Peak.Mass - IonMasses[Peak.IonType])
                        NewSkew = abs(Peak.Mass - Mass)
                        if (OldSkew < NewSkew + 0.1):
                            continue
                        # Ok, we're stealing the peak away from the old guy:
                        Peaks[Peak.IonType].remove(Peak)
                        RawIntensity[Peak.IonType] -= Peak.Intensity
                    Peak.IonType = ThisIon
                    if not Peaks.has_key(ThisIon):
                        Peaks[ThisIon] = []
                    Peaks[ThisIon].append(Peak)
                    RawIntensity[ThisIon] = RawIntensity.get(ThisIon, 0) + Peak.Intensity
        ############################################################
        # We are now ready to write out our features.  We want to err on the side of spewing out
        # too much stuff.
        SectorEdges = []
        SectorEdges.append(PM * 0.33)
        SectorEdges.append(PM * 0.66)
        for BreakIndex in range(MinBreakPoint, MaxBreakPoint + 1):
            Mass = Peptide.Masses[BreakIndex]
            Str = "%s\t%s\t%s\t"%(Peptide.GetModdedName(), FilePath, BreakIndex)
            Str += "%s\t%s\t"%(Peptide.Aminos[BreakIndex - 1], Peptide.Aminos[BreakIndex])
            NBase = 0
            NAcid = 0
            for Amino in Peptide.Aminos[:BreakIndex]:
                if Amino == "K" or Amino == "R":
                    NBase += 1
                if Amino == "D" or Amino == "E":
                    NAcid += 1
            CBase = 0
            CAcid = 0
            for Amino in Peptide.Aminos[BreakIndex:]:
                if Amino == "K" or Amino == "R":
                    CBase += 1
                if Amino == "D" or Amino == "E":
                    CAcid += 1
            Str += "%s\t%s\t%s\t%s\t"%(NBase, NAcid, CBase, CAcid)
            Sector = Mass / PM # report sector as a float, we can tweak sector edges later if we like
            Str += "%s\t"%Sector
            Str += "%s\t"%Spectrum.Charge
            for IonName in IonTypeNames:
                ThisIon = (IonName, BreakIndex)
                Intensity = RawIntensity.get(ThisIon, 0)
                Ratio = Intensity / self.AverageGrass
                ClaimedPeaks = Peaks.get(ThisIon,[])
                if not ClaimedPeaks:
                    Skew = 0
                else:
                    MeanMZ = 0
                    TotalIntensity = 0.1 # avoid division-by-zero on retard spectra 
                    for Peak in ClaimedPeaks:
                        MeanMZ += Peak.Intensity * Peak.Mass
                        TotalIntensity += Peak.Intensity
                    MeanMZ = MeanMZ / float(TotalIntensity)
                    Skew = IonMasses[ThisIon] - MeanMZ
                Str += "%.1f\t%.1f\t%.1f\t"%(IonMasses[ThisIon], Ratio, Skew)
            print Str
            self.OutputFile.write(Str+"\n")
    def GatherFeatures(self, AnnotationFileName, FileNameField, AnnotationField, SpectrumBaseDir, FilePosField = None):
        self.OpenOutputFile()
        File = open(AnnotationFileName, "r")
        MinFieldCount = max(FileNameField, AnnotationField) + 1
        LineNumber = 0
        # Iterate over the lines from AnnotationFileName.  For each one, compute features for every
        # covered breakpoint by calling Gather
        for FileLine in File.xreadlines():
            LineNumber += 1
            Bits = FileLine.strip().split("\t")
            if len(Bits) < MinFieldCount:
                continue
            FileName = Bits[FileNameField].replace("/","\\")
            FilePath = os.path.join(SpectrumBaseDir, FileName)
            if not os.path.exists(FilePath):
                print "Skip nonexistent file '%s' from line %s of %s"%(FilePath, LineNumber, AnnotationFileName)
                continue
            if FilePosField != None:
                FilePath = "%s:%s"%(FilePath, Bits[FilePosField])
            Annotation = Bits[AnnotationField]
            try:
                Peptide = GetPeptideFromModdedName(Annotation)
            except:
                print "Unable to parse peptide annotation '%s' from line %s of %s"%(Annotation, LineNumber, AnnotationFileName)
                continue
            Spectrum = MSSpectrum.SpectrumClass()
            try:
                Spectrum.ReadPeaks(FilePath)
            except:
                traceback.print_exc()
                print "Unable to parse spectrum file '%s' from line %s of %s"%(FilePath, LineNumber, AnnotationFileName)
                continue
            # We have the peptide and the spectrum; pass them to Gather()
            self.Gather(FilePath, Peptide, Spectrum)
        File.close()
    def ComputeFTableProbability(self, FeatureFileName, FeatureNumbers, FeatureMungers, ChildFeatureNumber, ChildMunger,
                                 PossibleChildValues):
        File = open(FeatureFileName, "r")
        File.readline() # discard
        File.readline() # discard
        TotalScore = 0
        TotalCorrectGuess = 0
        TotalIncorrectGuess = 0
        LineNumber = 0
        for FileLine in File.xreadlines():
            LineNumber += 1
            Bits = FileLine.strip().split("\t")
            # Get the child feature value:
            if ChildMunger:
                ChildValue = ChildMunger(Bits[ChildFeatureNumber])
            else:
                ChildValue = int(Bits[ChildFeatureNumber])
            # Get the parent feature values:
            ParentValues = []
            for Index in range(len(FeatureNumbers)):
                if FeatureMungers[Index]:
                    ParentValue = FeatureMungers[Index](Bits[FeatureNumbers[Index]])
                else:
                    ParentValue = int(Bits[FeatureNumbers[Index]])
                ParentValues.append(ParentValue)
            # Add to the count:
            BestProb = 0
            BestChild = None
            Values = ParentValues[:]
            Values.append(ChildValue)
            for Possible in PossibleChildValues:
                Values[-1] = Possible
                Probability = self.ProbabilityTable[tuple(Values)]
                if Probability  > BestProb:
                    BestProb = Probability
                    BestChild = Possible
            if BestChild == ChildValue:
                TotalCorrectGuess += 1
            else:
                TotalIncorrectGuess += 1
            Values[-1] = ChildValue
            Prob = self.ProbabilityTable.get(tuple(Values), 0.5)
            TotalScore += math.log(Prob)
        print "# Total score over %d training lines: %.4f"%(LineNumber, TotalScore)
        print "# Overall, %d correct and %d incorrect guesses (%.2f%%)"%(TotalCorrectGuess, TotalIncorrectGuess,
                                                                         100*TotalCorrectGuess / float(TotalCorrectGuess+TotalIncorrectGuess))
        print "\n\n"
    def ConstructProbabilityTable(self, FeatureFileName, FeatureNumbers, FeatureMungers, ChildFeatureNumber, ChildMunger,
        OutputFileName):
        """
        Given a list of parent features (and possibly a list of functions to process the parent features),
        generate a probability table for the child feature.
        """
        print "# Probability table %s"%OutputFileName
        print "# Feature file %s"%FeatureFileName
        print "# Feature numbers %s"%str(FeatureNumbers)
        print "# Child feature %s"%ChildFeatureNumber
        OutputFile = open(OutputFileName, "w")
        File = open(FeatureFileName, "r")
        Header = File.readline()
        Header = File.readline()
        HeaderBits = Header.split("\t")
        ParentCountTable = {}
        CountTable = {}
        PossibleValues = {} # column number -> dict of possible values
        for Index in FeatureNumbers:
            PossibleValues[Index] = {}
        PossibleValues[ChildFeatureNumber] = {}
        TotalLines = 0
        LineNumber = 0
        for FileLine in File.xreadlines():
            LineNumber += 1
            Bits = FileLine.strip().split("\t")
            # Get the child feature value:
            if ChildMunger:
                ChildValue = ChildMunger(Bits[ChildFeatureNumber])
            else:
                ChildValue = int(Bits[ChildFeatureNumber])
            PossibleValues[ChildFeatureNumber][ChildValue] = 1 # Flag a possible value
            # Get the parent feature values:
            ParentValues = []
            for Index in range(len(FeatureNumbers)):
                if FeatureMungers[Index]:
                    ParentValue = FeatureMungers[Index](Bits[FeatureNumbers[Index]])
                else:
                    ParentValue = int(Bits[FeatureNumbers[Index]])
                PossibleValues[FeatureNumbers[Index]][ParentValue] = 1 # Flag a possible value
                ParentValues.append(ParentValue)
            # Add to the count:
            Values = ParentValues[:]
            ParentCountTable[tuple(Values)] = ParentCountTable.get(tuple(Values), 0) + 1
            Values.append(ChildValue)
            Values = tuple(Values)
            CountTable[Values] = CountTable.get(Values, 0) + 1
            #print LineNumber, Values
            TotalLines += 1
        print "# Read in %d lines of features."%LineNumber
        for Key in PossibleValues.keys():
            List = PossibleValues[Key].keys()
            List.sort()
            PossibleValues[Key] = List
        self.ProbabilityTable = {}
        # Let's iterate over all possible tuples of the form (Parent1Value, Parent2Value, ..., ParentNValue, FeatureValue)
        # Iteration rules:
        # If the last value-index is too large, decrease list size by one.  If that takes you to length 0 stop, otherwise increment the
        #  precending list entry and carry on. 
        # If the list is short, extend it with the first value for the next parent.
        # If the list is full-length, increment the last parent value.  
        CurrentParent = 0
        CurrentValues = []
        while (1):
            ValueLen = len(CurrentValues)
            #print "CurrentValues:", CurrentValues
            if ValueLen and CurrentValues[-1] >= len(PossibleValues[FeatureNumbers[ValueLen-1]]):
                CurrentValues = CurrentValues[:-1]
                if not CurrentValues:
                    break
                CurrentValues[-1] += 1
                continue
            if ValueLen < len(FeatureNumbers):
                CurrentValues.append(0)
                continue
            BaseStr = ""
            ValueList = []
            for Index in range(ValueLen):
                PValue = PossibleValues[FeatureNumbers[Index]][CurrentValues[Index]]
                BaseStr += "%s\t"%PValue
                ValueList.append(PValue)
            ValueTuple = tuple(ValueList)
            # The probability for (A, B, C, X) is equal to the total (A, B, C, X)
            # divided by the total entries of the form (A, B, C, Z) for all Z.
            for PossibleCValue in PossibleValues[ChildFeatureNumber]:
                FullValueList = ValueList[:]
                FullValueList.append(PossibleCValue)
                Count = CountTable.get(tuple(FullValueList), 0)
                Probability = Count / float(ParentCountTable.get(ValueTuple, 1))
                Str = BaseStr + "%s\t%s"%(PossibleCValue, Probability)
                self.ProbabilityTable[tuple(FullValueList)] = min(0.99, max(0.01, Probability))
                print Str
                OutputFile.write(Str+"\n")
            # Special case: If no parents, stop now!
            if len(CurrentValues) == 0:
                break
            CurrentValues[-1] += 1
        ##################
        File.close()
        self.ComputeFTableProbability(FeatureFileName, FeatureNumbers, FeatureMungers, ChildFeatureNumber, ChildMunger, PossibleValues[ChildFeatureNumber])
    def ParseCommandLine(self, Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:w:i:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-r":
                self.ParseFilePath = Value
            elif Option == "-w":
                self.WriteOutputPath = Value
            elif Option == "-i":
                self.IntensityScheme = int(Value)
    def Main(self):
        if not self.ParseFilePath:
            print "* Error: Please specify a file (or directory) to parse"
            print UsageInfo
            return
        
def IntensityLeveler(Str):
    IntensityRatio = float(Str)
    if IntensityRatio < 0.1:
        return 0
    elif IntensityRatio < 2:
        return 1
    elif IntensityRatio < 10:
        return 2
    return 3

def SectorLeveler(Str):
    Sector = float(Str)
    if Sector > 0.66:
        return 2
    if Sector > 0.33:
        return 1
    return 0

def AcidBaseLeveler(Str):
    return min(1, int(Str))

def FeatureExtractionFunkyTown():
    Bob = ScorpionFeatureExtractor()
    try:
        os.makedirs("Scorpion")
    except:
        continue
    #FeatureFileName = "Temp.txt"
    FeatureFileName = "ScorpionBreakFeatures.txt"
    Bob.ConstructProbabilityTable(FeatureFileName, [], [], 24, IntensityLeveler, "Scorpion\\None.txt")
    Bob.ConstructProbabilityTable(FeatureFileName, [12], [IntensityLeveler], 24, IntensityLeveler, "Scorpion\\y.txt")
    Bob.ConstructProbabilityTable(FeatureFileName, [15], [IntensityLeveler], 24, IntensityLeveler, "Scorpion\\b.txt")
    Bob.ConstructProbabilityTable(FeatureFileName, [9], [SectorLeveler], 24, IntensityLeveler, "Scorpion\\Sector.txt")
    Bob.ConstructProbabilityTable(FeatureFileName, [10], [None], 24, IntensityLeveler, "Scorpion\\Charge.txt")
    Bob.ConstructProbabilityTable(FeatureFileName, [9, 10], [SectorLeveler, None], 24, IntensityLeveler, "Scorpion\\SectorCharge.txt")
    Bob.ConstructProbabilityTable(FeatureFileName, [5], [None], 24, IntensityLeveler, "Scorpion\\NBase.txt")
    Bob.ConstructProbabilityTable(FeatureFileName, [6], [None], 24, IntensityLeveler, "Scorpion\\NAcid.txt")
    Bob.ConstructProbabilityTable(FeatureFileName, [7], [None], 24, IntensityLeveler, "Scorpion\\CBase.txt")
    Bob.ConstructProbabilityTable(FeatureFileName, [8], [None], 24, IntensityLeveler, "Scorpion\\CAcid.txt")
    Bob.ConstructProbabilityTable(FeatureFileName, [12, 9, 10], [IntensityLeveler, SectorLeveler, None], 24, IntensityLeveler, "Scorpion\\ySectorCharge.txt")
    Bob.ConstructProbabilityTable(FeatureFileName, [15, 9, 10], [IntensityLeveler, SectorLeveler, None], 24, IntensityLeveler, "Scorpion\\bSectorCharge.txt")
    Bob.ConstructProbabilityTable(FeatureFileName, [12, 9, 10,5], [IntensityLeveler, SectorLeveler, None, AcidBaseLeveler], 24, IntensityLeveler, "Scorpion\\ySectorChargeNBase.txt")
    Bob.ConstructProbabilityTable(FeatureFileName, [12, 9, 10,6], [IntensityLeveler, SectorLeveler, None, AcidBaseLeveler], 24, IntensityLeveler, "Scorpion\\ySectorChargeNAcid.txt")
    Bob.ConstructProbabilityTable(FeatureFileName, [12, 9, 10,7], [IntensityLeveler, SectorLeveler, None, AcidBaseLeveler], 24, IntensityLeveler, "Scorpion\\ySectorChargeCBase.txt")
    Bob.ConstructProbabilityTable(FeatureFileName, [12, 9, 10,8], [IntensityLeveler, SectorLeveler, None, AcidBaseLeveler], 24, IntensityLeveler, "Scorpion\\ySectorChargeCAcid.txt")
    Bob.ConstructProbabilityTable(FeatureFileName, [5,12], [None,IntensityLeveler], 24, IntensityLeveler, "Scorpion\\yNBase.txt")
    Bob.ConstructProbabilityTable(FeatureFileName, [6,12], [None,IntensityLeveler], 24, IntensityLeveler, "Scorpion\\yNAcid.txt")
    Bob.ConstructProbabilityTable(FeatureFileName, [7,12], [None,IntensityLeveler], 24, IntensityLeveler, "Scorpion\\yCBase.txt")
    Bob.ConstructProbabilityTable(FeatureFileName, [8,12], [None,IntensityLeveler], 24, IntensityLeveler, "Scorpion\\yCAcid.txt")




UsageInfo = """
Scorpion: Compute features for Ion-based Scoring
Arguments:
-r [FileName]: Read annotations from this file (or directory)
-w [FileName]: Write features to this file (or directory)
-i [IntensityScheme]: Use the specified scheme for getting intensity levels
"""

if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "<no psyco>"
    Bob = ScorpionFeatureExtractor()
    Bob.ParseCommandLine(sys.argv[1:])
    Bob.Main()
##
##    
##    print "Commencing feature extraction..."
##    Global.FixedMods = {} 
##    Bob.GatherFeatures("OMICS2002.pv", 0, 2, "", 15)
##    Global.FixedMods = {"C": 57.0518} #%%% USUALLY, but not always...
##    Bob.GatherFeatures("OMICS2004.pv", 0, 2, "", 15)
##    Bob.GatherFeatures("HEK.pv", 0, 2, "", 15)
##    ##Bob.GatherFeatures("AriAnnotations.IKKB.txt", 0, 1, "") # ikkb_z2 directory is already there
##    ##Bob.GatherFeatures("AriAnnotations.Lens.Sprot.txt", 0, 1, "c:\\GoodLensData")
##    ##Global.FixedMods = {}
##    ##Bob.GatherFeatures("AriAnnotations.ISB.txt", 0, 1, "c:\\source\\bafna\\AllSpectra")
##    print "Feature extraction complete."
