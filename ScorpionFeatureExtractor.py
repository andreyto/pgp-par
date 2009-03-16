"""
Scorpion Breakpoint Feature Extractor
Given a list of annotations of the form (spectrum file, peptide), write out - for training purposes - the
features for each 'covered' breakpoint.  A breakpoint is 'covered' if its b *or* y peak lie between the
spectrum's lowest and highest peak.

Features are gathered in this way:
#######Ions are processed in precedence order (y, y2, b, b2, y-h, b-h, y-n, b-n, a)
The neutral loss ions that are important in phosphorylation are:
                    b-phos, b-water, b-phos-water
                    b2-water, b2-phos
                    y-water, y-phos
                    y2-water, y2-phos
                    a ions, and phosphoric acid loss are unused, as they appear less than once per spectrum
                    water and ammonia loss are groups together in my mind.
                    
For each mass m in the ion ladder, accumulate all the intensity in the interval (m-delta, m+delta), and flag
those peaks as annotated.  Exception: if a peak is already annotated by a higher-priority ion type, then you
aren't allowed to "take it away" and re-annotate it (stealing some intensity from the original), unless you
improved the mass skew by at least 0.1

Once all ion types are processed:
Then, each ion type for each breakpoint has an intensity LEVEL assigned (from cutoffs based on the overall intensity)
Then, the feature vectors are written out for each breakpoint.

The features should hopefully not need to be recalculated much, since they have only a few parameters:
- Initial window filter (start: 50 daltons, keep top 7)
- Delta for peaks (start: 0.5)
- Intensity level-ing scheme (start: Take bottom third of peaks, and call the average intensity X; your intensity
is based on your total intensity divided by X; cutoffs are 0.1, 2, 10)

Once the features are computed we can decide on the specifics of the graph connectivity.  And once the graph
is hooked up and the probability tables saved, we can generate a big stack of features for training the
regression model / SVM / magic 8-ball.

At the moment, there are two files output, one for gross verbosity, and the other for training.  The training
file will only output the ratio for each peak (along with meta data for the break).  Other peak related data
is only output in the verbose file.  Additionally the values for Training will be classifier [0,1,2,3] and not
a continuous float value.
"""

UsageInfo = """ ScropionFeatureExtractorPhosphorylation.py
Makes feature vectors about break point in a group of spectra

Required options:
 -r [FileName] - The name of the results file to parse.  If a directory is
    specified, then all .txt files within the directory will be combined into
    one report
 -s [DirName] - The root directory of the spectra.

Additional options:
 -c [Charge] - charge for all spectra in the training set. Default value 2
 -p [Column] - column containing the path information
 -a [Column] - column containing the annotation information
 -P          - Special case for phosphorylation spectra
"""
import os
import sys
import getopt
import math
import traceback
import MSSpectrum
import ScorpionUtils
from Utils import *
Initialize()

class ScorpionFeatureExtractor(ScorpionUtils.ScorpionParser):
    def __init__(self, OutputFileName = "ScorpionBreakFeatures.txt"):
        self.OutputFileName = OutputFileName
        self.OutputHandle = None
        self.PathColumn = 0
        self.AnnotationColumn =1
        self.AnnotationFileName = None
        self.SpectraDirectory = None
        self.SpectraCharge = 2 #default
        self.Delta= 0.5
        ScorpionUtils.ScorpionParser.__init__(self)
        ## See the inherited class for the columns used in the output!!
    def OpenOutputFile(self):
        if self.OutputHandle:
            return
        self.OutputHandle = open(self.OutputFileName, "w")
        #see the ScorpionUtils for the columns used
        Header = self.BreakInfoColumns.String #inherited
        Header += self.FeaturedIons.String #inherited
        self.OutputHandle.write(Header+"\n")
        Str = ""
        for Index in range(len(Header.split("\t"))):
            Str += "%s\t"%Index
        self.OutputHandle.write(Str + "\n")
        
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
            print "\n** WARNING: Parent mass is off by %.2f! %s\n"%(BestDiff, Peptide.GetFullModdedName())
        Spectrum.Charge = BestCharge
    def SetAverageGrass(self, Spectrum):
        "Sets a baseline intensities for noise peaks"
        SortedList = []
        for Peak in Spectrum.Peaks:
            SortedList.append(Peak.Intensity)
        SortedList.sort()
        Total = 0
        GrassCount = len(SortedList) / 3
        for Intensity in SortedList[:GrassCount]:
            Total += Intensity
        AverageGrass = Total / float(GrassCount)
        self.AverageGrass = AverageGrass
    def Gather(self, FilePath, Peptide, Spectrum):
        print FilePath #%%%
        PM = Peptide.Masses[-1] + 19
        (Mean,Stdev) = Spectrum.GetLogMeanStdev()
        #print "For a spectrum M:S %.4f\t%.4f"%(Mean,Stdev)
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
        # Accumulate the dictionary RawIntensity (iontype, breakindex) -> total intensity found
        # And dictionary Peaks (iontype, breakindex) -> list of all peaks scooped up by that theopeak
        IonTypeNames = self.FeaturedIons.IonTypeNames #inherited from ScorpionUtils
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
            if self.SpecialCase == "Phosphorylation":
                ModIndexArray = Peptide.Modifications.keys()
                PhosIndex = ModIndexArray[0] #HARD CODED ASSUMPTION.  only one phos per peptide.  no other modifications
                if BreakIndex > PhosIndex: 
                    PhosOnN = 1
                else:
                    PhosOnN = 0
            else:
                PhosOnN = 0 #setting to zero always should mean it contains zero information for nonphos searches
            Str += "%s\t"%PhosOnN
            for IonName in IonTypeNames:
                ThisIon = (IonName, BreakIndex)
                Intensity = RawIntensity.get(ThisIon, 0)
                #print "RawIntensity is %.4f.  AverageGrass is %.4f.%d.%s"%(Intensity,self.AverageGrass,BreakIndex,IonName)
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
                    try:                        
                        MeanMZ = MeanMZ / float(TotalIntensity)
                    except:
                        print "wtf hax0r"
                        print TotalIntensity
                        print ClaimedPeaks
                        for Peak in ClaimedPeaks:
                            print Peak.Mass, Peak.Intensity
                        sys.stdin.readline()
                    Skew = IonMasses[ThisIon] - MeanMZ
                #Hacking
                #Ratio = IntensityLevelerStdev(Mean, Stdev, Intensity)
                Str += "%.1f\t%.1f\t%.1f\t"%(IonMasses[ThisIon], Ratio, Skew)
            #print Str
            self.OutputHandle.write(Str+"\n")
        

    def GatherFeatures(self):
        # read the input file, and get features about each annotation
        self.OpenOutputFile()
        File = open(self.AnnotationFileName, "r")
        MinFieldCount = max(self.PathColumn, self.AnnotationColumn) + 1
        LineNumber = 0
        for FileLine in File.xreadlines():
            LineNumber += 1
            Bits = FileLine.strip().split("\t")
            if len(Bits) < MinFieldCount:
                continue
            FileName = Bits[self.PathColumn].replace("/","\\")
            (FileName,ByteOffset) = FileName.split(":")
            ByteOffset = int(ByteOffset)
            FilePath = os.path.join(self.SpectraDirectory, FileName)
            if not os.path.exists(FilePath):
                print "Skip nonexistent file '%s' from line %s of %s"%(FilePath, LineNumber, self.AnnotationFileName)
                continue
            Annotation = Bits[self.AnnotationColumn]
            try:
                Peptide = GetPeptideFromModdedName(Annotation)
            except:
                print "Unable to parse peptide annotation '%s' from line %s of %s"%(Annotation, LineNumber, self.AnnotationFileName)
                continue
            Spectrum = MSSpectrum.SpectrumClass()
            try:
                Handle = open(FilePath, "rb")
                Handle.seek(ByteOffset)
                Spectrum.ReadPeaksFromFile(Handle, FilePath)
                Spectrum.FilterPeaks(100,7) # not in Stephen's original code, but in inspect.
                if self.SpecialCase == "Phosphorylation":
                    self.RemovePMPeaks(Spectrum)
                self.CorrectParentMass(Peptide, Spectrum)
                self.SetAverageGrass(Spectrum)
                Handle.close()
            except:
                traceback.print_exc()
                print "Unable to parse spectrum file '%s' from line %s of %s"%(FilePath, LineNumber, self.AnnotationFileName)
                continue
            # We have the peptide and the spectrum; pass them to Gather()
            self.Gather(FilePath, Peptide, Spectrum)
        File.close()

    def RemovePMPeaks(self,Spectrum):
        # This function removes peaks related to the parent mass, that might confuse the program later
        # they don't correspond to any b/y breakages, so don't gather them.
        PM = Spectrum.ParentMass
        Charge = Spectrum.Charge
        PhosphateLoss = 98
        PhosphateAndWaterLoss = 98 + 18
        PhosphateLossPeak = (PM + Charge -1 -PhosphateLoss) / Charge
        #print "This is my m/z %f and this is my loss of 49 %f"%(Spectrum.PrecursorMZ, PhosphateLossPeak)
        Peak = Spectrum.GetPeak(PhosphateLossPeak)
        #print "I am going to remove a peak with intensity %s"%Peak.Intensity
        if Peak:
            Spectrum.Peaks.remove(Peak)
        PhosphateAndWaterLossPeak = (PM + Charge -1 -PhosphateAndWaterLoss) / Charge
        Peak = Spectrum.GetPeak(PhosphateAndWaterLossPeak)
        if Peak:
            Spectrum.Peaks.remove(Peak)

    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:s:p:a:c:P")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-r":
                # -r results file(s)
                if not os.path.exists(Value):
                    print "** Error: couldn't find results file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.AnnotationFileName = Value
            if Option == "-s":
                if not os.path.exists(Value):
                    print "** Error: couldn't find spectra dir file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.SpectraDirectory = Value
            if Option == "-p":
                #column containing the path
                self.PathColumn = int (Value)
            if Option == "-a":
                self.AnnotationColumn = int (Value)
            if Option == "-c":
                self.SpectraCharge = int (Value)
            if Option == "-P": #phosphorylation flag
                ScorpionUtils.ScorpionParser.SetSpecialCase(self,"Phosphorylation")
        if not OptionsSeen.has_key("-r") or not OptionsSeen.has_key("-s"):
            print UsageInfo
            sys.exit(1)



def IntensityLevelerStdev(Mean,Stdev,Intensity):
    if Intensity < 0.001: #avoid zero
        Intensity = 0.001
    LogIntensity = math.log(Intensity)
    if LogIntensity < Mean - Stdev:
        return 0
    elif LogIntensity < Mean:
        return 1
    elif LogIntensity < Mean + Stdev:
        return 2
    elif LogIntensity < Mean + Stdev + Stdev:
        return 3
    return 4

if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "<no psyco>"
    Bob = ScorpionFeatureExtractor()
    Bob.ParseCommandLine(sys.argv[1:])
    print "Commencing feature extraction..."
    Global.FixedMods = {"C":57.0518} #%%% USUALLY, but not always...
    Bob.GatherFeatures()
    print "Feature extraction complete."
