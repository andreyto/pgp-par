"""
BuildConsensusSpectrum:
- Given many similar spectra which we consider to be the "same" (either due to
Inspect annotations or due to clustering), we'd like to generate a CONSENSUS
SPECTRUM.  The consensus spectrum should contain less noise than the individual
spectra, and mass errors should be decreased.  Also, the computation shouldn't
require too much time; for efficiency, it should require just one I/O pass through
the spectra.

Current pseudocode:

# Accumulate intensity and peak count:
For spectrum S:
  Read in peaks from disk
  Intensity1 = max(max peak intensity, grass peak intensity * 20)
  For peak P in spectrum S:
    ScaledIntensity = (Intensity(P) / Intensity1)
    MassBin = mass of P, rounded to nearest 0.1Da
    TotalIntensity[MassBin] += ScaledIntensity
    PeakCount[MassBin]++
    PeakCount[MassBin-1]++
    PeakCount[MassBin+1]++
    
# Process intensity into a peak list:
Iterate over peaks P from high to low peak-count:
  if Assimilated[P], continue
  Peak P assimilates neighboring peaks if their total intensity is lower
  new spectrum receives this assimilated peak
"""
import PyInspect
import MSSpectrum
import os
import sys
import traceback
import cPickle

USE_COUNT_FLAG = 1

# Scaling factors, for levels of peak presence from 0% to 100%.
ScalingFactors = {}
for X in range(0, 101):
    #ScalingFactors[X] = 0.95 + 0.05 * (1.0 + X / 100.0)**5
    ScalingFactors[X] = 0.95 + 0.05 * (1.0 + X / 100.0)**5
    #print X, ScalingFactors[X]


class ConsensusBuilder:
    def __init__(self, Charge = None):
        self.SpectrumCount = 0
        self.TotalMZ = 0
        self.Intensity = {}
        self.PeakCount = {}
        self.Charge = Charge
        self.SignalPeakCount = 0
    def DebugPrint(self):
        Bins = self.Intensity.keys()
        Bins.sort()
        MinBin = min(Bins)
        MaxBin = max(Bins)
        for Bin in range(MinBin, MaxBin + 1):
            Str = "%s\t%s\t%s\t"%(Bin, self.Intensity.get(Bin, 0), self.PeakCount.get(Bin, 0))
            print Str
        
    def AddSpectrum(self, Spectrum):
        """
        Add one more spectrum to the binned intensity and peak counts
        """
        self.TotalMZ += Spectrum.PrecursorMZ
        self.SpectrumCount += 1
        if not self.Charge:
            self.Charge = Spectrum.Charge # ASSUMED: all spectra have same charge!
        if not self.SignalPeakCount:
            # Expect to see this many "signal" peaks:
            self.SignalPeakCount = int(round((Spectrum.ParentMass / 100.0) * 4))
        MaxIntensity = 0
        IntensityList = []
        for Peak in Spectrum.Peaks:
            MaxIntensity = max(MaxIntensity, Peak.Intensity)
            IntensityList.append(Peak.Intensity)
        IntensityList.sort()
        # Spectra with high signal-to-noise are weighted a bit more heavily:
        if len(IntensityList) > self.SignalPeakCount:
            MinimumPresenceIntensity = IntensityList[-self.SignalPeakCount]
            ScalingFactor = min(20.0 / MaxIntensity, 1.0 / MinimumPresenceIntensity)
        else:
            ScalingFactor = 20.0 / MaxIntensity
            MinimumPresenceIntensity = 0
        #print "%s peaks; Spectrum has scaling factor %s (grass %s, max peak %s)"%(len(IntensityList), ScalingFactor, MinimumPresenceIntensity, MaxIntensity)
        for Peak in Spectrum.Peaks:
            MZBin = int(round(Peak.Mass * 10))
            ScaledIntensity = Peak.Intensity * ScalingFactor
            self.Intensity[MZBin] = self.Intensity.get(MZBin, 0) + ScaledIntensity
            if Peak.Intensity > MinimumPresenceIntensity:
                self.PeakCount[MZBin] = self.PeakCount.get(MZBin, 0) + 1
                self.PeakCount[MZBin-1] = self.PeakCount.get(MZBin-1, 0) + 1
                self.PeakCount[MZBin-2] = self.PeakCount.get(MZBin-2, 0) + 1
                self.PeakCount[MZBin+1] = self.PeakCount.get(MZBin+1, 0) + 1
                self.PeakCount[MZBin+2] = self.PeakCount.get(MZBin+2, 0) + 1
            #self.Intensity[MZBin] = self.Intensity.get(MZBin, 0) + ScaledIntensity
    def PickleCluster(self, PicklePath):
        """
        Serialize the information we've read from many spectra.  This method is used
        if we want to reserve the option to add to the cluster later.
        """
        File = open(PicklePath, "wb")
        cPickle.dump(self.Charge, File)
        cPickle.dump(self.TotalMZ, File)
        cPickle.dump(self.SpectrumCount, File)
        cPickle.dump(self.Intensity, File)
        cPickle.dump(self.PeakCount, File)
        File.close()
    def UnpickleCluster(self, PicklePath):
        """
        Sister method to PickleCluster - load a cluster from disk.
        """
        File = open(PicklePath, "rb")
        self.Charge = cPickle.load(File)
        self.TotalMZ = cPickle.load(File)
        self.SpectrumCount = cPickle.load(File)
        self.Intensity = cPickle.load(File)
        self.PeakCount = cPickle.load(File)
        File.close()
        
    def ProduceConsensusSpectrum(self):
        Spectrum = MSSpectrum.SpectrumClass()
        Spectrum.Charge = self.Charge
        Spectrum.PrecursorMZ = self.TotalMZ / float(max(self.SpectrumCount, 1))
        Spectrum.ParentMass = (Spectrum.PrecursorMZ * Spectrum.Charge) - (Spectrum.Charge - 1)*1.0078
        # Iterate over intensity entries:
        self.AssimilationFlag = {}
        SortedList = []
        for (Bin, Score) in self.Intensity.items():
            SortedList.append((Score, Bin))
        SortedList.sort()
        SortedList.reverse()
        Spectrum.Peaks = []
        for (Intensity, Bin) in SortedList:
            if self.AssimilationFlag.get(Bin, 0):
                continue
            Peak = MSSpectrum.PeakClass(Bin / 10.0, Intensity)
            for NearBin in (Bin-1, Bin-2, Bin-3, Bin+1, Bin+2, Bin+3):
                if self.AssimilationFlag.get(NearBin, 0):
                    continue
                self.AssimilationFlag[NearBin] = 1
                Peak.Intensity += self.Intensity.get(NearBin, 0)
            # Scale the intensity by the peak count, IF the cluster is large:
            if USE_COUNT_FLAG and self.SpectrumCount > 4:
                FractionPresent = self.PeakCount.get(Bin,0) / float(self.SpectrumCount)
                Bin = min(100, int(round(FractionPresent * 100)))
                #print FractionPresent, Bin, ScalingFactors[Bin]
                Peak.Intensity *= ScalingFactors[Bin]
            Spectrum.Peaks.append(Peak)
        Spectrum.Peaks.sort()
        return Spectrum
    def AssimilateCluster(self, OtherCluster):
        """
        Assimilate the other consensus information into our consensus spectrum.
        """
        self.SpectrumCount += OtherCluster.SpectrumCount
        self.TotalMZ += OtherCluster.TotalMZ
        for (Key, Value) in OtherCluster.Intensity.items():
            self.Intensity[Key] = self.Intensity.get(Key, 0) + Value
        for (Key, Value) in OtherCluster.PeakCount.items():
            self.PeakCount[Key] = self.PeakCount.get(Key, 0) + Value
def TestConsensus(AnnotationsFile, Annotation, Charge = 2):
    """
    Build a consensus spectrum for a collection of spectra, and verify
    that it looks good.
    """
    #AnnotationsFile = "ConsensusTest.txt"
    #Annotation = "M.D+173VTIQHPWFK.R"
    SpectrumDir = "e:\\ms\\lens\\spectra"
    TestOutputFile = "ConsensusTest.dta"
    Builder = ConsensusBuilder()
    InputFile = open(AnnotationsFile, "rb")
    BestScores = []
    #MGFFile = open("TestTest.mgf", "wb") #%TEMP
    #MGFFile.close() #%TEMP
    for FileLine in InputFile.xreadlines():
        Bits = FileLine.split("\t")
        try:
            FilePath = Bits[0]
            FilePos = int(Bits[15])
            SpectrumCharge = int(Bits[4])
            MQScore = float(Bits[5])
        except:
            continue
        if SpectrumCharge != Charge:
            continue
        FilePath = os.path.join(SpectrumDir, FilePath.replace("/", "\\").split("\\")[-1])
        BestScores.append((MQScore, FilePath, FilePos))
        Spectrum = MSSpectrum.SpectrumClass()
        Spectrum.ReadPeaks("%s:%s"%(FilePath, FilePos))
        Spectrum.SetCharge(Charge)
        Builder.AddSpectrum(Spectrum)
        # Try ari's thingy:
        #Spectrum.WriteMGFPeaks("TestTest.mgf") #%TEMP
    InputFile.close()
    print "Build consensus spectrum for %s members."%Builder.SpectrumCount
    Consensus = Builder.ProduceConsensusSpectrum()
    Consensus.WritePeaks(TestOutputFile)
    #Command = "MakeConsensus -mgf TestTest.mgf > %s"%(TestOutputFile) #%TEMP
    #print Command #%TEMP
    #os.system(Command) #%TEMP
    PySpectrum = PyInspect.Spectrum(TestOutputFile, 0)
    Results = PySpectrum.ScorePeptideDetailed(Annotation)
    ConsensusScore = Results[0]
    print "Consensus spectrum score: %.2f (%.2f, %.2f, %.2f, %d)"%(ConsensusScore, Results[1], Results[2], Results[3], Results[4])
    # Compare the CONSENSUS score to the average for the BEST FIVE:
    BestScores.sort()
    TopHits = 0
    for (Score, Path, Pos) in BestScores[-5:]:
        PySpectrum = PyInspect.Spectrum(Path, Pos)
        Score = PySpectrum.ScorePeptide(Annotation)
        TopHits += Score
    BestFiveAverage = TopHits / 5.0
    ScoreGain = ConsensusScore - BestFiveAverage
    print "Consensus %s vs top-5 average %s (%s)"%(ConsensusScore, BestFiveAverage, ScoreGain)
    return ScoreGain

def TestMain():
    TestCases = [("Consensus.GNTIEIQGDDAPSLWVYGFSDR.txt", "K.GNTIEIQGDDAPSLWVYGFSDR.V", 2),
                 ("ConsensusTest.txt", "-.M+42DVTIQHPWFK.R", 1),
                 ("ConsensusTest.txt", "-.M+42DVTIQHPWFK.R", 2),
                 ("Consensus.R.QD-17DHGYISR.E.txt", "R.Q-17DDHGYISR.E", 1),
                 ("Consensus.R.QD-17DHGYISR.E.txt", "R.Q-17DDHGYISR.E", 2),
                 ]
    ResultTotal = 0
    ResultCount = 0
    #TestCases = TestCases[0:1] # TEMP%!
    for (AnnotationFile, Annotation, Charge) in TestCases:
        ResultTotal += TestConsensus(AnnotationFile, Annotation, Charge)
        ResultCount += 1
    print "OVERALL RESULTS: Average MQGain is %s"%(ResultTotal / float(max(1, ResultCount)))
    
if __name__ == "__main__":
    # Given the filename of a cluster, print verbose info:
    FileName = sys.argv[1]
    Bob = ConsensusBuilder()
    Bob.UnpickleCluster(FileName)
    Bob.DebugPrint()