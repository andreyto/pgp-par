"""
This is an auxillary module for ComputePTMFeatures.py.  It measurs the similarity between
two spectra.  The simplest case is a comparison of two spectra with the same annotation.
We are also able to compare two spectra which differ by a PTM (provided the PTM position
is known).  
"""

import MSSpectrum
import Label
import math
from Utils import *
Initialize()

class SpectralSimilarity:
    """
    Container class that holds the spectra, and measures similarity
    ASSUMPTION: Either AnnotationA and AnnotationB are the same,
    or AnnotationA is modified and AnnotationB is NOT.
    """
    def __init__(self, SpectrumA, SpectrumB, AnnotationA, AnnotationB):
        """
        SpectrumA and SpectrumB are SpectrumClass objects from MSSpectrum.py
        AnnotationA and AnnotationB are strings of the Inspect annotation of the
        spectrum (if the annotation was R.THISPEPK.M  this method would get
        passed "THISPEPK")
        """
        self.BinMultiplier = 1.0 # default
        # Spectra can be either MSSpectrum objects or file paths:
        if isinstance(SpectrumA, MSSpectrum.SpectrumClass):
            self.SpectrumA = SpectrumA
        else:
            Spectrum = MSSpectrum.SpectrumClass()
            Spectrum.ReadPeaks(SpectrumA)
            self.SpectrumA = Spectrum
            self.SpectrumA.FilterPeaks() 
            self.SpectrumA.RankPeaksByIntensity()
        if isinstance(SpectrumB, MSSpectrum.SpectrumClass):
            self.SpectrumB = SpectrumB
        else:
            Spectrum = MSSpectrum.SpectrumClass()
            Spectrum.ReadPeaks(SpectrumB)
            self.SpectrumB = Spectrum
            self.SpectrumB.FilterPeaks()
            self.SpectrumB.RankPeaksByIntensity()            
        # Annotations can be either a string, or a peptide object:
        if isinstance(AnnotationA, PeptideClass):
            self.PeptideA = AnnotationA
        else:
            self.PeptideA = GetPeptideFromModdedName(AnnotationA)
        if isinstance(AnnotationB, PeptideClass):
            self.PeptideB = AnnotationB
        else:
            self.PeptideB = GetPeptideFromModdedName(AnnotationB)
    def ComputeSimilarity(self):
        """
        This method determines how similar two spectra are.  It can use a variety of methods,
        and honestly this may turn into a big switch box
        """
        DPSimScore = self.DotProduct()
        print "The Scaled Dot Product of these two is %f"%DPSimScore
    def LabelPeaks(self, PeakTolerance = 0.5):
        """
        Should be called once, *if* the annotations differ, so that we can match
        corresponding peaks which have different masses due to PTM attachment.
        """
        # Label the spectra so that I can know which peaks belong to what name eg. b7
        Label.LabelSpectrum(self.SpectrumA, self.PeptideA, PeakTolerance)
        Label.LabelSpectrum(self.SpectrumB, self.PeptideB, PeakTolerance)
    def DotProductSignal(self, BinMultiplier = 1.0, MaxIntensityRank = 50, EnableShift = 1, VerboseFlag = 0, HashByRank = 0):
        """
        Variant of Dot Product, incorporating a correction factor introduced in [Parag and Mallick, 2006ish]
        """
        self.MaxIntensityRank = MaxIntensityRank
        self.BinMultiplier = BinMultiplier
        self.HashByRank = HashByRank
        # set up hashes
        HashA = {}
        HashB = {}
        # Populate HashA and HashB:
        if EnableShift and len(self.PeptideA.Modifications) > 0:
            self.HashPeaksWithShift(HashA, self.SpectrumA, self.PeptideA, VerboseFlag)
        else:
            self.HashPeaks(HashA, self.SpectrumA)
        if EnableShift and len(self.PeptideB.Modifications) > 0:
            self.HashPeaksWithShift(HashB, self.SpectrumB, self.PeptideB, VerboseFlag)
        else:
            self.HashPeaks(HashB, self.SpectrumB)
        #Do Dot Product
        MaxBins = max(HashA.keys())
        MaxBins = max(MaxBins, max(HashB.keys()))
        DotProduct = 0
        SumSqA = 0
        SumSqB = 0
        TotalIntensityA = 0
        TotalIntensityB = 0
        for I in range(MaxBins):
            A = HashA.get(I, 0)
            B = HashB.get(I, 0)
            TotalIntensityA += A
            TotalIntensityB += B
            if VerboseFlag and (A or B):
                print "%s\t%s\t%s\t%s\t"%(I, A, B, A*B)
            DotProduct += A * B
            SumSqA += A * A
            SumSqB += B * B
        #print "Dot Product is %f"%DotProduct
        #print "sqrt thing is %f"%math.sqrt(SumSqA * SumSqB)
        OddsCollision = 1.0 / (max(self.SpectrumA.ParentMass, self.SpectrumB.ParentMass) * BinMultiplier)
        DotProduct -= TotalIntensityA * TotalIntensityB * OddsCollision
        return DotProduct / math.sqrt(SumSqA * SumSqB)
    def DotProduct(self, BinMultiplier = 1.0, MaxIntensityRank = 50, EnableShift = 1, VerboseFlag = 0, HashByRank = 0):
        """
        This method measures similarity between spectra by calculating their dot product.
        It is written to work on spectra might be PTM shifted versions of each other.
        e.g. SAMMY and SAM+16MY.  If the peptide annotation has a PTM, then the peaks
        are shifted back.
        Variables:
        Rank = indicates that you want to use a rank based dotproduct, instead of intensity based
        """
        self.MaxIntensityRank = MaxIntensityRank
        self.BinMultiplier = BinMultiplier
        self.HashByRank = HashByRank
        # set up hashes
        HashA = {}
        HashB = {}
        # Populate HashA and HashB:
        if EnableShift and len(self.PeptideA.Modifications) > 0:
            self.HashPeaksWithShift(HashA, self.SpectrumA, self.PeptideA, VerboseFlag)
        else:
            self.HashPeaks(HashA, self.SpectrumA)
        if EnableShift and len(self.PeptideB.Modifications) > 0:
            self.HashPeaksWithShift(HashB, self.SpectrumB, self.PeptideB, VerboseFlag)
        else:
            self.HashPeaks(HashB, self.SpectrumB)
        #Do Dot Product
        MaxBins = max(HashA.keys())
        MaxBins = max(MaxBins, max(HashB.keys()))
        DotProduct = 0
        SumSqA = 0
        SumSqB = 0
        for I in range(MaxBins):
            A = HashA.get(I, 0)
            B = HashB.get(I, 0)
            if VerboseFlag and (A or B):
                print "%s\t%s\t%s\t%s\t%s\t"%(I, I/self.BinMultiplier, A, B, A*B)
            DotProduct += A * B
            SumSqA += A * A
            SumSqB += B * B
        #print "Dot Product is %f"%DotProduct
        #print "sqrt thing is %f"%math.sqrt(SumSqA * SumSqB)
        return DotProduct / math.sqrt(SumSqA * SumSqB)
    def HashPeaksWithShift(self, Hash, Spectrum, Peptide, VerboseFlag = 0):
        """
        Takes a peptide, spectrum, and hashtable
        and puts values into the hash so that it can be dotproducted (now there's a verb
        from a noun!)
        The caveat is that we shift peaks by the mass of the PTM.
        WARNING: Currently written for only ONE PTM per peptide
        """
        ModIndices = Peptide.Modifications.keys()
        ModList = Peptide.Modifications[ModIndices[0]]
        FirstMod = ModList[0]
        ModMass = FirstMod.Mass
        ModIndex = ModIndices[0] + 1
        PeptideLength = len(Peptide.Aminos)
        for Peak in Spectrum.Peaks:
            if self.MaxIntensityRank != None and Peak.IntensityRank > self.MaxIntensityRank:
                continue
            Bin = int(round(Peak.Mass * self.BinMultiplier)) # default
            if Peak.IonType:
                if Peak.IonType.Name in ("b", "b-h2o", "b-nh3", "b-h2o-h2o", "b-h2o-nh3", "a"):
                    if Peak.PeptideIndex >= ModIndex:
                        Bin = int(round((Peak.Mass - ModMass) * self.BinMultiplier))
                        if VerboseFlag:
                            print "Peak at %s is %s %s, shift left to %s"%(Peak.Mass, Peak.IonType.Name, Peak.PeptideIndex, Bin)
                if Peak.IonType.Name in ("b2",):
                    if Peak.PeptideIndex >= ModIndex:
                        Bin = int(round((Peak.Mass - ModMass/2.0) * self.BinMultiplier))
                        if VerboseFlag:
                            print "Peak at %s is %s %s, shift halfleft to %s"%(Peak.Mass, Peak.IonType.Name, Peak.PeptideIndex, Bin)
                if Peak.IonType.Name in ("y", "y-h2o", "y-nh3", "y-h2o-nh3", "y-h2o-h2o"):
                    if (PeptideLength - Peak.PeptideIndex) < ModIndex:
                        Bin = int(round((Peak.Mass - ModMass) * self.BinMultiplier))
                        if VerboseFlag:
                            print "Peak at %s is %s %s, shift right to %s"%(Peak.Mass, Peak.IonType.Name, Peak.PeptideIndex, Bin)
                if Peak.IonType.Name in ("y2",):
                    if (PeptideLength - Peak.PeptideIndex) < ModIndex:
                        Bin = int(round((Peak.Mass - ModMass/2.0) * self.BinMultiplier))
                        if VerboseFlag:
                            print "Peak at %s is %s %s, shift halfright to %s"%(Peak.Mass, Peak.IonType.Name, Peak.PeptideIndex, Bin)
            Value = Peak.Intensity
            if self.HashByRank:
                Value = 10.0 / (10 + Peak.IntensityRank)                        
            Hash[Bin] = Hash.get(Bin, 0) + Value
    def HashPeaks(self, Hash, Spectrum):
        """
        Hashes Peak intensities into integer bins
        """
        for Peak in Spectrum.Peaks:
            if self.MaxIntensityRank != None and Peak.IntensityRank > self.MaxIntensityRank:
                continue # only deal with the good peaks
            Bin = int(round(Peak.Mass * self.BinMultiplier))
            Value = Peak.Intensity
            if self.HashByRank:
                Value = 10.0 / (10 + Peak.IntensityRank)
            Hash[Bin] = Hash.get(Bin, 0) + Value
    def GetSharedPeakCount(self, PeakWeight, RankWeight, SkewMultiplier = 1.0,
            PeakCountDivisor = 40, EnableShift = 1, VerboseFlag = 0):
        """
        Measure the shared peak count between two spectra.
        Iterate over the top N peaks of spectrum A (where N = (ParentMass / 100)*5).
        If the peak is present (modulo epsilon) as one of the top N peaks in spectrum B,
        then receive score PeakWeight + RankWeight / IntensityRank.  The sum of these
        scores is then scaled relative to the maximum attainable score.
        """
        SkewMultipliers = []
        for X in range(5):
            SkewMultipliers.append(SkewMultiplier ** X)
        self.BinMultiplier = 1.0
        SortedPeaksA = []
        for Peak in self.SpectrumA.Peaks:
            SortedPeaksA.append((Peak.Intensity, Peak))
        SortedPeaksA.sort()
        SortedPeaksA.reverse()
        N = int(round(self.SpectrumA.ParentMass / float(PeakCountDivisor)))
        #print "N = %s/%s = %s"%(self.SpectrumA.ParentMass, PeakCountDivisor, N)
        self.SpectrumB.RankPeaksByIntensity()
        # Populate HashedPeaks[Bin] = list of peaks in SpectrumB
        # which fall into Bin.  Only peaks
        # with rank <= N are hashed.
        HashedPeaks = {}
        PeptideLength = len(self.PeptideB.Aminos)
        ModIndices = self.PeptideB.Modifications.keys()
        if len(ModIndices):
            ModIndex = ModIndices[0] + 1
            ModMass = self.PeptideB.Modifications[ModIndices[0]][0].Mass
        else:
            ModIndex = None
        for Peak in self.SpectrumB.Peaks:
            if Peak.IntensityRank > N:
                continue
            Peak.ShiftedMass = Peak.Mass
            if EnableShift and ModIndex != None and Peak.IonType:
                if Peak.IonType.Name in ("b", "b-h2o", "b-nh3", "b-h2o-h2o", "b-h2o-nh3", "a"):
                    if Peak.PeptideIndex >= ModIndex:
                        Peak.ShiftedMass -= ModMass
                if Peak.IonType.Name in ("b2",):
                    if Peak.PeptideIndex >= ModIndex:
                        Peak.ShiftedMass -= (ModMass / 2.0)
                if Peak.IonType.Name in ("y", "y-h2o", "y-nh3", "y-h2o-nh3", "y-h2o-h2o"):
                    if (PeptideLength - Peak.PeptideIndex) < ModIndex:
                        Peak.ShiftedMass -= ModMass
                if Peak.IonType.Name in ("y2",):
                    if (PeptideLength - Peak.PeptideIndex) < ModIndex:
                        Peak.ShiftedMass -= (ModMass / 2.0)
            Bin = int(round(Peak.ShiftedMass * self.BinMultiplier))
            if not HashedPeaks.has_key(Bin):
                HashedPeaks[Bin] = []
            HashedPeaks[Bin].append(Peak)
##        ########################
##        # Debug peak hashing:
##        Keys = HashedPeaks.keys()
##        Keys.sort()
##        for Key in Keys:
##            Str = "%s: "%Key
##            for Peak in HashedPeaks[Key]:
##                Str += "(#%d %.1f, %.2f)"%(Peak.IntensityRank, Peak.Mass, Peak.Intensity)
##            print Str
##        ########################        
        OverallScore = 0
        MaxScore = 0
        PeptideLength = len(self.PeptideA.Aminos)
        ModIndices = self.PeptideA.Modifications.keys()
        if len(ModIndices):
            ModIndex = ModIndices[0] + 1
            ModMass = self.PeptideA.Modifications[ModIndices[0]][0].Mass
        else:
            ModIndex = None
        for PeakIndex in range(min(N, len(SortedPeaksA))):
            Peak = SortedPeaksA[PeakIndex][1]
            Peak.ShiftedMass = Peak.Mass
            if EnableShift and ModIndex != None and Peak.IonType:
                if Peak.IonType.Name in ("b", "b-h2o", "b-nh3", "b-h2o-h2o", "b-h2o-nh3", "a"):
                    if Peak.PeptideIndex >= ModIndex:
                        Peak.ShiftedMass -= ModMass
                if Peak.IonType.Name in ("b2",):
                    if Peak.PeptideIndex >= ModIndex:
                        Peak.ShiftedMass -= (ModMass / 2.0)
                if Peak.IonType.Name in ("y", "y-h2o", "y-nh3", "y-h2o-nh3", "y-h2o-h2o"):
                    if (PeptideLength - Peak.PeptideIndex) < ModIndex:
                        Peak.ShiftedMass -= ModMass
                if Peak.IonType.Name in ("y2",):
                    if (PeptideLength - Peak.PeptideIndex) < ModIndex:
                        Peak.ShiftedMass -= (ModMass / 2.0)
            Bin = int(round(Peak.ShiftedMass * self.BinMultiplier))
            BestPeak = None
            BestScore = 0
            for NearBin in (Bin - 1, Bin, Bin + 1):
                PeakList = HashedPeaks.get(NearBin, [])
                for PeakB in PeakList:
                    Skew = abs(Peak.ShiftedMass - PeakB.ShiftedMass)
                    SkewDeciDaltons = int(round(abs(Skew) * 10))
                    if SkewDeciDaltons >= 5:
                        continue
                    #Score = PeakWeight / float(RankWeight + Peak.IntensityRank + PeakB.IntensityRank)
                    Score = PeakWeight*10 + RankWeight*10 / (10.0 + Peak.IntensityRank)
                    Score *= SkewMultipliers[SkewDeciDaltons]
                    if Score > BestScore:
                        BestScore = Score
                        BestPeak = PeakB
                        BestPeakScoreMultiplier = SkewMultipliers[SkewDeciDaltons]
            if VerboseFlag:
                Str = "PeakA #%d %.1f (bin %d):"%(Peak.IntensityRank, Peak.Mass, Bin)
                if Peak.IonType:
                    Str += " (%s %s)"%(Peak.IonType.Name, Peak.PeptideIndex)
                print Str
                if BestPeak:
                    Str = "  Best near peak #%d %.1f ==> %s"%(BestPeak.IntensityRank, BestPeak.Mass, BestScore)
                    if BestPeak.IonType:
                        Str += " (%s %s)"%(BestPeak.IonType.Name, BestPeak.PeptideIndex)
                    print Str
            OverallScore += BestScore
            MaxScore += PeakWeight*10 + RankWeight*10 / (10.0 + Peak.IntensityRank)
            #MaxScore += PeakWeight / float(RankWeight + Peak.IntensityRank + Peak.IntensityRank)
        return OverallScore / float(MaxScore)
    def ComputeCorrelationCoefficient(self, BinMultiplier = 1.0, MaxIntensityRank = 50, EnableShift = 1, VerboseFlag = 0, HashByRank = 0):
        """
        Compute similarity between two spectra by computing the
        correlation coefficient of the binned intensities.
        """
        self.BinMultiplier = BinMultiplier
        self.MaxIntensityRank = MaxIntensityRank
        self.HashByRank = HashByRank
        # set up hashes
        HashA = {}
        HashB = {}
        # Populate HashA and HashB:
        if EnableShift and len(self.PeptideA.Modifications) > 0:
            self.HashPeaksWithShift(HashA, self.SpectrumA, self.PeptideA, VerboseFlag)
        else:
            self.HashPeaks(HashA, self.SpectrumA)
        if EnableShift and len(self.PeptideB.Modifications) > 0:
            self.HashPeaksWithShift(HashB, self.SpectrumB, self.PeptideB, VerboseFlag)
        else:
            self.HashPeaks(HashB, self.SpectrumB)
        MinBin = min(HashA.keys())
        MinBin = min(MinBin, min(HashB.keys()))
        MaxBin = max(HashA.keys())
        MaxBin = max(MaxBin, max(HashB.keys()))
        TotalA = 0
        TotalB = 0
        BinCount = MaxBin - MinBin + 1
        for Bin in range(MinBin, MaxBin + 1):
            A = HashA.get(Bin, 0)
            B = HashB.get(Bin, 0)
            TotalA += A
            TotalB += B
        MeanA = TotalA / float(BinCount)
        MeanB = TotalB / float(BinCount)
        if VerboseFlag:
            print "MeanA %s over %s bins"%(MeanA, BinCount)
            print "MeanB %s over %s bins"%(MeanB, BinCount)
        SigmaSumA = 0
        SigmaSumB = 0
        for Bin in range(MinBin, MaxBin + 1):
            A = HashA.get(Bin, 0)
            B = HashB.get(Bin, 0)
            SigmaSumA += (A - MeanA)**2
            SigmaSumB += (B - MeanB)**2
        VarianceA = SigmaSumA / float(BinCount)
        SigmaA = math.sqrt(VarianceA)
        VarianceB = SigmaSumB / float(BinCount)
        SigmaB = math.sqrt(VarianceB)
        if VerboseFlag:
            print "A has variance %s stddev %s"%(VarianceA, SigmaA)
            print "B has variance %s stddev %s"%(VarianceB, SigmaB)
        CovarianceSum = 0
        for Bin in range(MinBin, MaxBin + 1):
            A = HashA.get(Bin, 0)
            B = HashB.get(Bin, 0)
            CovarianceSum += (A - MeanA) * (B - MeanB)
        Covariance = CovarianceSum / float(BinCount - 1)
        CorrelationCoefficient = Covariance / (SigmaA * SigmaB)
        if VerboseFlag:
            print "Covariance %s, corr.coeff %s"%(Covariance, CorrelationCoefficient)
        return CorrelationCoefficient

def Test():
    FileName = "..\mzxml\Dicty-HeavyCells-11.mzXML"
    FileHandle = open(FileName,"rb")
    FileHandle.seek(60145805)
    S1 = MSSpectrum.SpectrumClass()
    S1.ReadPeaksFromFile(FileHandle,FileName) #it also sets the PrecursorMZ
    Annotation1 = "MKRKLLK"

    FileName = "..\mzxml\Dicty-HeavyCells-13.mzXML"
    FileHandle = open(FileName,"rb")
    FileHandle.seek(59307113)
    S2 = MSSpectrum.SpectrumClass()
    S2.ReadPeaksFromFile(FileHandle,FileName) #it also sets the PrecursorMZ
    Annotation2 = "MKRKLLK"

    FileName = "..\mzxml\Dicty-HeavyCells-12.mzXML"
    FileHandle = open(FileName,"rb")
    FileHandle.seek(87683754)
    S3 = MSSpectrum.SpectrumClass()
    S3.ReadPeaksFromFile(FileHandle,FileName) #it also sets the PrecursorMZ
    Annotation3 = "MKIFIIK"

    FileName = "..\mzxml\Dicty-HeavyCells-05.mzXML"
    FileHandle = open(FileName,"rb")
    FileHandle.seek(102201432)
    S5 = MSSpectrum.SpectrumClass()
    S5.ReadPeaksFromFile(FileHandle,FileName) #it also sets the PrecursorMZ
    #S5.FilterPeaks()
    Annotation5 = "NWNGQPVGVPQGQYANMNYAR"   

    FileName = "..\mzxml\Dicty-HeavyCells-12.mzXML"
    FileHandle = open(FileName,"rb")
    FileHandle.seek(112303085)
    S6 = MSSpectrum.SpectrumClass()
    S6.ReadPeaksFromFile(FileHandle,FileName) #it also sets the PrecursorMZ
    Annotation6 = "NWNGQPVGVPQGQYANMNYAR+14"   

    FileName = "..\mzxml\Dicty-HeavyCells-12.mzXML"
    FileHandle = open(FileName,"rb")
    FileHandle.seek(111782847)
    S7 = MSSpectrum.SpectrumClass()
    S7.ReadPeaksFromFile(FileHandle,FileName) #it also sets the PrecursorMZ
    Annotation7 = "NWNGQPVGVPQGQYANMNYAR"   

    Simm = SpectralSimilarity(S5,S5,Annotation5,Annotation5)
    Simm.SpectralAlignment()

if __name__ == "__main__":
    # Command-line arguments: Two spectra, then two annotations.
    if len(sys.argv)<5:
        print "Not enough arguments:", sys.argv
        sys.exit(-1)
    Comparator = SpectralSimilarity(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    Comparator.LabelPeaks()
    Result = Comparator.GetSharedPeakCount(1, 0, 0.9, PeakCountDivisor = 20, VerboseFlag = 1)
    print "Shared 1 0 0.8:", Result
    print "\n\n"
    Result = Comparator.GetSharedPeakCount(0, 1, 0.9, PeakCountDivisor = 20, VerboseFlag = 1)
    print "Shared 0 1 0.8:", Result
    print "\n\n"
    
##    Result = Comparator.GetSharedPeakCount(0, 1, PeakCountDivisor = 20, VerboseFlag = 1)
##    print "Shared 0 1 1.0:", Result
##    sys.exit(0)
##    print "Shared 0 1 0.66:", Comparator.GetSharedPeakCount(0, 1, 0.66, PeakCountDivisor = 5)
##    print "Shared 0 1 0.66:", Comparator.GetSharedPeakCount(0, 1, 0.66, PeakCountDivisor = 50)
##    
##    Result = Comparator.DotProduct(VerboseFlag = 1)
##    print "Dot product similarity score 1.0:", Result
##    Result = Comparator.DotProduct(2.0, VerboseFlag = 1)
##    print "Dot product similarity score 2.0:", Result
    Result = Comparator.DotProduct(0.5, VerboseFlag = 1, HashByRank = 1)
    print "Dot product similarity score 0.5:", Result
##    print "Shared 1 0:", Comparator.GetSharedPeakCount(1, 0)
##    print "Shared 0 1:", Comparator.GetSharedPeakCount(0, 1)
##    
##    print "Shared 1 1:", Comparator.GetSharedPeakCount(1, 1)
##    print "Cov/corr:", Comparator.ComputeCorrelationCoefficient()
    Command = "label.py \"%s\" %s"%(sys.argv[1], sys.argv[3])
    os.system(Command)
    Command = "label.py \"%s\" %s"%(sys.argv[2], sys.argv[4])
    os.system(Command)
