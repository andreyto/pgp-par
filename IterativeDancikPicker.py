"""
Iterative Dancik offset-frequency-function maker.  Perform this iteration relative to b series,
y series, b2 series, y2 series; do it separately for charge 2 and for charge 3:
- Start with an emtpy set of offsets X
- Iterate over a corpus of annotated spectra.  Compute a histogram of peak offsets from points in the
series, omitting any peaks accounted for by an offset in the set X.
- Find the most abundant (by count or intensity) offset. Add this offset to the set X.

The goals of this script are:
- Refine the set of ions used by ScorpionTrainer / Scorpion
- Generate a nice graph for the dissertation
"""
import os
import sys
import traceback
import math
import random
import getopt
import MSSpectrum
from Utils import *
Initialize()
HYDROGEN = 1.0078
BIN_DIVISOR = 10.0

class OffsetClass:
    def __init__(self):
        self.Mass = None
        self.Score = None
        self.PeakCount = None
        self.TotalIntensity = None
        self.Type = None # "prefix", "prefix2", "suffix", "suffix2"
    def __str__(self):
        return "<Offset %s %s>"%(self.Type, self.Mass)
OffsetTypes = ["prefix", "prefix2", "suffix", "suffix2"]

class FrequencyFunctionMaker:
    def __init__(self):
        self.MaxOffsetCount = 50
        self.CountHistogram = {}
        self.IntensityHistogram = {}
        self.OffsetDict = {}
        for OffsetType in OffsetTypes:
            self.OffsetDict[OffsetType] = {}
        self.OffsetList = []
        self.TotalPeakCount = 0
        self.TotalPeakIntensity = 0
        self.Charge = None
        self.CorpusFileName = None
        self.QuickParseFlag = 0
        self.UnfilteredTotalIntensity = None
        self.UnfilteredTotalCount = None
    def ParseCommandLine(self, Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:c:s:Q")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-r":
                self.CorpusFileName = Value
            elif Option == "-c":
                self.Charge = int(Value)
            elif Option == "-Q":
                self.QuickParseFlag = 1
                
    def Main(self):
        if not self.Charge or not self.CorpusFileName:
            print UsageInfo
            sys.exit(-1)
        # OffsetList is a list of OffsetClass instances, ordered from best to worst:
        self.OffsetList = []
        # OffsetDict[OffsetType] maps a mass bin to the best OffsetClass within 0.3Da:
        self.OffsetDict = {}
        for OffsetType in OffsetTypes:
            self.OffsetDict[OffsetType] = {}
        self.OverallOutputFile = open("IDPOutput.%s.txt"%(self.Charge), "wb")
        while len(self.OffsetList) < self.MaxOffsetCount:
            self.CountHistogram = {}
            for OffsetType in OffsetTypes:
                self.CountHistogram[OffsetType] = {}
            self.IntensityHistogram = {}
            for OffsetType in OffsetTypes:
                self.IntensityHistogram[OffsetType] = {}
            self.TotalPeakCount = 0
            self.TotalPeakIntensity = 0
            self.SpectraProcessed = 0
            self.CountUnexplainedPeaks()
            print "Unexplained peaks at offset 0:"
            for OffsetType in OffsetTypes:
                Count = self.CountHistogram[OffsetType].get(0, 0)
                Intensity = self.IntensityHistogram[OffsetType].get(0, 0)
                print "%s: Count %s intensity %s"%(OffsetType, Count, Intensity)
            OffsetPlotFileName = "Offsets.%s.%s.txt"%(self.Charge, len(self.OffsetList))
            self.PlotOffsets(OffsetPlotFileName)
            self.SelectNewOffset()
    def PlotOffsets(self, OffsetPlotFileName):
        File = open(OffsetPlotFileName, "wb")
        HeaderStr = "Bin\tMass\t"
        for OffsetType in OffsetTypes:
            HeaderStr += "%sCount\t%sIntensity\t"%(OffsetType, OffsetType)
        File.write(HeaderStr + "\n")
        AllKeys = {}
        for OffsetType in OffsetTypes:
            for Key in self.CountHistogram[OffsetType].keys():
                AllKeys[Key] = 1
        Keys = AllKeys.keys()
        Keys.sort()
        for Key in Keys:
            Str = "%s\t%s\t"%(Key, Key / BIN_DIVISOR)
            for OffsetType in OffsetTypes:
                Str += "%s\t%s\t"%(self.CountHistogram[OffsetType].get(Key, 0),
                                   self.IntensityHistogram[OffsetType].get(Key, 0))
            File.write(Str + "\n")
        File.close()
    def AddOffset(self, BestOffsetBin, BestOffsetType):
        print "\n\n>>>NEW OFFSET!"
        Str = "Best new offset type %s found at %s"%(BestOffsetType, BestOffsetBin / BIN_DIVISOR)
        self.OverallOutputFile.write(Str + "\n")
        print Str
        Count = self.CountHistogram[BestOffsetType].get(BestOffsetBin, 0)
        Percent = 100 * Count / float(max(1, self.TotalPeakCount))
        Str = "This offset explains %s (%.2f%%) of remaining peaks"%(Count, Percent)
        self.OverallOutputFile.write(Str + "\n")
        print Str
        Intensity = self.IntensityHistogram[BestOffsetType].get(BestOffsetBin, 0)
        Percent = 100 * Intensity / float(max(1, self.TotalPeakIntensity))
        Str = "This offset explains %s (%.2f%%) of remaining intensity"%(Intensity, Percent)
        self.OverallOutputFile.write(Str + "\n")
        print Str
        CuteStr = "Offset\t%s\t%s\t%s\t"%(len(self.OffsetList), BestOffsetType, BestOffsetBin)
        CuteStr += "%s\t%s\t"%(Count, Intensity)
        if self.UnfilteredTotalCount == None:
            self.UnfilteredTotalCount = self.TotalPeakCount
        if self.UnfilteredTotalIntensity == None:
            self.UnfilteredTotalIntensity = self.TotalPeakIntensity            
        CuteStr += "%s\t%s\t%s\t"%(BestOffsetBin / BIN_DIVISOR, Count / float(self.UnfilteredTotalCount), Intensity / self.UnfilteredTotalIntensity)
        self.OverallOutputFile.write(CuteStr + "\n")
        Offset = OffsetClass()
        Offset.Type = BestOffsetType
        Offset.Mass = BestOffsetBin / BIN_DIVISOR
        Offset.Bin = BestOffsetBin
        Offset.PeakCount = Count
        Offset.TotalIntensity = Intensity
        self.OffsetList.append(Offset)
        # Update self.OffsetDict:
        for NearBin in range(BestOffsetBin - 2, BestOffsetBin + 3):
            if not self.OffsetDict[BestOffsetType].has_key(NearBin):
                self.OffsetDict[BestOffsetType][NearBin] = Offset
        self.OverallOutputFile.flush()        
    def SelectNewOffset(self):
        BestScore = None
        BestOffsetBin = None
        BestOffsetType = None
        # Compute a score based on total intensity and on total count.
        for OffsetType in OffsetTypes:
            Bins = self.CountHistogram[OffsetType].keys()
            for Bin in Bins:
                if self.OffsetDict[OffsetType].has_key(Bin):
                    continue # Too close to an already-picked bin
                CountFraction = (self.CountHistogram[OffsetType][Bin] / float(max(1, self.TotalPeakCount)))
                IntensityFraction = (self.IntensityHistogram[OffsetType][Bin] / float(max(1, self.TotalPeakIntensity)))
                Score = CountFraction + IntensityFraction
                if Score > BestScore:
                    BestOffsetBin = Bin
                    BestScore = Score
                    BestOffsetType = OffsetType
                    print "Best so far: %s for %s %s"%(BestScore, BestOffsetBin, BestOffsetType)
        self.AddOffset(BestOffsetBin, BestOffsetType)
    def CountUnexplainedPeaks(self):
        """
        Iterate over file(s) in the corpus, and count the peaks.
        """
        if os.path.isdir(self.CorpusFileName):
            for FileName in os.listdir(self.CorpusFileName):
                Path = os.path.join(self.CorpusFileName, FileName)
                self.CountUnexplainedPeaksFromFile(Path)
        else:
            self.CountUnexplainedPeaksFromFile(self.CorpusFileName)
    def CountUnexplainedPeaksFromSpectrum(self, SpectrumPath, FilePos, Peptide):
        self.SpectraProcessed += 1
        # Construct the spectrum:
        Spectrum = MSSpectrum.SpectrumClass()
        Spectrum.ReadPeaks("%s:%s"%(SpectrumPath, FilePos))
        Spectrum.FilterPeaks()
        # We'll scale peak intensities relative to the spectrum's total intensity:
        TotalIntensity = 0
        for Peak in Spectrum.Peaks:
            TotalIntensity += Peak.Intensity
        # Get the anchor points for the series:
        AnchorPoints = {}
        PrefixPoints = {}
        SuffixPoints = {}
        for OffsetType in OffsetTypes:
            AnchorPoints[OffsetType] = self.GetAnchorPoints(Peptide, OffsetType)
        # Reckon the offsets of peaks relative to these anchor points.
        for Peak in Spectrum.Peaks:
            # Determine whether this peak is explained by any offset:
            ExplainedFlag = 0
            for OffsetType in OffsetTypes:
                PrefixPoints[OffsetType] = None
                SuffixPoints[OffsetType] = None
                for AnchorPoint in AnchorPoints[OffsetType]:
                    if AnchorPoint < Peak.Mass and (PrefixPoints[OffsetType] == None or PrefixPoints[OffsetType] < AnchorPoint):
                        PrefixPoints[OffsetType] = AnchorPoint
                    if AnchorPoint > Peak.Mass and (SuffixPoints[OffsetType] == None or SuffixPoints[OffsetType] > AnchorPoint):
                        SuffixPoints[OffsetType] = AnchorPoint
                if PrefixPoints[OffsetType] != None:
                    OffsetPrefix = Peak.Mass - PrefixPoints[OffsetType]
                    BinPrefix = int(round(OffsetPrefix * BIN_DIVISOR))
                    if self.OffsetDict[OffsetType].has_key(BinPrefix):
                        #Offset = self.OffsetDict[OffsetType][BinPrefix]
                        #print "Peak %s in %s explained as %s from %s"%(Peak.Mass, SpectrumPath, Offset, Peptide)
                        ExplainedFlag = 1
                        break
                if SuffixPoints[OffsetType] != None:
                    OffsetSuffix = Peak.Mass - SuffixPoints[OffsetType]
                    BinSuffix = int(round(OffsetSuffix * BIN_DIVISOR))
                    if self.OffsetDict[OffsetType].has_key(BinSuffix):
                        #Offset = self.OffsetDict[OffsetType][BinSuffix]
                        #print "Peak %s in %s explained as %s from %s"%(Peak.Mass, SpectrumPath, Offset, Peptide)
                        ExplainedFlag = 1
                        break
            if ExplainedFlag:
                continue
            ###########################################################################
            # Unexplained peak: Add to the histograms!
            ScaledIntensity = Peak.Intensity / float(max(0.0001, TotalIntensity))
            for OffsetType in OffsetTypes:
                if PrefixPoints[OffsetType] != None:
                    OffsetPrefix = Peak.Mass - PrefixPoints[OffsetType]
                    BinPrefix = int(round(OffsetPrefix * BIN_DIVISOR))
                    for NearBin in range(BinPrefix - 2, BinPrefix + 3):
                        self.CountHistogram[OffsetType][NearBin] = self.CountHistogram[OffsetType].get(NearBin, 0) + 1
                        self.IntensityHistogram[OffsetType][NearBin] = self.IntensityHistogram[OffsetType].get(NearBin, 0) + ScaledIntensity
                if SuffixPoints[OffsetType] != None:
                    OffsetSuffix = Peak.Mass - SuffixPoints[OffsetType]
                    BinSuffix = int(round(OffsetSuffix * BIN_DIVISOR))
                    for NearBin in range(BinSuffix - 2, BinSuffix + 3):
                        if NearBin == BinSuffix:
                            Multiplier = 1.0
                        else:
                            Multiplier = 0.0
                        self.CountHistogram[OffsetType][NearBin] = self.CountHistogram[OffsetType].get(NearBin, 0) + 1
                        self.IntensityHistogram[OffsetType][NearBin] = self.IntensityHistogram[OffsetType].get(NearBin, 0) + ScaledIntensity
            self.TotalPeakIntensity += ScaledIntensity
            self.TotalPeakCount += 1
        #print "Processed %s spectra, %s intensity"%(self.SpectraProcessed, self.TotalPeakIntensity)
    def CountUnexplainedPeaksFromFile(self, FileName):
        """
        Parse spectra (and annotations) from this file.  Filter peaks, and count where the many
        not-yet-annotated peaks fall.
        """
        File = open(FileName, "rb")
        LineNumber = 0
        for FileLine in File.xreadlines():
            LineNumber += 1
            if LineNumber % 100 == 0:
                print "%s Line %s..."%(FileName, LineNumber)
                if self.QuickParseFlag:
                    break
            Bits = FileLine.split("\t")
            if FileLine[0] == "#":
                continue # header line
            try:
                Charge = int(Bits[4])
                Peptide = GetPeptideFromModdedName(Bits[2])
                SpectrumPath = self.FixFilePath(Bits[0])
                FilePos = int(Bits[15])
            except:
                print "* Can't parse line %s of file %s"%(LineNumber, FileName)
                traceback.print_exc()
                return # BAIL OUT!
            if Charge != self.Charge:
                continue
            self.CountUnexplainedPeaksFromSpectrum(SpectrumPath, FilePos, Peptide)
        File.close()
    def GetAnchorPoints(self, Peptide, OffsetType):
        """
        Get the masses of the anchor points (b series, y series, b2 series, or y2 series)
        """
        MassList = []
        ParentMass = Peptide.Masses[-1] + 19
        for PrefixResidueMass in Peptide.Masses:
            if OffsetType == "prefix":
                Mass = PrefixResidueMass + HYDROGEN
            elif OffsetType == "prefix2":
                Mass = (PrefixResidueMass + HYDROGEN + HYDROGEN) / 2.0
            elif OffsetType == "suffix":
                Mass = (ParentMass - PrefixResidueMass)
            elif OffsetType == "suffix2":
                Mass = (ParentMass - PrefixResidueMass + HYDROGEN) / 2.0
            MassList.append(Mass)
        return MassList
    def FixFilePath(self, FilePath):
        return FilePath
        #return r"F:\ftproot\briggs\hek293\H293b-total-try-2nd-digest-c-500ug-2D34-121905-LTQ1\H293b-total-try-2nd-digest-c-500ug-2D34-121905-LTQ1-17.mzXML"
    
UsageInfo = """
IterativeDancikPicker arguments:
-r [FileName]: Corpus of annotations
-c [Charge]: Charge state to handle

Example:
 IterativeDancikPicker.py -r TrainingCorpus -c 2
  
"""
        

if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(no psyco)"
    Dancik = FrequencyFunctionMaker()
    Dancik.ParseCommandLine(sys.argv[1:])
    Dancik.Main()
        