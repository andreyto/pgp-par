"""
This script is designed to discover the characteristics of TRUE TAGS.  A *tag* is a set of
four PRMs, separated by amino acid masses, which can be interpreted as a tripeptide.  True tags
are those which correspond to the actual protein being examined by mass-spec.  We use
the PeptideClass is to represent a tag.

Here's the overall plan:
- Sort all spectra by total ion current.  Limit attention to the most intense half.
- For each of these spectra, filter the peaks (using a 50/6 sliding window filter).  Iterate
over true tags, and look for each one where all four nodes are represented by any b or y peak.
For each such "capturable tag", here's the sort of information we want to know:
  -- Probability that peak #n (sorted by intensity) is a b or a y ion in such a tag
  -- Frequency of different forms of ion support (b+y, or b+y+y-h2o, or whatever)
  -- Histogram of node-to-node skew (rounded to nearest hundredth of an amu)
  -- Histogram of overall skew
  -- Odds that the same (with 1.0) peak supplies both a b and a y peak to this tag
  -- Odds of breakage before and after amino acids (for now, just P and G)
"""
from Utils import *
import traceback
import types
import cPickle
import math
import MSSpectrum
import struct
import Label
import MakeImage

# List of node
NodeTypes = []
for Charge in range(1, 4):
    for Tri in range(3):
        NodeType = (Charge, Tri)
        NodeTypes.append(NodeType)


# It's assumed that SpectrumDir is present and includes a large-enough training set.
if IS_WINDOWS:
    #SpectrumDir = r"c:\source\bafna\spectra"
    #SpectrumDir = r"d:\research\bafna\spectra"
    SpectrumDir = "Spectra"
    HappyDir = r"d:\research\bafna\HappySpectraPMC"
else:
    SpectrumDir = "Spectra"
    HappyDir = "HappySpectraPMC"

def FindIntenseSpectra():
    """
    Print a list of the spectra in SpectrumDir, ordered by charge and by intensity.
    """
    Listing = []
    for FileName in os.listdir(SpectrumDir):
        FilePath = os.path.join(SpectrumDir, FileName)
        try:
            Spectrum = MSSpectrum.SpectrumClass()
            Spectrum.ReadPeaks(FilePath)
            #Spectrum.RankWithFilter(50, 6)
            Intensity = 0
            for Peak in Spectrum.Peaks:
                Intensity += Peak.Intensity
            Listing.append((Spectrum.Charge, Intensity, FileName))
        except:
            traceback.print_exc()
    Listing.sort()
    for Tuple in Listing:
        print "%s\t%s\t%s\t"%Tuple

class StatMeister:
    "Main singleton class for accumulating parameters from the training set."
    MinimumWitnessSampleSize = 10
    def __init__(self, PeptideOracle):
        self.TotalCount = {}
        self.TotalTrue = {}
        self.WitnessCount = {}
        self.WitnessTrue = {}
        self.BRankCount = {}
        self.BRankTrue = {}
        self.YRankCount = {}
        self.YRankTrue = {}
        self.IsotopeCount = {}
        self.IsotopeTrue = {}
        #
        self.PeptideOracle = PeptideOracle
        self.Epsilon = 0.5
        ##
        self.TotalTags = 0 
        self.TotalSpectra = 0
        ###
        self.SkewHistogram = {}
        self.TotalSkew = 0
        self.TotalAbsSkew = 0
        self.SkewCount = 0
        self.TotalSkewHistogram = {}
        self.TotalTagSkew = 0
        self.TotalAbsTagSkew = 0
        self.TagSkewCount = 0
        ###
        self.BIonAssignments = {} # rank -> count
        self.YIonAssignments = {} # rank -> count
        self.BYIonAssignments = {} # rank -> count; peak works as a B and as a Y ion
        self.TotalIonAssignments = {} # rank -> overall count
        ####
        self.WitnessSetsB = {}
        self.WitnessSetsBTrue = {}
        self.WitnessSetsStrictB = {}
        self.WitnessSetsStrictBTrue = {}
        self.TriWitnessCountB = {"Hi":0, "Med":0, "Lo":0}
        self.TriWitnessSetsB = {"Hi":{}, "Med":{}, "Lo":{}}
        self.TriWitnessSetsBTrue = {"Hi":{}, "Med":{}, "Lo":{}}
        self.TriWitnessSetsStrictB = {"Hi":{}, "Med":{}, "Lo":{}}
        self.TriWitnessSetsStrictBTrue = {"Hi":{}, "Med":{}, "Lo":{}}
        self.WitnessSetsY = {}
        self.WitnessSetsYTrue = {}
        self.WitnessSetsStrictY = {}
        self.WitnessSetsStrictYTrue = {}
        self.TriWitnessCountY = {"Hi":0, "Med":0, "Lo":0}
        self.TriWitnessSetsY = {"Hi":{}, "Med":{}, "Lo":{}}
        self.TriWitnessSetsYTrue = {"Hi":{}, "Med":{}, "Lo":{}}
        self.TriWitnessSetsStrictY = {"Hi":{}, "Med":{}, "Lo":{}}
        self.TriWitnessSetsStrictYTrue = {"Hi":{}, "Med":{}, "Lo":{}}
        
        self.WitnessSetCount = 0
        ##########
        self.IsotopicB = {}
        self.IsotopicY = {}
        self.IsotopicAll = {}
        self.IsotopicEvil = {}
        ##########
        self.YPeakFound = {} # key: position, or None for overall; value = count of all the times a y peak was present
        self.YPeakCount = {} # incremented every time we look for a YPeakFound element
        self.BPeakFound = {}
        self.BPeakCount = {}
        self.PreProlineYPeakFound = {}
        self.PreProlineYPeakCount = {}
        self.PostProlineYPeakFound = {}
        self.PostProlineYPeakCount = {}
        self.PreGlycineYPeakFound = {}
        self.PreGlycineYPeakCount = {}
        self.PostGlycineYPeakFound = {}
        self.PostGlycineYPeakCount = {}
        self.PreProlineBPeakFound = {}
        self.PreProlineBPeakCount = {}
        self.PostProlineBPeakFound = {}
        self.PostProlineBPeakCount = {}
        self.PreGlycineBPeakFound = {}
        self.PreGlycineBPeakCount = {}
        self.PostGlycineBPeakFound = {}
        self.PostGlycineBPeakCount = {}
        self.B2TrueCount = {}
        self.B2TotalTrueCount = {}
        self.Y2TrueCount = {}
        self.Y2TotalTrueCount = {}
        self.ResetStats()
    def ReportProlineStatsHelper(self, Title, FoundDict, CountDict):
        print Title
        Keys = CountDict.keys()
        Keys.sort()
        for Key in Keys:
            Count = CountDict[Key]
            Found = FoundDict.get(Key, 0)
            print "%s: %d of %d (%.2f%%)"%(Key, Found, Count, 100*Found/float(Count))
        
    def ReportProlineStats(self):
        self.ReportProlineStatsHelper("B peaks", self.BPeakFound, self.BPeakCount)
        self.ReportProlineStatsHelper("Y peaks", self.YPeakFound, self.YPeakCount)
        self.ReportProlineStatsHelper("B peaks pre proline", self.PreProlineBPeakFound, self.PreProlineBPeakCount)
        self.ReportProlineStatsHelper("B peaks post proline", self.PostProlineBPeakFound, self.PostProlineBPeakCount)
        self.ReportProlineStatsHelper("Y peaks pre proline", self.PreProlineYPeakFound, self.PreProlineYPeakCount)
        self.ReportProlineStatsHelper("Y peaks post proline", self.PostProlineYPeakFound, self.PostProlineYPeakCount)
        self.ReportProlineStatsHelper("B peaks pre Glycine", self.PreGlycineBPeakFound, self.PreGlycineBPeakCount)
        self.ReportProlineStatsHelper("B peaks post Glycine", self.PostGlycineBPeakFound, self.PostGlycineBPeakCount)
        self.ReportProlineStatsHelper("Y peaks pre Glycine", self.PreGlycineYPeakFound, self.PreGlycineYPeakCount)
        self.ReportProlineStatsHelper("Y peaks post Glycine", self.PostGlycineYPeakFound, self.PostGlycineYPeakCount)
        
    def AccumulateProlineIons(self, Spectrum, Peptide):
        Mass = 0
        for AminoIndex in range(len(Peptide.Aminos)):
            YAminoIndex = len(Peptide.Aminos) - AminoIndex - 1
            # Look for b/y peaks for this cut point
            Mass += Global.AminoMass[Peptide.Aminos[AminoIndex]]
            Peak = Spectrum.GetPeak(Mass + 1.0078)
            self.BPeakCount[None] = self.BPeakCount.get(None, 0) + 1
            self.BPeakCount[AminoIndex] = self.BPeakCount.get(AminoIndex, 0) + 1
            if Peak:
                self.BPeakFound[None] = self.BPeakFound.get(None, 0) + 1
                self.BPeakFound[AminoIndex] = self.BPeakFound.get(AminoIndex, 0) + 1
            if Peptide.Aminos[AminoIndex] == "P":
                self.PostProlineBPeakCount[None] = self.PostProlineBPeakCount.get(None, 0) + 1
                self.PostProlineBPeakCount[AminoIndex] = self.PostProlineBPeakCount.get(AminoIndex, 0) + 1
                if Peak:
                    self.PostProlineBPeakFound[None] = self.PostProlineBPeakFound.get(None, 0) + 1
                    self.PostProlineBPeakFound[AminoIndex] = self.PostProlineBPeakFound.get(AminoIndex, 0) + 1
            if AminoIndex<len(Peptide.Aminos)-1 and Peptide.Aminos[AminoIndex+1] == "P":
                self.PreProlineBPeakCount[None] = self.PreProlineBPeakCount.get(None, 0) + 1
                self.PreProlineBPeakCount[AminoIndex] = self.PreProlineBPeakCount.get(AminoIndex, 0) + 1
                if Peak:
                    self.PreProlineBPeakFound[None] = self.PreProlineBPeakFound.get(None, 0) + 1
                    self.PreProlineBPeakFound[AminoIndex] = self.PreProlineBPeakFound.get(AminoIndex, 0) + 1
            if Peptide.Aminos[AminoIndex] == "G":
                self.PostGlycineBPeakCount[None] = self.PostGlycineBPeakCount.get(None, 0) + 1
                self.PostGlycineBPeakCount[AminoIndex] = self.PostGlycineBPeakCount.get(AminoIndex, 0) + 1
                if Peak:
                    self.PostGlycineBPeakFound[None] = self.PostGlycineBPeakFound.get(None, 0) + 1
                    self.PostGlycineBPeakFound[AminoIndex] = self.PostGlycineBPeakFound.get(AminoIndex, 0) + 1
            if AminoIndex<len(Peptide.Aminos)-1 and Peptide.Aminos[AminoIndex+1] == "G":
                self.PreGlycineBPeakCount[None] = self.PreGlycineBPeakCount.get(None, 0) + 1
                self.PreGlycineBPeakCount[AminoIndex] = self.PreGlycineBPeakCount.get(AminoIndex, 0) + 1
                if Peak:
                    self.PreGlycineBPeakFound[None] = self.PreGlycineBPeakFound.get(None, 0) + 1
                    self.PreGlycineBPeakFound[AminoIndex] = self.PreGlycineBPeakFound.get(AminoIndex, 0) + 1
                    
            Peak = Spectrum.GetPeak(Spectrum.ParentMass - Mass)
            self.YPeakCount[None] = self.YPeakCount.get(None, 0) + 1
            self.YPeakCount[YAminoIndex] = self.YPeakCount.get(YAminoIndex, 0) + 1
            if Peak:
                self.YPeakFound[None] = self.YPeakFound.get(None, 0) + 1
                self.YPeakFound[YAminoIndex] = self.YPeakFound.get(YAminoIndex, 0) + 1
            if Peptide.Aminos[AminoIndex] == "P":
                self.PostProlineYPeakCount[None] = self.PostProlineYPeakCount.get(None, 0) + 1
                self.PostProlineYPeakCount[YAminoIndex] = self.PostProlineYPeakCount.get(YAminoIndex, 0) + 1
                if Peak:
                    self.PostProlineYPeakFound[None] = self.PostProlineYPeakFound.get(None, 0) + 1
                    self.PostProlineYPeakFound[YAminoIndex] = self.PostProlineYPeakFound.get(YAminoIndex, 0) + 1
            if AminoIndex<len(Peptide.Aminos)-1 and Peptide.Aminos[AminoIndex+1] == "P":
                self.PreProlineYPeakCount[None] = self.PreProlineYPeakCount.get(None, 0) + 1
                self.PreProlineYPeakCount[YAminoIndex] = self.PreProlineYPeakCount.get(YAminoIndex, 0) + 1
                if Peak:
                    self.PreProlineYPeakFound[None] = self.PreProlineYPeakFound.get(None, 0) + 1
                    self.PreProlineYPeakFound[YAminoIndex] = self.PreProlineYPeakFound.get(YAminoIndex, 0) + 1
            if Peptide.Aminos[AminoIndex] == "G":
                self.PostGlycineYPeakCount[None] = self.PostGlycineYPeakCount.get(None, 0) + 1
                self.PostGlycineYPeakCount[YAminoIndex] = self.PostGlycineYPeakCount.get(YAminoIndex, 0) + 1
                if Peak:
                    self.PostGlycineYPeakFound[None] = self.PostGlycineYPeakFound.get(None, 0) + 1
                    self.PostGlycineYPeakFound[YAminoIndex] = self.PostGlycineYPeakFound.get(YAminoIndex, 0) + 1
            if AminoIndex<len(Peptide.Aminos)-1 and Peptide.Aminos[AminoIndex+1] == "G":
                self.PreGlycineYPeakCount[None] = self.PreGlycineYPeakCount.get(None, 0) + 1
                self.PreGlycineYPeakCount[YAminoIndex] = self.PreGlycineYPeakCount.get(YAminoIndex, 0) + 1
                if Peak:
                    self.PreGlycineYPeakFound[None] = self.PreGlycineYPeakFound.get(None, 0) + 1
                    self.PreGlycineYPeakFound[YAminoIndex] = self.PreGlycineYPeakFound.get(YAminoIndex, 0) + 1
                
    def ResetStats(self):
        for NodeType in NodeTypes:
            self.TotalCount[NodeType] = 0
            self.TotalTrue[NodeType] = 0
            self.WitnessCount[NodeType] = {}
            self.WitnessTrue[NodeType] = {}
            self.BRankCount[NodeType] = {}
            self.BRankTrue[NodeType] = {}
            self.YRankCount[NodeType] = {}
            self.YRankTrue[NodeType] = {}
            self.IsotopeCount[NodeType] = {}
            self.IsotopeTrue[NodeType] = {}
            self.B2TrueCount[NodeType] = 0
            self.B2TotalTrueCount[NodeType] = 0
            self.Y2TrueCount[NodeType] = 0
            self.Y2TotalTrueCount[NodeType] = 0
    def AccumulateStats(self, SpectrumListFileName, FileNameBit, DirName, Limit = None):
        "Main method for training: Iterate over the spectra on this list, and accumulate stats for each."
        SpectrumListFile = open(SpectrumListFileName, "r")
        Count = 0
        # Iterate over all spectra.  Get stats for each:
        for FileLine in SpectrumListFile.xreadlines():
            Bits = FileLine.strip().split()
            FileName = Bits[FileNameBit]
            if not FileName:
                continue
            if not self.PeptideOracle.has_key(FileName):
                print "SKIP spectrum file with unknown peptide:", FileName
                continue
            Peptide = self.PeptideOracle[FileName]
            if type(Peptide)==type(""):
                Peptide = PeptideClass()
                Peptide.Aminos = self.PeptideOracle[FileName]
                Peptide.ComputeMasses()
            FilePath = os.path.join(DirName, FileName)
            print FileName, Peptide.Aminos
            try:
                self.AccumulateStatsFromFile(FilePath, Peptide)
            except:
                print "*"*70
                print FilePath
                traceback.print_exc()
            # Stop after n spectra, if a limit was provided:
            Count += 1
            if Limit!=None and Count>Limit:
                break
            # Reporting stats along the way:
##            if Count%5 == 0:
##                self.ReportStats("TaggingModelTemp.dat")
        # Summarize:
        SpectrumListFile.close()
        self.ReportStats("TaggingModel.dat")
        print "ALL STATS GENERATED."

    def SummarizeIonTypes(self, WitnessSets, WitnessSetsTrue, Title):
        print "Witness sets (%s)"%Title
        Keys = WitnessSets.keys()
        Keys.sort()
        for Key in Keys:
            True = WitnessSetsTrue.get(Key, 0)
            Total = WitnessSets[Key]
            Percent = True/float(Total)
            print "%s\t%s\t%s\t%s\t"%(Key, True, Total, Percent)
        print             
    def SummarizeSkewHistogram(self, Histo, Title):
        print "\n\nSummary of skew histogram %s:\n"%Title
        Keys = Histo.keys()
        Keys.sort()
        if not len(Keys):
            return
        StartKey = Keys[0]
        EndKey = Keys[-1] + 1
        for Key in range(StartKey, EndKey):
            print "%s\t%s\t"%(Key, Histo.get(Key, 0))
    def SummarizeIonAssignments(self, Assignments, AssignmentCount, Name):
        print "\n\nSummary of ion-type assignments for ion type %s:\n"%Name
        Keys = AssignmentCount.keys()
        Keys.sort()
        MaxIndex = Keys[-3] # skip left, right
        Ranks = list(range(MaxIndex))
        Ranks.append("Left")
        Ranks.append("Right")
        for RankIndex in Ranks: #range(MaxIndex):
            Assigned = Assignments.get(RankIndex, 0)
            Total = AssignmentCount.get(RankIndex, 0)
            if Total == 0:
                Percent = ""
            else:
                Percent = Assigned / float(Total)
            print "%s\t%s\t%s\t%s\t"%(RankIndex, Assigned, Total, Percent)
    def GetProbStruct(self, DictTrue, DictCount, NodeType, Key):
        True = DictTrue[NodeType].get(Key, 0)
        Count = DictCount[NodeType].get(Key, 0)
        if Count == 0:
            # No data was found at all.  So, the odds ratio is 1, and its log is 0:
            return struct.pack("<f", 0)
        if True == 0:
            if Count > 100:
                Score = -6.9 # Odds ratio is zero-ish, but we can't take ln(0), so take ln(.001)
            elif Count > 10:
                Score = -3.0 # Odds ratio is zero-ish, but we can't take ln(0), so take ln(.05)
            else:
                Score = -0.69 # Odds ratio is zero-ish, but we can't take ln(0), so take ln(.1)
        else:
            print "TC:", True, Count, self.TotalTrue[NodeType], self.TotalCount[NodeType]
            Score = math.log(True / float(Count)) - math.log(self.TotalTrue[NodeType]/float(self.TotalCount[NodeType]))
            #Score = math.log(True / float(self.TotalTrue[NodeType])) - math.log(Count / float(self.TotalCount[NodeType]))
        return struct.pack("<f", Score)
    def ReportStats(self, BinaryFileName):
        File = open(BinaryFileName, "wb")
        for NodeType in NodeTypes:
            print NodeType
            Ratio = self.TotalTrue[NodeType] / float(max(self.TotalCount[NodeType], 1))
            if Ratio > 0:
                LogRatio = math.log(Ratio)
            else:
                LogRatio = -5
            Str = struct.pack("<f", LogRatio)
            File.write(Str)
##            print "B TRUE:"
##            print self.BRankTrue[NodeType]
##            print "B COUNT:"
##            print self.BRankCount[NodeType]
            print "B ranks for %s:"%str(NodeType)
            for BRank in range(22):
                True =self.BRankTrue[NodeType].get(BRank, 0) 
                Count = self.BRankCount[NodeType].get(BRank, 0)
                print BRank, True, Count, 100*True/float(max(Count, 1))
                Str = self.GetProbStruct(self.BRankTrue, self.BRankCount, NodeType, BRank)
                print "B rank %s:"%BRank, Str
                File.write(Str)
            print "Y ranks for %s:"%str(NodeType)
            for YRank in range(22):
                True = self.YRankTrue[NodeType].get(YRank, 0) 
                Count = self.YRankCount[NodeType].get(YRank, 0)
                print YRank, True, Count, 100*True/float(max(Count, 1))
                Str = self.GetProbStruct(self.YRankTrue, self.YRankCount, NodeType, YRank)
                #print "Y rank %s:"%YRank, Str
                File.write(Str)
            print "Isotope flags for %s:"%str(NodeType)
            for IsoFlag in range(16):
                Count = self.IsotopeCount[NodeType].get(IsoFlag, 0)
                True = self.IsotopeTrue[NodeType].get(IsoFlag, 0)
                print IsoFlag, True, Count, 100 * True/float(max(Count, 1))
                Str = self.GetProbStruct(self.IsotopeTrue, self.IsotopeCount, NodeType, IsoFlag)
                File.write(Str)
            CuteResults = []
            for IonFlag in range(512):
                FlagName = self.GetFlagName(IonFlag)
                TrueCount = self.WitnessTrue[NodeType].get(IonFlag, 0)
                TotalCount = self.WitnessCount[NodeType].get(IonFlag, 0)
                #print IonFlag, self.GetFlagName(IonFlag), TrueCount, TotalCount
                Plentiful = 0
                if TotalCount < self.MinimumWitnessSampleSize:
                    (TrueCount, TotalCount) = self.GetFewerIonCount(IonFlag, self.WitnessCount[NodeType], self.WitnessTrue[NodeType])
                else:
                    Plentiful = 1
                if TrueCount == 0:
                    if TotalCount == 0:
                        Score = -1.5
                    else:
                        Score = -6.9 # log(.01)
                else:
                    Score = math.log(TrueCount / float(self.TotalTrue[NodeType])) - math.log(TotalCount / float(self.TotalCount[NodeType]))
                if Plentiful:
                    CuteResults.append((Score, IonFlag, FlagName, TrueCount, TotalCount, "%.2f"%(100*TrueCount/float(TotalCount))))                    
                #print "->",IonFlag, self.GetFlagName(IonFlag), TrueCount, TotalCount, Score
                File.write(struct.pack("<f", Score))
                
            CuteResults.sort()
            CuteResults.reverse()
            print "Charge %s sector %s %.2f%%"%(NodeType[0], NodeType[-1], 100*float(self.TotalTrue[NodeType]) / float(max(1,self.TotalCount[NodeType])))
            for Cute in CuteResults[:10]:
                print Cute
            print "B2 odds in this region: %s/%s (%.2f%%)"%(self.B2TrueCount[NodeType], self.B2TotalTrueCount[NodeType],
                                                            self.B2TrueCount[NodeType] / float(max(1, self.B2TotalTrueCount[NodeType])))
            print "Y2 odds in this region: %s/%s (%.2f%%)"%(self.Y2TrueCount[NodeType], self.Y2TotalTrueCount[NodeType],
                                                            self.Y2TrueCount[NodeType] / float(max(1, self.Y2TotalTrueCount[NodeType])))
            
        # And we're done:
        File.close()
    def AccumulateNeoStats(self, Spectrum, PRM, Peptide, BPeak = None, YPeak = None):
        # Determine the NodeType of this PRM:
        if (PRM > Spectrum.ParentMass * 0.66):
            Tri = 2
        elif (PRM > Spectrum.ParentMass * 0.33):
            Tri = 1
        else:
            Tri = 0
        NodeType = (min(3, Spectrum.Charge), Tri)
        #NodeType = (Spectrum.Charge,)
        # Classify the PRM as true or false:
        TruePRMFlag = 0
        for Mass in Peptide.Masses:
            if abs(PRM - Mass) < 0.5:
                TruePRMFlag = 1
        # Overall totals for this node type.  p(true)
        self.TotalCount[NodeType] = self.TotalCount.get(NodeType, 0) + 1
        self.TotalTrue[NodeType] = self.TotalTrue.get(NodeType, 0) + TruePRMFlag
        # Find all the witnesses:
        IonTypes = []
        # b:
        if not BPeak:
            BPeak = self.AccumulateNeoStatsHelper(None, Spectrum, IonTypes, PRM + 1.0078, "b", None)
        else:
            IonTypes.append("b")
        if (BPeak):
            BIntensity = BPeak.Intensity
        else:
            BIntensity = 0
        # y:
        if not YPeak:
            YPeak = self.AccumulateNeoStatsHelper(None, Spectrum, IonTypes, Spectrum.ParentMass - PRM, "y", None)
        else:
            IonTypes.append("y")
        if (YPeak):
            YIntensity = YPeak.Intensity
        else:
            YIntensity = 0
        self.AccumulateNeoStatsHelper(BIntensity, Spectrum, IonTypes, PRM - 16.026549128, "b-nh3", "b")        
        self.AccumulateNeoStatsHelper(BIntensity, Spectrum, IonTypes, PRM - 26.994914640, "a", "b")
        self.AccumulateNeoStatsHelper(BIntensity, Spectrum, IonTypes, (PRM + 2*1.0078)/2.0, "b2", "b")
        self.AccumulateNeoStatsHelper(YIntensity, Spectrum, IonTypes, Spectrum.ParentMass - PRM - 17.010564720, "y-nh3", "y")
        self.AccumulateNeoStatsHelper(YIntensity, Spectrum, IonTypes, Spectrum.ParentMass - PRM - 18.010564720, "y-h2o", "y")
        self.AccumulateNeoStatsHelper(YIntensity, Spectrum, IonTypes, (Spectrum.ParentMass - PRM + 1.0078)/2.0, "y2", "y")
        if (TruePRMFlag):
            # Track how many true PRMs have b2 and y2 ions:
            self.B2TotalTrueCount[NodeType] += 1
            if "b2" in IonTypes:
                self.B2TrueCount[NodeType] += 1
            self.Y2TotalTrueCount[NodeType] += 1
            if "y2" in IonTypes:
                self.Y2TrueCount[NodeType] +=1 
        IonTypeFlag = self.GetIonFlag(IonTypes)
        IonTypes.sort()
        ##print IonTypes
        self.WitnessCount[NodeType][IonTypeFlag] = self.WitnessCount[NodeType].get(IonTypeFlag, 0) + 1
        self.WitnessTrue[NodeType][IonTypeFlag] = self.WitnessTrue[NodeType].get(IonTypeFlag, 0) + TruePRMFlag
        ######################################################
        # B RANK
        RankLookup = [0,1,2,3,4,5,6,7,8,9,
                      10,10,10,10,10,
                      11,11,11,11,11,
                      12,12,12,12,12,
                      13,13,13,13,13,
                      14,14,14,14,14,
                      15,15,15,15,15,
                      16,16,16,16,16,16,16,16,16,16,
                      17,17,17,17,17,17,17,17,17,17,
                      18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,
                      19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,]
        if not BPeak:
            BRank = 21
        elif BPeak.IntensityRank >= len(RankLookup):
            BRank = 20
        else:
            BRank = RankLookup[BPeak.IntensityRank]
        self.BRankCount[NodeType][BRank] = self.BRankCount[NodeType].get(BRank, 0) + 1
        self.BRankTrue[NodeType][BRank] = self.BRankTrue[NodeType].get(BRank, 0) + TruePRMFlag
        ######################################################
        # Y RANK
        if not YPeak:
            YRank = 21
        elif YPeak.IntensityRank >= len(RankLookup):
            YRank = 20
        else:
            YRank = RankLookup[YPeak.IntensityRank]
        self.YRankCount[NodeType][YRank] = self.YRankCount[NodeType].get(YRank, 0) + 1
        self.YRankTrue[NodeType][YRank] = self.YRankTrue[NodeType].get(YRank, 0) + TruePRMFlag
        ######################################################
        # Isotopes:
        if not BPeak:
            BHas = 0
            BIs = 0
        else:
            BHas = BPeak.HasPlausibleIsotopicPeak
            BIs = BPeak.IsPlausibleIsotopicPeak
        if not YPeak:
            YHas = 0
            YIs = 0
        else:
            YHas = YPeak.HasPlausibleIsotopicPeak
            YIs = YPeak.IsPlausibleIsotopicPeak
        IsotopeFlag = (YIs * 1) + (YHas * 2) + (BIs * 4) + (BHas * 8)
        self.IsotopeCount[NodeType][IsotopeFlag] = self.IsotopeCount[NodeType].get(IsotopeFlag, 0) + 1
        self.IsotopeTrue[NodeType][IsotopeFlag] = self.IsotopeTrue[NodeType].get(IsotopeFlag, 0) + TruePRMFlag
    def AccumulateWitnessSetStats(self, Spectrum, Peptide):
        "Called once per spectrum.  Get stats on what sets of witness peaks are most common."
        for Peak in Spectrum.Peaks:
            # NOTE: For now, we don't consider double-losses, because there isn't enough training data
            # to draw strong conclusions from.
            ################################################################################
            # Consider this peak as a b peak.  What other ions support this interpretation?
            PRM = Peak.Mass - 1.008
            self.AccumulateNeoStats(Spectrum, PRM, Peptide, Peak, None)
            PRM = Spectrum.ParentMass - Peak.Mass
            self.AccumulateNeoStats(Spectrum, PRM, Peptide, None, Peak)
    def GetFlagName(self, IonFlag):
        Str = ""
        if IonFlag & 1:
            Str += "b, "
        if IonFlag & 16:
            Str += "y, "
        if IonFlag & 2:
            Str += "b-h2o, "
        if IonFlag & 4:
            Str += "b-nh3, "
        if IonFlag & 8:
            Str += "a, "
        if IonFlag & 32:
            Str += "y-h2o, "
        if IonFlag & 64:
            Str += "y-nh3, "
        if IonFlag & 128:
            Str += "b2, "
        if IonFlag & 256:
            Str += "y2, "
        return Str
    def GetIonFlag(self, IonTypes):
        Flag = 0
        if "b" in IonTypes:
            Flag += 1
        if "b-h2o" in IonTypes:
            Flag += 2
        if "b-nh3" in IonTypes:
            Flag += 4
        if "a" in IonTypes:
            Flag += 8
        if "y" in IonTypes:
            Flag += 16
        if "y-h20" in IonTypes:
            Flag += 32
        if "y-nh3" in IonTypes:
            Flag += 64
        if "b2" in IonTypes:
            Flag += 128
        if "y2" in IonTypes:
            Flag += 256
        return Flag
        
    def GetFewerIonCount(self, KeyFlag, Dict, TrueDict):
        # If our ion set is empty or b or y, then we can't drop anything:
        if (KeyFlag == 0x001 or KeyFlag == 0 or KeyFlag == 0x0010):
            return (TrueDict.get(KeyFlag, 0), Dict.get(KeyFlag, 1))
        Count = Dict.get(KeyFlag, 0)
        # Stop dropping flags, if we dropped enough already:
        if Count > self.MinimumWitnessSampleSize:
            return (TrueDict.get(KeyFlag, 0), Count)
        BestSubOdds = 0
        # Assume that we can take the best sub-flag and extend it.  (More ion types
        # is never WORSE than fewer ion types)
        BestTrueCount = 0
        BestTotalCount = 0
        for (SecondaryIonType, Flag) in [("b-h2o", 0x2),
                                ("b-nh3", 0x4),
                                ("a", 0x8),
                                ("y-h2o", 0x20),
                                ("y-nh3", 0x40),
                                ("b2", 0x80),
                                ("y2", 0x100)]:
            if not (KeyFlag & Flag):
                continue
            SubKeyFlag = KeyFlag - Flag
            (TrueCount, TotalCount) = self.GetFewerIonCount(SubKeyFlag, Dict, TrueDict)
            Odds = TrueCount / float(max(TotalCount, 1))
            if Odds >= BestSubOdds:
                BestSubOdds = Odds
                BestTrueCount = TrueCount
                BestTotalCount = TotalCount
        return (TrueDict.get(KeyFlag, 0) + BestTrueCount, Dict.get(KeyFlag, 0) + BestTotalCount)
    def WriteBinaryWitnessScores(self, Dict, TrueDict, File):
        ScoresByFlag = {}
        for Flag in range(512):
            Key = [] 
##        for Key in Dict.keys():
##            Flag = 0
##            if ("b" in Key):
##                Flag += 0x0001
##            if ("b-h2o" in Key):
##                Flag += 0x0002
##            if ("b-nh3" in Key):
##                Flag += 0x0004
##            if ("a" in Key):
##                Flag += 0x0008
##            if ("y" in Key):
##                Flag += 0x0010
##            if ("y-h2o" in Key):
##                Flag += 0x0020
##            if ("y-nh3" in Key):
##                Flag += 0x0040
##            if ("b2" in Key):
##                Flag += 0x0080
##            if ("y2" in Key):
##                Flag += 0x0100
            if Flag & 0x1:
                Key.append("b")
            if Flag & 0x2:
                Key.append("b-h2o")
            if Flag & 0x4:
                Key.append("b-nh3")
            if Flag & 0x8:
                Key.append("a")
            if Flag & 0x10:
                Key.append("y")
            if Flag & 0x20:
                Key.append("y-h2o")
            if Flag & 0x40:
                Key.append("y-nh3")
            if Flag & 0x80:
                Key.append("b2")
            if Flag & 0x100:
                Key.append("y2")
            Key.sort()
            Key = tuple(Key)
            TrueCount = TrueDict.get(Key, 0)
            TotalCount = Dict.get(Key, 0)
            Odds = TrueCount / float(max(TotalCount, 1))
            print "%s: %s / %s = %.2f%%"%(Key, TrueCount, TotalCount, Odds*100)
            if TotalCount < self.MinimumWitnessSampleSize:
                (SubKey, TrueCount, TotalCount) = self.GetFewerIonCount(Key, Flag, Dict, TrueDict)
                Odds = TrueCount / float(max(TotalCount, 1))
                print "  -> %s: %s / %s = %.2f%%"%(SubKey, TrueCount, TotalCount, Odds*100)
                
            if Odds == 0:
                LogOdds = -3.0
            else:
                LogOdds = math.log(Odds)
            ScoresByFlag[Flag] = LogOdds # convert probability to log score
            
        Flags = ScoresByFlag.keys()
        for Flag in range(512):
            File.write(struct.pack("<f", ScoresByFlag.get(Flag, -3.0)))
            #print Flag, ScoresByFlag.get(Flag, -3.0)

    def AccumulateNeoStatsHelper(self, MaxAllowedIntensity, Spectrum, IonTypes, Mass, IonType, RequireIonType):
        if RequireIonType != None and RequireIonType not in IonTypes:
            return # missing a prereq!
        LeftIndex = 0
        RightIndex = len(Spectrum.Peaks)
        MinMass = Mass - 0.5
        MaxMass = Mass + 0.5
        BestIntensity = 0
        BestPeak = None
        while (LeftIndex + 1 < RightIndex):
            MiddleIndex = (LeftIndex + RightIndex) / 2
            if (Spectrum.Peaks[MiddleIndex].Mass > MinMass):
                RightIndex = MiddleIndex
            else:
                LeftIndex = MiddleIndex
        for Index in range(LeftIndex, len(Spectrum.Peaks)):
            Peak = Spectrum.Peaks[Index]
            if (Peak.Mass < MinMass):
                continue
            if (Peak.Mass > MaxMass):
                break
            if (MaxAllowedIntensity==None or Peak.Intensity < MaxAllowedIntensity):
                BestIntensity = Peak.Intensity
                BestPeak = Peak
        if BestIntensity != 0:
            IonTypes.append(IonType)
        return BestPeak
    def AccumulateStatsFromFile(self, FilePath, Peptide):
        # Read in the spectrum.  Filter peaks by window, and
        # assign the intensity-based RankIndex of each peak:
        Spectrum = MSSpectrum.SpectrumClass()
        Spectrum.ReadPeaks(FilePath)
        Spectrum.ParentMass = 19 + Peptide.Masses[-1]
        Spectrum.RankWithFilter(50, 6)
        Spectrum.FindIsotopicPeaks()
        Spectrum.RankPeaksByIntensity()
        Spectrum.CorrectPeptide = Peptide
        self.AccumulateWitnessSetStats(Spectrum, Peptide)
        self.TotalSpectra += 1

def GetHistogramOddsTable(Histogram):
    Keys = Histogram.keys()
    Keys.sort()
    Extreme = Keys[-1]
    Extreme = max(Extreme, abs(Keys[0]))
    ExtremeCounts = {}
    OverallCount = 0
    for Value in range(Extreme, 0, -1):
        ExtremeCounts[Value] = ExtremeCounts.get(Value+1, 0)
        ExtremeCounts[Value] += Histogram.get(Value, 0)
        ExtremeCounts[Value] += Histogram.get(-Value, 0)
        OverallCount += Histogram.get(Value, 0)
        OverallCount += Histogram.get(-Value, 0)
    Histo = {}
    Histo[0] = math.log(1.0)
    for Value in range(1, Extreme+1):
        Odds = ExtremeCounts[Value] / float(OverallCount)
        print "%s\t%s"%(Value, Odds)
        if Odds == 0:
            Histo[Value] = -999999
        else:
            Histo[Value] = math.log(Odds)
    return Histo

def LoadHistogram(FileName, IsFloat = 0):
    File = open(FileName, "r")
    Histogram = {}
    for FileLine in File.xreadlines():
        Bits = FileLine.strip().split()
        if len(Bits)!=2:
            continue
        if IsFloat:
            Histogram[int(Bits[0])] = float(Bits[1])
        else:
            Histogram[int(Bits[0])] = int(Bits[1])
    File.close()
    return Histogram

def PrepareSkewTables():
    HistoStep = LoadHistogram("SingleStepSkew.txt")
    HistoStep = GetHistogramOddsTable(HistoStep)
    HistoFull = LoadHistogram("FullTagSkew.txt")
    HistoFull = GetHistogramOddsTable(HistoFull)
    PickleFile = open("SkewOdds.dat", "wb")
    cPickle.dump(HistoStep, PickleFile)
    cPickle.dump(HistoFull, PickleFile)
    PickleFile.close()

def PrepareIonTypeTables():
    BHisto = LoadHistogram("BIonOdds.txt", 1)
    for (Key, Value) in BHisto.items():
        if Value == 0:
            BHisto[Key] = -999
        else:
            BHisto[Key] = math.log(Value)
    YHisto = LoadHistogram("YIonOdds.txt", 1)
    for (Key, Value) in YHisto.items():
        if Value == 0:
            YHisto[Key] = -999
        else:
            YHisto[Key] = math.log(Value)
    PickleFile = open("IntensityRankIonOdds.dat", "wb")
    cPickle.dump(BHisto, PickleFile)
    cPickle.dump(YHisto, PickleFile)
    PickleFile.close()

def MungeWitnessSets():
    File = open("StatsMungeMe.txt", "r")
    for FileLine in File.xreadlines():
        Bits = FileLine.split()
        if len(Bits)<4:
            continue
        del Bits[-3]
        Bits[-2]=":"
        try:
            Bits[-1] = str(math.log(float(Bits[-1])))
        except:
            pass
        print "%s,"%string.join(Bits, " ") 
        

def Main():
    Initialize()
    ######################
    #FindIntenseSpectra()
    #####################
    GetSergeiPeptides()
    Oracle = Global.SergeiStuff 
    Stats = StatMeister(Oracle,)
    Stats.AccumulateStats("IntenseISB.txt", -1, SpectrumDir)
    #Stats.AccumulateStats(SpectrumListFileName)

def MergeBinaryWitnessFiles(FileNameList, OutputFileName = "IonWitnessScores.dat"):
    OutputFile = open(OutputFileName, "wb")
    OutputFile.write(struct.pack("<i", 512))
    OutputFile.write(struct.pack("<i", len(FileNameList)))
    ExpectLen = 8
    for FileName in FileNameList:
        print "Pos:", OutputFile.tell()
        File = open(FileName, "rb")
        Block = File.read()
        File.close()
        OutputFile.write(Block)
        print "Pos:", OutputFile.tell()
    OutputFile.close()

if __name__ == "__main__":
    Main()
    #MergeBinaryWitnessFiles(["ISBCharge1.dat", "ISBCharge2.dat","ISBCharge3.dat"])

