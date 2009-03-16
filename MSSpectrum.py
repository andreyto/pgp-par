"""
Classes representing an MS/MS spectrum and its peaks.
"""
import traceback
import sys
import os
import types
import string
import re
import struct
import math
import base64 # for mzxml parsing
import ParseXML
from Utils import *
Initialize()

# Some systems have old versions of base64, without the new decode interface.
# As a workaround, set "B64Decode" to the decoding function.
if hasattr(base64, "b64decode"):
    B64Decode = base64.b64decode
    B64Encode = base64.b64encode
else:
    B64Decode = base64.decodestring
    B64Encode = base64.encodestring

class PeakClass:
    """
    One peak from an ms/ms spectrum.  Mostly just a mass, but we track the intensity and
    (perhaps putative) ion-type, as well.  Note that this may be a spectral peak (in
    which case the mass is a spectral mass) or a PRM peak (in which case the mass is
    a prefix residue mass).  Each spectral peak gives rise to n PRM peaks, where n
    is the number of ion types available.
    """
    def __init__(self, Mass, Intensity):
        self.Mass = Mass
        self.Intensity = Intensity
        self.IonType = None # Assigned only for PRM peaks
        self.Score = 0 # Score, based on current filtering scheme.
        self.PeptideIndex = None
        self.FilterRank = None
        self.IntensityRank = None 
        self.IsPlausibleIsotopicPeak = 0
        self.HasPlausibleIsotopicPeak = 0
    def __cmp__(self, Other):
        "Sort two peak objects - compare the masses"
        if not isinstance(Other, PeakClass):
            return 1
        if (self.Mass < Other.Mass):
            return -1
        if (self.Mass > Other.Mass):
            return 1
        return 0
    def __str__(self):
        return "<peak %s>"%self.Mass
    def GetPeakMass(self, ParentMass):
        "Get the corresponding spectral mass from this PRM peak"
        return self.IonType.GetPeakMass(self.Mass, ParentMass)
    def GetPRMMass(self, ParentMass):
        "Get the corresponding PRM mass from this spectral peak"
        return self.IonType.GetPRMMass(self.Mass, ParentMass)
    def PrintMe(self):
        print "Printing information for a PeakClass object"
        print "Mass %f Intensity %f "%(self.Mass,self.Intensity)
        print "IonType %s PeptideIndex %s IntensityRank %d"%(self.IonType, self.PeptideIndex, self.IntensityRank)

class SpectrumClass:
    """
    Mass-spec data, and some functions to filter peaks and such.  
    """
    InstanceCount = 0
    def __del__(self):
        SpectrumClass.InstanceCount -= 1
    def __init__(self, Scoring = None):
        SpectrumClass.InstanceCount += 1
        # Init our attributes:
        self.Name = None
        self.ParentMass = None
        self.PrecursorMZ = None
        self.PrecursorIntensity = None
        self.Charge = 1 # default
        self.Peaks = None # list of PeakClass instances
        self.PRMPeaks = None # list of PeakClass instances
        # The actual parent peptide (instance of PeptideClass), if known:
        self.CorrectPeptide = None
        self.Scoring = Scoring
    def GetSignalToNoise(self):
        "Return signal-to-noise ratio for this spectrum"
        Intensities = []
        for Peak in self.Peaks:
            Intensities.append(Peak.Intensity)
        Intensities.sort()
        IntenseCount = min(len(Intensities), 5)
        if not IntenseCount:
            return 0
        Signal = Intensities[-IntenseCount/2]
        Noise = Intensities[len(Intensities)/2]
        return Signal / float(Noise)
    def GetTotalIntensity(self):
        Intensity = 0
        for Peak in self.Peaks:
            Intensity += Peak.Intensity
        return Intensity 
    def SetCharge(self, NewCharge):
        self.Charge = NewCharge
        self.ParentMass = self.PrecursorMZ * NewCharge - (NewCharge - 1)*1.0078
    def ReadPeaksMGF(self, File):
        self.Peaks = []
        for FileLine in File.xreadlines():
            if FileLine[:8] == "END IONS":
                break
            if FileLine[:6] == "CHARGE":
                # Note: "2+ and 3+" is NOT supported.  Use the MultiCharge option,
                # or include two scans in the input file.
                Charge = int(FileLine[7:9].replace("+",""))
                continue
            if FileLine[:7] == "PEPMASS":
                #self.ParentMass = float(FileLine[8:])
                self.PrecursorMZ = float(FileLine[8:].split()[0])
                continue
            Bits = FileLine.split()
            try:
                Mass = float(Bits[0])
                Intensity = float(Bits[1])
            except:
                continue # some other header line we don't eat.
            Peak = PeakClass(Mass, Intensity)
            self.Peaks.append(Peak)
        if self.Charge == 0:
            # Guess!
            self.Charge = 2
        self.ParentMass = (self.PrecursorMZ * self.Charge) - (self.Charge-1)*1.0078
        #print "PrecursorMZ %s charge %s"%(self.PrecursorMZ, self.Charge)
    def ReadPeaksDTA(self, File):
        "Read a spectrum from a file, assuming .dta or .pkl or .ms2 format."
        HeaderLine = File.readline()
        if not HeaderLine.strip():
            HeaderLine = File.readline()
            if not HeaderLine.strip():
                HeaderLine = File.readline()
                if not HeaderLine.strip():
                    HeaderLine = File.readline()
                    if not HeaderLine.strip():
                        HeaderLine = File.readline()
        #print "HeaderLine: '%s'"%HeaderLine.strip()
        Bits = HeaderLine.strip().split()
        if HeaderLine[:7]=="CHARGE=":
            self.Charge = int(HeaderLine[7])
            HeaderLine = File.readline()
            self.PrecursorMZ = float(HeaderLine[8:])
            self.ParentMass = (self.PrecursorMZ * self.Charge) - (self.Charge-1)*1.0078
        elif Bits[0] == "Z": # MS2 format:
            self.Charge = int(Bits[1])
            self.ParentMass = float(Bits[2])
        elif HeaderLine[0] == ":": # MS2 colon format:
            HeaderLine = File.readline()
            Bits = HeaderLine.strip().split()
            self.ParentMass = float(Bits[0])
            self.Charge = int(Bits[1]) # always an integer!
        elif len(Bits) == 3: #PKL format:
            self.PrecursorMZ = float(Bits[0])
            self.Charge = int(Bits[2])
            self.ParentMass = (self.PrecursorMZ * self.Charge) - (self.Charge-1)*1.0078
        else:
            self.ParentMass = float(Bits[0])
            self.Charge = int(Bits[1]) # always an integer!
        if self.Charge == 0:
            # Guess!
            self.Charge = 2
            self.PrecursorMZ = self.ParentMass
            self.ParentMass = (self.PrecursorMZ * self.Charge) - (self.Charge-1)*1.0078
            #print "Prescursor MZ is %.2f, so guess a parent mass of %.2f"%(self.PrecursorMZ, self.ParentMass)
        else:
            self.PrecursorMZ = (self.ParentMass + (self.Charge-1)*1.0078) / self.Charge
        self.Peaks = []
        for FileLine in File.xreadlines():
            Bits = FileLine.split()
            # Skip comments:
            if FileLine[0] == "#":
                continue
            Bits = FileLine.split()
            if not Bits:
                break
            if Bits[0] == "Z":
                continue # special for ms2: ignore.
            if len(Bits) > 2:
                break # no more!
            try:
                Mass = float(Bits[0])
                Intensity = float(Bits[1])
            except:
                # It's over, over, over.
                break
            Peak = PeakClass(Mass, Intensity)
            self.Peaks.append(Peak)
        #File.close()
        self.Peaks.sort() # sort by mass
    def ReadPeakDTALine(self, FileLine):
        FileLine = FileLine.strip()
        Bits = FileLine.split()
        if len(Bits) < 2:
            return # blank (or broken) line, skip
        Peak = PeakClass(float(Bits[0]), float(Bits[1]))
        # If this is a labeled .dta file, read the ion types and the peptide indices:
        if len(Bits)>2:
            Peak.IonType = Global.AllIonDict.get(Bits[2], None)
        if len(Bits)>3:
            try:
                Peak.PeptideIndex = int(Bits[3])
            except:
                pass # silent failure (a novel)
        self.Peaks.append(Peak)
    def RankPeaksByIntensity(self):
        "Set Peak.IntensityRank for each of our peaks."
        PeaksSortedByIntensity = []
        for Peak in self.Peaks:
            PeaksSortedByIntensity.append((Peak.Intensity, Peak))
        PeaksSortedByIntensity.sort()
        PeaksSortedByIntensity.reverse()
        for Index in range(len(PeaksSortedByIntensity)):
            PeaksSortedByIntensity[Index][1].IntensityRank = Index
    def ReadPeaksPKL(self, File):
        "Read peaks from a file in .pkl format"
        HeaderLine = File.readline()
        Bits = HeaderLine.split()
        if len(Bits)!=3:
            # .pkl files should have precursor m/z, precursor peak intensity, and
            # guessed charge.  If we don't have three pieces, then this isn't a
            # valid .pkl file...
            raise ValueError, "Invalid input file: Header line '%s' not a .pkl header."%HeaderLine
        self.PrecursorMZ = float(Bits[0])
        self.PrecursorIntensity = float(Bits[1])
        self.Charge = int(Bits[2])
        # We hope to be called with an actual charge.  If we didn't get one at all, then guess 2.
        if not self.Charge:
            self.Charge = 2
        self.ParentMass = (self.PrecursorMZ * self.Charge) - (1.0078 * (self.Charge - 1))
        ##print "Prec %.2f times charge %s gives pm %s"%(self.PrecursorMZ, self.Charge, self.ParentMass)
        # All subsequent lines: Mass and intensity
        self.Peaks = []
        for FileLine in File.xreadlines():
            Bits = FileLine.split()
            if len(Bits) > 2:
                break
            self.ReadPeakDTALine(FileLine)
        File.close()
        self.Peaks.sort() # sort by mass
    def ReadPeaksFromFile(self, File, FileName):
        # Given a file 'blah123.dta', name the spectrum 'blah123'.
        self.Name = os.path.split(FileName)[1]
        (self.Name, FileExtension) = os.path.splitext(self.Name)
        FileExtension = FileExtension.lower()
        # Strip ".mzxml:444279" to ".mzxml":
        if FileExtension.find(":")!=-1:
            FileExtension = FileExtension[:FileExtension.find(":")]
        # Use the appropriate parser
        if FileExtension == ".pkl":
            self.ReadPeaksPKL(File)
        elif FileExtension == ".mzxml":
            self.ReadPeaksMZXML(File)
        elif FileExtension == ".mzdata":
            self.ReadPeaksMZData(File)
        elif FileExtension == ".mgf":
            self.ReadPeaksMGF(File)
        else:
            # default case: DTA
            self.ReadPeaksDTA(File)
    def ReadPeaks(self, FileName, FilePos = None):
        """
        Instantiator - Read a spectrum from a file.  Sets ParentMass, Charge,
        and Peaks list.  Doesn't filter, yet.
        """
        if FilePos == None:
            try:
                ColonBits = FileName.split(":")
                FilePos = int(ColonBits[-1])
                FileName = string.join(ColonBits[:-1], ":")
            except:
                FilePos = 0
        try:
            File = open(FileName, "rb")
        except:
            print "Error in ReadPeaks(): File '%s' couldn't be opened."%FileName
            traceback.print_exc()
            return
        File.seek(FilePos)
        self.ReadPeaksFromFile(File, FileName)
        self.FilePath = FileName
        self.FilePos = FilePos
        File.close()
    def ReadPeaksMZData(self, File):
        """
        Parse peaks from an .mzdata format file.  This format is slightly inferior
        to .mzxml, and not as commonly used.
        """
        ParseXML.GetSpectrumPeaksMZData(self, File)
        self.Charge = 2 # guess!
        self.ParentMass = (self.PrecursorMZ * self.Charge) - (self.Charge - 1) * 1.0078
        return 
        
    def ReadPeaksMZXML(self, File):
        """
        Parse peaks from an .mzXML format file.  Assumes we've already scanned to the
        desired file offset.
        """
        ParseXML.GetSpectrumPeaksMZXML(self, File)
        self.Charge = 2 # guess!
        self.ParentMass = (self.PrecursorMZ * self.Charge) - (self.Charge - 1) * 1.0078
        return 
    def DebugPrint(self, ShowPeaks = 0):
        """
        Print information on our spectrum, for debugging.
        """
        print "Spectrum '%s' has parent mass %f,\n   charge %f, and %d peaks"%(self.Name, self.ParentMass,
            self.Charge, len(self.Peaks))
        if self.CorrectPeptide:
            print " True parent peptide is: %s"%self.CorrectPeptide
        if ShowPeaks:
            for Peak in self.Peaks:
                if Peak.IonType:
                    print "  %f\t%f\t%s"%(Peak.Mass, Peak.Intensity, Peak.IonType.Name)
                else:
                    print "  %f\t%f"%(Peak.Mass, Peak.Intensity)
    def GetBestPeak(self, Mass, MaxIntensity = None, Epsilon = 1.0):
        "Used in labeling.  Find the best nearby peak whose intensity doesn't exceed our limit."
        if MaxIntensity == 0:
            return (None, None)
        BestPeak = None
        BestPeakError = None
        ClosestError = None
        for Peak in self.Peaks:
            Error = Peak.Mass - Mass
            if Error < -Epsilon:
                continue
            if Error > Epsilon:
                break
            if MaxIntensity and Peak.Intensity > MaxIntensity:
                continue # forbid neutral losses which are taller than the original            
            if (BestPeak == None or BestPeak.Intensity < Peak.Intensity):
                BestPeak = Peak
                BestPeakError = Error
        return (BestPeak, BestPeakError)
        
    def GetPeak(self, Mass, Epsilon = 1.0):
        """
        Get the closest peak to the specified mass, with a maximum error of Epsilon.
        """
        ClosestPeak = None
        ClosestError = None
        for Peak in self.Peaks:
            Error = abs(Peak.Mass - Mass)
            if Error < Epsilon:
                if (ClosestPeak == None or ClosestError > Error):
                    ClosestPeak = Peak
                    ClosestError = Error
            if Peak.Mass > Mass:
                break
        return ClosestPeak
    def GetPRMPeak(self, Mass, Epsilon = 1.0):
        """
        Get the closest peak to the specified mass, with a maximum error of Epsilon.
        """
        ClosestPeak = None
        ClosestError = None
        for Peak in self.PRMPeaks:
            Error = abs(Peak.Mass - Mass)
            if Error < Epsilon:
                if (ClosestPeak == None or ClosestError > Error):
                    ClosestPeak = Peak
                    ClosestError = Error
            if Peak.Mass > Mass:
                break
        return ClosestPeak
    def GetPRMPeaks(self, Mass, Epsilon = 1.0):
        """
        Get all peaks within Epsilon of Mass
        """
        Peaks = []
        for Peak in self.PRMPeaks:
            Error = abs(Peak.Mass - Mass)
            if Error < Epsilon:
                Peaks.append(Peak)
            if Peak.Mass > Mass:
                break
        return Peaks
    def AssignIonTypesFromPeptide(self):
        """
        Assign ion types to our peaks, based on the CorrectPeptide.
        """
        # The true PRMPeaks are sums of peptide masses.  Iterate
        # over the length of the peptide:
        LeftMass = 0
        for Index in range(0, len(self.CorrectPeptide.Aminos)):
            LeftMass += Global.AminoMass[self.CorrectPeptide.Aminos[Index]]
            LeftMass += Global.FixedMods.get(self.CorrectPeptide.Aminos[Index], 0)
            # For this PRMPeak, look for all the possible spectral peaks
            # corresponding to the various ion types:
            for IonType in AllIons:
                Mass = IonType.GetPeakMass(LeftMass, self.ParentMass)
                Peak = self.GetPeak(Mass, 1.0)
                if Peak:
                    Peak.IonType = IonType
                    Peak.Pep = self.CorrectPeptide[:Index+1]
    def ApplyWindowFilter(self, RegionCutoffs, WindowSizes, MaxRankInclusive):
        """
        Apply this window-fiter to our peaks.  RegionCutoffs describe the edges
        of "early", "medium" and "late" spectral portions; WindowSizes are the
        sizes (in AMUs) of windows for these portions.  MaxRankInclusive is the
        worst rank to keep.  
        """
        #print "Apply window:", WindowSizes, RegionCutoffs, MaxRankInclusive
        GoodPeaks = []
        # List of region-edges:
        Borders = []
        for Cutoff in RegionCutoffs:
            Borders.append(self.ParentMass * Cutoff)
        NextBorderIndex = 0
        LastBorderIndex = len(RegionCutoffs)
        WindowIndex = 0
        BadPeakIntensityList = []
        for Peak in self.Peaks:
            while (NextBorderIndex < LastBorderIndex and Peak.Mass > Borders[NextBorderIndex]):
                NextBorderIndex += 1
                WindowIndex += 1
            WindowSize = WindowSizes[WindowIndex]
            MinMass = Peak.Mass - WindowSize/2
            MaxMass = Peak.Mass + WindowSize/2
            List = []
            for OtherPeak in self.Peaks:
                if OtherPeak.Mass > MaxMass:
                    break
                if OtherPeak.Mass > MinMass:
                    List.append((OtherPeak.Intensity, OtherPeak))
            List.sort()
            List.reverse() # best to worst
            if (len(List) < MaxRankInclusive+1) or (Peak.Intensity >= List[MaxRankInclusive][0]):
                GoodPeaks.append(Peak)
            else:
                BadPeakIntensityList.append(Peak.Intensity)
        #print "Kept %d of %d original peaks."%(len(GoodPeaks), len(self.Peaks))
        self.Peaks = GoodPeaks
        if len(BadPeakIntensityList):
            BadPeakIntensityList.sort() 
            return BadPeakIntensityList[len(BadPeakIntensityList)/2]
        else:
            return -1
    def FilterPeaks(self, WindowSize = 50, PeakCount = 6):
        self.ApplyWindowFilter([], (WindowSize,), PeakCount - 1)       
    def WritePeaks(self, FilePath):
        "Write out a .dta file."
        File = open(FilePath, "w")
        File.write("%f\t%d\n"%(self.ParentMass, self.Charge))
        for Peak in self.Peaks:
            File.write("%f\t%f\n"%(Peak.Mass, Peak.Intensity))
        File.close()
    def WritePKLPeaks(self,FilePath):
        """"
        Append to the end of a .pkl file.  Note this APPENDS a file.
        if no precursor intensity is known, then we say zero.  I hope that does not break things.
        """
        FileHandle = open(FilePath, "a")
        if self.PrecursorIntensity:
            FileHandle.write("%s %s %s\n"%(self.PrecursorMZ,self.PrecursorIntensity,self.Charge))
        else:
            FileHandle.write("%s 0.0 %s\n"%(self.PrecursorMZ,self.Charge))
        for Peak in self.Peaks:
            FileHandle.write("%f\t%f\n"%(Peak.Mass, Peak.Intensity))
        FileHandle.write("\n") #need a blank line to separate different scans
        FileHandle.close
    def WriteMGFPeaks(self, TheFile, Title = "Spectrum", ScanNumber = None):
        """
        Append to the end of an mgf file.  Pass in an open file, or
        (as a string) the path of a file to be APPENDED to.
        """
        if type(TheFile) == type(""):
            File = open(FilePath, "a")
        else:
            File = TheFile
        File.write("BEGIN IONS\n")
        File.write("TITLE=%s\n"%Title)
        if ScanNumber != None:
            File.write("SCAN=%s\n"%ScanNumber)
        File.write("CHARGE=%d\n"%self.Charge)
        File.write("PEPMASS=%f\n"%self.PrecursorMZ)
        for Peak in self.Peaks:
            File.write("%f\t%f\n"%(Peak.Mass, Peak.Intensity))
        File.write("END IONS\n")
        if type(TheFile) == type(""):
            File.close()
    def WriteMZXMLPeaks(self, File, ScanNumber):
        PeakCount = len(self.Peaks)
        Str = """<scan num="%s" msLevel="2" peaksCount="%s" polarity="+" scanType="Full" lowMz="125" highMz="2000" """%(ScanNumber, PeakCount)
        Str += """\n<precursorMz """
        if self.PrecursorIntensity:
            Str += """ precursorIntensity = "%.2f" """%self.PrecursorIntensity
        Str += ">%.5f</precursorMz>\n"%self.PrecursorMZ
        Str += """\n<peaks precision="32" byteOrder="network" pairOrder="m/z-int">"""
        PeakString = ""
        for Peak in self.Peaks:
            PeakString += struct.pack(">ff", Peak.Mass, Peak.Intensity)
        PeakString = B64Encode(PeakString)
        Str += PeakString
        Str += "</peaks>\n</scan>\n"
        File.write(Str + "\n")
        
    def GetTopUnexplainedPeak(self):
        TopUXRank = len(self.Peaks)
        for Peak in self.Peaks:
            if Peak.IonType == None and Peak.IntensityRank < TopUXRank:
                TopUXRank = Peak.IntensityRank
        return TopUXRank
    def FindIsotopicPeaks(self):
        for PeakIndex in range(len(self.Peaks)):
            Peak = self.Peaks[PeakIndex]
            RoundMass = int(round(Peak.Mass))
            ExpectedFraction = Global.IsotopeWeights.get(RoundMass, None)
            if ExpectedFraction==None:
                continue
            for IsotopePeakIndex in range(PeakIndex+1, len(self.Peaks)):
                OtherPeak = self.Peaks[IsotopePeakIndex]
                if OtherPeak.Mass < Peak.Mass + 0.8:
                    continue
                if OtherPeak.Mass > Peak.Mass + 1.2:
                    break
                Fraction = OtherPeak.Intensity / Peak.Intensity
                # magic numbers ahoy:
                if abs(Fraction - ExpectedFraction) < 0.5 or (abs(Fraction - ExpectedFraction) < 0.8 and OtherPeak.Mass > Peak.Mass + 0.9 and OtherPeak.Mass < Peak.Mass + 1.1): 
                    OtherPeak.IsPlausibleIsotopicPeak = 1
                    Peak.HasPlausibleIsotopicPeak = 1
    def GetExplainedIntensity(self):
        """
        Callable *after* the spectrum has been labeled.  Returns the percentage
        of total spectral intensity that has been explained by labels.  All
        things being equal, a candidate peptide with a higher explained
        intensity is BETTER.
        """
        TotalIntensity = 0
        ExplainedIntensity = 0
        for Peak in self.Peaks:
            TotalIntensity += Peak.Intensity
            if Peak.IonType != None:
                ExplainedIntensity += Peak.Intensity
            #print "%s\t%s\t%s\t%s\t%s\t"%(Peak.Mass, Peak.Intensity, Peak.IonType, Peak.PeptideIndex, Peak.RescueFlag)
        return ExplainedIntensity / float(max(1, TotalIntensity))
    def GetExplainedIons(self, Peptide, DynamicRangeMin = 150, DynamicRangeMax = 2000):
        "Return the percentage of b and y peaks present."
        Annotated = {}
        
        PhosphorylationFlag = 0
        PhosB = [0]*40
        PhosY = [0]*40
        for (Pos, ModList) in Peptide.Modifications.items():
            for Mod in ModList:
                if Mod.Name == "Phosphorylation":
                    PhosphorylationFlag = 1
        for Peak in self.Peaks:
            if Peak.IonType:
                Annotated[(Peak.IonType.Name, Peak.PeptideIndex)] = 1
        Count = 0
        Present = 0
        TotalCutPresent = 0
        PM = 19 + Peptide.Masses[-1]
        for Index in range(len(Peptide.Masses)):
            CutPresent = 0
            BMass = Peptide.Masses[Index] + 1.0078
            if BMass > DynamicRangeMin and BMass < DynamicRangeMax:
                Count += 1
                BPresent = Annotated.get(("b", Index),0)
                BPresent |= Annotated.get(("b-p", Index),0)
                Present += BPresent
                CutPresent |= BPresent
            YMass = PM - Peptide.Masses[Index]
            if YMass > DynamicRangeMin and YMass < DynamicRangeMax:
                Count += 1
                YPresent = Annotated.get(("y", len(Peptide.Aminos) - Index), 0)
                YPresent |= Annotated.get(("y-p", Index),0)
                Present += YPresent
                CutPresent |= YPresent
            # Count the CUT POINTS that are witnessed:
            if (Index and Index<len(Peptide.Masses)-1) and CutPresent:
                TotalCutPresent += 1
        CutCount = len(Peptide.Masses) - 1
        #print Peptide.Aminos, "%s cut points of %s"%(TotalCutPresent, CutCount)
        return (Present, Count, Present / max(1, float(Count)),
                TotalCutPresent, CutCount, TotalCutPresent/max(1, float(CutCount)))
            
        
    def GetExplainedPeaks(self, MaxRank = 24):
        """
        Returns the percentage of the top n peaks that have been explained
        by peak labeling, where n = MaxRank.  The output of GetExplainedPeaks()
        should be high for a good candidate peptide.
        """
        TotalGoodPeaks = 0
        ExplainedGoodPeaks = 0
        for Peak in self.Peaks:
            if Peak.IntensityRank <= MaxRank:
                TotalGoodPeaks += 1
                #print Peak.IntensityRank, Peak.Mass, Peak.IonType
                if Peak.IonType != None:
                    ExplainedGoodPeaks += 1
        return ExplainedGoodPeaks / float(max(1, TotalGoodPeaks))
    def GetLogMeanStdev(self):
        """computes the mean and standard deviation of the peak intensities
        This can be done at any time before or after filtering.  just
        make sure that you know what it means.
        This computes things based on the LOG intensity values
        """
        IntensitySum = 0.0
        NumPeaks = len(self.Peaks)
        for Peak in self.Peaks:
            IntensitySum += math.log(Peak.Intensity)
        Mean = IntensitySum / NumPeaks
        VarSum = 0.0
        for Peak in self.Peaks:
            Diff = math.log(Peak.Intensity) - Mean
            VarSum += Diff*Diff
        Variance = VarSum / NumPeaks
        Stdev = math.sqrt(Variance)
        return (Mean,Stdev)
    