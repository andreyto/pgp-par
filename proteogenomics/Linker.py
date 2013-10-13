"""
Analyze cross-linked peptides from the Briggs lab.
"""
import os
import sys
import MSSpectrum
import cStringIO
from Utils import *
Initialize()

class Hyrule:
    def __init__(self):
        pass
    def LoadExpectedPeptides(self):
        """
        Parse regular search results, to get a list of peptides that we know
        to be present in the sample.
        """
        self.KnownPeptides = {}
        self.SpectrumPValues = {}
        LabeledPeptides = {}
        File = open("DTSSP.fixed.out", "rb")
        for FileLine in File.xreadlines():
            Bits = FileLine.split("\t")
            Annotation = Bits[2][2:-2]
            Peptide = GetPeptideFromModdedName(Annotation)
            if self.KnownPeptides.has_key(Annotation):
                (OldPeptide, OldCount) = self.KnownPeptides[Annotation]
                self.KnownPeptides[Annotation] = (OldPeptide, OldCount + 1)
            else:
                self.KnownPeptides[Annotation] = (Peptide, 1)
            if Annotation.find("+")!=-1:
                LabeledPeptides[Annotation] = LabeledPeptides.get(Annotation, 0) + 1
        print "Found a total of %s expected peptides."%len(self.KnownPeptides.keys())
        #####################
        print "Labeled peptides"
        SortedList = []
        for (Key, Count) in LabeledPeptides.items():
            SortedList.append((Count, Key))
        SortedList.sort()
        SortedList.reverse()
        for (Count, Key) in SortedList:
            print "%s: %s"%(Key, Count)
        return 
        Dir = r"E:\ms\Briggs\ResultsLinkerFixed"
        for FileName in os.listdir(Dir):
            Path = os.path.join(Dir, FileName)
            print "Parse results from %s..."%FileName
            File = open(Path, "rb")
            for FileLine in File.xreadlines():
                Bits = FileLine.split("\t")
                try:
                    PValue = float(Bits[10])
                    Peptide = GetPeptideFromModdedName(Bits[2][2:-2])
                    Spectrum = (FileName, int(Bits[1]))
                except:
                    continue
                if PValue <= 0.05:
                    self.KnownPeptides[Peptide.Aminos] = self.KnownPeptides.get(Peptide.Aminos, 0) + 1
                self.SpectrumPValues[Spectrum] = PValue
            File.close()
    def GetExpectedMasses(self, LinkerMassDelta):
        print "Compute expected masses of cross-linked peptides."
        self.ExpectedMasses = {} # mass (in daltons) -> list of peptide pairs
        Keys = self.KnownPeptides.keys()
        for PeptideIndexA in range(len(Keys)):
            KeyA = Keys[PeptideIndexA]
            (PeptideA, CountA) = self.KnownPeptides[KeyA]
            #MassA = PeptideA.Masses[-1] + 19
            MassA = GetMass(PeptideA.Aminos) + 19
            for PeptideIndexB in range(PeptideIndexA, len(Keys)):
                KeyB = Keys[PeptideIndexB]
                (PeptideB, CountB) = self.KnownPeptides[KeyB]
                #MassB = PeptideB.Masses[-1] + 19
                MassB = GetMass(PeptideB.Aminos) + 19
                Abundance = min(CountA, CountB)
                LinkedMass = MassA + MassB + LinkerMassDelta
                Bin = int(round(LinkedMass))
                if not self.ExpectedMasses.has_key(Bin):
                    self.ExpectedMasses[Bin] = []
                self.ExpectedMasses[Bin].append((Abundance, KeyA, KeyB))
        for List in self.ExpectedMasses.values():
            List.sort()
            List.reverse()
    def FindUnexplainedScans(self):
        Dir = r"E:\ms\Briggs\CrossLink\Uncleaved"
        for FileName in os.listdir(Dir):
            (Stub, Extension) = os.path.splitext(FileName)
            if Extension.lower() != ".mzxml":
                continue
            print "Find unexplained scans in %s..."%FileName
            Path = os.path.join(Dir, FileName)
            self.FindUnexplainedScansFromFile(Path, FileName)
    def FindUnexplainedScansFromFile(self, Path, FileName):
        File = open(Path, "rb")
        MZXML = File.read()
        MZXMLLen = len(MZXML)
        File.close()
        Pos = -1
        print "ScanNumber\tPeakCount\tSignalToNoise\tMZ\tOldPValue\tch2PM\tPeptides2\tch3PM\tPeptides3\t"
        while (1):
            ScanPos = MZXML.find("<scan", Pos + 1)
            if ScanPos == -1:
                break
            NextScanPos = MZXML.find("<scan", ScanPos + 1)
            if NextScanPos == -1:
                NextScanPos = MZXMLLen
            Pos = ScanPos
            ScanXML = MZXML[ScanPos:NextScanPos]
            if ScanXML.find("msLevel=\"2\"") == -1:
                continue # skip ms1 scans!
            QuotePos1 = ScanXML.find('"')
            QuotePos2 = ScanXML.find('"', QuotePos1 + 1)
            ScanNumber = int(ScanXML[QuotePos1+1:QuotePos2])
            #print "Scan %s from %s-%s"%(ScanNumber, ScanPos, NextScanPos)
            Spectrum = MSSpectrum.SpectrumClass()
            FakeFile = cStringIO.StringIO(ScanXML)
            Spectrum.ReadPeaksMZXML(FakeFile)
            PeakCount = len(Spectrum.Peaks)
            if PeakCount < 11:
                continue
            IntensityList = []
            for Peak in Spectrum.Peaks:
                IntensityList.append(Peak.Intensity)
            IntensityList.sort()
            BigIntensity = IntensityList[-5]
            SmallIntensity = 0
            for Val in IntensityList[:-10]:
                SmallIntensity += Val
            SmallIntensity /= float(len(IntensityList) - 10)
            SignalToNoise = BigIntensity / float(SmallIntensity)
            Str = "%s\t%s\t%s\t%s\t"%(ScanNumber, PeakCount, SignalToNoise, Spectrum.PrecursorMZ)
            OldPValueKey = (FileName, ScanNumber)
            Str += "%s\t"%self.SpectrumPValues.get(OldPValueKey, "")
            #######################
            Spectrum.SetCharge(2)
            Str += "%s\t"%Spectrum.ParentMass
            Bin = int(round(Spectrum.ParentMass))
            Charge2Peptides = ""
            for NearBin in range(Bin - 1, Bin + 2):
                PepList = self.ExpectedMasses.get(NearBin, [])
                for PepPair in PepList[:5]:
                    Charge2Peptides += "%s,"%str(PepPair)
            Str += "%s\t"%Charge2Peptides
            #######################
            Spectrum.SetCharge(3)
            Str += "%s\t"%Spectrum.ParentMass
            Bin = int(round(Spectrum.ParentMass))
            Charge3Peptides = ""
            for NearBin in range(Bin - 1, Bin + 2):
                PepList = self.ExpectedMasses.get(NearBin, [])
                for PepPair in PepList[:5]:
                    Charge3Peptides += "%s,"%str(PepPair)
            Str += "%s\t"%Charge3Peptides
            print Str

if __name__ == "__main__":
    LinkFinder = Hyrule()
    LinkFinder.LoadExpectedPeptides()
    # Uncleavable MassDelta: -2H + 8C + 2 Ox + 12H + 2C
    LinkFinder.GetExpectedMasses(138)
    #LinkFinder.GetExpectedMasses(160)
    LinkFinder.FindUnexplainedScans()