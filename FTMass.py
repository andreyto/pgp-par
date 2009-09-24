"""
If we have several spectra available for a PTM, then perhaps we can determine its mass
to high accuracy, particularly if some spectra were acquired on FT instruments.
Plan:
- Parse the feature table output by TrainPTMFeatures.py
- Parse annotations-files, produced by SelectSites.py
- Read the parent masses from the spectra annotated with PTMs.  Get a list of all
the mass-deviations.  Get the mean and standard deviation.
"""
import sys
import os
import traceback
import BasicStats
import getopt
import MSSpectrum
from Utils import *
Initialize()
from TrainPTMFeatures import FormatBits

class PTMSite:
    "Simple class for examining a particular site."
    def __init__(self):
        self.Bits = None
        self.Mass = None
        self.Masses = [] # all spectra
        self.FTMasses = [] # 

class FTMassMaster:
    def __init__(self):
        self.SiteDict = {} # (DBPos, ModMass) -> PTMSite instance
        self.HeaderLines = []
        self.DBPosCache = {}
        # Some hard-coded paths:
        self.DBPath = "Database\\ShewHalf.trie"
        self.FeaturePath = None #"PTMScore\\FeaturesShew10.txt"
        #self.PTMAnnotationPath = "ptmscore\shewptm10\PTMAnnotations.txt"
        self.SpectrumDir = r"f:\ftproot\shewanella\MSData"
        self.QuickParseFlag = 0 # for debugging!
        self.OutputFileName = "PMCorrectedFeatures.txt"
    def LoadFTScanInfo(self):
        """
        Populate self.FTFlagDict, which will tell us whether each run was
        performed on an FT instrument.
        """
        self.FTFlagDict = {}
        Path = r"ShewExperimentInformation.txt"
        File = open(Path, "rb")
        for FileLine in File.xreadlines():
            Bits = FileLine.split("\t")
            Stub = Bits[1]
            if Bits[4] == "LTQ_FT":
                self.FTFlagDict[Stub] = 1
                print "FT run:", Stub
        File.close()
    def LoadDB(self):
        File = open(self.DBPath, "rb")
        self.DB = File.read()
        File.close()
    def LoadSites(self):
        """
        Populate the dictionary self.SiteDict, where
        self.SiteDict[(DBPos, ModMass)] is a PTMSite instance
        """
        File = open(self.FeaturePath, "rb")
        LineNumber = 0
        for FileLine in File.xreadlines():
            LineNumber += 1
            if LineNumber % 500 == 0:
                print "PTMSite line %s..."%LineNumber
            if FileLine[0] == "#":
                self.HeaderLines.append(FileLine)
                continue
            Bits = list(FileLine.split("\t"))
            Bits[-1] = Bits[-1].replace("\r","").replace("\n","")
            try:
                Peptide = GetPeptideFromModdedName(Bits[FormatBits.OriginalAnnotation])
            except:
                self.HeaderLines.append(FileLine)
                continue # header line, skip it
            DBPos = self.DB.find(Peptide.Aminos)
            if not Peptide.Modifications.keys():
                print "** Warning: Unmodified peptide in line %s of %s"%(LineNumber, self.FeaturePath)
                print Bits
                continue
            ModifiedResidueNumber = Peptide.Modifications.keys()[0]
            ModMass = int(round(Peptide.Modifications[ModifiedResidueNumber][0].Mass))
            DBPos += ModifiedResidueNumber
            SiteKey = (DBPos, ModMass)
            if self.SiteDict.has_key(SiteKey):
                print "** Error: Duplicate site!", SiteKey
            #print SiteKey, Bits[FormatBits.Peptide]
            Site = PTMSite()
            Site.DBPos = DBPos
            Site.Mass = ModMass
            Site.Bits = Bits
            self.SiteDict[SiteKey] = Site
        File.close()
    def ParsePTMResultsFromFile(self, FilePath):
        """
        We've already parsed a list of the PTMs.  Now we're ready to parse the annotations
        of all spectra, and along the way we'll get a distribution of masses for
        each site.
        """
        if os.path.isdir(FilePath):
            FileNames = os.listdir(FilePath)
            FileNames.sort()
            if self.QuickParseFlag:
                FileNames = FileNames[:10]
            for FileNameIndex in range(len(FileNames)):
                SubFileName = FileNames[FileNameIndex]
                print "(%s/%s) %s (so far: %s:%s)"%(FileNameIndex, len(FileNames), SubFileName, self.AnnotationCount, self.FTAnnotationCount)
                SubPath = os.path.join(FilePath, SubFileName)
                self.ParsePTMResultsFromFile(SubPath)
            return
        # Check the file extension.  For now, we consider ONLY .txt files.
        (Stub, Extension) = os.path.splitext(FilePath)
        if Extension != ".txt":
            return 
        File = open(FilePath, "rb")
        LineNumber = 0
        FTFlag = None
        for FileLine in File.xreadlines():
            LineNumber += 1
            if LineNumber % 1000 == 0:
                print "Line %s..."%LineNumber
            if self.QuickParseFlag and LineNumber > 5000:
                break
            Bits = FileLine.split("\t")
            try:
                Annotation = Bits[2]
                Peptide = GetPeptideFromModdedName(Bits[2])
                Charge = int(Bits[4])
            except:
                continue # header or other garbage line
            #######################################################
            # On the first line of the file, we determine whether this
            # is an FT run, and set the value of FTFlag.
            if FTFlag == None:
                FileName = Bits[0].replace("/", "\\").split("\\")[-1]
                Stub = FileName.split(".")[0]
                FTFlag = self.FTFlagDict.get(Stub, 0)
                if Stub[-4:] == "_dta":
                    FTFlag = FTFlag or self.FTFlagDict.get(Stub[:-4], 0)
                print "%s: %s"%(Stub, FTFlag)
            self.AnnotationCount += 1
            if FTFlag:
                self.FTAnnotationCount += 1
            if not len(Peptide.Modifications.keys()):
                continue # for now, ignore all unmodified peptides.
            ModIndex = Peptide.Modifications.keys()[0]
            ############################################################
            # Look up DBPos, the index (within self.DB) of the modified residue.
            if self.DBPosCache.has_key(Annotation):
                DBPos = self.DBPosCache[Annotation]
            else:
                DBPos = self.DB.find(Peptide.Aminos)
                if DBPos == -1:
                    print "** Warning: Didn't find peptide '%s' from line %s in the database"%(Peptide.Aminos, LineNumber)
                    continue
                DBPos += ModIndex
                self.DBPosCache[Annotation] = DBPos
            Mass = int(round(Peptide.Modifications[ModIndex][0].Mass))
            SiteKey = (DBPos, Mass)
            if not self.SiteDict.has_key(SiteKey):
                print "** Warning: Didn't find modification '%s' from line %s"%(Bits[2], LineNumber)
                print SiteKey
                continue
            ############################################################
            # Read the spectrum from the file, and save the mass.
            FilePath = os.path.join(self.SpectrumDir, FileName)
            FilePos = int(Bits[15])
            Site = self.SiteDict[SiteKey]
            try:
                SpectrumFile = open(FilePath, "rb")
                SpectrumFile.seek(FilePos)
                Spectrum = MSSpectrum.SpectrumClass()
                Spectrum.ReadPeaksFromFile(SpectrumFile, FileName)
                #print FilePath, FilePos, Spectrum.ParentMass, GetMass(Peptide.Aminos), Peptide.Aminos
                ParentMass = Spectrum.PrecursorMZ * Charge - (Charge - 1)*1.0078
                MassDifference = ParentMass - (GetMass(Peptide.Aminos) + 19)
                Site.Masses.append(MassDifference)
                if FTFlag:
                    Site.FTMasses.append(MassDifference)
##                
##                if FileLine[:8] != "PEPMASS=":
##                    FileLine = SpectrumFile.readline() # (M+H)+ and charge
##                if FileLine[:8] != "PEPMASS=":
##                    FileLine = SpectrumFile.readline() # (M+H)+ and charge
##                ParentMass = float(FileLine.split()[0])
##                MassDifference = ParentMass - (GetMass(Peptide.Aminos) + 19)
##                Site.Masses.append(MassDifference)
##                if FTFlag:
##                    Site.FTMasses.append(MassDifference)
            except:
                traceback.print_exc()
                print "** Error reading spectrum from %s:%s (line %s)"%(FilePath, FilePos, LineNumber)
    def OutputSitesVerbose(self):
        """
        Output EVERY MASS SEEN for every site.  This is a rather long-winded report!
        """
        Keys = self.SiteDict.keys()
        Keys.sort()
        OutputFile = open("MassDetails.txt", "wb")
        OutputFileFT = open("MassDetailsFT.txt", "wb")        
        for SiteKey in Keys:
            Site = self.SiteDict[SiteKey]
            DBPos = Site.Bits[1]
            ModMass = Site.Bits[2]
            AA = Site.Bits[3]
            if not AA.strip():
                continue
            BaseStr = "%s\t%s\t%s\t%s\t"%(Site.Bits[0], AA, ModMass, DBPos)
            for Mass in Site.Masses:
                Str = BaseStr + "%s\t"%Mass
                OutputFile.write(Str + "\n")
            for Mass in Site.FTMasses:
                Str = BaseStr + "%s\t"%Mass
                OutputFileFT.write(Str + "\n")
        OutputFile.close()
        OutputFileFT.close()
    def OutputSites(self):
        """
        Now that we've parsed in the sites, and parsed the mass offsets from
        the original spectra, let's output information on the mass offsets.
        """
        Keys = self.SiteDict.keys()
        Keys.sort()
        OutputFile = open(self.OutputFileName, "wb")
        for HeaderLine in self.HeaderLines:
            OutputFile.write(HeaderLine)
        for SiteKey in Keys:
            (DBPos, ModMass) = SiteKey
            Site = self.SiteDict[SiteKey]
            # Pad site-bits if necessary:
            while len(Site.Bits) < 26:
                Site.Bits.append("")
##            #########################################
##            # Verbose output:
##            Site.Masses.sort()
##            print 
##            print Site.Bits
##            for Mass in Site.Masses:
##                print "%s %.3f"%(int(round(Mass)), Mass)
##            ##################################
            # Find the most common mass (rounded to nearest integer):
            MonoisotopicDict = {}
            for Mass in Site.Masses:
                IntMass = int(round(Mass))
                MonoisotopicDict[IntMass] = MonoisotopicDict.get(IntMass, 0) + 1
            SortedDict = []
            for (Key, Count) in MonoisotopicDict.items():
                SortedDict.append((Count, Key))
            SortedDict.sort()
            SortedDict.reverse()
            if SortedDict:                
                MonoisotopicMass = SortedDict[0][1]
                NearMasses = []
                for Mass in Site.Masses:
                    if abs(Mass - MonoisotopicMass) < 0.75:
                        NearMasses.append(Mass)
                (Mean, StdDev) = BasicStats.GetMeanStdDev(NearMasses)
                Site.Bits.append("%.6f"%Mean)
                Site.Bits.append("%.6f"%StdDev)                
            else:
                MonoisotopicMass = ModMass
                Site.Bits.append("")
                Site.Bits.append("")
                StdDev = ""
            # Add bits for the average mass delta:
            FTScanCount = len(Site.FTMasses)
            Site.Bits.append(str(FTScanCount))
            if FTScanCount:
                NearMasses = []
                for Mass in Site.FTMasses:
                    if abs(Mass - MonoisotopicMass) < 0.75:
                        NearMasses.append(Mass)
                (Mean, StdDev) = BasicStats.GetMeanStdDev(NearMasses)
                Site.Bits.append("%.6f"%Mean)
                Site.Bits.append("%.6f"%StdDev)
                #Site.Bits.append(Str)
            else:
                Site.Bits.append("")
                Site.Bits.append("")
            Str = string.join(Site.Bits, "\t")
            OutputFile.write(Str + "\n")
        OutputFile.close()
    def Main(self):
        print "\nLoad DB and FT scan list..."
        self.LoadDB()
        self.LoadFTScanInfo()
        print "\nLoad PTM sites..."
        self.LoadSites()
        print "\nParse PTMAnnotations..."
        self.AnnotationCount = 0
        self.FTAnnotationCount = 0
        self.ParsePTMResultsFromFile(self.ParsePTMPath)
        print "\nOutput mass details..."
        self.OutputSitesVerbose()
        self.OutputSites()


UsageInfo = """
FTMass.py - measure accurate masses for PTMs.

Arguments:
 -r [ResultsFile]: File (or directory) to parse annotated sites from
 -f [FeatureFile]: Feature file (from TrainPTMFeatures) to read from
 -s [SpectrumDir]: Directory where .mzxml (or .mgf) files are stored
 -w [OutputFile]: File to write annotations to.
"""
if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(psyco not found - running in non-optimized mode)"
    Master = FTMassMaster()
    (Options, Args) = getopt.getopt(sys.argv[1:], "r:f:w:s:Qd:")
    OptionsSeen = {}
    for (Option, Value) in Options:
        OptionsSeen[Option] = 1
        if Option == "-r":
            # -r results file(s)
            if not os.path.exists(Value):
                print "** Error: couldn't find results file '%s'\n\n"%Value
                print UsageInfo
                sys.exit(1)
            Master.ParsePTMPath = Value
        elif Option == "-Q":
            Master.QuickParseFlag = 1
        elif Option == "-f":
            # -f feature file(s)
            if not os.path.exists(Value):
                print "** Error: couldn't find results file '%s'\n\n"%Value
                print UsageInfo
                sys.exit(1)
            Master.FeaturePath = Value
        elif Option == "-d":
            # -d database file(s)
            if not os.path.exists(Value):
                print "** Error: couldn't find database file '%s'\n\n"%Value
                print UsageInfo
                sys.exit(1)
            Master.DBPath = Value
        elif Option == "-s":
            # -s spectrum dir
            if not os.path.exists(Value):
                print "** Error: couldn't find database file '%s'\n\n"%Value
                print UsageInfo
                sys.exit(1)
            Master.SpectrumDir = Value
        elif Option == "-w":
            # -w output filen
            Master.OutputFileName = Value
    Master.Main()