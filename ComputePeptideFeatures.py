"""
Compute a collection of features describing *peptides* obtained from
a database search.  We'll determine which combination of features
can best separate the valid peptides from the invalid ones.
"""
import os
import sys
import struct
import traceback
import MSSpectrum
import PyInspect
import shutil
import math
import time
import ResultsParser
import getopt
import BuildConsensusSpectrum
from Utils import *
Initialize()

MAX_MODLESS_CLUSTER_SIZE = 50

class FormatBits:
    ProteinName = 1
    Peptide = 3
    Charge = 4
    TrueProteinFlag = 5
    FirstFeature = 7
    BestFScore = 9
    LastFeature = 21
    ModelScore = 22
    ModelPValue = 23
    FeatureCount = LastFeature - FirstFeature + 1

# helper classes are COPY-PASTA from ComputePTMFeatures
class PeptideSpecies:
    InstanceCount = 0
    # Peptides[(Annotation, Charge)] -> Peptide species
    def __init__(self):
        self.HitCount = 0
        self.ModifiedFlag = 0
        self.ModMass = 0
        self.DBPos = 0
        self.ModDBPos = 0
        self.Spectra = []
        self.Peptide = None
        self.SpectrumCount = 0
        PeptideSpecies.InstanceCount += 1
    def __del__(self):
        if PeptideSpecies:
            PeptideSpecies.InstanceCount -= 1
    def __str__(self):
        return self.Annotation

class PeptideFeatureComputer(ResultsParser.ResultsParser, ResultsParser.SpectrumOracleMixin):
    def __init__(self):
        self.PValueCutoff = 0.25 # default
        self.OutputPath = "PeptideFeatures.txt" # default
        self.TempFileDir = "PeptideScoring"
        self.CachedFixedFilePaths = []
        self.CachedFilePaths = []
        self.Peptides = {}
        self.QuickParseFlag = 0
        self.MZXMLOracle = {}
        self.SpectrumDir = None
        self.AllLineCount = 0
        self.MinDBPos = None
        self.MaxDBPos = None
        self.DBPath = None
        ResultsParser.ResultsParser.__init__(self)
        ResultsParser.SpectrumOracleMixin.__init__(self)
    def GetValidProteinFlag(self, Species):
        # Normally we prepend "xxx" to the bogus names:
        if Species.ProteinName[:3] == "XXX":
            return 0
        return 1
    def RememberString(self, StringList, NewString):
        """
        Return the index of NewString within StringList, adding to the list if necessary.
        We keep a list of mzxml file names and store indexes into the list, to avoid
        the memory hit required to store each occurrence of the name.
        """
        try:
            Index = StringList.index(NewString)
            return Index
        except:
            StringList.append(NewString)
            return len(StringList) - 1
    def ParseCommandLine(self):
        (Options, Args) = getopt.getopt(sys.argv[1:], "d:r:w:s:M:t:lp:X:Y:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-r":
                # -r results file(s)
                self.ResultsFileName = Value
            elif Option == "-d":
                self.DBPath = Value
            elif Option == "-w":
                self.OutputPath = Value
            elif Option == "-l":
                self.QuickParseFlag = 1
            elif Option == "-s":
                self.SpectrumDir = Value
            elif Option == "-M":
                self.PopulateSpectrumOracle(Value)
            elif Option == "-t":
                self.TempFileDir = Value
            elif Option == "-p":
                self.PValueCutoff = float(Value)
            elif Option == "-X":
                self.MinDBPos = int(Value)
            elif Option == "-Y":
                self.MaxDBPos = int(Value)
            else:
                print "* Error: Unrecognized option %s"%Option
        if not self.DBPath:
            return 0
        if not self.ResultsFileName:
            return 0
        return 1
    def WipeDir(self, Dir):
        try:
            shutil.rmtree(Dir)
        except:
            pass
    def LoadDB(self):
        """
        Load the database searched.  For future reference, we want the protein names as well.
        """
        # Populate self.DB with the contents of the .trie file
        File = open(self.DBPath, "rb")
        self.DB = File.read()
        File.close()
        # Populate self.ProteinNames and self.ProteinPositions by parsing the index file:
        print "Populate protein names..."
        self.ProteinNames = []
        self.ProteinPositions = []
        IndexPath = os.path.splitext(self.DBPath)[0] + ".index"
        File = open(IndexPath, "rb")
        BlockSize = struct.calcsize("<qi80s")
        while 1:
            Block = File.read(BlockSize)
            if not Block:
                break
            Tuple = struct.unpack("<qi80s", Block)
            Name = Tuple[-1]
            NullPos = Name.find("\0")
            if NullPos != -1:
                Name = Name[:NullPos]
            self.ProteinNames.append(Name)
            self.ProteinPositions.append(Tuple[1])
        File.close()
        print "Protein names populated."
    def ParseResultsFile(self, FilePath):
        if os.path.isdir(FilePath):
            print "NOTE: Skipping results sub-directory '%s'"%FilePath
            return
        LineNumber = 0
        OldSpectrum = None
        File = open(FilePath, "rb")
        for FileLine in File.xreadlines():
            self.AllLineCount += 1
            LineNumber += 1
            if LineNumber % 5000 == 0:
                print "  Line %s..."%LineNumber
                if self.QuickParseFlag:
                    break
            if FileLine[0] == "#":
                continue
            Bits = FileLine.strip().split("\t")
            if len(Bits) < self.Columns.PValue:
                continue # not valid!
            Spectrum = (Bits[0], Bits[1])
            if Spectrum == OldSpectrum:
                continue
            OldSpectrum = Spectrum
            PValue = float(Bits[self.Columns.PValue])
            Annotation = Bits[self.Columns.Annotation]
            Charge = int(Bits[self.Columns.Charge])
            AnnotationKey = (Annotation, Charge)
            ApproxDBPos = min(len(self.DB), int(Bits[self.Columns.DBPos]) + 200)
            DBPos = self.DB.rfind(Annotation[2:-2], 0, ApproxDBPos)
            if DBPos == -1:
                print "* warning: peptide '%s' isn't near %s (%s)"%(Annotation, Bits[self.Columns.DBPos], ApproxDBPos)
                DBPos = self.DB.rfind(Annotation[2:-2])
            if self.MaxDBPos != None:
                if DBPos > self.MaxDBPos or DBPos < self.MinDBPos:
                    continue
            ##############################################################
            # If we've never seen this annotation before, then create a PeptideSpecies object
            # and record it in self.Peptides
            Species = self.Peptides.get(AnnotationKey, None)
            if not Species:
                Species = PeptideSpecies()
                Species.Peptide = GetPeptideFromModdedName(Annotation)
                Species.ProteinName = Bits[self.Columns.ProteinName]
                self.Peptides[AnnotationKey] = Species
                # Get the database position of the peptide:
                ApproxDBPos = max(0, int(Bits[self.Columns.DBPos]) - 1000)
                Species.DBPos = self.DB.find(Species.Peptide.Aminos, ApproxDBPos)
                # Get the residue-number of the peptide:
                StarPos = self.DB.rfind("*", 0, Species.DBPos)
                if StarPos == -1:
                    Species.ResidueNumber = Species.DBPos
                else:
                    Species.ResidueNumber = Species.DBPos - StarPos
                Species.Annotation = Annotation
                Species.Charge = Charge
            if Species.DBPos == -1:
                print "* skipping unknown peptide: %s"%Annotation
                continue
            MQScore = float(Bits[self.Columns.MQScore])
            FScore = float(Bits[self.Columns.FScore])
            self.AnnotationCount += 1
            # Store TUPLES for the best spectra:
            FileName = Bits[0].replace("/","\\").split("/")[-1]
            FileNameIndex = self.RememberString(self.CachedFilePaths, FileName)
            DeltaScore = float(Bits[self.Columns.DeltaScore])
            ByteOffset = int(Bits[self.Columns.FileOffset])
            Tuple = (-FScore, -MQScore, -DeltaScore, FileNameIndex, ByteOffset)
            Species.Spectra.append(Tuple)
            Species.SpectrumCount += 1
            if len(Species.Spectra) > MAX_MODLESS_CLUSTER_SIZE:
                Species.Spectra.sort()
                Species.Spectra = Species.Spectra[:MAX_MODLESS_CLUSTER_SIZE]
            continue
        File.close()
        print "...After parsing %s, we know of %s peptides"%(FilePath, len(self.Peptides.keys()))
    def FixSpectrumPath(self, Path):
        """
        Given a spectrum file-path, return the correct path to the .mzXML file
        on local disk.
        """
        # (Useful if you run searches on the grid in temp-folders)
        Bits = Path.replace("/", "\\").split("\\")
        FileName = Bits[-1]
        # handle -s option:
        if self.SpectrumDir:
            return os.path.join(self.SpectrumDir, FileName)
        # handle -M option:
        Stub = os.path.splitext(FileName)[0]
        SpectrumPath = self.SpectrumOracle.get(Stub, None)
        if not SpectrumPath:
            print "** Warning: Couldn't figure out spectrum path for '%s'"%Path
        return SpectrumPath
    def PrepareClusterDirectories(self):
        self.ConsensusClusterDir = os.path.join(self.TempFileDir, "Clusters")
        self.ConsensusSpectrumDir = os.path.join(self.TempFileDir, "Spectra")
        print "Prepare cluster directories..."
        if self.MaxDBPos == None or self.MinDBPos == 0:
            self.WipeDir(self.ConsensusClusterDir)
            self.WipeDir(self.ConsensusSpectrumDir)
        MakeDirectory(self.ConsensusClusterDir)
        MakeDirectory(self.ConsensusSpectrumDir)
        for AA in "ACDEFGHIKLMNOPQRSTUVWY":
            PathA = os.path.join(self.ConsensusClusterDir, AA)
            PathB = os.path.join(self.ConsensusSpectrumDir, AA)
            for Path in (PathA, PathB):
                try:
                    os.makedirs(Path)
                except:
                    pass
        print "...cluster directories prepared."
    def ComputeFeaturesMain(self):
        self.PrepareClusterDirectories()
        print "Load database..."
        self.LoadDB()
        print "Parse annotations..."
        self.AnnotationCount = 0
        self.BestPeptideHits = {}
        self.ProcessResultsFiles(self.ResultsFileName, self.ParseResultsFile)
        # Fix file paths:
        for FilePath in self.CachedFilePaths:
            FixedPath = self.FixSpectrumPath(FilePath)
            self.CachedFixedFilePaths.append(FixedPath)
        self.ProduceConsensusSpectra()
        # Compute and output features:
        self.ComputeFeaturesAllPeptides()
    def OutputInfoHeader(self):
        Header = "#DBPos\tProtein\tResidue\tPeptide\tCharge\tValidFlag\tBestSpectrum\tLogSpectrumCount\tLogLength\t"
        Header += "BestFScore\tBestMQScore\tBestDeltaScore\tConsensusMQ\t"
        Header += "Length\tTotalPRMScore\tMedianPRMScore\tFractionY\tFractionB\tIntensity\tNTT\t"
        Header += "LogSpecCount\tLogLength\t" 
        self.OutputFile.write(Header + "\n")
        Header = "#0\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\t12\t13\t14\t15\t16\t17\t18\t19\t20\t21\t22\t23"
        self.OutputFile.write(Header + "\n")
        Header = "#Feature\t\t\t\t\t\t\t0\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\t12\t13"
        self.OutputFile.write(Header + "\n")
    def ComputeFeaturesAllPeptides(self):
        """
        Calculate, and report, all the features for each peptide
        """
        self.OutputFile = open(self.OutputPath, "wb")
        self.OutputInfoHeader()
        Keys = self.Peptides.keys()
        Keys.sort()
        for KeyIndex in range(len(Keys)):
            Key = Keys[KeyIndex]
            Species = self.Peptides[Key]
            print "(%s/%s) %s"%(KeyIndex, len(Keys), Species.Peptide.Aminos)
            try:
                Features = self.ComputeFeatures(Species)
            except:
                traceback.print_exc()
                print "** Error: Unable to compute features for %s"%Species
                continue
            Str = "%s\t%s\t%s\t"%(Species.DBPos, Species.ProteinName, Species.ResidueNumber) 
            Str += "%s\t"%Species.Annotation
            Str += "%s\t"%Species.Charge
            for Feature in Features:
                Str += "%s\t"%Feature
            print Str
            self.OutputFile.write(Str + "\n")
            # We're done with this PTM now, so let's forget about it:
            del self.Peptides[Key]
        self.OutputFile.close()
    def ComputeFeatures(self, Species):
        """
        Compute scoring-features for this peptide species, return them as a list
        """
        Features = []
        # Feature: Is the PTM from a valid protein?  (Note: this feature is not INPUT for the
        # model, it's our desired output)
        Feature = self.GetValidProteinFlag(Species)
        Features.append(Feature)
        # The best spectrum observed (meta-data, not a scoring feature)
        BestFScore = -999
        BestMQScore = -999
        BestDeltaScore = None
        for Tuple in Species.Spectra:
            FScore = -Tuple[0]
            MQScore = -Tuple[1]
            if FScore > BestFScore:
                BestFScore = FScore
                BestDeltaScore = -Tuple[2]
                FilePath = self.CachedFixedFilePaths[Tuple[3]]
                BestFSpectrum = "%s:%s"%(FilePath, Tuple[4])
            if MQScore > BestMQScore:
                BestMQScore = -Tuple[1]
        Features.append(BestFSpectrum)
        # TRUE FEATURES BEGIN NOW:
        Features.append(math.log(Species.SpectrumCount))
        Features.append(math.log(len(Species.Peptide.Aminos)))
        Features.append(BestFScore)
        Features.append(BestMQScore)
        Features.append(BestDeltaScore)
        # Feature: Consensus annotation score (and score-features) for this peptide
        Species.ConsensusScore = None
        self.GetConsensusMQScore(Species, Features)
        Features.append(math.log(Species.SpectrumCount + 1))
        Features.append(math.log(len(Species.Peptide.Aminos) - 5))
        return Features
    def GetConsensusMQScore(self, Species, Features):
        # Load in the consensus spectrum, and score the peptide annotation:
        try:
            PySpectrum = PyInspect.Spectrum(Species.ConsensusPath, 0)
            Species.PySpectrum = PySpectrum
            ScoreList = PySpectrum.ScorePeptideDetailed(Species.Annotation)
            Species.ConsensusScore = ScoreList[0]
            for ScoreItem in ScoreList[:8]:
                Features.append(ScoreItem)
            print "PyInspect score %s -> %s"%(Species.Annotation, ScoreList[0])
        except:
            traceback.print_exc()
            for DummyIndex in range(8):
                Features.append(0)
    def ProduceConsensusSpectra(self):
        """
        Build a consensus spectrum for each peptide species, using its cluster
        """
        Keys = self.Peptides.keys()
        for PeptideIndex in range(len(Keys)):
            if (PeptideIndex % 100 == 0):
                print "For peptide %s/%s..."%(PeptideIndex, len(Keys))
            AnnotationKey = Keys[PeptideIndex]
            (Annotation, Charge) = AnnotationKey
            Species = self.Peptides[AnnotationKey]
            Species.ClusterPath = os.path.join(self.ConsensusClusterDir, Annotation[2], "%s.%s.cls"%(Annotation.replace("*", "-"), Charge))
            Species.ConsensusPath = os.path.join(self.ConsensusSpectrumDir, Annotation[2], "%s.%s.dta"%(Annotation.replace("*", "-"), Charge))
            Builder = BuildConsensusSpectrum.ConsensusBuilder(Species.Charge)
            MeanMQ = 0
            for Tuple in Species.Spectra:
                MeanMQ -= Tuple[0]
            MeanMQ /= float(len(Species.Spectra))
            for Tuple in Species.Spectra:
                # Omit from the consensus spectra with very poor scores:
                MQScore = -Tuple[0]
                if MQScore < MeanMQ - 3.0:
                    continue
                SpectrumFilePath = self.CachedFixedFilePaths[Tuple[3]]
                Spectrum = MSSpectrum.SpectrumClass()
                SpectrumFile = open(SpectrumFilePath, "rb")
                SpectrumFile.seek(Tuple[4])
                try:
                    Spectrum.ReadPeaksFromFile(SpectrumFile, SpectrumFilePath)
                except:
                    print "* Error parsing %s:%s"%(SpectrumFilePath, Tuple[4])
                    raise
                Spectrum.SetCharge(Charge)
                SpectrumFile.close()
                Builder.AddSpectrum(Spectrum)
                # Special (and easy) case: If we only saw one spectrum, then write it
                # out without changing it!
                if len(Species.Spectra) == 1:
                    Spectrum.WritePeaks(Species.ConsensusPath)
            if len(Species.Spectra) > 1:
                Spectrum = Builder.ProduceConsensusSpectrum()
                Spectrum.WritePeaks(Species.ConsensusPath)
    def PopulateMZXMLOracle(self, RootDirectory):
        """
        Used when mzxml files are spread over multiple subdirectories.
        MZXMLOracle[Stub] = full path to the corresponding MZXML file
        Used with -M option (not with -s option)
        """
        if not os.path.exists(RootDirectory):
            return
        for SubFileName in os.listdir(RootDir):
            SubFilePath = os.path.join(RootDir, SubFileName)
            if os.path.isdir(SubFilePath):
                self.PopulateMZXMLOracle(SubFilePath)
                continue
            (Stub, Extension) = os.path.splitext(SubFileName)
            if Extension.lower() == ".mzxml":
                self.MZXMLOracle[Stub] = os.path.join(RootDirectory, SubFileName)

UsageInfo = """
ComputePeptideFeatures: Generate feature values for PTMs observed on a data-set.

Arguments:
 -r [ResultsFile]: Name of the results file (or directory)
 -d [DBPath]: Path to the .trie file searched
 -w [OutputFile]: Output file name.  Defaults to PeptideFeatures.txt
 -s [Directory]: Spectrum directory.  (Use either -s or -M)
 -M [RootDir]: Root directory for mzXML files. 
 -t [TempFileDir]: Write consensus spectra, etc. to the specified directory
 -X DBStart
 -Y DBEnd
"""        

if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(psyco not installed - no optimization for you!)"
    Trainer = PeptideFeatureComputer()
    Result = Trainer.ParseCommandLine()
    if not Result:
        print UsageInfo
    else:
        Trainer.ComputeFeaturesMain()