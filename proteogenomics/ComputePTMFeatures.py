"""
Plan:
Output a large collection of features for each post-translational modification accepted on a
search of a part-bogus database.  All modifications on the bogus proteins are incorrect.
An equivalent number of modifications on the non-bogus proteins are incorrect.  Let's compute
a variety of features for the PTMs observed.

Input:
A collection of annotated spectra, output by SelectSites.py
Output:
A file listing all the observed modification sites, with various features computed.  

Then, we train a model to distinguish between good (correct DB) and bad (incorrect DB)
modifications.  Model types: LDA, logistic regression, SVM, etc.

(Another possible experiment: Search unmodified spectra against a mutated database,
judge correct precisely those modifications which "undo" the mutations)
"""

import os
import sys
import struct
import traceback
import getopt
import MSSpectrum
import PyInspect
import random
import shutil
import time
import math
import cPickle
import BasicStats
import ResultsParser
import BuildConsensusSpectrum
import SpectralSimilarity
import StripPTM
random.seed(1)
from Utils import *
from TrainPTMFeatures import FormatBits
Initialize()

AMINO_ACIDS = "ACDEFGHIKLMNOPQRSTUVWY" # O and U are included, for now.
INVALID_MASS = 99999

# Retain at most this many spectra for an unmodified peptide:
MAX_MODLESS_CLUSTER_SIZE = 100

# For running the Python profiler:
PROFILE_FLAG = 0

class PeptideSpecies:
    """
    Peptides[(Annotation, Charge)] -> PeptideSpecies instance
    The PeptideSpecies remembers a list of spectra, the modification position, the modification mass.
    """
    InstanceCount = 0
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

class SpectrumInfoClass(ResultsParser.ResultsParser):
    """
    Information about a single scan.  We remember only the info we'll need later.
    """
    InstanceCount = 0
    def __init__(self, Bits, Trainer):
        ResultsParser.ResultsParser.__init__(self)
        self.FileNameIndex = Trainer.RememberString(Trainer.CachedFilePaths, Bits[0])
        self.MQScore = float(Bits[self.Columns.MQScore])
        self.DeltaScore = float(Bits[self.Columns.DeltaScore])
        self.ByteOffset = int(Bits[self.Columns.FileOffset])
        SpectrumInfoClass.InstanceCount += 1
    def __cmp__(self, Other):
        """
        Sort from BEST to WORST match.
        """
        if self.MQScore > Other.MQScore:
            return -1
        if self.MQScore < Other.MQScore:
            return 1
        return 0
    def __del__(self):
        if SpectrumInfoClass:
            SpectrumInfoClass.InstanceCount -= 1
        
class PTMFeatureComputer(ResultsParser.ResultsParser, ResultsParser.SpectrumOracleMixin):
    def __init__(self):
        self.PValueCutoff = 0.1 # default
        self.ResultsFileName = None
        self.DBPath = None
        self.OutputDir = "PTM"
        self.OutputPath = os.path.join(self.OutputDir, "PTMFeatures.txt")
        self.ConsensusClusterDir = None
        self.ConsensusSpectrumDir = None
        # Peptides keeps a list of SpectrumInfo objects for each peptide species
        # we've observed.  The keys have the form (Annotation, Charge) and the values
        # are lists SpectrumInfo instance.
        self.Peptides = {}
        self.PTMs = {} # keys of the form (DBPos, Mass)
        self.CoverageThreshold = 2 # at least this many spectra to consider a residue 'covered'.
        self.QuickParseFlag = 0 # if true, then parse only the first n lines
        self.PoolFlag = 0
        self.ModelType = None
        self.SisterProteins = {} # protein index -> sister protein's index
        self.MZXMLOracle = {}
        self.ModelTrainFilePath = "PTMFeatures.All.txt"
        self.ModelTestFilePath = None
        # Dictionary of unmodified peptides, for computing the coverage level:
        self.UnmodifiedPeptides = {}
        self.FeatureSelectionFlag = None
        self.CachedFilePaths = []
        self.CachedFixedFilePaths = []
        self.StartOutputDBPos = 0
        self.RequiredFileNameChunk = None
        ResultsParser.ResultsParser.__init__(self)
        ResultsParser.SpectrumOracleMixin.__init__(self)
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
    def LoadDB(self):
        """
        Load the database searched.  For future reference, we want the protein names as well.
        """
        # Populate self.DB with the contents of the .trie file
        File = open(self.DBPath, "rb")
        self.DB = File.read()
        File.close()
        # Populate self.ProteinNames and self.ProteinPositions by parsing the index file:
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
        # Initialize our coverage arrays:
        self.Coverage = [0] * len(self.DB)
        self.ModCoverage = [0] * len(self.DB)
        self.PeptideCoverage = [0] * len(self.DB)
        self.ModPeptideCoverage = [0] * len(self.DB)
        # Find sister-proteins.
        # A shuffled protein name should be the standard protein's name
        # with the characters "XXX" prepended.  (If the standard protein name
        # is very long, the last characters may "slide off the edge")
        for IndexA in range(len(self.ProteinNames)):
            Name = self.ProteinNames[IndexA]
            SisterName = None
            if Name[:3] == "XXX":
                SisterName = Name[3:]
            if SisterName:
                for IndexB in range(len(self.ProteinNames)):
                    Name = self.ProteinNames[IndexB][:77]
                    if Name == SisterName:
                        self.SisterProteins[IndexA] = IndexB
                        self.SisterProteins[IndexB] = IndexA
    def ComputeProteinCoverage(self):
        """
        Compute residue-level coverage.  And, compute
        what fraction of each protein is covered. 
        """
        for Species in self.Peptides.values():
            DBPos = Species.DBPos
            if Species.ModifiedFlag:
                for Pos in range(DBPos, DBPos + len(Species.Peptide.Aminos)):
                    self.ModCoverage[Pos] += Species.SpectrumCount # count modified spectra
                    self.ModPeptideCoverage[Pos] += 1 # count distinct peptide species
            else:
                for Pos in range(DBPos, DBPos + len(Species.Peptide.Aminos)):
                    self.Coverage[Pos] += Species.SpectrumCount # count unmodified spectra
                    self.PeptideCoverage[Pos] += 1 # count distinct unmodified peptide species
        #########################################################
        # Compute percentage of each protein that's covered:
        self.ProteinCoverageLevels = []
        for ProteinIndex in range(len(self.ProteinPositions)):
            StartPos = self.ProteinPositions[ProteinIndex]
            if ProteinIndex < len(self.ProteinPositions) - 1:
                EndPos = self.ProteinPositions[ProteinIndex + 1]
            else:
                EndPos = len(self.DB)
            #print "Protein %s (%s) from %s-%s"%(ProteinIndex, self.ProteinNames[ProteinIndex], StartPos, EndPos)
            CoverFlags = 0
            for Pos in range(StartPos, EndPos):
                if self.Coverage[Pos] >= self.CoverageThreshold:
                    CoverFlags += 1
            ProteinLength = EndPos - StartPos
            #print "  -> Coverage %s/%s = %s"%(CoverFlags, ProteinLength, CoverFlags / float(ProteinLength))
            self.ProteinCoverageLevels.append(CoverFlags / float(ProteinLength))
        # SAVE protein coverage levels:
        CoveragePath = os.path.join(self.OutputDir, "Coverage.dat")
        CoverageFile = open(CoveragePath, "wb")
        for DBPos in range(len(self.DB)):
            Str = struct.pack("<II", self.Coverage[DBPos], self.ModCoverage[DBPos])
            CoverageFile.write(Str)
        CoverageFile.close()
        # Boost protein coverage levels based upon sister proteins:
        for Index in range(len(self.ProteinNames)):
            SisterIndex = self.SisterProteins.get(Index, None)
            if SisterIndex == None:
                continue
            self.ProteinCoverageLevels[Index] = max(self.ProteinCoverageLevels[Index], self.ProteinCoverageLevels[SisterIndex])
            #print "%s and %s are sisters, with coverage %s"%(self.ProteinNames[Index], self.ProteinNames[SisterIndex], self.ProteinCoverageLevels[Index])
    def FixPeptideSpecies(self):
        """
        Iterate over all the peptide species we observed.  Strip "obviously unnecessary" PTMs.
        """
        Keys = self.Peptides.keys()
        for Key in Keys:
            Annotation = Key[0]
            Species = self.Peptides[Key]
            Result = StripPTM.StripNeedlessModifications(self.DB, Annotation)
            if not Result:
                continue
            (DBPos, FixedAnnotation) = Result
            # If the annotation wasn't changed, then continue.
            if FixedAnnotation == Annotation:
                continue
            Species.Peptide = GetPeptideFromModdedName(FixedAnnotation)
            Species.Annotation = FixedAnnotation
            Species.DBPos = DBPos
            ModKeys = Species.Peptide.Modifications.keys()
            if len(ModKeys):
                ModIndex = ModKeys[0]
                Species.ModifiedFlag = 1
                Species.ModMass = int(round(Species.Peptide.Modifications[ModIndex][0].Mass))
                Species.ModDBPos = DBPos + ModIndex
            else:
                Species.ModifiedFlag = 0
                Species.ModMass = 0
                Species.ModDBPos = None
            del self.Peptides[Key]
            # Either merge into the existing species with this fixed annotation,
            # or move into the empty pigeonhole:
            FixedKey = (FixedAnnotation, Key[1])
            OldSpecies = self.Peptides.get(FixedKey, None)
            if OldSpecies:
                OldSpecies.Spectra.extend(Species.Spectra)
                OldSpecies.SpectrumCount += Species.SpectrumCount
                if len(OldSpecies.Spectra) > MAX_MODLESS_CLUSTER_SIZE:
                    OldSpecies.Spectra.sort()
                    OldSpecies.Spectra = Species.Spectra[:MAX_MODLESS_CLUSTER_SIZE]
            else:
                self.Peptides[FixedKey] = Species
                if len(Species.Spectra) > MAX_MODLESS_CLUSTER_SIZE:
                    Species.Spectra.sort()
                    Species.Spectra = Species.Spectra[:MAX_MODLESS_CLUSTER_SIZE]
    def ParsePTMsFromResultsFile(self, FilePath):
        """
        Callback for parsing one Inspect results-file.  Our job here
        is to populate self.BestSpectra, and self.PTMs.
        Note: It's POSSIBLE that we'll spot some modified annotations
        which can be "trivially fixed" to produce unmodified annotations.
        Examples: T.T+101LAPTTVPITSAK.A, Y.E+163NPNFTGK.K
        Because of this, we keep a dictionary self.FixedAnnotation,
        where keys are raw annotations and values are fixed-up annotations.
        """
        if self.RequiredFileNameChunk:
            Pos = FilePath.find(self.RequiredFileNameChunk)
            if Pos == -1:
                return
        if os.path.isdir(FilePath):
            print "NOTE: Skipping results sub-directory '%s'"%FilePath
            return
        try:
            File = open(FilePath, "rb")
        except:
            print "** Unable to open results file '%s'"%FilePath
            return
        LineNumber = 0
        OldSpectrum = None
        for FileLine in File.xreadlines():
            LineNumber += 1
            if LineNumber % 1000 == 0:
                print "  Line %s..."%LineNumber
                if self.QuickParseFlag:
                    break
            if FileLine[0] == "#":
                continue
            Bits = FileLine.strip().split("\t")
            if len(Bits) < 15:
                continue # not valid!
            Spectrum = (Bits[0], Bits[1])
            if Spectrum == OldSpectrum:
                continue
            OldSpectrum = Spectrum
            PValue = float(Bits[self.Columns.PValue])
            if PValue > self.PValueCutoff:
                continue
            Annotation = Bits[self.Columns.Annotation]
            Charge = int(Bits[self.Columns.Charge])
            AnnotationKey = (Annotation, Charge)
            ##############################################################
            # If we've never seen this annotation before, then create a PeptideSpecies object
            # and record it in self.Peptides
            Species = self.Peptides.get(AnnotationKey, None)
            if not Species:
                Species = PeptideSpecies()
                Species.Peptide = GetPeptideFromModdedName(Annotation)
                Mods = []
                for (Index, List) in Species.Peptide.Modifications.items():
                    for Mod in List:
                        Mods.append((Index, Mod))
                if len(Mods):
                    Species.ModifiedFlag = 1
                    Species.ModMass = int(Mods[0][1].Mass)
                    Species.ModAA = Species.Peptide.Aminos[Mods[0][0]]
                else:
                    Species.ModifiedFlag = 0
                Species.ProteinName = Bits[self.Columns.ProteinName]
                self.Peptides[AnnotationKey] = Species
                # Get the database position of the peptide:
                Species.DBPos = self.DB.find(Species.Peptide.Aminos)
                if len(Species.Peptide.Modifications.keys()):
                    ModIndex = Species.Peptide.Modifications.keys()[0]
                    Species.ModDBPos = Species.DBPos + ModIndex
                else:
                    Species.ModDBPos = None
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
                del self.Peptides[AnnotationKey] # remove the Species that was just created!
                continue
            MQScore = float(Bits[self.Columns.MQScore])
            self.AnnotationCount += 1
            ##############################################################
            # Populate Species.Spectra:
            try:
                Info = SpectrumInfoClass(Bits, self)
            except:
                print "** Error: Couldn't parse spectrum info from line %s of file %s"%(LineNumber, FilePath)
                traceback.print_exc()
                continue
            Species.Spectra.append(Info)
            Species.SpectrumCount += 1
            if not Species.ModifiedFlag:
                if len(Species.Spectra) > MAX_MODLESS_CLUSTER_SIZE:
                    Species.Spectra.sort()
                    Species.Spectra = Species.Spectra[:MAX_MODLESS_CLUSTER_SIZE]
            else:
                pass
        File.close()
    def WipeDir(self, Dir):
        try:
            shutil.rmtree(Dir)
        except:
            pass
    def ComputeFeaturesMain(self):
        """
        Main method:
        - Load the database searched
        - Iterate over the results-file, to get a list of PTMs
        - Iterate over the PTMs, and write out features for each one.
        """
        self.ConsensusClusterDir = os.path.join(self.OutputDir, "Clusters")
        self.ClusterScanListDir = os.path.join(self.OutputDir, "ClusterMembers")
        self.ConsensusSpectrumDir = os.path.join(self.OutputDir, "Spectra")
        if not self.StartOutputDBPos:
            # Make sure necessary directories exist, and clean up any OLD output:
            print "Prepare cluster directories..."
            self.WipeDir(self.ConsensusClusterDir)
            self.WipeDir(self.ConsensusSpectrumDir)
            self.WipeDir(self.ClusterScanListDir)
            MakeDirectory(self.ConsensusClusterDir)
            MakeDirectory(self.ConsensusSpectrumDir)
            MakeDirectory(self.ClusterScanListDir)
            for AA in AMINO_ACIDS:
                PathA = os.path.join(self.ConsensusClusterDir, AA)
                PathB = os.path.join(self.ConsensusSpectrumDir, AA)
                PathC = os.path.join(self.ClusterScanListDir, AA)
                for Path in (PathA, PathB, PathC):
                    MakeDirectory(Path)
        else:
            print "CONTINUING ComputePTMFeatures from DBPosition %s"%self.StartOutputDBPos
        print "Load database..."
        self.LoadDB()
        print "Parse annotations..."
        self.AnnotationCount = 0
        self.PTMAnnotationCount = 0
        self.BestModlessHits = {}
        self.ProcessResultsFiles(self.ResultsFileName, self.ParsePTMsFromResultsFile)
        # Fix annotations:
        self.FixPeptideSpecies()
        # Fix file paths:
        for FilePath in self.CachedFilePaths:
            FixedPath = self.FixSpectrumPath(FilePath)
            self.CachedFixedFilePaths.append(FixedPath)
        self.PairModifiedUnmodifiedPeptides()
        print "Produce CONSENSUS SPECTRA for modified and unmodified petpides..."
        StartTime = time.clock()
        self.ProduceConsensusSpectra()
        EndTime = time.clock()
        print "Elapsed time: %s"%(EndTime - StartTime)
        print "Compute protein coverage..."
        self.ComputeProteinCoverage()
        print "Count spectra (and sites) by PTM type..."
        self.ComputeTotalSpectraForModType()
        print "Generate non-redundant PTM list..."
        self.ListDistinctPTMs()
        print "Compute features and output PTM info..."
        self.ComputeFeaturesAllPTMs()
    def PairModifiedUnmodifiedPeptides(self):
        for Species in self.Peptides.values():
            if not Species.ModifiedFlag:
                continue
            ModlessAnnotation = "%s.%s.%s"%(Species.Peptide.Prefix, Species.Peptide.Aminos, Species.Peptide.Suffix)
            ModlessKey = (ModlessAnnotation, Species.Charge)
            ModlessSpecies = self.Peptides.get(ModlessKey, None)
            Species.Modless = ModlessSpecies
    def ListDistinctPTMs(self):
        """
        Populate self.PTMs; keys are (DBPos, Mass) and values are simple objects
        with lists of peptide species.
        """
        for Species in self.Peptides.values():
            if not Species.ModifiedFlag:
                continue
            Index = Species.Peptide.Modifications.keys()[0]
            ModifiedPos = Species.DBPos + Index
            Key = (ModifiedPos, Species.ModMass)
            if not self.PTMs.has_key(Key):
                PTM = Bag()
                PTM.SpeciesList = []
                self.PTMs[Key] = PTM
                ModIndex = Species.Peptide.Modifications.keys()[0]
                PTM.DBPos = Species.ModDBPos
            PTM.SpeciesList.append(Species)
            Species.PTM = PTM
    def ComputeTotalSpectraForModType(self):
        """
        Populate a dictionary of the form (AA, Mass) -> SpectrumCount.  If a modificaiton is seen
        at multiple sites, it is more likely to be valid.
        """
        self.ModTypeSpectrumCount = {}
        self.ModTypeSiteCount = {}
        for Species in self.Peptides.values():
            if Species.ModifiedFlag:
                ModTypeKey = (Species.ModAA, Species.ModMass)
                self.ModTypeSpectrumCount[ModTypeKey] = self.ModTypeSpectrumCount.get(ModTypeKey, 0) + Species.SpectrumCount
                self.ModTypeSiteCount[ModTypeKey] = self.ModTypeSiteCount.get(ModTypeKey, 0) + 1
        pass
    def OutputPTMInfoHeader(self):
        """
        Output column headers, plus some general-purpose information such as the number
        of spectra parsed and the database size.
        """
        Header = "#Group\tDBPosition\tMass\tAminoAcid\tProtein\tResidueNumber\t"
        Header += "Peptide\tCharge\tValidProteinFlag\tFacultativeFlag\tBestSpectrum\t"
        Header += "BestModlessSpectrum\tBestModlessMQScore\tBigDBAnn\tBigDBScore\tSpectra\tModlessSpectra\tBestMQScore\t"
        Header += "BestDeltaScore\tPeptideCount\tConsensusMQScore\tPeptideLength\tCutScoreTotal\t"
        Header += "MedianCutScore\tYPresent\tBPresent\tBYIntensity\tNTT\tModdedFraction\tProteinCoverage\t"
        Header += "SpectraThisModType\tSitesThisModType\tUnmodifiedPeptideCount\tDot0.5\tShared01\tShared11\t"
        Header += "Correlation\tLogSpectrumCount\tLogPeptideLength\tLogSpecThisType\tLogSitesThisType\t"
        Header += "DeltaVsBigDB\tModelScore\tModelPValue\tSitePValue\tKnownModType\tKnownModAnnotation\t"
        Header += "KnownModScore\tKnownModSitePValue\t"
        self.OutputFile.write(Header + "\n")
        # Two more header lines, for feature-numbers and column-numbers:
        Header = "#0\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\t12\t13\t14\t15\t16\t17\t18\t19\t20\t21\t22\t23\t24\t25\t26\t27\t28\t29\t30\t31\t32\t33\t34\t35\t36\t37\t38\t39\t40\t41\t42\t43\t44\t45\t46\t47\t48\t49\t"
        self.OutputFile.write(Header + "\n")
        Header = "#Feature\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t0\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\t12\t13\t14\t15\t16\t17\t18\t19\t20\t21\t22\t23\t24\t25\t26\t"
        self.OutputFile.write(Header + "\n")
        ProteinRecordCount = self.DB.count("*")
        DBResidueSize = len(self.DB) - ProteinRecordCount
        self.OutputFile.write("#DatabaseSize\t%s\t\n"%DBResidueSize)
        self.OutputFile.write("#AnnotationCount\t%s\t\n"%self.AnnotationCount)
        SiteCount = len(self.PTMs.values())
        self.OutputFile.write("#SiteCount\t%s\t\n"%SiteCount)
    def ProduceConsensusSpectra(self):
        """
        We adopt a brute-force strategy: Output a consensus spectrum for each
        modified peptide species.  And, if the equivalent unmodified peptide
        was observed, then output a consensus spectrum for the unmodified
        peptide, too.
        We write a consensus CLUSTER of each modified peptide species.  
        Later on in processing, we may try MERGING two of these clusters
        (e.g. EAM+16APK, EAMA+16PK)
        """
        ClustersBuilt = {} # keep track of clusters that have ALREADY been built
        Keys = self.Peptides.keys()
        for PeptideIndex in range(len(Keys)):
            if (PeptideIndex % 100 == 0):
                print "For peptide %s/%s..."%(PeptideIndex, len(Keys))
            AnnotationKey = Keys[PeptideIndex]
            (Annotation, Charge) = AnnotationKey
            Species = self.Peptides[AnnotationKey]
            if not Species.ModifiedFlag:
                continue
            Species.ClusterPath = os.path.join(self.ConsensusClusterDir, Annotation[2], "%s.%s.cls"%(Annotation.replace("*", "-"), Charge))
            Species.ConsensusPath = os.path.join(self.ConsensusSpectrumDir, Annotation[2], "%s.%s.dta"%(Annotation.replace("*", "-"), Charge))
            ClusterContentPath = os.path.join(self.ClusterScanListDir, Annotation[2], "%s.%s.txt"%(Annotation.replace("*", "-"), Charge))
            ClusterContentFile = open(ClusterContentPath, "wb")
            Builder = BuildConsensusSpectrum.ConsensusBuilder(Species.Charge)
            MeanMQ = 0
            for Info in Species.Spectra:
                MeanMQ += Info.MQScore
            MeanMQ /= float(len(Species.Spectra))
            ValidSpectra = 0
            for Info in Species.Spectra:
                # Omit from the consensus spectra with very poor scores:
                if Info.MQScore < MeanMQ - 3.0:
                    continue
                SpectrumFilePath = self.CachedFixedFilePaths[Info.FileNameIndex]
                # Keep track of where these scans came from:
                ClusterContentFile.write("%s\t%s\t\n"%(SpectrumFilePath, Info.ByteOffset))
                Spectrum = MSSpectrum.SpectrumClass()
                SpectrumFile = open(SpectrumFilePath, "rb")
                SpectrumFile.seek(Info.ByteOffset)
                Spectrum.ReadPeaksFromFile(SpectrumFile, SpectrumFilePath)
                if not Spectrum.PrecursorMZ:
                    print "* Error: Unable to read spectrum from '%s:%s'"%(SpectrumFilePath, Info.ByteOffset)
                    continue
                ValidSpectra += 1
                Spectrum.SetCharge(Charge)
                SpectrumFile.close()
                Builder.AddSpectrum(Spectrum)
                # Special (and easy) case: If we only saw one spectrum, then write it
                # out without changing it!
                if len(Species.Spectra) == 1:
                    Spectrum.WritePeaks(Species.ConsensusPath)
            # Write the modded cluster to disk, since we may try to augment it later:
            Builder.PickleCluster(Species.ClusterPath)
            if len(Species.Spectra) > 1:
                Spectrum = Builder.ProduceConsensusSpectrum()
                Spectrum.WritePeaks(Species.ConsensusPath)
            ClusterContentFile.close()
            # If we have unmodified peptides for this species, build their cluster:
            if Species.Modless:
                Species.Modless.ConsensusPath = os.path.join(self.ConsensusSpectrumDir, Species.Modless.Annotation[2], "%s.%s.dta"%(Species.Modless.Annotation.replace("*", "-"), Charge))
                Species.Modless.ClusterPath = os.path.join(self.ConsensusClusterDir, Species.Modless.Annotation[2], "%s.%s.cls"%(Species.Modless.Annotation.replace("*", "-"), Charge))
                ModlessKey = (Species.Modless.Annotation, Species.Modless.Charge)
                if ClustersBuilt.has_key(ModlessKey):
                    pass
                else:
                    ModlessMeanMQ = 0
                    for Info in Species.Modless.Spectra:
                        ModlessMeanMQ += Info.MQScore
                    ModlessMeanMQ /= float(len(Species.Modless.Spectra))
                    Builder = BuildConsensusSpectrum.ConsensusBuilder(Species.Charge)
                    for Info in Species.Modless.Spectra:
                        # Omit from the consensus spectra with very poor scores:
                        if Info.MQScore < ModlessMeanMQ - 3.0:
                            continue
                        SpectrumFilePath = self.CachedFixedFilePaths[Info.FileNameIndex]
                        Spectrum = MSSpectrum.SpectrumClass()
                        SpectrumFile = open(SpectrumFilePath, "rb")
                        SpectrumFile.seek(Info.ByteOffset)
                        Spectrum.ReadPeaksFromFile(SpectrumFile, SpectrumFilePath)
                        Spectrum.SetCharge(Charge)
                        SpectrumFile.close()
                        Builder.AddSpectrum(Spectrum)
                    Spectrum = Builder.ProduceConsensusSpectrum()
                    Spectrum.WritePeaks(Species.Modless.ConsensusPath)
                    Builder.PickleCluster(Species.Modless.ClusterPath)
                    ClustersBuilt[ModlessKey] = 1
    def ComputeFeaturesAllPTMs(self):
        """
        Compute, and output, features for each modification site.
        """
        self.OutputFile = open(self.OutputPath, "wb")
        self.OutputPTMInfoHeader()
        # Use self.ConsensusCreatedFlags to flag which unmodified peptides we have
        # already generated consensus spectra for, so that we don't do the same one twice and waste time:
        self.ConsensusCreatedFlags = {}
        # Order the peptides by (ModDBPos, Annotation, Charge).  It's important to order things
        # in this way so that, when we combine the *large* output files for HEK293, we can keep
        # consistent 'cursors' moving through each of our input files.
        print "Sorting %s peptides..."%len(self.Peptides.values())
        SortedKeys = []
        for Peptide in self.Peptides.values():
            Key = (Peptide.ModDBPos, Peptide.Annotation, Peptide.Charge)
            SortedKeys.append(Key)
        SortedKeys.sort()
        for KeyIndex in range(len(SortedKeys)):
            (ModDBPos, Annotation, Charge) = SortedKeys[KeyIndex]
            Key = (Annotation, Charge)
            Species = self.Peptides[Key]
            if not Species.ModifiedFlag:
                continue
            if Species.DBPos < self.StartOutputDBPos:
                continue
            print "(%s/%s) PTM: %+d on db residue %d"%(KeyIndex, len(SortedKeys), Species.ModMass, Species.DBPos + Species.Peptide.Modifications.keys()[0])
            sys.stdout.flush()
            try:
                Features = self.ComputePTMFeatures(Species)
            except:
                traceback.print_exc()
                print "** Error: Unable to compute PTM features for %s"%Species
                continue
            Str = "%s\t"%self.OutputPath
            Str += "%s\t%+d\t%s\t"%(Species.ModDBPos, Species.ModMass, Species.ModAA) 
            Str += "%s\t"%Species.ProteinName
            Str += "%s\t"%(Species.ResidueNumber + Species.Peptide.Modifications.keys()[0])
            Str += "%s\t"%Species.Annotation
            Str += "%s\t"%Species.Charge
            for Feature in Features:
                Str += "%s\t"%Feature
            print Str
            self.OutputFile.write(Str + "\n")
            # We're done with this PTM now, so let's forget about it:
            del self.Peptides[Key]
        self.OutputFile.close()
    def ComputePTMFeatures(self, Species):
        """
        Compute scoring-features for this peptide species, return them as a list
        """
        Features = []
        # Feature: Is the PTM from a valid protein?  (Note: this feature is not INPUT for the
        # model, it's our desired output)
        Feature = self.GetValidProteinFlag(Species)
        Features.append(Feature)
        # Important question: Is this PTM constituitive, or facultative?  (In other words:
        # is there a spectra annotated with an UNMODIFIED peptide for this PTM type?)
        # Set flag to "1" if the PTM is facultative:
        if Species.Modless:
            Features.append("1")
        else:
            Features.append("")
        # The best spectrum observed (meta-data, not a scoring feature)
        BestMQScore = -999
        BestDeltaScore = None
        for Info in Species.Spectra:
            if Info.MQScore > BestMQScore:
                BestMQScore = Info.MQScore
                BestDeltaScore = Info.DeltaScore
                FilePath = self.CachedFixedFilePaths[Info.FileNameIndex]
                BestMQSpectrum = ("%s:%s"%(FilePath, Info.ByteOffset))
        Features.append(BestMQSpectrum)
        # The best MODLESS spectrum observed (meta-data, not a scoring feature)
        if Species.Modless:
            BestMQScoreModless = -999
            for Info in Species.Modless.Spectra:
                if Info.MQScore > BestMQScoreModless:
                    BestMQScoreModless = Info.MQScore
                    FilePath = self.CachedFixedFilePaths[Info.FileNameIndex]
                    BestMQSpectrum = ("%s:%s"%(FilePath, Info.ByteOffset))
            Features.append(BestMQSpectrum)
            Features.append(str(BestMQScoreModless))
        else:
            Features.append("")
            Features.append("")
        # Feature: Annotation, and MQScore, from a search versus big-DB.  (This feature
        # will be spiked in later)
        Features.append("")
        Features.append("")
        # Feature: Number of spectra annotated with this PTM
        Features.append(Species.SpectrumCount)
        # Feature: Number of spectra for the *unmodified* peptide version:
        if Species.Modless:
            Feature = Species.Modless.SpectrumCount
            Features.append(Feature)
        else:
            Features.append(0)
        # Feature: Best MQScore for this PTM, and the best delta-score for that scan:
        Features.append(BestMQScore)
        Features.append(BestDeltaScore)
        # Feature: Number of peptide species observed for this PTM on this residue
        Features.append(len(Species.PTM.SpeciesList))
        # Feature: Consensus annotation score (and score-features) for this peptide
        Species.ConsensusScore = None
        self.GetConsensusMQScore(Species, Features)
        # Feature: Presence of unmodified peptides covering the residue of interest
        ModlessCount = self.Coverage[Species.PTM.DBPos]
        ModdedSpectrumCount = self.ModCoverage[Species.PTM.DBPos]
        ModdedFraction = ModdedSpectrumCount / float(ModlessCount + ModdedSpectrumCount)
        Features.append(ModdedFraction)
        # Feature: Coverage of the protein of interest (ONLY FOR FACULTATIVE!)
        ProteinIndex = self.GetProteinIndex(Species.DBPos)
        ProteinCoverage = self.ProteinCoverageLevels[ProteinIndex]
        Features.append(ProteinCoverage)
        # Feature: Number of annotations using this modification-type
        ModTypeKey = (Species.ModAA, Species.ModMass)
        ModTypeSpectrumCount = self.ModTypeSpectrumCount.get(ModTypeKey, 0)
        Features.append(ModTypeSpectrumCount)
        # Feature: Number of sites using this modification-type
        ModTypeSiteCount = self.ModTypeSiteCount.get(ModTypeKey, 0)
        Features.append(ModTypeSiteCount)
        # Feature: Number of unmodified peptide species for this site
        ModlessPeptides = self.PeptideCoverage[Species.DBPos]
        Features.append(ModlessPeptides)
        # Features for FACULTATIVE PTMs only:
        # These features have been commented out, since we no longer pursue a
        # special model for facultative PTMs.
        if 0: #Species.Modless:
            Comparator = SpectralSimilarity.SpectralSimilarity(Species.ConsensusPath,
               Species.Modless.ConsensusPath, Species.Annotation, Species.Modless.Annotation)
            Comparator.LabelPeaks(0.5)
            Similarity = Comparator.DotProduct(0.5, HashByRank = 1)
            Features.append(Similarity)
            Similarity = Comparator.GetSharedPeakCount(0, 1)
            Features.append(Similarity)
            Similarity = Comparator.GetSharedPeakCount(1, 1)
            Features.append(Similarity)
            CorrelationCoefficient = Comparator.ComputeCorrelationCoefficient(1.0)
            Features.append(CorrelationCoefficient)
            del Comparator
        else:
            # This PTM is constitutive, so omit the spectrum-comparison features:
            Features.append("") # dot
            Features.append("") # shared-peaks
            Features.append("") # shared-peaks
            Features.append("") # correlation
        # Feature: Log of spectrum-count
        Features.append(math.log(1.0 + Species.SpectrumCount))
        # Feature: Log of peptide-length
        Features.append(math.log(len(Species.Peptide.Aminos)))
        # Feature: Log of same-modtype-spectrum-count
        Features.append(math.log(1.0 + ModTypeSpectrumCount))
        # Feature: Log of same-modyupe-site-count
        Features.append(math.log(1.0 + ModTypeSiteCount))
        # Feature: Delta-score versus big-db search result.  To be spiked in later!
        Features.append("")
        # Free the PySpectrum object now:
        Species.PySpectrum = None
        return Features
    def GetValidProteinFlag(self, PTM):
        # Normally we prepend "xxx" to the bogus names:
        if PTM.ProteinName[:3] == "XXX":
            return 0
        return 1
    def GetProteinIndex(self, DBPos):
        for ProteinIndex in range(len(self.ProteinPositions)):
            Pos = self.ProteinPositions[ProteinIndex]
            if Pos > DBPos:
                return ProteinIndex - 1
        return len(self.ProteinPositions) - 1
    def PreComputeAminosForMasses(self):
        """
        PepNovo often gives us partial interpretations - e.g. a
        peptide that starts at 250Da.  We "fill in" the prefix and
        suffix to generate a (not necessarily optimal) full-length
        peptide.
        """
        Aminos = "ACDEFGHILMNPQSTVWYRK" # PREFER ending in R or K.
        self.AAStrings = {}
        TotalMass = 0
        for AA1 in Aminos:
            Mass1 = Global.AminoMass[AA1]
            TotalMass = int(round(Mass1))
            self.AAStrings[TotalMass] = "%c"%(AA1)
            for AA2 in Aminos:
                Mass2 = Global.AminoMass[AA2]
                TotalMass = int(round(Mass1 + Mass2))
                self.AAStrings[TotalMass] = "%c%c"%(AA1, AA2)
                for AA3 in Aminos:
                    Mass3 = Global.AminoMass[AA3]
                    TotalMass = int(round(Mass1 + Mass2 + Mass3))
                    self.AAStrings[TotalMass] = "%c%c%c"%(AA1, AA2, AA3)
                    for AA4 in Aminos:
                        Mass4 = Global.AminoMass[AA4]
                        TotalMass = int(round(Mass1 + Mass2 + Mass3 + Mass4))
                        self.AAStrings[TotalMass] = "%c%c%c%c"%(AA1, AA2, AA3, AA4)
    def AddSpectrumToCluster(self, InputFilePath, InputFilePos, ClusterFile, Charge):
        """
        Append the specified scan to an ever-growing .mgf file
        Returns 1 if successful, 0 if failed
        """
        try:
            SpectrumFile = open(InputFilePath, "rb")
        except:
            print "** Error: couldn't open spectrum data file %s"%InputFilePath
            return 0
        SpectrumFile.seek(InputFilePos)
        Spectrum = MSSpectrum.SpectrumClass()
        try:
            Spectrum.ReadPeaksFromFile(SpectrumFile, InputFilePath)
        except:
            traceback.print_exc()
            print "***Can't parse:", InputFilePath, FileOffset
            return 0
        SpectrumFile.close()
        ParentMass = Spectrum.PrecursorMZ * Charge - (Charge - 1)*1.0078 #Peptide.Masses[-1] + 19
        #MZ = (ParentMass + (Info.Charge - 1)*1.0078) / Info.Charge
        # Now write out this spectrum to the cluster:
        self.ClusterScanNumber += 1 # ASSUMED: The caller set this to 0 at the start of the cluster!
        ClusterFile.write("BEGIN IONS\n")
        ClusterFile.write("TITLE=%s:%s\n"%(InputFilePath, InputFilePos))
        ClusterFile.write("SCAN=%s\n"%self.ClusterScanNumber)
        ClusterFile.write("CHARGE=%s\n"%Charge)
        ClusterFile.write("PEPMASS=%s\n"%ParentMass)
        for Peak in Spectrum.Peaks:
            ClusterFile.write("%s %s\n"%(Peak.Mass, Peak.Intensity))
        ClusterFile.write("END IONS\n")
        #ClusterFile.close()
        return 1
    def GetConsensusMQScore(self, Species, Features):
        """
        Feature: MQScore of the consensus spectrum.
        - Write specra to a cluster (done by ProduceConsensusSpectra)
        - Generate a consensus-spectrum for the cluster (done by ProduceConsensusSpectra)
        - Load the consensus-spectrum
        - Score the spectrum
        """
        # Load in the consensus spectrum, and score the peptide annotation:
        try:
            print ">>PyConsensus spectrum:", Species.ConsensusPath
            PySpectrum = PyInspect.Spectrum(Species.ConsensusPath, 0)
            Species.PySpectrum = PySpectrum
            print ">>ScorePeptideDetailed(%s)"%Species.Annotation
            ScoreList = PySpectrum.ScorePeptideDetailed(Species.Annotation)
            Species.ConsensusScore = ScoreList[0]
            for ScoreItem in ScoreList:
                Features.append(ScoreItem)
            print "PyInspect score %s -> %s"%(Species.Annotation, ScoreList[0])
        except:
            traceback.print_exc()
            for X in range(8):
                Features.append(0)
    def ParseCommandLine(self, Arguments):
        (Options, Args) = getopt.getopt(Arguments, "d:r:w:s:M:lp:c:Z:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-r":
                # -r results file(s)
                self.ResultsFileName = Value
            elif Option == "-c":
                self.RequiredFileNameChunk = Value
            elif Option == "-d":
                self.DBPath = Value
            elif Option == "-M":
                self.PopulateSpectrumOracle(Value)
            elif Option == "-w":
                self.OutputDir = Value
                self.OutputPath = os.path.join(self.OutputDir, "PTMFeatures.txt")
            elif Option == "-l":
                self.QuickParseFlag = 1
            elif Option == "-s" or Option == "-M":
                self.PopulateSpectrumOracle(Value)
                #self.SpectrumDir = Value
            elif Option == "-p":
                self.PValueCutoff = float(Value)
            elif Option == "-Z": # secret debugging option: Start output from DB position
                self.StartOutputDBPos = int(Value)
            else:
                print "* Error: Unrecognized option %s"%Option
            
UsageInfo = """
ComputePTMFeatures: Generate feature values for PTMs observed on a data-set.
Run this AFTER running SelectSites, and BEFORE running TrainPTMFeatures.

Arguments:
 -r [ResultsFile]: Name of the results file (or directory)
 -d [DBPath]: Path to the .trie file searched
 -w [OutputDir]: Output file directory.  Features are written to
    PTMFeatures.txt within this directory.  Clusters and other info
    is written in (or below) this directory.
 -M [RootDir]: Root directory for mzXML files. 
"""
           
if __name__ == "__main__":
    if not PROFILE_FLAG:
        try:
            import psyco
            psyco.full()
        except:
            print "(psyco not installed; running unoptimized)"
    Trainer = PTMFeatureComputer()
    Trainer.ParseCommandLine(sys.argv[1:])
    if not Trainer.ResultsFileName or not Trainer.DBPath:
        print UsageInfo
        sys.exit(-1)
    if PROFILE_FLAG:
        import profile
        profile.run("Trainer.ComputeFeaturesMain()")
    else:
        Trainer.ComputeFeaturesMain()
    
    