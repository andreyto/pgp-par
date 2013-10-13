"""
CombinePTMFeatures:
- Parse the output of several runs of ComputePTMFeatures.
- Accumulate consensus clusters, and accumulate coverage
- Write the results to a merged directory
"""
import sys
import os
import math
import getopt
import string
import struct
import traceback
import shutil
import BuildConsensusSpectrum
import ResultsParser
import PyInspect
from Utils import *
from TrainPTMFeatures import FormatBits

UsageInfo = """
-r [PATH]: Directory containing the results-directories to merge
-w [PATH]: Output file
-d [PATH]: Database path

Optional:
-M [DIR]: Directory subtree where mzxml files *really* live.
    Spectrum paths will be corrected to use these paths.
-x: If set, prepare output directories.  (Should be set for the first
    run of a batch, and NOT for any others)
-s [POS]: Start DBPosition
-e [POS]: End DBPosition
-q: Quick-parse flag
-c [STRING]: Required filename chunk
"""

class PTMFeatureMerger(ResultsParser.SpectrumOracleMixin):
    def __init__(self):
        self.OutputDir = None
        self.Peptides = {} # (annotation, charge) -> peptide species
        self.HeaderLines = []
        self.HeaderLinesParsed = 0
        self.DBStart = None
        self.DBEnd = None
        self.QuickParseFlag = 0
        self.TotalSpectrumCount = 0
        self.SpectrumRoot = None
        ResultsParser.SpectrumOracleMixin.__init__(self)
    def WipeDir(self, Dir):
        try:
            shutil.rmtree(Dir)
        except:
            pass
    def ParseCommandLine(self, Arguments):
        PrepareDirsFlag = 0
        (Options, Args) = getopt.getopt(Arguments, "d:r:w:s:e:qxM:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-r":
                self.PTMFeatureDirectory = Value
            elif Option == "-d":
                self.DBPath = Value
            elif Option == "-w":
                self.OutputPath = Value
            elif Option == "-s":
                self.DBStart = int(Value)
            elif Option == "-e":
                self.DBEnd = int(Value)
            elif Option == "-q":
                self.QuickParseFlag = 1
            elif Option == "-x":
                PrepareDirsFlag = 1
            elif Option == "-M":
                self.SpectrumRoot = Value
            else:
                print "* Error: Unrecognized option %s"%Option
        if not self.OutputPath:
            print "* Please specify an output file (-w)"
            print UsageInfo
            sys.exit(-1)
        self.OutputDir = os.path.split(self.OutputPath)[0]
        #self.OutputPath = os.path.join(self.OutputDir, "PTMFeatures.txt")
        self.ClusterDir = os.path.join(self.OutputDir, "Clusters")
        self.SpectrumDir = os.path.join(self.OutputDir, "Spectra")
        self.ClusterMemberDir = os.path.join(self.OutputDir, "ClusterMembers")
        print "Prepare directories..."
        if PrepareDirsFlag:
            self.WipeDir(self.OutputDir)
            MakeDirectory(self.OutputDir)
            for Dir in (self.ClusterDir, self.SpectrumDir, self.ClusterMemberDir):
                MakeDirectory(Dir)
                for AA in "ACDEFGHIKLMNPQRSTVWY":
                    MakeDirectory(os.path.join(Dir, AA))
        return 1 # success
    def LoadDB(self):
        # Populate self.DB with the contents of the .trie file
        File = open(self.DBPath, "rb")
        self.DB = File.read()
        File.close()
        self.Coverage = [0] * len(self.DB)
        self.ModCoverage = [0] * len(self.DB)
    def OutputCoverage(self):
        CoveragePath = os.path.join(self.OutputDir, "Coverage.dat")
        CoverageFile = open(CoveragePath, "wb")
        for DBPos in range(len(self.DB)):
            Str = struct.pack("<II", self.Coverage[DBPos], self.ModCoverage[DBPos])
            CoverageFile.write(Str)
        CoverageFile.close()
    def OutputPTMFeatures(self):
        File = open(self.OutputPath, "wb")
        for FileLine in self.HeaderLines:
            File.write(FileLine)
        for Peptide in self.Peptides.values():
            String = string.join(Peptide.Bits, "\t")
            File.write(String + "\n")
    def GrabClusterMembers(self, Directory, Peptide):
        InputPath = os.path.join(Directory, "ClusterMembers", Peptide.Annotation[2], "%s.%s.txt"%(Peptide.Annotation, Peptide.Charge))
        OutputPath = os.path.join(self.ClusterMemberDir, Peptide.Annotation[2], "%s.%s.txt"%(Peptide.Annotation, Peptide.Charge))
        InputFile = open(InputPath, "rb")
        OutputFile = open(OutputPath, "a+b")
        OutputFile.write(InputFile.read())
        InputFile.close()
        OutputFile.close()
    def BuildNewPeptide(self, Cursor):
        """
        We've parsed the first of the cursors that contains this peptide species.  Build a peptide
        and start accumulating clusters.
        """
        Bits = Cursor.Bits
        
    def AddNewPeptide(self, Directory, Bits):
        """
        We're parsing a peptide (from the ComputePTMFeatures output) which we haven't
        seen before.  Create a Peptide object and populate it.
        """
        Peptide = Bag()
        Peptide.Bits = list(Bits)
        Peptide.Charge = int(Bits[FormatBits.Charge])
        Peptide.Annotation = Bits[FormatBits.Peptide]
        SpectrumCount = int(Bits[FormatBits.SpectrumCount])
        self.TotalSpectrumCount += SpectrumCount
        Peptide.Peptide = GetPeptideFromModdedName(Peptide.Annotation)
        Peptide.ModlessAnnotation = "%s.%s.%s"%(Peptide.Peptide.Prefix, Peptide.Peptide.Aminos, Peptide.Peptide.Suffix)
        # Grab the cluster members right away:
        self.GrabClusterMembers(Directory, Peptide)
        # Add a cluster for the peptide:
        Peptide.Cluster = BuildConsensusSpectrum.ConsensusBuilder(Peptide.Charge)
        ClusterPath = os.path.join(Directory, "Clusters", Peptide.Annotation[2], "%s.%s.cls"%(Peptide.Annotation, Peptide.Charge))
        Peptide.Cluster.UnpickleCluster(ClusterPath)
        # Add a cluster for the modless peptide, if the cluster-file exists:
        ModlessClusterPath = os.path.join(Directory, "Clusters", Peptide.Annotation[2], "%s.%s.cls"%(Peptide.ModlessAnnotation, Peptide.Charge))
        if os.path.exists(ModlessClusterPath):
            Peptide.ModlessCluster = BuildConsensusSpectrum.ConsensusBuilder(Peptide.Charge)
            Peptide.ModlessCluster.UnpickleCluster(ModlessClusterPath)
        else:
            Peptide.ModlessCluster = None
        # Add the peptide to our dictionary:
        Key = (Peptide.Annotation, Peptide.Charge)
        self.Peptides[Key] = Peptide
    def AssimilatePeptide(self, Dir, Peptide, Bits):
        """
        We're parsing a peptide (from the ComputePTMFeatures output) which we 
        have already seen.  Adjust the features of our Peptide object - accumulate
        spectrum counts, et cetera.
        """
        # Best modless spectrum and MQScore.  These may be empty for the new
        # file bits, for the existing peptide, or both.
        ScoreStr = Peptide.Bits[FormatBits.BestModlessMQScore]
        if ScoreStr:
            OldModlessMQScore = float(ScoreStr)
        else:
            OldModlessMQScore = None
        ScoreStr = Bits[FormatBits.BestModlessMQScore]
        if ScoreStr:
            ModlessMQScore = float(ScoreStr)
        else:
            ModlessMQScore = None
        if ModlessMQScore > OldModlessMQScore:
            Peptide.Bits[FormatBits.BestModlessMQScore] = Bits[FormatBits.BestModlessMQScore]
            Peptide.Bits[FormatBits.BestModlessSpectrumPath] = Bits[FormatBits.BestModlessSpectrumPath]
        # Best modded spectrum, mqscore, delta-score:
        OldBestMQScore = float(Peptide.Bits[FormatBits.BestMQScore])
        BestMQScore = float(Bits[FormatBits.BestMQScore])
        if BestMQScore > OldBestMQScore:
            Peptide.Bits[FormatBits.BestMQScore] = Bits[FormatBits.BestMQScore]
            Peptide.Bits[FormatBits.BestDeltaScore] = Bits[FormatBits.BestDeltaScore]
            Peptide.Bits[FormatBits.BestSpectrumPath] = Bits[FormatBits.BestSpectrumPath]
        # Spectra:
        CurrentSpectra = int(Peptide.Bits[FormatBits.SpectrumCount])
        NewBlockSpectra = int(Bits[FormatBits.SpectrumCount])
        TotalSpectra = CurrentSpectra + NewBlockSpectra
        self.TotalSpectrumCount += NewBlockSpectra
        Peptide.Bits[FormatBits.SpectrumCount] = str(Spectra)
        # Modless spectra
        Spectra = int(Peptide.Bits[FormatBits.ModlessSpectrumCount])
        Spectra += int(Bits[FormatBits.ModlessSpectrumCount])
        Peptide.Bits[FormatBits.ModlessSpectrumCount] = str(Spectra)
        # Accumulate modded spectra into the cluster:
        ClusterPath = os.path.join(Dir, "Clusters", Peptide.Annotation[2], "%s.%s.cls"%(Peptide.Annotation, Peptide.Charge))
        TempCluster = BuildConsensusSpectrum.ConsensusBuilder(Peptide.Charge)
        TempCluster.UnpickleCluster(ClusterPath)
        Peptide.Cluster.AssimilateCluster(TempCluster)
        # Accumulate modless spectra into the modless cluster:
        ClusterPath = os.path.join(Dir, "Clusters", Peptide.Annotation[2], "%s.%s.cls"%(Peptide.ModlessAnnotation, Peptide.Charge))
        if os.path.exists(ClusterPath):
            TempCluster = BuildConsensusSpectrum.ConsensusBuilder(Peptide.Charge)
            TempCluster.UnpickleCluster(ClusterPath)
            if Peptide.ModlessCluster:
                Peptide.ModlessCluster.AssimilateCluster(TempCluster)
            else:
                Peptide.ModlessCluster = TempCluster
        # Consensus MQScore handled at the end
        # Spectra/sites this mod type handled at the *VERY* end, possibly after
        #   multiple runs of CombinePTMFeatures!
        # Log spectrum-count handled at the end
    def MergeResultsFromFile(self, Path):
        Dir = os.path.split(Path)[0]
        File = open(Path, "rb")
        LineNumber = 0
        for FileLine in File.xreadlines():
            LineNumber += 1
            if LineNumber % 100 == 0:
                print "%s line %s (%s peptides, %s spectra)"%(Path, LineNumber, len(self.Peptides.keys()), self.TotalSpectrumCount)
                if self.QuickParseFlag:
                    break
            if FileLine[0] == "#":
                if not self.HeaderLinesParsed:
                    self.HeaderLines.append(FileLine)
                continue
            Bits = FileLine.strip().split("\t")
            # Skip any blank lines:
            if len(Bits) < 2:
                continue
            try:
                DBPos = int(Bits[FormatBits.DBPos])
                Annotation = Bits[FormatBits.Peptide]
                Charge = int(Bits[FormatBits.Charge])
            except:
                print "* Warning: Line %s of %s isn't valid!"%(LineNumber, Path)
                traceback.print_exc()
                continue
            # Ignore any peptides which don't fall within our database region of interest:
            if self.DBStart != None and DBPos < self.DBStart:
                continue
            if self.DBEnd != None and DBPos >= self.DBEnd:
                continue
            Key = (Annotation, Charge)
            Peptide = self.Peptides.get(Key, None)
            if Peptide:
                self.AssimilatePeptide(Dir, Peptide, Bits)
            else:
                self.AddNewPeptide(Dir, Bits)
        File.close()
        self.HeaderLinesParsed = 1
    def FinalizePTMFeatures(self):
        """
        Some PTM feature processing, such as building and scoring a consensus spectrum,
        should happen just once.  Those steps happen here, *after* each input file
        has been parsed.
        """
        for Peptide in self.Peptides.values():
            # Update log-spectrum-count:
            Spectra = int(Peptide.Bits[FormatBits.SpectrumCount])
            Peptide.Bits[FormatBits.LogSpectrumCount] = str(math.log(Spectra))
            # Write out the consensus MODLESS cluster:
            if Peptide.ModlessCluster:
                Path = os.path.join(self.ClusterDir, Peptide.Annotation[2], "%s.%s.cls"%(Peptide.ModlessAnnotation, Peptide.Charge))
                Peptide.ModlessCluster.PickleCluster(Path)
            # Write out the consensus MODLESS spectrum:
            if Peptide.ModlessCluster:
                Path = os.path.join(self.SpectrumDir, Peptide.Annotation[2], "%s.%s.dta"%(Peptide.ModlessAnnotation, Peptide.Charge))
                Spectrum = Peptide.ModlessCluster.ProduceConsensusSpectrum()
                Spectrum.WritePeaks(Path)
            # Write out the CLUSTER:
            Path = os.path.join(self.ClusterDir, Peptide.Annotation[2], "%s.%s.cls"%(Peptide.Annotation, Peptide.Charge))
            Peptide.Cluster.PickleCluster(Path)
            # Write out the consensus SPECTRUM:
            ConsensusSpectrumPath = os.path.join(self.SpectrumDir, Peptide.Annotation[2], "%s.%s.dta"%(Peptide.Annotation, Peptide.Charge))
            Spectrum = Peptide.Cluster.ProduceConsensusSpectrum()
            Spectrum.WritePeaks(ConsensusSpectrumPath)
            # Compute consensus spectrum features:
            PySpectrum = PyInspect.Spectrum(ConsensusSpectrumPath, 0)
            ScoreList = PySpectrum.ScorePeptideDetailed(Peptide.Annotation)
            Peptide.Bits[FormatBits.ConsensusMQScore] = str(ScoreList[0])
            Peptide.Bits[FormatBits.ConsensusPRMScore] = str(ScoreList[1])
            Peptide.Bits[FormatBits.ConsensusBYPresent] = str(ScoreList[2])
            Peptide.Bits[FormatBits.ConsensusTopPeaks] = str(ScoreList[3])
            Peptide.Bits[FormatBits.NTT] = str(ScoreList[4])
            # Compute comparison features:
            if Peptide.ModlessCluster:
                Peptide.Bits[FormatBits.SisterAnnotationFlag] = "1"
                pass #NOTE: skip these features, since we don't really use them!
            
    def AssimilateDatabaseCoverage(self, CoverageFilePath):
        """
        Read Coverage.dat from one of the ComputePTMFeatures runs.
        Accumulate total coverage.
        """
        StructSize = struct.calcsize("<II")
        File = open(CoverageFilePath, "rb")
        for DBPos in range(len(self.DB)):
            Block = File.read(StructSize)
            (Coverage, ModCoverage) = struct.unpack("<II", Block)
            self.Coverage[DBPos] += Coverage
            self.ModCoverage[DBPos] += ModCoverage
            if DBPos % 10000 == 0:
                print "%s/%s..."%(DBPos, len(self.DB))
        File.close()
    def WriteSingletonPeptide(self, Species, Cursor):
        """
        Here's a peptide that was found in only one of the input files.  That makes our job very
        easy; all we need to do is write out the fileline, and copy over: cluster members,
        cluster file, spectrum file, modless cluster file (if it exists), modless spectrum file (if
        it exists)
        """
        # Write the file line:
        OutputLine = string.join(Cursor.Bits, "\t")
        self.OutputFile.write(OutputLine + "\n")
        # Copy the needful files:
        Annotation = Species.Annotation
        Charge = Species.Charge
        AA = Species.AA
        SourcePath = os.path.join(Cursor.Directory, "ClusterMembers", AA, "%s.%s.txt"%(Annotation, Charge))
        TargetPath = os.path.join(self.OutputDir, "ClusterMembers", AA, "%s.%s.txt"%(Annotation, Charge))
        shutil.copyfile(SourcePath, TargetPath)
        SourcePath = os.path.join(Cursor.Directory, "Clusters", AA, "%s.%s.cls"%(Annotation, Charge))
        TargetPath = os.path.join(self.OutputDir, "Clusters", AA, "%s.%s.cls"%(Annotation, Charge))
        shutil.copyfile(SourcePath, TargetPath)
        SourcePath = os.path.join(Cursor.Directory, "Spectra", AA, "%s.%s.dta"%(Annotation, Charge))
        TargetPath = os.path.join(self.OutputDir, "Spectra", AA, "%s.%s.dta"%(Annotation, Charge))
        shutil.copyfile(SourcePath, TargetPath)
        SourcePath = os.path.join(Cursor.Directory, "Clusters", AA, "%s.%s.cls"%(Species.ModlessAnnotation, Charge))
        TargetPath = os.path.join(self.OutputDir, "Clusters", AA, "%s.%s.cls"%(Species.ModlessAnnotation, Charge))
        if os.path.exists(SourcePath):
            shutil.copyfile(SourcePath, TargetPath)
        SourcePath = os.path.join(Cursor.Directory, "Spectra", AA, "%s.%s.dta"%(Species.ModlessAnnotation, Charge))
        TargetPath = os.path.join(self.OutputDir, "Spectra", AA, "%s.%s.dta"%(Species.ModlessAnnotation, Charge))
        if os.path.exists(SourcePath):
            shutil.copyfile(SourcePath, TargetPath)
    def BuildPeptideSpecies(self, Cursor):
        Species = Bag()
        Bits = Cursor.Bits
        Species.Annotation = Bits[FormatBits.Peptide]
        Species.AA = Species.Annotation[2]
        Species.Peptide = GetPeptideFromModdedName(Species.Annotation)
        Species.ModlessAnnotation = "%s.%s.%s"%(Species.Peptide.Prefix, Species.Peptide.Aminos, Species.Peptide.Suffix)
        Species.Charge = int(Bits[FormatBits.Charge])
        Species.Bits = Bits[:]
        Species.Cluster = None # ConsensusBuilder, instantaited later
        Species.ModlessCluster = None # ConsensusBuilder, instantaited later
        # Fix the spectrum PATHS:
        Species.Bits[FormatBits.BestSpectrumPath] = self.FixSpectrumPath(Species.Bits[FormatBits.BestSpectrumPath])
        if Species.Bits[FormatBits.BestModlessSpectrumPath]:
            Species.Bits[FormatBits.BestModlessSpectrumPath] = self.FixSpectrumPath(Species.Bits[FormatBits.BestModlessSpectrumPath])
        return Species
    def AssimilatePeptideSpectra(self, Species, Cursor):
        """
        Accumulate total spectra for this species.
        """
        ###############################
        # Adjust features - best spectrum, number of spectra, etc.
        # Best modless spectrum and MQScore.  These may be empty for the new
        # file bits, for the existing peptide, or both.
        ScoreStr = Species.Bits[FormatBits.BestModlessMQScore]
        if ScoreStr:
            OldModlessMQScore = float(ScoreStr)
        else:
            OldModlessMQScore = None
        ScoreStr = Cursor.Bits[FormatBits.BestModlessMQScore]
        if ScoreStr:
            ModlessMQScore = float(ScoreStr)
        else:
            ModlessMQScore = None
        if ModlessMQScore > OldModlessMQScore:
            Species.Bits[FormatBits.BestModlessMQScore] = Cursor.Bits[FormatBits.BestModlessMQScore]
            Species.Bits[FormatBits.BestModlessSpectrumPath] = Cursor.Bits[FormatBits.BestModlessSpectrumPath]
        # Best modded spectrum, mqscore, delta-score:
        OldBestMQScore = float(Species.Bits[FormatBits.BestMQScore])
        BestMQScore = float(Cursor.Bits[FormatBits.BestMQScore])
        if BestMQScore > OldBestMQScore:
            Species.Bits[FormatBits.BestMQScore] = Cursor.Bits[FormatBits.BestMQScore]
            Species.Bits[FormatBits.BestDeltaScore] = Cursor.Bits[FormatBits.BestDeltaScore]
            Species.Bits[FormatBits.BestSpectrumPath] = Cursor.Bits[FormatBits.BestSpectrumPath]
        # Spectra:
        CurrentSpectra = int(Species.Bits[FormatBits.SpectrumCount])
        NewBlockSpectra = int(Cursor.Bits[FormatBits.SpectrumCount])
        TotalSpectra = CurrentSpectra + NewBlockSpectra
        self.TotalSpectrumCount += NewBlockSpectra
        Species.Bits[FormatBits.SpectrumCount] = str(TotalSpectra)
        # Log of spectrum-count:
        Species.Bits[FormatBits.LogSpectrumCount] = str(math.log(TotalSpectra))
        # Modless spectra
        Spectra = int(Species.Bits[FormatBits.ModlessSpectrumCount])
        Spectra += int(Cursor.Bits[FormatBits.ModlessSpectrumCount])
        Species.Bits[FormatBits.ModlessSpectrumCount] = str(Spectra)
        ###############################
        # Accumulate a list of cluster members:
        ClusterMemberFileName = "%s.%s.txt"%(Species.Annotation, Species.Charge)
        ClusterMemberPath = os.path.join(self.OutputDir, "ClusterMembers", Species.AA, ClusterMemberFileName)
        ClusterMemberFile = open(ClusterMemberPath, "a+b")
        CursorMemberPath = os.path.join(Cursor.Directory, "ClusterMembers", Species.AA, ClusterMemberFileName)
        CursorMemberFile = open(CursorMemberPath, "rb") 
        Text = CursorMemberFile.read()
        ClusterMemberFile.write(Text)
        CursorMemberFile.close()
        ClusterMemberFile.close()
        ###############################
        # Accumulate members of the modded cluster:
        ClusterPath = os.path.join(Cursor.Directory, "Clusters", Species.AA, "%s.%s.cls"%(Species.Annotation, Species.Charge))
        CursorCluster = BuildConsensusSpectrum.ConsensusBuilder(Species.Charge)
        CursorCluster.UnpickleCluster(ClusterPath)
        if not Species.Cluster:
            Species.Cluster = CursorCluster 
        else:
            Species.Cluster.AssimilateCluster(CursorCluster)
        ###############################
        # Accumulate members of the modless cluster:
        ClusterPath = os.path.join(Cursor.Directory, "Clusters", Species.AA, "%s.%s.cls"%(Species.ModlessAnnotation, Species.Charge))
        if os.path.exists(ClusterPath):
            CursorCluster = BuildConsensusSpectrum.ConsensusBuilder(Species.Charge)
            CursorCluster.UnpickleCluster(ClusterPath)
            if not Species.ModlessCluster:
                Species.ModlessCluster = CursorCluster 
            else:
                Species.ModlessCluster.AssimilateCluster(CursorCluster)
    def WriteCompletedPeptide(self, Species):
        """
        We've read data for this peptide species from TWO OR MORE cursors.  Now we'll write out
        one line to the output file, and output our consensus spectrum.
        """
        # Write the file line:
        OutputLine = string.join(Species.Bits, "\t")
        self.OutputFile.write(OutputLine + "\n")
        # Cluster:
        ClusterOutputPath = os.path.join(self.OutputDir, "Clusters", Species.AA, "%s.%s.cls"%(Species.Annotation, Species.Charge))
        Species.Cluster.PickleCluster(ClusterOutputPath)
        # Consensus spectrum:
        ConsensusSpectrum = Species.Cluster.ProduceConsensusSpectrum()
        ConsensusSpectrumPath = os.path.join(self.OutputDir, "Spectra", Species.AA, "%s.%s.dta"%(Species.Annotation, Species.Charge))
        ConsensusSpectrum.WritePeaks(ConsensusSpectrumPath)
        if Species.ModlessCluster:
            # Modless cluster:
            ClusterOutputPath = os.path.join(self.OutputDir, "Clusters", Species.AA, "%s.%s.cls"%(Species.ModlessAnnotation, Species.Charge))
            Species.ModlessCluster.PickleCluster(ClusterOutputPath)
            # Modless consensus spectrum:
            ConsensusSpectrum = Species.ModlessCluster.ProduceConsensusSpectrum()
            ConsensusSpectrumPath = os.path.join(self.OutputDir, "Spectra", Species.AA, "%s.%s.dta"%(Species.ModlessAnnotation, Species.Charge))
            ConsensusSpectrum.WritePeaks(ConsensusSpectrumPath)
    def MergeResults(self):
        print "Load db..."
        self.LoadDB()
        # Measure combined db coverage:
        print "Combine database coverage..."
        self.CombineDatabaseCoverage()
        print "Populate MZXML oracle..."
        self.PopulateMZXMLOracle(self.SpectrumRoot)
        self.OutputFile = open(self.OutputPath, "wb")
        for FileLine in self.HeaderLines:
            self.OutputFile.write(FileLine)        
        class FeatureCursor:
            """
            Wrapper for a feature file - tracks the "next" peptide.  
            """
            def __init__(self, Path):
                self.Path = Path
                self.File = open(Path, "rb")
                self.NextKey = None
                self.Bits = None
                self.Directory = os.path.split(Path)[0]
                self.HeaderLines = []
            def Close(self):
                self.File.close()
                self.NextKey = None
            def GetNextPeptide(self):
                # Read one or more lines (skipping header or invalid lines), and remember
                # the next peptide to be processed
                while (1):
                    FileLine = self.File.readline()
                    if not FileLine:
                        self.NextKey = None
                        return None # EOF
                    # Skip blank or comment lines:
                    if FileLine[0] == "#":
                        self.HeaderLines.append(FileLine)
                        continue
                    if not FileLine.strip():
                        continue
                    # Attempt to parse the line:
                    Bits = FileLine.strip().split("\t")
                    try:
                        ModDBPos = int(Bits[FormatBits.DBPos])
                        ModMass = int(Bits[FormatBits.ModificationMass])
                        Annotation = Bits[FormatBits.Peptide]
                        Charge = int(Bits[FormatBits.Charge])
                    except:
                        traceback.print_exc()
                        continue # skip invalid line
                    # We know our next key, so stop now:
                    self.NextKey = (ModDBPos, ModMass, Annotation, Charge)
                    self.Bits = Bits
                    break
        self.FeatureCursors = []
        # List the directories that need to be parsed, and open the files:
        for SubDirectory in os.listdir(self.PTMFeatureDirectory):
            Dir = os.path.join(self.PTMFeatureDirectory, SubDirectory)
            if not os.path.isdir(Dir):
                continue
            FeatureFilePath = os.path.join(Dir, "PTMFeatures.txt")
            if not os.path.exists(FeatureFilePath):
                print "* Warning: Subdirectory %s doesn't contain a feature file!"%Dir
                continue
            CoverageFilePath = os.path.join(Dir, "Coverage.dat")
            if not os.path.exists(CoverageFilePath):
                print "* Warning: Subdirectory %d doesn't contain a coverage file!"%Dir
                continue
            Cursor = FeatureCursor(FeatureFilePath)
            Cursor.GetNextPeptide()
            self.FeatureCursors.append(Cursor)
        # Output the header lines from an (arbitrary) cursor:
        for HeaderLine in self.FeatureCursors[0].HeaderLines:
            self.OutputFile.write(HeaderLine)
        # Loop through the peptides, until all cursors hit EOF.  At each stage, process
        # the peptide with the first key - and process it from any and all cursors which
        # are now pointing at it.
        while (1):
            CursorsForThisKey = 0
            CursorThisKey = None
            FirstKey = None
            print 
            for Cursor in self.FeatureCursors:
                if Cursor.NextKey != None and (FirstKey == None or FirstKey > Cursor.NextKey):
                    FirstKey = Cursor.NextKey
                    CursorThisKey = Cursor
                    CursorsForThisKey = 0
                if Cursor.NextKey == FirstKey:
                    CursorsForThisKey += 1
                if Cursor.NextKey == FirstKey:
                    print "*Cursor %s: %s"%(Cursor.Path, Cursor.NextKey) 
                else:
                    print "Cursor %s: %s"%(Cursor.Path, Cursor.NextKey) 
            if FirstKey == None:
                break
            Species = self.BuildPeptideSpecies(CursorThisKey)
            print "Next species is '%s', found in %s files."%(Species.Annotation, CursorsForThisKey)
            ########################################################################
            # Shortcut: If CursorsForThisKey == 1, then we don't need to re-build the
            # cluster or the consensus spectrum!
            if CursorsForThisKey == 1:
                self.WriteSingletonPeptide(Species, CursorThisKey)
                CursorThisKey.GetNextPeptide()
                continue
            ########################################################################
            # Standard case: Two or more cursors have the same peptide.  Assimilate the
            # total spectra, and cluster members, from each one!
            for Cursor in self.FeatureCursors:
                if Cursor.NextKey == FirstKey:
                    self.AssimilatePeptideSpectra(Species, Cursor)
                    Cursor.GetNextPeptide()
            # We've assimilated all spectra into the cluster; now write out the
            # totaled spectrum count and the consensus spectra!
            self.WriteCompletedPeptide(Species)
        ########################################################
        # All peptides have been written out; every cursor is at EOF.
        # Loop over cursors and finish up:
        for Cursor in self.FeatureCursors:
            Cursor.Close()
    def CombineDatabaseCoverage(self):
        """
        Iterate over directories, and compute the total coverage for each database residue.
        """
        for SubDirectory in os.listdir(self.PTMFeatureDirectory):
            Dir = os.path.join(self.PTMFeatureDirectory, SubDirectory)
            if not os.path.isdir(Dir):
                continue
            CoverageFilePath = os.path.join(Dir, "Coverage.dat")
            if not os.path.exists(CoverageFilePath):
                print "* Warning: Subdirectory %s doesn't contain a coverage file!"%Dir
                continue
            print "Assimilate from %s..."%CoverageFilePath
            self.AssimilateDatabaseCoverage(CoverageFilePath)
        self.OutputCoverage()

if __name__ == "__main__":
    try:
        import psyco
    except:
        print "(psyco not available - no optimization for you)"
    Merger = PTMFeatureMerger()
    Result = Merger.ParseCommandLine(sys.argv[1:])
    if not Result:
        print UsageInfo
        sys.exit(-1)
    Merger.MergeResults()
    