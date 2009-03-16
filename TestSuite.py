"""
The nightly build and test suite is an important tool for releasing robust software.
These tests run automatically every night on the latest code snapshot.  If the
build fails, or if Inspect doesn't work properly, the suite will complain.  Loudly.

The standard way to run the test suite is on a development machine which has a C compiler,
expat, Subversion, and libsvm installed.
- The unattended build procedure works on Windows with Visual Studio 2003 (tested version 7.1.3008),
  and works on a Linux machine with recent Python (2.3+) and gcc.  Other operating systems
  and configurations may require minor changes here or in distutils
- If libsvm is NOT installed, it must be downloaded and installed in ~/libsvm (unix)
  or hard-coded path (windows)
- If subversion is NOT installed, you can still run the test suite by using the "-r" option.
  Make sure that the TestSuite directory (and its contents) are present in the working directory.

"""
import os
import getopt
import shutil
import traceback
import distutils
import distutils.dir_util
import ResultsParser
from Utils import *
Initialize()

SVN_URL = "svn://bioinfo2.ucsd.edu/Inspect"
if sys.platform == "win32":
    InspectExecutableName = "inspect.exe"
    InspectCommand = "inspect"
    CopyCommand = "copy"
else:
    InspectExecutableName = "inspect"
    InspectCommand = "./inspect"
    CopyCommand = "cp"

UsageInfo = """
TestSuite.py: Build inspect, and run system tests.

By default, we copy files from the WORKING DIRECTORY into the testbed.
Use the "-c" option to check out the current version from the repository
instead.  To skip copying altogether, use the "-f" option.

Required arguments:
 -d [DIRECTORY]: Specifies the test bed, where the build will be run.  ANY
   FILES CURRENLY IN THE TEST BED WILL BE REMOVED, unless you specify
   the -c flag!
   The standard place to put this is a sister directory to your Inspect
   directory.
 
Options:
 -f: Don't copy files into the testbed (Use this to skip file copying
     when re-running tests)
 -c: Run tests using the current checked-in files
 -b: Skip build (assumes a current build is in TestBed already)
 -t: Skip tests
 -w [FILENAME]: File to write test results to.  stdout and stderr will be piped
   to the file.

Example (normal run):
TestSuite -w TestResults.txt

Example (tests only):
TestSuite.py -f -b
"""

def NukeDirectory(FileName):
    """
    Remove a directory (and its files and its subdirectories), ignoring errors.
    This is stronger than shutil.rmtree(), which wimps out as soon as it finds
    a read-only file.  (Subversion generates many read-only files)
    """
    if sys.platform == "win32":
        if not os.path.exists(FileName):
            return
        if os.path.isdir(FileName):
            for SubFileName in os.listdir(FileName):
                Path = os.path.join(FileName, SubFileName)
                NukeDirectory(Path)
            os.system("rmdir \"%s\""%FileName)
            return
        Command = "del /F \"%s\" 2> nul"%FileName
        #print Command
        os.system(Command)
    else:
        Command = "rm -rf \"%s\""%FileName
        os.system(Command)

class AnnotationParser(ResultsParser.ResultsParser):
    def ParseFile(self, FileName):
        """
        Build a dictionary, mapping scan numbers to annotation strings.
        """
        self.Annotations = {}
        File = open(FileName, "rb")
        for FileLine in File.xreadlines():
            Bits = FileLine.split("\t")
            try:
                Scan = int(Bits[self.Columns.ScanNumber])
                Annotation = Bits[self.Columns.Annotation]
                ProteinName = Bits[self.Columns.ProteinName]
            except:
                continue
            if self.Annotations.has_key(Scan):
                continue
            AnnotationInfo = Bag()
            self.Annotations[Scan] = AnnotationInfo 
            AnnotationInfo.Annotation = Annotation
            AnnotationInfo.ProteinName = ProteinName
        File.close()

class InspectProctor:
    def __init__(self):
        self.CopyFilesFlag = 1
        self.CheckoutFlag = 0
        self.BuildFlag = 1
        self.TestFlag = 1
        self.ErrorCount = 0
        self.TestBedRoot = None
        self.TestBedDirectory = None
        self.TestFilesDirectory = None
        self.WriteResultsPath = None #os.path.abspath("TestSuiteResults.txt") # default
        self.ErrorFile = None
        self.CopyReleaseFilesFlag = 0
    def ParseCommandLine(self):
        (Options, Args) = getopt.getopt(sys.argv[1:], "cfbtd:w:r")
        for (Option, Value) in Options:
            if Option == "-c":
                self.CheckoutFlag = 1
            elif Option == "-b":
                self.BuildFlag = 0
            elif Option == "-f":
                self.CopyFilesFlag = 0
            elif Option == "-t":
                self.TestFlag = 0
            elif Option == "-d":
                self.TestBedRoot = Value
                self.TestBedDirectory = os.path.join(self.TestBedRoot, "Inspect")
                self.TestFilesDirectory = os.path.join(self.TestBedDirectory, "TestSuite")
            elif Option == "-w":
                self.WriteResultsPath = Value
            elif Option == "-r":
                self.CopyReleaseFilesFlag = 1
        if self.CheckoutFlag and not self.CopyFilesFlag:
            print "Use -c or -f, not both!"
            print UsageInfo
            sys.exit(-1)
        if not self.TestBedRoot:
            print UsageInfo
            sys.exit(-1)
    def CopyReleaseFiles(self):
        """
        Copy, to the test bed, all the files listed in ReleaseFiles.txt
        """
        File = open("ReleaseFiles.txt", "rb")
        for FileLine in File.xreadlines():
            if FileLine[0] == "#":
                continue
            FileName = FileLine.strip()
            if not FileName:
                continue
            if sys.platform == "win32":
                SourceFileName = FileName.replace("/", "\\")
            else:
                SourceFileName = FileName
            DestinationPath = os.path.join(self.TestBedDirectory, SourceFileName)
            # Create subdirectories, if necessary:
            TargetDir = os.path.split(DestinationPath)[0]
            try:
                os.makedirs(TargetDir)
            except:
                pass
            # Copy the file itself:
            Command = "%s %s %s"%(CopyCommand, SourceFileName, DestinationPath)
            try:
                print Command
                os.system(Command)
            except:
                self.ReportError("Unable to copy '%s'"%FileName)
        File.close()
        ###################
        # Also, copy over expat:
        if sys.platform == "win32":
            SourceDir = os.path.join("expat", "lib")
            TargetDir = os.path.join(self.TestBedDirectory, "expat", "lib")
            self.CopyTree(SourceDir, TargetDir)
        ###################
        # Also, copy over testsuite files:
        SourceDir = "TestSuite"
        TargetDir = os.path.join(self.TestBedDirectory, "TestSuite")
        self.CopyTree(SourceDir, TargetDir)
    def CopyTree(self, SourceDir, TargetDir):
        try:
            os.makedirs(TargetDir)
        except:
            pass
        for FileName in os.listdir(SourceDir):
            SourcePath = os.path.join(SourceDir, FileName)
            TargetPath = os.path.join(TargetDir, FileName)
            if os.path.isdir(SourcePath):
                self.CopyTree(SourcePath, TargetPath)
                continue
            Command = "%s \"%s\" \"%s\""%(CopyCommand, SourcePath, TargetPath)
            try:
                print Command
                os.system(Command)
            except:
                self.ReportError("Unable to copy '%s' to '%s'"%(SourcePath, TargetPath))
    def CheckOutCode(self):
        """
        DELETE anything currently in the test bed, and then check out the latest
        stuff from subversion. 
        """
        print "Stomp directory '%s'..."%self.TestBedRoot
        try:
            #shutil.rmtree(self.TestBedRoot, 1)
            NukeDirectory(self.TestBedRoot)
        except:
            traceback.print_exc()
            pass
        MakeDirectory(self.TestBedRoot)
        try:
            os.chdir(self.TestBedRoot)
            Command = "svn checkout %s"%SVN_URL            
            os.system(Command)
        except:
            traceback.print_exc()
            self.ReportError()
            print "* Error checking out code to test bed '%s'"%self.TestBedDirectory
    def BuildCode(self):
        """
        Remove the old inspect executable.  Build a new inspect executable.
        """
        os.chdir(self.TestBedDirectory)
        self.BuildInspectExecutable()
        self.BuildPyInspect()
        self.BuildPySVM()
    def BuildInspectExecutable(self):
        import BuildInspect # delayed import until now, to import the new stuff!
        if os.path.exists(InspectExecutableName):
            try:
                os.remove(InspectExecutableName)
            except:
                traceback.print_exc()
                self.ReportError("Unable to delete inspect.exe!")
        try:
            BuildInspect.BuildInspect(1)
        except:
            traceback.print_exc()
            self.ReportError("Unable to build inspect.exe")
        if not os.path.exists(InspectExecutableName):
            self.ReportError("Inspect executable not present after the build!")
    def BuildPySVM(self):
        BuildDir = "build"
        NukeDirectory(BuildDir)
        # Note: distutils tries to cache paths that it "knows" are present.  We
        # play a bit of a dirty trick on poor distutils by nuking directories
        # ourselves.  The caching isn't needed, so we stop it, by resetting _path_created.
        # (It's a bit naughty to touch _path_created, but if we didn't, the SECOND
        # build - of PySVM - would fail)
        distutils.dir_util._path_created = {}
        try:
            import ReleasePySVM
            ReleasePySVM.Main(["build"])
        except:
            traceback.print_exc()
            self.ReportError("Unable to build PySVM")
            return
        if not os.path.exists(BuildDir):
            self.ReportError("Unable to build PySVM")
            return
        if sys.platform == "win32":
            DesiredFile = "PySVM.pyd"
        else:
            DesiredFile = "PySVM.so"
        if os.path.exists(DesiredFile):
            try:
                os.remove(DesiredFile)
            except:
                self.ReportError("Unable to clear and rebuild %s"%DesiredFile)
                return
        for SubDir in os.listdir(BuildDir):
            SubDirPath = os.path.join(BuildDir, SubDir)
            for SubFile in os.listdir(SubDirPath):
                if SubFile == DesiredFile:
                    shutil.copy(os.path.join(SubDirPath, SubFile), DesiredFile)
                    return
        self.ReportError("Desired library %s not found after build attempt"%DesiredFile)
    def BuildPyInspect(self):
        BuildDir = "build"
        NukeDirectory(BuildDir)
        distutils.dir_util._path_created = {}
        try:
            import ReleasePyInspect
            ReleasePyInspect.Main(["build"])
        except:
            traceback.print_exc()
            self.ReportError("Unable to build PyInspect")
            return
        if not os.path.exists(BuildDir):
            self.ReportError("Unable to build PyInspect")
            return
        if sys.platform == "win32":
            DesiredFile = "PyInspect.pyd"
        else:
            DesiredFile = "PyInspect.so"
        if os.path.exists(DesiredFile):
            try:
                os.remove(DesiredFile)
            except:
                self.ReportError("Unable to clear and rebuild %s"%DesiredFile)
                return
        for SubDir in os.listdir(BuildDir):
            SubDirPath = os.path.join(BuildDir, SubDir)
            for SubFile in os.listdir(SubDirPath):
                if SubFile == DesiredFile:
                    shutil.copy(os.path.join(SubDirPath, SubFile), DesiredFile)
                    return
        self.ReportError("Desired library %s not found after build attempt"%DesiredFile)
    def RunTests(self):
        import SystemTest
        os.chdir(self.TestBedDirectory)
        if not os.path.exists(InspectExecutableName):
            self.ReportError("** No inspect executable - tests aborted!")
            return
        # Run the system tests, which are distributed with Inspect:
        Runner = SystemTest.InspectRunner()
        Runner.RunTests()
        if Runner.ErrorCount:
            self.ReportError("SystemTest failed", Runner.ErrorCount)
        # Run some more demanding tests:
        self.RunOMICS04ModlessTests()
        self.RunOMICS04ModdedTests()
        self.TestMZXMLBlindSearch()
        self.TestPhosphopeptideSearch()
    def DeleteFileIfPresent(self, Path):
        """
        Delete a file (if it currently exists).  Increment ErrorCount if deletion fails.
        """
        if os.path.exists(Path):
            try:
                os.remove(Path)
            except:
                pass
        if os.path.exists(Path):
            self.ReportError("Unable to remove '%s'"%Path)
            return 0
        return 1
    def VerifyFilePresent(self, Path):
        """
        Check to be sure a file is present.  If it's absent, increment ErrorCount.
        """
        if not os.path.exists(Path):
            self.ReportError("Expected file '%s' doesn't exist!"%Path)
            return 0
        ## check to see if the file is 0 kb, that would be bad
        FileSize = os.path.getsize(Path)
        if FileSize == 0:
            self.ReportError("Expected file '%s' is 0 kb"%Path)
            return 0                             
        return 1
    def VerifyFileAbsent(self, Path):
        """
        Check to be sure a file is absent.  If it's present, increment ErrorCount.
        """
        if os.path.exists(Path):
            self.ReportError("Undesirable file '%s' found!"%Path)
            return 0
        return 1
    def RunCommand(self, Command):
        print
        print ">>", Command
        try:
            os.system(Command)
        except:
            traceback.print_exc()
            return 0
        return 1
    def RunOMICS04ModdedTests(self):
        """
        Test database search and PTMFinder functionality.
        - Run an unrestrictive modification search
        - Identify modification sites        
        """
        os.chdir(self.TestBedDirectory)
        # Ensure that the .trie file exists:
        FullTriePath = self.PrepOMICS04Database()
        # Now, run an MS-Alignment search on an .mgf file.
        InspectOutputPath = os.path.join(self.TestBedDirectory, "OMICS04Mod.out")
        InspectInputPath = os.path.join(self.TestBedDirectory, "OMICS04Mod.in")
        InspectErrorPath = os.path.join(self.TestBedDirectory, "OMICS04Mod.err")
        SpectrumPath = os.path.join(self.TestFilesDirectory, "0400-0710mz_02_dta.mgf")
        # Wipe old output:
        self.DeleteFileIfPresent(InspectOutputPath)
        self.DeleteFileIfPresent(InspectErrorPath)
        # Write input params:
        InspectInputFile = open(InspectInputPath, "wb")
        InspectInputFile.write("protease,trypsin\n")
        InspectInputFile.write("db,OMICS04Full.trie\n")
        InspectInputFile.write("mod,+57,C,fix\n")
        InspectInputFile.write("mods,1\n")
        InspectInputFile.write("blind,1,1\n")
        InspectInputFile.write("spectra,%s\n"%SpectrumPath)
        InspectInputFile.close()
        # Run inspect:
        Command = "%s -i OMICS04Mod.in -o OMICS04Mod.out -e OMICS04Mod.err"%(InspectCommand)
        print Command
        if not os.path.exists(InspectOutputPath):
            self.RunCommand(Command)
        # Verify that we have output, and no errors:
        self.VerifyFilePresent(InspectOutputPath)
        self.AssertNoInspectErrors(InspectErrorPath)
        # Prepare a "known mods" file for the PTMAnalysis script
        KnownModPath = "OMICS04KnownMods.txt"
        File = open(KnownModPath, "wb")
        File.write("16\toxidation\tMW\t\n")
        File.write("17\tpyroglutamate\tQ\tN\n")
        File.close()
        # Now, run the PTMAnalysis script.
        if sys.platform == "win32":
            Command = "PTMAnalysis.py"
        else:
            Command = "python PTMAnalysis.py"
        BigDBPath = os.path.abspath(os.path.join("TestSuite", "DictyShort.trie"))
        PTMAnalysisOutputPath = "OMICS04PTMTable.txt"
        Command += " -r %s -d %s"%(InspectOutputPath, FullTriePath)
        Command += " -S 0.5 -s %s -w %s"%(self.TestFilesDirectory, PTMAnalysisOutputPath)
        Command += " -B %s -k %s"%(BigDBPath, KnownModPath)
        Command += " -m lda"
        print Command
        try:
            os.system(Command)
        except:
            traceback.print_exc()
            self.ReportError()
        Result = self.VerifyFilePresent(PTMAnalysisOutputPath)
        if not Result:
            return
        # Parse a couple of sites:
        FoundSiteA = 0
        FoundSiteB = 0
        from TrainPTMFeatures import FormatBits
        File = open(PTMAnalysisOutputPath, "rb")
        for FileLine in File.xreadlines():
            Bits = FileLine.split("\t")
            try:
                PValue = float(Bits[44])
                DBPos = int(Bits[1])
                Mass = int(Bits[2])
            except:
                continue
            # oxidation:
            if PValue < 0.75 and abs(DBPos - 8109) <= 2 and abs(Mass - 16) <= 2:
                FoundSiteA = 1
            # G->D substitution:
            # Note - PValues for this search are not particularly good, mainly because
            # it's hard to fit a mixture model when there are so few data points.
            # That's ok, the false discovery rate is good, and that's what
            # really counts.
            if PValue < 0.75 and abs(DBPos - 2710) <= 2 and abs(Mass - 58) <= 2:
                FoundSiteB = 1
        if not FoundSiteA:
            self.ReportError("PTMAnalysis warning: Didn't find expected site A (K.AFKDEDTQAM+16PFR.V)")
        if not FoundSiteB:
            self.ReportError("PTMAnalysis warning: Didn't find expected site B (K.WENG+58ECAQK.K)")
    def PrepOMICS04Database(self):
        # First, convert from FASTA to TRIE format:
        DatabaseInputPath = os.path.join(self.TestFilesDirectory, "OMICS04.fasta")
        TempTriePath = os.path.join(self.TestFilesDirectory, "OMICS04.trie")
        TempIndexPath = os.path.join(self.TestFilesDirectory, "OMICS04.index")
        TriePath = os.path.join(self.TestBedDirectory, "Database", "OMICS04.trie")
        IndexPath = os.path.join(self.TestBedDirectory, "Database", "OMICS04.index")
        FullTriePath = os.path.join(self.TestBedDirectory, "Database", "OMICS04Full.trie")
        if os.path.exists(FullTriePath):
            return FullTriePath
        self.DeleteFileIfPresent(TriePath)
        if sys.platform == "win32":
            Command = "PrepDB.py"
        else:
            Command = "python PrepDB.py"
        Command += " FASTA \"%s\""%DatabaseInputPath
        self.RunCommand(Command)
        self.VerifyFilePresent(TempTriePath)
        shutil.move(TempTriePath, TriePath)
        self.VerifyFilePresent(TriePath)
        shutil.move(TempIndexPath, IndexPath)
        self.VerifyFilePresent(IndexPath)
        ################################################################################
        # Next, run ShuffleDB to generate a full-length database:
        self.DeleteFileIfPresent(FullTriePath)
        if sys.platform == "win32":
            Command = "ShuffleDB.py"
        else:
            Command = "python ShuffleDB.py"
        Command += " -r \"%s\" -w \"%s\" -p"%(TriePath, FullTriePath)
        self.RunCommand(Command)
        self.VerifyFilePresent(FullTriePath)
        return FullTriePath        
    def TestPhosphopeptideSearch(self):
        """
        Search an mgf file containing some phosphopeptide spectra
        """
        InspectOutputPath = os.path.join(self.TestBedDirectory, "Phospho.out")
        InspectInputPath = os.path.join(self.TestBedDirectory, "Phospho.in")
        InspectErrorPath = os.path.join(self.TestBedDirectory, "Phospho.err")
        SpectrumPath = os.path.join(self.TestFilesDirectory, "Phosphopeptides.mgf")
        DBPath = os.path.abspath(os.path.join(self.TestFilesDirectory, "PhosphoDB.trie"))
        # Wipe old output:
        self.DeleteFileIfPresent(InspectOutputPath)
        self.DeleteFileIfPresent(InspectErrorPath)
        # Write input params:
        InspectInputFile = open(InspectInputPath, "wb")
        InspectInputFile.write("protease,trypsin\n")
        InspectInputFile.write("db,%s\n"%DBPath)
        InspectInputFile.write("mod,+57,C,fix\n")
        InspectInputFile.write("mods,1\n")
        InspectInputFile.write("mod,+80,STY,opt,phosphorylation\n")
        InspectInputFile.write("spectra,%s\n"%SpectrumPath)
        InspectInputFile.close()
        # Run inspect:
        Command = "%s -i Phospho.in -o Phospho.out -e Phospho.err"%(InspectCommand)
        self.RunCommand(Command)
        # Verify that we have output, and no errors:
        self.VerifyFilePresent(InspectOutputPath)
        self.AssertNoInspectErrors(InspectErrorPath)
        # Check the results:
        Parser = AnnotationParser()
        Parser.ProcessResultsFiles(InspectOutputPath, Parser.ParseFile)
        SpotChecks = [(0, "KQNSphosPTTYGASR"),
                      (1, "YDSphosSTDLQAGR"),
                      (2, "VSDDSphosESESGDKEATAPLIQR"),
                      (3, "VHSphosYTDLAYR"),
                      (4, "RGSphosVYHVPLNIVQADAVR"),
                      (5, "LGSphosLVGQDSGYVGGLPK")
                      ]
        for (ScanNumber, DesiredAnnotation) in SpotChecks:
            Match = Parser.Annotations.get(ScanNumber, None)
            if not Match:
                self.ReportError("* Error: phospho search has no match for scan %s (should be %s)"%(ScanNumber, DesiredAnnotation))
                continue
            DesiredPeptide = GetPeptideFromModdedName(DesiredAnnotation)
            ActualPeptide = GetPeptideFromModdedName(Match.Annotation)
            # Be forgiving if one amino acid is missing due to misplaced endpoints:
            Pos = ActualPeptide.Aminos.find(DesiredPeptide.Aminos[2:-2])
            if Pos == -1:
                self.ReportError("* Error: phospho search has '%s' for scan %s (should be %s)"%(Match.Annotation, ScanNumber, DesiredAnnotation))
                continue
            print "Scan %s ok (%s ~= %s)"%(ScanNumber, Match.Annotation, DesiredAnnotation)
    def RunOMICS04ModlessTests(self):
        """
        Test some database search functionality.
        - Prepare a database from a FASTA file (PrepDB.py)
        - Search the database
        - Summarize the proteins present in the database (Summary.py)
        """
        os.chdir(self.TestBedDirectory)
        FullTriePath = self.PrepOMICS04Database()
        ################################################################################
        # Now, run an inspect search on an .mgf file.
        InspectOutputPath = os.path.join(self.TestBedDirectory, "OMICS04.out")
        InspectInputPath = os.path.join(self.TestBedDirectory, "OMICS04.in")
        InspectErrorPath = os.path.join(self.TestBedDirectory, "OMICS04.err")
        SpectrumPath = os.path.join(self.TestFilesDirectory, "0400-0710mz_02_dta.mgf")
        # Wipe old output:
        self.DeleteFileIfPresent(InspectOutputPath)
        self.DeleteFileIfPresent(InspectErrorPath)
        # Write input params:
        InspectInputFile = open(InspectInputPath, "wb")
        InspectInputFile.write("protease,trypsin\n")
        InspectInputFile.write("db,OMICS04Full.trie\n")
        InspectInputFile.write("mod,+57,C,fix\n")
        InspectInputFile.write("spectra,%s\n"%SpectrumPath)
        InspectInputFile.close()
        # Run inspect:
        Command = "%s -i OMICS04.in -o OMICS04.out -e OMICS04.err"%(InspectCommand)
        self.RunCommand(Command)
        
        # Verify that we have output, and no errors:
        self.VerifyFilePresent(InspectOutputPath)
        self.AssertNoInspectErrors(InspectErrorPath)
        DBPath = os.path.join("Database", "OMICS04Full.trie")
        ################################################################################
        # Run Summary.py to list the proteins found.
        if sys.platform == "win32":
            Command = "Summary.py"
        else:
            Command = "python Summary.py"
        SummaryHTMLPath = "OMICSSummary.html"
        self.DeleteFileIfPresent(SummaryHTMLPath)
        Command += " -r \"%s\" -d \"%s\" -w \"%s\" "%(InspectOutputPath, DBPath, SummaryHTMLPath)
        self.RunCommand(Command)
        self.VerifyFilePresent(SummaryHTMLPath)
        SubDBPath = os.path.join("Database", "OMICSSubDB.trie")
        self.DeleteFileIfPresent(SubDBPath)
        Command += " -r \"%s\" -d \"%s\" -b \"%s\" "%(InspectOutputPath, DBPath, SubDBPath)
        self.RunCommand(Command)
        self.VerifyFilePresent(SubDBPath)
        
        ################################################################################
        # Now, run PValue.py to get a collection of annotations for a given
        # false discovery rate.
        PValueOutputPath = os.path.join(self.TestBedDirectory, "OMICS04.5.txt")
        self.DeleteFileIfPresent(PValueOutputPath)
        if sys.platform == "win32":
            Command = "PValue.py"
        else:
            Command = "python PValue.py"
        
        Command += " -r \"%s\" -w \"%s\" -S 0.5 -d \"%s\" -a "%(InspectOutputPath, PValueOutputPath, DBPath)
        self.RunCommand(Command)
        self.VerifyFilePresent(PValueOutputPath)
        # Finally, let's spot-check a few scan numbers to ensure they were correctly annotated.
        Parser = AnnotationParser()
        Parser.ProcessResultsFiles(PValueOutputPath, Parser.ParseFile)
        SpotChecks = [(357, "A.DFAEDKDVCK.N"),
                      (561, "L.ALVYGEATSR.R"),
                      (641, "K.AEFVEVTK.L")]
        for (ScanNumber, DesiredAnnotation) in SpotChecks:
            Match = Parser.Annotations.get(ScanNumber, None)
            if Match == None:
                print "* Erorr: In results file '%s', scan %s has no match (should be %s)"%(PValueOutputPath, ScanNumber, DesiredAnnotation)
            elif Match.Annotation != DesiredAnnotation:
                self.ReportError("* Error: In results file '%s', the annotation for scan %s is %s != %s"%(PValueOutputPath, ScanNumber, Match.Annotation, DesiredAnnotation))
    def TestMZXMLBlindSearch(self):
        """
        Test that we can PARSE an mzXML file, and SEARCH a multi-spectrum file against
        multiple databases in unrestrictive-modifications mode.
        """
        ##K.FQELAQGSTNNDLTSIN-17GLSK.F
        ##M.S+42ITVANQQELNEALATFK.N
        ##R.VETGIIKPGM-48VVTFAPAGLSTEVK.S
        ##K.CD-17LSHVTEYDQLPIQAK.N
        ##K.Y+57LNEQIAQYNQYLQDCR.L
        ##R.QI-17NGQQADYTVNPEIVAQAK.L
        ##K.QI-17DYLLASMPAVEGGFDYNSFAEK.L
        ##K.AM+16AADTTVMDGPDSEGDMFER.K
        os.chdir(self.TestBedDirectory)
        ################################################################################
        # Now, run an inspect search on an .mgf file.
        # Wipe old output:
        InspectOutputPath = os.path.join(self.TestBedDirectory, "DictySmall.out")
        InspectInputPath = os.path.join(self.TestBedDirectory, "DictySmall.in")
        InspectErrorPath = os.path.join(self.TestBedDirectory, "DictySmall.err")
        SpectrumPath = os.path.join(self.TestFilesDirectory, "DictySmall.mzXML")
        self.DeleteFileIfPresent(InspectOutputPath)
        self.DeleteFileIfPresent(InspectErrorPath)
        # Write input params:
        InspectInputFile = open(InspectInputPath, "wb")
        InspectInputFile.write("protease,trypsin\n")
        InspectInputFile.write("sequencefile,%s\n"%os.path.join("Database", "CommonContaminants.fasta"))
        AbsolutePath = os.path.abspath(os.path.join(self.TestFilesDirectory, "DictyShort.trie"))
        InspectInputFile.write("db,%s\n"%AbsolutePath)
        InspectInputFile.write("mod,+57,C,fix\n")
        InspectInputFile.write("blind,1\n")
        InspectInputFile.write("spectra,%s\n"%SpectrumPath)
        InspectInputFile.close()
        # Run inspect:
        Command = "%s -i DictySmall.in -o DictySmall.out -e DictySmall.err"%(InspectCommand)
        self.RunCommand(Command)
        # Verify that we have output, and no errors:
        self.VerifyFilePresent(InspectOutputPath)
        self.AssertNoInspectErrors(InspectErrorPath)
        # Verify that the search results are correct (or 'delta-correct'!)
        Parser = AnnotationParser()
        Parser.ProcessResultsFiles(InspectOutputPath, Parser.ParseFile)
        SpotChecks = [(1, "K.FQELAQGSTNNDLTSIN-17GLSK.F", "DDB0191175 |Protein| gene: cadA"),
                      (3, "R.VETGIIKPGM-48VVTFAPAGLSTEVK.S", "DDB0191134 |Protein| gene: efaAII"),
                      (8, "K.AM+16AADTTVMDGPDSEGDMFER.K", "DDB0184465 |Protein| gene: BEC6V2_0_00694"),
                      ]
        for (ScanNumber, DesiredAnnotation, ProteinName) in SpotChecks:
            Match = Parser.Annotations.get(ScanNumber, None)
            if not Match:
                self.ReportError("* Error: BlindMZXML search has no match for scan %s (should be %s)"%(ScanNumber, DesiredAnnotation))
                continue
            DesiredPeptide = GetPeptideFromModdedName(DesiredAnnotation)
            ActualPeptide = GetPeptideFromModdedName(Match.Annotation)
            # Be forgiving if one amino acid is missing due to misplaced endpoints:
            Pos = ActualPeptide.Aminos.find(DesiredPeptide.Aminos[2:-2])
            if Pos == -1:
                self.ReportError("* Error: BlindMZXML search has '%s' for scan %s (should be %s)"%(Match.Annotation, ScanNumber, DesiredAnnotation))
                continue
            if Match.ProteinName != ProteinName :
                self.ReportError("* Error: BlindMZXML search has protein '%s' for scan %s (should be %s)"%(Match.Annotation, Match.ProteinName, ProteinName))
                continue
            print "Scan %s ok (%s ~= %s)"%(ScanNumber, Match.Annotation, DesiredAnnotation)
    def CopySVNFiles(self):
        if sys.platform == "win32":
            Command = "SubversionCopy.py"
        else:
            Command = "python SubversionCopy.py"
        Command += " -w \"%s\""%self.TestBedDirectory
        self.RunCommand(Command)
    def Main(self):
        """
        Main entry point for the nightly checkout/build/test.
        """
        if self.WriteResultsPath:
            self.OutputFile = open(self.WriteResultsPath, "wb")
            print " Piping stdout and stderr to '%s'"%self.WriteResultsPath
            sys.stderr = self.OutputFile
            sys.stdout = self.OutputFile
        if not self.TestBedDirectory:
            print UsageInfo
            self.ReportError()
            return
        try:
            print "Make '%s'"%self.TestBedDirectory
            os.makedirs(self.TestBedDirectory)
        except:
            pass
        ErrorFilePath = os.path.join(self.TestBedDirectory, "TestSuiteErrors.txt")
        self.ErrorFile = open(ErrorFilePath, "wb")
        # Check out code from SVN, or copy code from the working directory:
        if self.CopyFilesFlag:
            if self.CheckoutFlag:
                self.CheckOutCode()
            elif self.CopyReleaseFilesFlag:
                self.CopyReleaseFiles()
            else:
                self.CopySVNFiles()
        if not os.path.exists(self.TestBedDirectory):
            self.ReportError("** Error: Test bed directory '%s' doesn't exist!"%self.TestBedDirectory)
            return
        # Build Inspect:
        if self.BuildFlag:
            self.BuildCode()
        # Run system tests:
        if self.TestFlag:
            self.RunTests()
    def AssertNoInspectErrors(self, Path):
        if not os.path.exists(Path):
            return
        File = open(Path, "rb")
        for FileLine in File.xreadlines():
            FileLine = FileLine.strip()
            if not FileLine:
                continue
            if FileLine[0] != "{":
                self.ReportError("Inspect error observed at %s:\n %s"%(Path, FileLine))
        File.close()
    def ReportError(self, Message = "", Count = 1):
        self.ErrorCount += Count
        if Message:
            Str = "*@* Error: %s"%Message
        else:
            Str = "*@* Error"
        print Str
        if self.ErrorFile:
            self.ErrorFile.write(Str + "\n")
        traceback.print_stack()
        if self.ErrorFile:
            for String in traceback.format_stack():
                self.ErrorFile.write(String)
            self.ErrorFile.write("\n")
            
if __name__ == "__main__":
    Proctor = InspectProctor()
    Proctor.ParseCommandLine()
    Proctor.Main()
    print "\n\n\n\n"
    if Proctor.ErrorCount == 0:
        print "---- All clear! ----"
        sys.exit(0)
    else:
        print "**** ENCOUNTERED %s ERROR(S) ****"%Proctor.ErrorCount
        sys.exit(-1)