"""
GridRunner - runs on the grid (mostly sleeping) to push through a lot of Inspect jobs.

Can run in BLIND mode or in STANDARD mode.  If we're running in BLIND mode, then each job
will search only a block of scans.  This means that we may submit several jobs for the same
mzxml file.  We do extra book-keeping to determine when the mzxml file is
completely searched.  This script is always run in combination with BatchRunner, which
copies in files to be searched, and tells which files (spectrum files and results files)
can be deleted

Main work loop:
- Read a list of files to erase: $SCRATCH/DeleteFiles.txt.  Erase them.
- Write a list of results-files with sizes to $SCRATCH/ResultsFiles.txt.  
- Write a list of mzxml-files with sizes to $SCRATCH/MZXMLFiles.txt
- Call the ScriptMongler from ClusterSub.py; it returns a list of job objects
- Check how many jobs are running currently
- Submit up to (JobLimit - CurrentJobCount) scripts
- Write the number of "pending mzxml files" to $SCRATCH/Status.txt (backlog size)
- Sleep for an hour
"""
from ClusterUtils import *
import os
import string
import sys
import getopt
import traceback
import time
import ClusterSub

TempFilePath = os.path.join(ScratchDir, "GridRunnerTemp.txt")

class Bag:
    "Stub class for passing simple bags of data around"
    pass

class GridRunner:
    def __init__(self):
        self.PendingJobs = []
        self.ScanCounts = {} 
        self.FileInfo = {} # stub -> file info object
        self.DBPath = None # For passing to ClusterSub
        self.BlindFlag = 0 
    def LoadScanCounts(self):
        """
        Called at startup, for a blind search.  Parse $SCRATCH/ScanCount.txt, to learn
        how many scans there are in each mzxml file.
        """
        FilePath = os.path.join(ScratchDir, "ScanCount.txt")
        File = open(FilePath, "rb")
        for FileLine in File.xreadlines():
            Bits = FileLine.split("\t")
            try:
                ScanCount = int(Bits[2])
                FileSize = int(Bits[4])
            except:
                continue
            Info = Bag()
            Info.Stub = Bits[1]
            Info.ScanCount = ScanCount
            Info.BlockCount = (ScanCount / BLIND_BLOCK_SIZE) + 1
            Info.FileSize = FileSize
            self.FileInfo[Info.Stub] = Info
            self.ScanCounts[Info.Stub] = ScanCount
        File.close()
    def CheckSpectrumFileSizes(self, Job):
        for FileName in Job.GetSpectrumFileNames():
            Stub = os.path.splitext(FileName)[0]
            Info = self.FileInfo.get(Stub, None)
            if not Info:
                print "UNKNOWN mzxml stub '%s' encountered in job"%Stub
                return 0
            DesiredSize = Info.FileSize
            Path = os.path.join(MZXMLDir, "%s.mzXML"%Stub)
            if not os.path.exists(Path):
                Path = os.path.join(MZXMLDir, "%s.mgf"%Stub)
            if not os.path.exists(Path):
                Path = os.path.join(MZXMLDir, "%s.ms2"%Stub)
            if not os.path.exists(Path):
                print "???? Can't find spectrum-file for pending job!", Stub, Path
                return 0
            Size = os.stat(Path).st_size
            if Size != Info.FileSize:
                print "???? Expected to see filesize %s actual %s for %s"%(Info.FileSize, Size, Path)
                return 0
        return 1 # success!
    def SubmitPendingJobs(self):
        """
        Check the number of jobs running, and then submit some jobs (if there are slots free,
        and if there are jobs waiting to be submitted).  Then write out status.
        """
        self.RunningJobCount = GetRunningJobCount()
        print "Currently running %s (of max %s) jobs"%(self.RunningJobCount, MAXIMUM_JOBS)
        os.chdir(JobDir)
        while self.RunningJobCount < MAXIMUM_JOBS and len(self.PendingJobs):
            Job = self.PendingJobs[0]
            # Check to make sure that this job's spectrum-file(s) has the expected size.
            # We perform this check because it's possible that BatchRunner.py is halfway
            # through pscp-ing the mzxml file to the grid; if that's the case, we mustn't
            # submit the job yet, or the search will be no good.
            Result = self.CheckSpectrumFileSizes(Job)
            if not Result:
                # Kick this job off the list, for now.  It'll be re-added,
                # and when it is, the spectrum file should be large enough!
                self.PendingJobs[1:]
                continue
            # Remove the job from the PENDING list:
            self.PendingJobs = self.PendingJobs[1:]
            # Flag the job as RUNNING:
            self.FlagRunningJob(Job)
            # Submit the job, at last:
            Command = "qsub %s"%Job.FileName
            print Command
            os.system(Command)
            self.RunningJobCount += 1
            time.sleep(1)
    def FlagRunningJob(self, Job):
        """
        Write a file to the Done directory to indicate that this
        job is qsub'd.  No new jobs should be created for these
        blocks.  (This "qsubbed" file semi-redundant with the
        ".running" and ".done" files)
        """
        JobStub = os.path.splitext(Job.FileName)[0]
        RunningPath = os.path.join(DoneDir, "%s.subbed"%JobStub)
        RunningFlagFile = open(RunningPath, "wb")
        Job.FlagRunningJob(JobStub, RunningFlagFile)
        RunningFlagFile.close()
        
    def WriteStatus(self):
        """
        Write our status to the file Status.txt.  The script BatchRunner will read it.
        """
        StatusPath = os.path.join(ScratchDir, "Status.txt")
        File = open(StatusPath, "wb")
        File.write("Running:%s\n"%(self.RunningJobCount))
        File.write("Pending:%s\n"%(len(self.PendingJobs)))
        File.close()
    def PrepareJobFiles(self):
        """
        Execute the ScriptMongler, from ClusterSub.py.
        Write out job scripts, and get a list of jobs which
        still need to be submitted.
        """
        os.chdir(JobDir)
        # Create a ScriptMongler and pass along all the command-line info:
        Mongler = ClusterSub.ScriptMongler()
        Mongler.DBPath = self.DBPath
        Mongler.ScanCounts = self.ScanCounts
        Mongler.BlindSearchFlag = self.BlindFlag
        print "INVOKE THE SCRIPT-MONGLER."
        # Run script-generation:
        try:
            PendingJobList = Mongler.BuildJobs()
        except:
            traceback.print_exc()
            PendingJobList = []
        self.PendingJobs = []
        # We only add jobs to our list if they correspond to one of the
        # spectrum files we want to search.  (There may be STRAY FILES
        # in the mzxml dir):
        for Job in PendingJobList:
            SpectrumFileNames = Job.GetSpectrumFileNames()
            for SpectrumFileName in SpectrumFileNames:
                Stub = os.path.splitext(SpectrumFileName)[0]
                if not self.ScanCounts.has_key(Stub):
                    print "* Skip job on irrelevant mzxml:", Stub
                    continue
            self.PendingJobs.append(Job)
    def ListResultsFiles(self):
        """
        List the files in the results directory to "ResultsFiles.txt".
        BatchRunner will read from ResultsFiles.txt
        """
        Extensions = [".txt"]
        SourceDirectory = os.path.join(ScratchDir, ResultsXDir)
        TargetFile = os.path.join(ScratchDir, "ResultsFiles.txt")
        self.ListFiles(SourceDirectory, Extensions, TargetFile)
    def ListMZXMLFiles(self):
        Extensions = [".mzxml", ".mgf", ".ms2"]
        SourceDirectory = os.path.join(ScratchDir, MZXMLDir)
        TargetFile = os.path.join(ScratchDir, "MZXMLFiles.txt")
        self.ListFiles(SourceDirectory, Extensions, TargetFile)
    def ListFiles(self, SourceDirectory, Extensions, TargetFile):
        """
        Iterate over file names in SourceDirectory.  If the file extension (case-insensitive)
        is one of the extensions listed, then write the file name (and size) to the
        output file.  Used for getting a list of MZXML files to be searched, and a list
        of results-files with their sizes.
        """
        OutputFile = open(TargetFile, "wb")
        for FileName in os.listdir(SourceDirectory):
            (Stub, Extension) = os.path.splitext(FileName)
            if Extension.lower() not in Extensions:
                continue
            FilePath = os.path.join(SourceDirectory, FileName)
            FileSize = os.stat(FilePath).st_size
            OutputFile.write("%s\t%s\t\n"%(FileName, FileSize))
        OutputFile.close()
    def DeleteOldFiles(self):
        """
        In $SCRATCH/DeleteFiles.txt is a list of files that we should delete.
        So, we'll do so.  The BatchRunner script decides what gets deleted.
        """
        # File format: One file-path per line; no special syntax.
        DelFilesPath = os.path.join(ScratchDir, "DeleteFiles.txt")
        try:
            File = open(DelFilesPath)
        except:
            print "DeleteFiles.txt not opened; not deleting anything."
            return
        FileLines = File.readlines()
        File.close()
        DelCount = 0
        for FileLine in FileLines:
            FilePath = FileLine.strip()
            if FilePath and os.path.exists(FilePath):
                try:
                    os.remove(FilePath)
                    DelCount += 1
                    print ">>> Deleted %s"%FilePath
                except:
                    print "** Error: Not able to delete file '%s'"%FilePath
        if DelCount:
            print ">>> Deleted %s old file(s)."%DelCount
    def PrintVerboseStatus(self):
        print "--=={{<<::>>}}==--"*4
        print ">>> Status:"
        Keys = self.FileInfo.keys()
        Keys.sort()
        for Stub in Keys:
            Info = self.FileInfo[Stub]
            PresenceFlag = 0
            Path = os.path.join(MZXMLDir, "%s.mzXML"%Stub)
            if os.path.exists(Path):
                PresenceFlag = 1
            Path = os.path.join(MZXMLDir, "%s.mgf"%Stub)
            if os.path.exists(Path):
                PresenceFlag = 1
            Path = os.path.join(MZXMLDir, "%s.ms2"%Stub)
            if os.path.exists(Path):
                PresenceFlag = 1
        print "Total of %s mzxml files known."%len(Keys)
    def CreateDoneFlagsFromResults(self):
        DoneFileName = os.path.join(DoneDir, "Results.done")
        DoneFile = open(DoneFileName, "wb")
        for FileName in os.listdir(ResultsXDir):
            Stub = os.path.splitext(FileName)[0]
            DoneFile.write("%s\t%s\t\n"%(Stub, Stub))
        DoneFile.close()
    def ParseCommandLine(self):
        "Parse command-line arguments"
        (Options, Args) = getopt.getopt(sys.argv[1:], "d:bx")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-d":
                # -d database filename
                DBFileName = Value
                # Add a .zip extension, if there's no extension:
                if DBFileName.find(".") == -1:
                    DBFileName += ".zip"
                DBPath = os.path.join(ScratchDir, DBFileName)
                if not os.path.exists(DBPath):
                    print "** Error: couldn't find database file '%s'\n\n"%DBPath
                    print UsageInfo
                    sys.exit(1)
                self.DBPath = DBPath
            elif Option == "-b":
                self.BlindFlag = 1
            elif Option == "-x":
                self.CreateDoneFlagsFromResults()
                sys.exit()
            else:
                print "** Unknown option:", Option
        if not self.DBPath or not os.path.exists(self.DBPath):
            print UsageInfo
            sys.exit(-1)
    def MainLoop(self):
        while (1):
            print ">>> Grid runner awake at %s"%time.asctime()
            self.DoStuff()
            time.sleep(ONE_HOUR)
    def DoStuff(self):
        """
        Main update method.
        """
        self.DeleteOldFiles()
        self.ListResultsFiles()
        self.ListMZXMLFiles()
        self.PrepareJobFiles()
        self.SubmitPendingJobs()
        self.WriteStatus()
        self.PrintVerboseStatus()
        
UsageInfo = """
>>> Grid runner - engine to launch grid jobs.
Arguments:
 -b: Activate blind search mode (MS-Alignment algorithm)
 -d: Database file name (assumed to live in the scratch directyory)
"""

if __name__ == "__main__":
    MakeGridDirectories()
    Runner = GridRunner()
    Runner.LoadScanCounts()
    Runner.ParseCommandLine()
    #################
    # Execute:
    Runner.MainLoop()
    