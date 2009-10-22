#!/usr/bin/env python
"""
Generate job scripts that'll search .mzxml files.  The individual searches
will be run by ClusterRunInspect.py.  This script assumes the existence of some
directories:
$SCRATCH/USER_NAME - root directory for searching stuff.  The database zipfile lives here.
$SCRATCH/USER_NAME/Inspect - Inspect should be installed here
$SCRATCH/USER_NAME/jobs - Run this script from this directory!  Job scripts live here.
$SCRATCH/USER_NAME/mzxml - All the mzxml files to be searched live here
$SCRATCH/USER_NAME/Done - A record is kept here of which .mzxml files have already been searched.
$SCRATCH/USER_NAME/output - stdout and stderr from jobs go here.
$SCRATCH/USER_NAME/CopyFlags - Usually empty.  Temporary files are written here as a crude
  semaphore system so that only N jobs perform heavy I/O at once on the NFS drive.
$SCRATCH/USER_NAME/ResultsX - Directory for results-files to be copied into.
"""
import os
import sys
import getopt
import shutil

from CountScans import CountScanBits
import ClusterUtils

# 50mb or so is plenty of stuff to search in one *unmodified* run:
MAX_MZMXML_PER_RUN = 50000000

class JobClass:
    """
    A job corresponds to a single shell-script to be submitted to the grid engine.
    
    TEMP: Each MASTER job has two MAIN sub-jobs, run in parallel, because nbcr
    hates serial jobs..
    
    In an unmodified search:
      Each MAIN job has a list of 1 or more sub jobs.  The subjobs have block number 0.
      The subjobs have mzxml file names attached.
    
    Example:
    Master
    --Main1
    ----FileA.mzxml,0
    ----FileB.mzxml,0
    --Main2
    ----FileC.mzxml,0
      
    In a blind search:
      Each MAIN job has a list of 2 subjobs.  The subjobs have associated mzxml filenames
      and block numbers.
    Example:
    Master
    --Main1: FileA.mzxml,4
    --Main2: FileA.mzxml,5
    """
    def __init__(self):
        self.SubJobs = None
        self.SpectrumFileName = None
        self.BlockNumber = None
        self.File = None
        self.FileName = None
        self.TotalScans = None

    def FlagRunningJob(self, Name, File):
        """
        Traverse the job tree (recursively).  Write out flags to indicate that
        these blocks are "searched" (or "submitted" or whatever)
        """
        if self.SpectrumFileName:
            File.write("%s\t%s\t%s\t\n"%(Name, self.SpectrumFileName, self.BlockNumber))
        if self.SubJobs:
            for SubJob in self.SubJobs:
                SubJob.FlagRunningJob(Name, File)
    def GetSpectrumFileNames(self):
        FileNameList = []
        if self.SubJobs:
            for SubJob in self.SubJobs:
                FileNameList.extend(SubJob.GetSpectrumFileNames())
        if self.SpectrumFileName:
            FileNameList.append(self.SpectrumFileName)
        return FileNameList
      
class ScriptMongler:
    def __init__(self,gridEnv):
        self.CheckDoneFlags = 1
        self.MZXMLListFile = None
        self.DBPath = None 
        self.BlindSearchFlag = 0
        # number of PTMs allowed per peptide; passed along to ClusterRunInspect
        self.PTMLimit = 0 
        self.ScanCounts = {} # stub -> number of scans
        self.IncludeCommonDB = 0
        self.gridEnv = gridEnv
        self.scanCountFile = None
        # Count of the number of master jobs, to be used as job array index
        self.master_job_count = 0

    def createMainJob(self, CurrentMasterJob, SpectrumFileName, BlockNumber):
        MainJob = JobClass()
        MainJob.SpectrumFileName = SpectrumFileName
        MainJob.BlockNumber = BlockNumber
        CurrentMasterJob.SubJobs.append(MainJob)

    def createMasterJob(self, BlockNumber):
        self.master_job_count += 1
        CurrentMasterJob = JobClass()
        CurrentMasterJob.FileName = "%s.in" % self.master_job_count
        JobScriptPath = os.path.join(self.gridEnv.JobDir, CurrentMasterJob.FileName)
        CurrentMasterJob.File = open(JobScriptPath, "w")
        CurrentMasterJob.SubJobs = []
        return CurrentMasterJob


    def closeMasterJob(self, PendingJobList, CurrentMasterJob):
        command = self.BuildInspectInputFile(CurrentMasterJob)
        CurrentMasterJob.File.write(command)
        CurrentMasterJob.File.close()
        PendingJobList.append(CurrentMasterJob)
        CurrentMasterJob = None
        return CurrentMasterJob

    def LoadScanCounts(self):
        """
        Parse $SCRATCH/ScanCount.txt, to learn
        how many scans there are in each mzxml file.
        """
        scanfile = open(self.scanCountFile, "r")
        for fileLine in scanfile.xreadlines():
            Bits = fileLine.split("\t")
            print "LoadScanCounts:", Bits
            try:
                ScanCount = int(Bits[CountScanBits.MaxScan]) 
            except:
                continue
            self.ScanCounts[Bits[CountScanBits.Stub]] = ScanCount
        scanfile.close()

    def BuildInspectInputFile(self,Job):
        """
        Create the body of an inspect input file specifying the spectra, mods, etc.
        Returns a string containing the desired contents of the file.
        """
        SpectraStr = ""
        for SpectrumFileName in Job.GetSpectrumFileNames():
            ScratchMZXMLPath = os.path.join(self.gridEnv.MZXMLDir, SpectrumFileName)
            if self.BlindSearchFlag:
                firstScanNumber = Job.BlockNumber * self.gridEnv.BLIND_BLOCK_SIZE
                lastScanNumber  = firstScanNumber + self.gridEnv.BLIND_BLOCK_SIZE
                SpectraStr += "spectra,%s,%s,%s\n"%(ScratchMZXMLPath, firstScanNumber, lastScanNumber)
            else:
                SpectraStr += "spectra,%s\n"%(ScratchMZXMLPath)
        Script = """
instrument,ESI-ION-TRAP
protease,Trypsin
mod,+57,C,fix
tagcount,25
PMTolerance,3.0
""" 
        Script += "\n"
        Script += SpectraStr
        Script += "\n"
        Script += "mods,%s\n"%self.PTMLimit 
        if self.BlindSearchFlag:
            Script += "blind,1\n"
        Script += "db,%s\n"%self.DBPath

        if self.IncludeCommonDB:
            Script += "db,%s\n"%self.CommonDBPath
        return Script

    def GetAlreadySearchedDict(self):
        """
        Return a dictionary of jobs which have already been searched.
        Keys have the form (MZXMLFileName, BlockNumber)
        """
        self.AlreadySearchedDict = {}
        if self.CheckDoneFlags:
            for DoneFileName in os.listdir(self.gridEnv.DoneDir):
                DoneFile = open(os.path.join(self.gridEnv.DoneDir, DoneFileName), "rb")
                for fileLine in DoneFile.xreadlines():
                    Bits = fileLine.strip().split("\t")
                    try:
                        BlockNumber = int(Bits[2])
                    except:
                        continue
                    Key = (Bits[1], BlockNumber)
                    self.AlreadySearchedDict[Key] = 1
                DoneFile.close()
##        # DEBUG: Print the done flags:
##        print "===ClusterSub done flags:"
##        Keys = self.AlreadySearchedDict.keys()
##        Keys.sort()
##        for Key in Keys:
##            print "  ",Key
##        print "===End of clusterSub done flags."
        return self.AlreadySearchedDict
    def GetSpectrumFileList(self):
        # If we were passed -m argument, look in that file for a list
        # of mzXML file names.  
        if self.MZXMLListFile:
            SpectrumFileNameList = []
            mzfile = open(self.MZXMLListFile)
            for fileLine in mzfile.xreadlines():
                fileLine = fileLine.strip()
                if not fileLine or fileLine[0] == '#':
                    continue
                basename = os.path.basename(fileLine)
                destFile = os.path.join( self.gridEnv.MZXMLDir, basename)
                if not os.path.exists( destFile ):
                    shutil.copy( fileLine, destFile )

                SpectrumFileNameList.append(basename)

            mzfile.close()
        else:
            # List all files from MZXMLDir with a valid extension:
            TempList = os.listdir(self.gridEnv.MZXMLDir)
            SpectrumFileNameList = []
            for Name in TempList:
                (Stub, Extension) = os.path.splitext(Name)
                if Extension.lower() in (".mzxml", ".mgf", ".ms2"):
                    SpectrumFileNameList.append(Name)
        SpectrumFileNameList.sort()
        return SpectrumFileNameList
    def BuildJobs(self):
        """
        Create input scripts (.sh files) to search all the spectra that
        still need searching.  Return a list of JOB objects.  
        """
        print "Mongler invoked: Blind flag %s"%self.BlindSearchFlag
        # Sanity-check our arguments:
        if not self.DBPath:
            print "** Error: Please specify a database path!"
            print UsageInfo
            return
        #JobScript = self.GetJobScript()
        self.GetAlreadySearchedDict()
        SpectrumFileNames = self.GetSpectrumFileList()
        if self.BlindSearchFlag:
            JobList = self.BuildJobsBlindSearch(SpectrumFileNames)
        else:
            JobList = self.BuildJobsStandardSearch(SpectrumFileNames)
        ########
        print "qsub -P %s -cwd -j y -o '$TASK_ID.log' -n inspect -t 1-%i" % (
            self.gridEnv.projectCode, len(JobList) )
        ########
        return JobList
    def BuildJobsStandardSearch(self, SpectrumFileNames):
        PendingJobList = []
        CurrentMasterJob = None
        for SpectrumFileName in SpectrumFileNames:
            (prefix, Extension) = os.path.splitext(SpectrumFileName)
            ScanCount = self.ScanCounts.get(prefix, None)
            if not ScanCount:
                # If we don't know the number of scans, let's assume it's few enough
                # to handle in a single job.
                ScanCount = 100
                print "CS: File %s has UNKNOWN scan-count, guess %s"%(prefix, ScanCount)
            else:
                print "CS: File %s has %s scans"%(prefix, ScanCount)
            BlockNumber = 0
            JobKey = (SpectrumFileName, BlockNumber)
            if self.AlreadySearchedDict.has_key(JobKey):
                print "SKIP already searched:", JobKey
                continue
            if not CurrentMasterJob:
                CurrentMasterJob = self.createMasterJob(BlockNumber)
                CurrentMasterJob.TotalScans = 0
            # Add a MAIN JOB, to search this file:
            self.createMainJob(CurrentMasterJob, SpectrumFileName, BlockNumber)
            CurrentMasterJob.TotalScans += ScanCount
            # If the Master job now has two Main jobs, then finish it off:
            if len(CurrentMasterJob.SubJobs) > 1:
                CurrentMasterJob = self.closeMasterJob(PendingJobList, CurrentMasterJob)
        # The loop over spectrum files and blocks is now complete.
        # If we're in the middle of handling a Master job, then finish it off.
        if CurrentMasterJob:
            CurrentMasterJob = self.closeMasterJob(PendingJobList, CurrentMasterJob)
        return PendingJobList
    def BuildJobsBlindSearch(self, SpectrumFileNames):
        """
        For a blind search, a job may be TOO LARGE and take a long time
        to run.  So, we split the search up into multiple jobs,
        one for each block of scans
        """
        PendingJobList = []
        CurrentMasterJob = None
        for SpectrumFileName in SpectrumFileNames:
            (Stub, Extension) = os.path.splitext(SpectrumFileName)
            ScanCount = self.ScanCounts.get(Stub, None)
            if not ScanCount:
                # If we don't know the number of scans, let's assume it's few enough
                # to handle in a single job.
                ScanCount = 100
                print "CS: File %s has UNKNOWN scan-count, guess %s"%(Stub, ScanCount)
            else:
                print "CS: File %s has %s scans"%(Stub, ScanCount)
            for FirstScanNumber in range(0, ScanCount, self.gridEnv.BLIND_BLOCK_SIZE):
                BlockNumber = FirstScanNumber / self.gridEnv.BLIND_BLOCK_SIZE
                JobKey = (SpectrumFileName, BlockNumber)
                if self.AlreadySearchedDict.has_key(JobKey):
                    print "SKIP already searched:", JobKey
                    continue
                # Create a new MASTER JOB, if necessary:
                if not CurrentMasterJob:
                    CurrentMasterJob = self.createMasterJob(BlockNumber)

                # Add a MAIN JOB, to search this file:
                self.createMainJob(CurrentMasterJob, SpectrumFileName, BlockNumber)
                # If the Master job now has two Main jobs, then finish it off:
                if len(CurrentMasterJob.SubJobs) > 1:
                    CurrentMasterJob = self.closeMasterJob(PendingJobList, CurrentMasterJob)
        # The loop over spectrum files and blocks is now complete.
        # If we're in the middle of handling a Master job, then finish it off.
        if CurrentMasterJob:
            CurrentMasterJob = self.closeMasterJob(PendingJobList, CurrentMasterJob)
        return PendingJobList
    def ParseCommandLine(self):
        """
        Parse command-line arguments.
        """
        (Options, Args) = getopt.getopt(sys.argv[1:], "d:am:bp:c:s:")
        OptionsSeen = {}
        if len(Options) == 0:
            print UsageInfo
            sys.exit(1)

        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-d":
                # -d database path
                if not os.path.exists(Value):
                    print "** Error: couldn't find database archive '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.DBPath = Value
            elif Option == "-c":
                if not os.path.exists(Value):
                    print "** Error: couldn't find database archive '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.IncludeCommonDB = 1
                self.CommonDBPath = Value               
            elif Option == "-m":
                # -m mzxml file name list
                if not os.path.exists(Value):
                    print "** Error: couldn't find mzxml-list file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.MZXMLListFile = Value
            elif Option == "-a":
                # -a Search all mzxml files; IGNORE any done-flags
                self.CheckDoneFlags = 0
            elif Option == "-b":
                # -b blind flag!
                self.BlindSearchFlag = 1
            elif Option == "-p":
                # -p PTM count
                self.PTMLimit = int(Value)
            elif Option == "-s":
                # -s ScanCount.txt file
                self.scanCountFile = Value
                destScanCountPath = os.path.join(self.gridEnv.ScratchDir,"ScanCount.txt")
                if not os.path.exists( destScanCountPath ):
                    shutil.copy( Value, destScanCountPath )
            else:
                print UsageInfo
                sys.exit(1)
                
        if not self.scanCountFile:
            print UsageInfo
            sys.exit(1)

        self.LoadScanCounts()
        # If blind search was requested, but PTM limit is zero, then set it to 1
        if self.BlindSearchFlag and self.PTMLimit == 0:
            print "* Setting PTM limit to 1, since blind search was requested!"
            self.PTMLimit = 1
        
UsageInfo = """
ClusterSub.py - Prepare scripts to search mzXML files over a grid.

Required Parameters:
 -d [DatabaseFile] - Give the full path to a database .zip file.  File named Foo.zip should
   contain Foo.trie and Foo.index (as produced by PrepDB.py)
 -s [FilePath] path to the ScanCount.txt file created by CountScans.py

Options:
 -c [CommonDBFile] - Give the full path to a database.zip file as above.  This is an
     additional database, expected to be the common contaminants.
 -a Search ALL mzxml files in the mzxml directory, ignore flags in the "done" directory
 -m [FileName] - Specify a text file listing the mzxml file names to search.
 -b Blind search - for PTMs!
 -p [Int] number of PTMs desired in the search

See the comments in this script for more details.
"""

if __name__ == "__main__":
    gridEnv = ClusterUtils.JCVIGridEnv()
    gridEnv.MakeGridDirectories()
    # Parse command-line options; see UsageInfo for explanations
    Mongler = ScriptMongler(gridEnv)
    Mongler.ParseCommandLine()
    Mongler.BuildJobs()
