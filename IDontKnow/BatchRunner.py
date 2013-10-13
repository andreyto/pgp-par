"""
BatchRunner - runs on bioinfo2 (or some other master PC), communicates with
GridRunner (which runs on the grid).  Keeps copying over files and launching jobs.
It's assumed that ONLY ONE batch runs on a grid for a user at one time.

THE PROCESS:
- Subclass BatchRunner in BatchRunnerXXX.py
- Generate ScanCount.txt using the -s argument
- Copy ScanCount.txt to $SCRATCH on the grid
- Copy XXXHalf.zip (database zipfile) to $SCRATCH on the grid
- On local PC, run BatchRunnerXXX.py
- On grid, from jobs directory, run nice python ../GridRunner.py -d PATH_TO_DATABASE -b

Update tasks:
- Copy ResultsFiles.txt from the grid, if available.  
- Copy MZXMLFiles.txt from the grid, if available
- Copy Status.txt from the grid, if available
- Copy down from the grid each results file we don't have (with correct size) yet.
  If we have the file with the correct size, add the results file and the mzxml file to the delete-list.
- Copy the delete-list to the grid
- If the MZXML backlog is 30 files or more, stop
- List all mzxml files which haven't been searched and aren't on the grid.  Copy up to 30-BACKLOGSIZE
  files to the grid.
- Sleep for an hour
- Also: Get "true status", parsing results files to get the scan numbers seen.
"""
import os
import time
import sys
import traceback
from ClusterUtils import *
import CountSpectraFromMZXML
import CountSpectraFromMGF

class Bag:
    pass

class BatchRunner:
    def __init__(self):
        self.MZXMLPaths = []
        self.FileInfo = {}
        self.ScanCounts = {}
        # ResultsDirectory is the full path to the directory where output files will be upt.
        # It must be set by the caller
        self.ResultsDirectory = ""
        self.GridRunningCount = None
        self.GridPendingCount = None
        # We use ReverseOrderFlag if we want TWO grid runners, one moving top-to-bottom
        # and one moving bottom-to-top.  They'll meet in the middle.
        self.ReverseOrderFlag = 0
        self.FilesToDelete = []
    def InitFileInfo(self):
        self.FileInfo = {}
        for MZXMLPath in self.MZXMLPaths:
            Info = Bag()
            (Dir, FileName) = os.path.split(MZXMLPath)
            (Stub, Extension) = os.path.splitext(FileName)
            Info.FileName = FileName
            Info.FilePath = MZXMLPath
            Info.Stub = Stub
            Info.AllDoneFlag = 0
            self.FileInfo[Stub] = Info
            Info.SearchScanCount = 0
            Info.SearchMinScanNumber = 9999999
            Info.SearchMaxScanNumber = 0
    def LoadScanCounts(self):
        """
        Called at startup, for a blind search.  Parse $SCRATCH/ScanCount.txt, to learn
        how many scans there are in each mzxml file.
        """
        FilePath = "ScanCount.txt"
        File = open(FilePath, "rb")
        for FileLine in File.xreadlines():
            Bits = FileLine.split("\t")
            # Name, Stub, Max, Count, FileSize
            try:
                #MinScan = int(Bits[1])
                MaxScan = int(Bits[2])
                ScanCount = int(Bits[3])
                FileSize = int(Bits[4])
            except:
                traceback.print_exc()
                continue
            Info = self.FileInfo.get(Bits[1], None)
            if not Info:
                #print "** Warning: Unknown filename in ScanCount.txt!", Bits
                # Temp: don't warn
                continue
            else:
                Info.MinScanNumber = 0
                Info.MaxScanNumber = MaxScan
                Info.ScanCount = ScanCount
                Info.FileSize = FileSize
                #self.ScanCounts[Bits[0]] = ScanCount
        File.close()
        ##################
        # Initialize RemoteResultsSize and ResultsSize lists:
        for (Key, Info) in self.FileInfo.items():
            if not hasattr(Info, "MaxScanNumber"):
                print "* Warning: Max scan number not known for %s %s"%(Key, Info.FilePath)
                continue
            BlockCount = (Info.MaxScanNumber / BLIND_BLOCK_SIZE) + 1
            Info.RemoteResultsSize = [None] * BlockCount
            Info.ResultsSize = [None] * BlockCount
    def SetGridNBCR(self):
        "Set some HARD-CODED paths for the grid."
        self.UserName = "stanner"
        self.ScratchDir = "/nas3/stanner"
        self.HostName = "kry"
        self.MZXMLDir = "/nas3/stanner/mzxml"
        self.GridResultsDir = "/nas3/stanner/ResultsX"
    def SetGridNitin(self):
        "Set some HARD-CODED paths for the grid."
        self.UserName = "ngupta"
        self.ScratchDir = "/nas/ngupta"
        self.HostName = "Nitin"
        self.MZXMLDir = "/nas/ngupta/mzxml"
        self.GridResultsDir = "/nas/ngupta/ResultsX"
    def SetGridFWNitin(self):
        "Set some HARD-CODED paths for the grid."
        self.UserName = "ngupta"
        self.ScratchDir = "/scratch/ngupta"
        self.HostName = "fwnitin"
        self.MZXMLDir = "/scratch/ngupta/mzxml"
        self.GridResultsDir = "/scratch/ngupta/ResultsX"
    def SetGridSam(self):
        "Set some HARD-CODED paths for the grid."
        self.UserName = "fwg.spayne"
        self.ScratchDir = "/scratch/spayne"
        self.HostName = "Sam"
        self.MZXMLDir = "/scratch/spayne/mzxml"
        self.GridResultsDir = "/scratch/spayne/ResultsX"
    def SetGridSamK(self):
        "Set some HARD-CODED paths for the grid."
        self.UserName = "spayne"
        self.ScratchDir = "/nas3/spayne"
        self.HostName = "samk"
        self.MZXMLDir = "/nas3/spayne/mzxml"
        self.GridResultsDir = "/nas3/spayne/ResultsX"
    def SetGridFW(self):
        "Set some HARD-CODED paths for the grid."
        self.UserName = "fwg.stanner"
        self.ScratchDir = "/scratch/stanner"
        self.HostName = "grid"
        self.MZXMLDir = "/scratch/stanner/mzxml"
        self.GridResultsDir = "/scratch/stanner/ResultsX"
    def CopyGridFile(self, FileName):
        RemotePath = "%s/%s"%(self.ScratchDir, FileName)
        Command = "pscp %s:%s ."%(self.HostName, RemotePath)
        try:
            os.system(Command)
        except:
            traceback.print_exc()
    def ParseResultsFileList(self):
        """
        We've copied ResultsFiles.txt down from the server.  Each line of ResultsFiles
        gives the name of a results file and its size.  We parse through these lines.
        If there are results files we don't have yet, we copy them down.  Soon we'll
        tell the server to erase those results-files which we already have.
        """
        ResultsFileName = "ResultsFiles.txt"
        if not os.path.exists(ResultsFileName):
            print "(No results file-list copied from grid yet)"
            return
        File = open(ResultsFileName, "rb")
        for FileLine in File.xreadlines():
            Bits = FileLine.split("\t")
            try:
                FileSize = int(Bits[1])
            except:
                continue
            FileName = Bits[0]
            RemotePath = "%s/%s"%(self.GridResultsDir, FileName)
            (Stub, Extension) = os.path.splitext(FileName)
            BlockNumber = 0
            Info = self.FileInfo.get(Stub, None)
            if not Info:
                # There's block info on this job, it seems:
                try:
                    (Stub, BlockNumber) = os.path.splitext(Stub)
                    BlockNumber = int(BlockNumber[1:]) # skip .
                    Info = self.FileInfo[Stub]
                except:
                    print "(Ignoring unknown results-file from grid: %s)"%FileName
                    continue
            Info.RemoteResultsSize[BlockNumber] = FileSize
            if Info.ResultsSize[BlockNumber] == None:
                # We haven't got a results-file for this mzxml file yet.  That means
                # this is a file we want to copy!
                pass
            elif Info.ResultsSize[BlockNumber] > Info.RemoteResultsSize[BlockNumber]:
                print "?? Strange: Local results file size %s > %s for file %s block %s"%(Info.ResultsSize[BlockNumber], Info.RemoteResultsSize[BlockNumber], Info.FilePath, BlockNumber)
                continue
            elif Info.ResultsSize[BlockNumber] == Info.RemoteResultsSize[BlockNumber]:
                # We've already got this one: Add it to the list of files to delete!
                self.FilesToDelete.append(RemotePath)
                continue
            elif Info.ResultsSize[BlockNumber] < Info.RemoteResultsSize[BlockNumber]:
                # We have a SMALLER file than the real results-file!
                pass
            #############################################################
            # If we reached this point in the loop, then we do want to copy down
            # the results-file from the grid:
            TargetPath = os.path.join(self.ResultsDirectory, FileName)
            Command = "pscp %s:%s %s"%(self.HostName, RemotePath, TargetPath)
            print Command
            try:
                os.system(Command)
            except:
                traceback.print_exc()
            # Now we've copied the file; verify the file size:
            if os.path.exists(TargetPath):
                Info.ResultsSize[BlockNumber] = os.stat(TargetPath).st_size
                if Info.ResultsSize[BlockNumber] == Info.RemoteResultsSize[BlockNumber]:
                    # We've got the results file, now it can be deleted from the server
                    self.FilesToDelete.append(RemotePath)
    def PrintVerboseStatus(self):
        print
        print "--=={{<<::>>}}==--"*4
        print ">>> Status:"
        FinishedCount = 0
        StartedCount = 0
        NotStartedCount = 0
        BlockCount = 0
        Keys = self.FileInfo.keys()
        Keys.sort()
        BatchRunnerStatusFile = open("BatchRunStatus.txt", "wb")
        for Key in Keys:
            Info = self.FileInfo[Key]
            # Decide what state this mzxml file is in.
            if Info.AllDoneFlag:
                FinishedCount += 1
                BatchRunnerStatusFile.write("  %s: Complete!\n"%(Key))
                continue
            Str = "|"
            StartedFlag = 0
            for Index in range(len(Info.RemoteResultsSize)):
                if Info.RemoteResultsSize[Index]:
                    StartedFlag = 1
                    if Info.RemoteResultsSize[Index] == Info.ResultsSize[Index]:
                        Str += " "
                    else:
                        Str += "."
                        BlockCount += 1
                else:
                    if Info.ResultsSize[Index]:
                        StartedFlag = 1
                        Str += " "
                    else:
                        Str += "#"
                        BlockCount += 1
            Str += "|"
            BatchRunnerStatusFile.write(" %s: %s\n"%(Key, Str))
            if StartedFlag:
                StartedCount += 1
            else:
                NotStartedCount += 1
        Str = "Summary: %s mzxml files.  %s complete, %s started, %s not started"%(len(Keys), FinishedCount, StartedCount, NotStartedCount)
        BatchRunnerStatusFile.write(Str + "\n")
        print Str
        Str = "%s blocks remain."%BlockCount
        BatchRunnerStatusFile.write(Str + "\n")
        print Str
        BatchRunnerStatusFile.close()
    def DoStuff(self):
        """
        Main update method.
        """
        # Get status from the grid:
        self.CopyGridFile("ResultsFiles.txt")
        self.CopyGridFile("MZXMLFiles.txt")
        self.CopyGridFile("Status.txt")
        self.ParseGridStatus()
        # Parse the list of results files, and copy down any files we don't already have:
        self.FilesToDelete = []
        self.ParseResultsFileList()
        self.CheckCompletedMZXMLFiles()
        self.DeleteUnwantedGridFiles()
        self.CopyMZXMLToGrid()
        self.PrintVerboseStatus()
        #self.CopyResultsFromGrid()
    def ParseGridStatus(self):
        """
        Parse Status.txt, where the grid job sends us general information such
        as the number of jobs running, the number of jobs awaiting submission.
        """
        try:
            File = open("Status.txt", "rb")
        except:
            print "* Status.txt not present; the grid may not be running yet."
            self.GridPendingCount = None
            self.GridPendingCount = None
            return
        for FileLine in File.xreadlines():
            Bits = FileLine.split(":")
            if len(Bits) != 2:
                continue
            if Bits[0] == "Running":
                self.GridRunningCount = int(Bits[1])
            elif Bits[0] == "Pending":
                self.GridPendingCount = int(Bits[1])
            else:
                print "* Warning: Skipping unknown line from status file:\n%s"%FileLine
        File.close()
    def CopyMZXMLToGrid(self):
        """
        Copy MZXML files to the grid, for searching.  We only copy mzxml files
        if the number of pending jobs is from 0 to 20.  That report comes
        from Status.txt (parsed in ParseGridStatus); so, we won't copy anything
        over until GridRunner gives us a report.
        """
        if self.GridPendingCount == None or self.GridPendingCount > 20:
            print "Not copying over any mzxml files.  (Pending count: %s)"%self.GridPendingCount
            return
        FilesCopied = 0
        # Search files in predictable (alphabetical) order
        Keys = self.FileInfo.keys()
        Keys.sort()
        if self.ReverseOrderFlag:
            Keys.reverse()
        for Key in Keys:
            Info = self.FileInfo[Key]
            if Info.AllDoneFlag:
                # This file has been searched already!
                continue
            if self.GridMZXMLFlags.get(Info.Stub, None):
                # This one's already on the grid!
                continue
            # SPECIAL CASE: If we have a pre-specified list of files to search,
            # then only select those files.
            Stub = Info.Stub
            if hasattr(self,"DesiredStubs") and not self.DesiredStubs.has_key(Stub):
                continue
            Command = "pscp %s %s:%s"%(Info.FilePath, self.HostName, self.MZXMLDir)
            print Command
            try:
                os.system(Command)
            except:
                traceback.print_exc()
            FilesCopied += 1
            if FilesCopied >= 10:
                break
        print "...Copied %s mzxml files to the grid for searching."%(FilesCopied)
    def DeleteUnwantedGridFiles(self):
        """
        In self.FilesToDelete, we have a list of (full paths to) files on the grid
        which can be deleted to save space.  (Key idea: We want to be able to launch
        a long line of searches such that storing all the mzxml files and results files
        on the grid at once is not feasible!)
        """
        File = open("DeleteFiles.txt", "wb")
        for FilePath in self.FilesToDelete:
            File.write("%s\n"%FilePath)
        File.close()
        print "There are %s file(s) awaiting destruction."%len(self.FilesToDelete)
        Command = "pscp DeleteFiles.txt %s:%s"%(self.HostName, self.ScratchDir)
        print Command
        os.system(Command)
    def CheckCompletedMZXMLFiles(self):
        """
        Check each mzxml file to see whether the searches are all done.
        Then parse MZXMLFiles.txt (a list of mzxml files on the grid).
        For each file on the
        list which is part of our job-batch and for which all results files are present and
        accounted for, add the mzxml file to the delete-list.
        """
        print ">>>CheckCompletedMZXMLFiles:"
        for Info in self.FileInfo.values():
            # Check whether we already have all the results for this guy.
            # Note that the grid may not STILL have the results files. 
            # We assume a results file is ok if its filesize is non-zero
            # and if there's not a larger results-file on the grid.
            AllResultsOK = 1
            #print "Check file %s..."%Info.FileName
            for Index in range(len(Info.ResultsSize)):
                if not Info.ResultsSize[Index]:
                    #print ":( No results for block %s!"%Index
                    AllResultsOK = 0
                if Info.ResultsSize[Index] < Info.RemoteResultsSize[Index]:
                    #print ":( Remote results results for block %s are large (%s vs %s)"%(Index, Info.ResultsSize[Index], Info.RemoteResultsSize[Index])
                    AllResultsOK = 0
            if AllResultsOK:
                Info.AllDoneFlag = 1
            else:
                Info.AllDoneFlag = 0
        ########################
        GridMZXMLFile = "MZXMLFiles.txt"
        if not os.path.exists(GridMZXMLFile):
            print "(Grid mzxml-file list not present)"
            return
        self.GridMZXMLFlags = {}
        File = open(GridMZXMLFile, "rb")
        for FileLine in File.xreadlines():
            Bits = FileLine.split("\t")
            if len(Bits) < 2:
                continue # blank line
            FileName = Bits[0]
            (Stub, Extension) = os.path.splitext(FileName)
            Info = self.FileInfo.get(Stub, None)
            if not Info:
                continue
            if Info.AllDoneFlag:
                RemotePath = "%s/%s"%(self.MZXMLDir, Info.FileName)
                self.FilesToDelete.append(RemotePath)
            self.GridMZXMLFlags[Stub] = 1
        File.close()
    def CheckCurrentResultsFiles(self):
        """
        At startup: Iterate over files in our results-directory, and check their sizes.
        This way, we know which files have already been searched.
        """
        print "-="*33
        print "Parse current results..."
        for FileName in os.listdir(self.ResultsDirectory):
            FilePath = os.path.join(self.ResultsDirectory, FileName)
            (Stub, Extension) = os.path.splitext(FileName)
            Info = self.FileInfo.get(Stub, None)
            BlockNumber = 0 
            if not Info:
                try:
                    (Stub, BlockNumber) = os.path.splitext(Stub)
                    BlockNumber = int(BlockNumber[1:]) # skip .
                    Info = self.FileInfo[Stub]
                except:
                    #traceback.print_exc()
                    #print "Ignoring unexpected file from results directory:", FileName
                    continue
            Size = os.stat(FilePath).st_size
            print Stub, BlockNumber, Size
            Info.ResultsSize[BlockNumber] = Size
        print "-="*33
    def Main(self):
        while (1):
            print ">>> Batch runner awake at %s"%time.asctime()
            self.InitFileInfo()
            self.LoadScanCounts()
##            for FilePath in self.MZXMLPaths:
##                FileName = os.path.split(FilePath)[1]
##                (Stub, Extension) = os.path.splitext(FileName)
##                FileInfo = Bag()
##                self.FileInfo[Stub] = FileInfo 
##                FileInfo.FilePath = FilePath
##                FileInfo.Stub = Stub
##                FileInfo.ResultsSize = None # unknown
##                FileInfo.RemoteResultsSize = None # unknown
            self.CheckCurrentResultsFiles()
            self.CheckCompletedMZXMLFiles()
            self.DoStuff()
            print "Sleeping for one hour."
            time.sleep(ONE_HOUR)
    def GetScanCounts(self):
        """
        A one-time parsing of all mzxml files to count the number of scans in each
        file.  (We need to know the number of scans so that we know how many
        search blocks will be done for that file; only 10,000 spectra are searched
        at one time!)
        """
        ScanCountFile = open("ScanCount.txt", "wb")
        # Write out the scan counts which we *do* know:
        for (Key, Info) in self.FileInfo.items():
            if hasattr(Info, "MaxScanNumber"):
                Str = "%s\t%s\t%s\t%s\t%s\t%s\t"%(Info.Stub,
                    Info.MinScanNumber, Info.MaxScanNumber, Info.ScanCount,
                    Info.FilePath, Info.FileSize)
                ScanCountFile.write(Str + "\n")
                print "Known:", Info.FileName
            else:
                print "NOT Known:", Info.FileName
        # Now look up the scan counts which we *don't* know:
        for (Key, Info) in self.FileInfo.items():
            if hasattr(Info, "MaxScanNumber"):
                continue
            print "Parse scan numbers from %s..."%(Info.FilePath)
            (Stub, Extension) = os.path.splitext(Info.FileName)
            if Extension.lower() == ".mgf":
                (MinScanNumber, MaxScanNumber, ScanCount) = CountSpectraFromMGF.CountScans(Info.FilePath)
            elif Extension.lower() == ".ms2":
                (MinScanNumber, MaxScanNumber, ScanCount) = CountSpectraFromMGF.CountScansMS2(Info.FilePath)
            else:
                (MinScanNumber, MaxScanNumber, ScanCount) = CountSpectraFromMZXML.CountScans(Info.FilePath)
            FullFileSize = os.stat(Info.FilePath).st_size
            Info.MinScanNumber = MinScanNumber
            Info.MaxScanNumber = MaxScanNumber
            Info.ScanCount = ScanCount
            Info.FileSize = FullFileSize
            Str = "%s\t%s\t%s\t%s\t%s\t%s\t"%(Stub, MinScanNumber, MaxScanNumber, ScanCount, Info.FilePath, FullFileSize)
            ScanCountFile.write(Str + "\n")
            ScanCountFile.flush()
        ScanCountFile.close()
    def AssessProgress(self):
        #FileInfo = {}
        #print self.FileInfo.keys()
        for MZXMLPath in self.MZXMLPaths:
            FileName = os.path.split(MZXMLPath)[-1]
            (Stub, Extension) = os.path.splitext(FileName)
        OutputFile = open("BatchRunnerProgress.txt", "wb")
        ResultsFileNames = os.listdir(self.ResultsDirectory)
        for ResultsFileName in ResultsFileNames: 
            (Stub, Extension) = os.path.splitext(ResultsFileName)
            if not self.FileInfo.has_key(Stub):
                (Stub2, Extension) = os.path.splitext(Stub)
                if self.FileInfo.has_key(Stub2):
                    Stub = Stub2
                else:
                    print "** Skipping unexpected results-dir file:", ResultsFileName
                    continue
            Info = self.FileInfo[Stub]
            Path = os.path.join(self.ResultsDirectory, ResultsFileName)
            File = open(Path, "rb")
            ScanCount = 0
            OldSpectrum = None
            MinScanNumber = 999999
            MaxScanNumber = -99999
            for FileLine in File.xreadlines():
                Bits = FileLine.split("\t")
                if len(Bits) < 3:
                    continue
                try:
                    ScanNumber = int(Bits[1])
                except:
                    continue
                Spectrum = (Bits[0], Bits[1])
                if Spectrum != OldSpectrum:
                    OldSpectrum = Spectrum
                    MinScanNumber = min(MinScanNumber, ScanNumber)
                    MaxScanNumber = max(MaxScanNumber, ScanNumber)
                    ScanCount += 1
            File.close()
            print "Parsed results from %s: %s scans searched"%(Path, ScanCount)
            Info.SearchMinScanNumber = min(MinScanNumber, Info.SearchMinScanNumber)
            Info.SearchMaxScanNumber = max(MaxScanNumber, Info.SearchMaxScanNumber)
            Info.SearchScanCount += ScanCount
        for MZXMLPath in self.MZXMLPaths:
            (FileDir, FileName) = os.path.split(MZXMLPath)
            (Stub, Extension) = os.path.splitext(FileName)
            if not self.FileInfo.has_key(Stub):
                print "???"
                continue # temp, for speedy-parsing
            Info = self.FileInfo[Stub]
            Str = "%s\t%s\t"%(FileDir, FileName)
            Str += "%s\t%s\t%s\t"%(Info.MinScanNumber, Info.MaxScanNumber, Info.ScanCount)
            Str += "%s\t%s\t%s\t"%(getattr(Info, "SearchMinScanNumber", ""),
                                   getattr(Info, "SearchMaxScanNumber", ""),
                                   getattr(Info, "SearchScanCount", ""))
            print Str
            OutputFile.write(Str + "\n")


    