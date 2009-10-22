"""
Main job-script for clusters.  

Don't run this program directly.  Use ClusterSub to generate job scripts, then submit
those scripts using the qsub command, then the job script will invoke this program.
And then this program runs Inspect.  (It sounds complicated at first, and it is.  But
such is the price of using the grid!)

Workflow:
- Set the "running" file in DoneFlags
- Make a scratch directory
- Copy spectrum file(s) to the scratch directory
- Copy database to the scratch directory
- Run inspect
- Copy the results back
- Set the "done" file in DoneFlags
- Clear the "running" file from DoneFlags

"""
import os
import sys
import traceback
import time
import getopt
import getpass

import ClusterUtils

# wait no longer than this for other jobs to finish copying:
MAX_WAIT_TIME = 3600

class JobRunner:
    def __init__(self,gridEnv):
        self.Stub = "XX"
        self.UniqueStub = "XXXXX"
        self.ScratchDir = None
        self.BlindSearchFlag = 0
        self.SpliceSearchFlag = 0
        self.SpectrumFileNames = []
        self.DBPath = None
        self.FirstScanNumber = None
        self.LastScanNumber = None
        self.BlockNumber = None
        self.PTMLimit = 0
        self.IncludeCommonDB =0
        self.CommonDBPath = None 
        self.gridEnv = gridEnv
    def CopyFilePolitely(self, SourcePath, DestinationPath):
        """
        Copy a file, respecting the limit on how many files may be
        copied at once.
        """
        # FlagFilePath must be UNIQUE!  Our stub is NOT unique!
        FlagFilePath = os.path.join(self.gridEnv.CopyFlagDir, self.UniqueStub)
        if os.path.exists(FlagFilePath):
            os.remove(FlagFilePath) # stomp old dead stuff
        StartTime = time.clock()
        while (1):
            FileNames = os.listdir(self.gridEnv.CopyFlagDir)
            if len(FileNames) >= self.gridEnv.MAX_SIMULTANEOUS_COPIES:
                # Copying operations are in progress already.
                if time.clock() - StartTime > MAX_WAIT_TIME:
                    print "** Tired of waiting to copy files; bail out!"
                    return 0
                # Wait for the other jobs to finish:
                time.sleep(120)
                continue
            # CLAIM a copying spot:
            File = open(FlagFilePath, "w")
            File.write("%s"%time.asctime())
            File.close()
            CopyStartTime = time.clock()
            try:
                FileSize = os.stat(SourcePath).st_size
            except:
                print "** Error: Can't get filesize for source file '%s'"%SourcePath
                FileSize = 0
            print "Copy %s (filesize %s)"%(SourcePath, FileSize)
            try:
                Command = "cp %s %s"%(SourcePath, DestinationPath)
                print Command            
                os.system(Command)
            finally:
                # RELINQUISH the copying slot:
                os.remove(FlagFilePath)
            try:
                DestFileSize = os.stat(DestinationPath).st_size
            except:
                print "** Error: Can't get filesize for destination file '%s'"%DestinationPath
                DestFileSize = 0
            if FileSize != DestFileSize:
                print "** Warning: Destination file %s has size %s != %s"%(DestinationPath, DestFileSize, FileSize)
            Now = time.clock()
            print "Copied %s to %s - took %s seconds in all, %s seconds in queue"%(SourcePath, DestinationPath,
                Now - StartTime, CopyStartTime - StartTime)
            sys.stdout.flush()
            return 1
    def GetScratchDirectory(self):
        TempFileName = "%s.temp"%self.UniqueStub
        if os.path.exists(TempFileName):
            os.remove(TempFileName)

        File = open(TempFileName, "rb")
        TempDirectory = File.read().strip()
        File.close()
        os.remove(TempFileName)
        return TempDirectory
    def CleanupScratchDir(self):
        if self.ScratchDir:
            Command = "rm -rf %s"%(self.ScratchDir)
            print Command
            try:
                os.system(Command)
            except:
                traceback.print_exc()
    def FlagJob(self, DoneFlag = 0):
        """
        Write a file to the Done directory, to indicate that a block is searching
        or complete.
        """
        if DoneFlag:
            TypeString = "done"
        else:
            TypeString = "running"
        FlagPath = os.path.join(self.gridEnv.DoneDir, "%s.%s"%(self.UniqueStub, TypeString))
        print ">>>Flag job:", FlagPath
        # Flag these mzxml files as "in progress":
        RunningFile = open(FlagPath, "wb")
        for SpectrumFileName in self.SpectrumFileNames:
            if self.FirstScanNumber:
                BlockNumber = self.FirstScanNumber / self.gridEnv.BLIND_BLOCK_SIZE
            else:
                BlockNumber = 0
            RunningFile.write("%s\t%s\t%s\t\n"%(self.UniqueStub, SpectrumFileName, BlockNumber))
        RunningFile.close()
        return FlagPath
    def CopyFilesToNode(self):
        """
        Copy over spectrum files and database to the compute node
        from the master scratch drive.  (Motivation: Only 2-3 jobs at a time
        can do heavy I/O on the master drive)
        """
        DBFileName = os.path.split(self.DBPath)[1]
        DBStub = os.path.splitext(DBFileName)[0]        
        # Copy over the MZXML files, for searching:
        for SpectrumFileName in self.SpectrumFileNames:
            MZXMLPath = os.path.join(self.gridEnv.MZXMLDir, SpectrumFileName)
            ScratchMZXMLPath = os.path.join(self.ScratchDir, SpectrumFileName)
            self.CopyFilePolitely(MZXMLPath, ScratchMZXMLPath)
            if not os.path.exists(ScratchMZXMLPath):
                print "** Error: MZXML file not found at '%s'"%ScratchMZXMLPath
                return
        # Copy over, and decompress, the database:
        TempZipPath = os.path.join(self.ScratchDir, DBFileName)
        self.CopyFilePolitely(self.DBPath, TempZipPath)
        os.chdir(self.ScratchDir)
        Command = "unzip %s"%TempZipPath
        print Command
        os.system(Command)
        self.UnzippedDBPath = os.path.join(self.ScratchDir, "%s.trie"%DBStub)
        if not os.path.exists(self.UnzippedDBPath):
            # No database?  That's bad!  Let's see whether it's in a subdirectory
            # with the same name as the stub:
            CheckPath2 = os.path.join(self.ScratchDir, DBStub, "%s.trie"%DBStub)
            if os.path.exists(CheckPath2):
                self.UnzippedDBPath = CheckPath2
            else:
                print "** Warning: Database file '%s' not found!"%self.UnzippedDBPath
        # now if desired copy over, and decompress the CommonDatabase
        if self.IncludeCommonDB:
            DBFileName = os.path.split(self.CommonDBPath)[1]
            DBStub = os.path.splitext(DBFileName)[0]        
            TempZipPath = os.path.join(self.ScratchDir,DBFileName)
            self.CopyFilePolitely(self.CommonDBPath, TempZipPath)
            os.chdir(self.ScratchDir)
            Command = "unzip %s"%TempZipPath
            print Command
            os.system(Command)
            self.UnzippedCommonDBPath = os.path.join(self.ScratchDir, "%s.trie"%DBStub)
            if not os.path.exists(self.UnzippedDBPath):
                # No database?  That's bad!  Let's see whether it's in a subdirectory
                # with the same name as the stub:
                CheckPath2 = os.path.join(self.ScratchDir, DBStub, "%s.trie"%DBStub)
                if os.path.exists(CheckPath2):
                    self.UnzippedCommonDBPath = CheckPath2
                else:
                    print "** Warning: Database file '%s' not found!"%self.UnzippedCommonDBPath            
            
        
    def Main(self):
        SpectrumFileName = self.SpectrumFileNames[0]
        self.Stub = os.path.splitext(SpectrumFileName)[0]
        if self.BlockNumber != None:
            self.Stub = "%s.%s"%(self.Stub, self.BlockNumber)
        self.UniqueStub = "%s.%s"%(self.Stub, os.getpid())
        OutputFileName = "%s.txt"%(self.Stub)
        # What node are we running on?  (Write to stdout)
        Command = "uname -a"
        os.system(Command)
        sys.stdout.flush()    
        # Build scratch directory:
        self.RunningFlagPath = self.FlagJob(0)
#        self.CopyFilesToNode()
        InspectInputPath = self.BuildInspectInputFile()
        OutputPath = os.path.join(self.gridEnv.ResultsXDir, OutputFileName)
        InspectStdout = os.path.join(self.ScratchDir, "%s.out"%self.UniqueStub)
        ##################################################################
        # Run inspect:
        Command = "%s/inspect -i %s -o %s -r %s > %s" % ( self.gridEnv.InspectDir,
                InspectInputPath, OutputPath, self.gridEnv.InspectDir, InspectStdout)
        print Command
        sys.stdout.flush()
        os.system(Command)
        #######################################
        # And, flag these XML files as done:
        try:
            os.remove(self.RunningFlagPath)
        except:
            pass
        self.FlagJob(1)
    def ParseCommandLine(self):
        """
        Parse command-line options:
        """
        (Options, Args) = getopt.getopt(sys.argv[1:], "d:x:bf:l:p:c:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-d":
                # -d database path
                if not os.path.exists(Value):
                    print "** Error: couldn't find database archive '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.DBPath = Value
            elif Option == "-b":
                # -b blind search flag
                self.BlindSearchFlag = 1
            elif Option == "-s":
                # -s splice-db search flag
                self.SpliceSearchFlag = 1
            elif Option == "-x":
                # Space-delimited list of XML filenames:
                for FileName in Value.split():
                    FileName = FileName.strip()
                    if len(FileName):
                        self.SpectrumFileNames.append(FileName)
            elif Option == "-f":
                # FIRST scan number:
                self.FirstScanNumber = int(Value)
                self.BlockNumber = self.FirstScanNumber / self.gridEnv.BLIND_BLOCK_SIZE
            elif Option == "-l":
                # LAST scan number:
                self.LastScanNumber = int(Value)
            elif Option == "-p":
                # Maximum PTM count:
                self.PTMLimit = int(Value)
            elif Option == "-c":
                if not os.path.exists(Value):
                    print "** Error: couldn't find database archive '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.IncludeCommonDB = 1
                self.CommonDBPath = Value
                
        if not(self.DBPath and self.SpectrumFileNames):
            print UsageInfo
            sys.exit(1) 

UsageInfo = """ClusterRunInspect.py <args>
-b [0|1] blind search flag
-d [dbPath] Path to database file to search
-x [mzXMLPath] Path to a mzXML file to search against
-p [int] Maximum PTM count
"""
        
if __name__ == "__main__":
    gridEnv = ClusterUtils.JCVIGridEnv()
    Runner = JobRunner(gridEnv)
    Runner.ParseCommandLine()
    try:
        Runner.Main()
    except:
        traceback.print_exc()
    Runner.CleanupScratchDir()