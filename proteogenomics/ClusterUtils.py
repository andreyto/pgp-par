"""
This script contains some miscellaneous constants and utilities that are helpful for
running jobs on the grid.  It's imported ON the grid.
"""
import getpass
import os


class ClusterEnv:
    # Temp: wait only a short while, to shake bugs out of the process:
    ONE_HOUR = 60 # * 60

    # Blind block size is 1000 spectra now, not 10000:  
    BLIND_BLOCK_SIZE = 1000

    def __init__(self):
        self.USER_NAME = getpass.getuser()
        self.HomeDir = "/home/%s"%self.USER_NAME
        self.CopyFlagDir = "/scratch/%s/CopyFlags"%self.USER_NAME
        self.MZXMLDir = "/scratch/%s/mzxml"%self.USER_NAME
        self.DoneDir = "/scratch/%s/Done"%self.USER_NAME
        self.InspectDir = "/scratch/%s/Inspect"%self.USER_NAME
        self.GenomeZipPath = "/scratch/%s/IPIIPI.zip"%self.USER_NAME
        self.ResultsXDir = "/scratch/%s/ResultsX"%self.USER_NAME
        self.OutputDir = "/scratch/%s/output"%self.USER_NAME
        self.JobDir = "/scratch/%s/jobs"%self.USER_NAME
        self.ScratchDir = "/scratch/%s"%self.USER_NAME
        self.MAX_SIMULTANEOUS_COPIES = 5
        self.MAXIMUM_JOBS = 63

    def MakeDirectory(self,Dir):
        try:
            os.makedirs(Dir)
        except:
            pass

    def MakeGridDirectories(self):
        """
        Create some directories for input files, results, etc.
        """
        self.MakeDirectory(self.ScratchDir)

        SubDirNames = [ "jobs", "mzxml", "Done", "output", "ResultsX",
                "Databases/Genomic","Databases/Predictions","Databases/Proteomic"]
        for SubDirName in SubDirNames:
            SubDir = os.path.join(self.ScratchDir, SubDirName)
            self.MakeDirectory(SubDir)

class JCVIGridEnv(ClusterEnv):
    def __init__(self, ScratchDir="/usr/local/depot/projects/PGP/run/%s/" % getpass.getuser(),
                 projectName=None):
        ClusterEnv.__init__(self)
        self.ScratchDir = ScratchDir
        if projectName:
            self.ScratchDir = os.path.join( self.ScratchDir,
                os.path.basename( projectName.rstrip(os.path.sep) ) )

        self.CopyFlagDir   = None
        self.MZXMLDir      = os.path.join(self.ScratchDir , "mzxml" )
        self.DoneDir       = os.path.join(self.ScratchDir , "Done" )
        self.InspectDir    = None
        self.GenomeZipPath = os.path.join(self.ScratchDir , "IPIIPI.zip" )
        self.ResultsXDir   = os.path.join(self.ScratchDir , "ResultsX" )
        self.OutputDir     = os.path.join(self.ScratchDir , "output" )
        self.JobDir        = os.path.join(self.ScratchDir , "jobs" )
        self.InspectDBDir  = os.path.join(self.ScratchDir , "Databases","Genomic")
        self.Contaminants  = os.path.join(self.ScratchDir , "Databases","Proteomic")
        self.projectCode   = '0438'

    def GetRunningJobCount(self):
        TempFilePath = "TempJobCount.txt"
        Command = "qstat | grep -c %s > %s"%(self.USER_NAME, TempFilePath)
        #print Command
        os.system(Command)
        File = open(TempFilePath, "rb")
        FileLine = File.readline()
        File.close()
        return int(FileLine)
