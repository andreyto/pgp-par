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

        SubDirNames = ["Inspect", "jobs", "mzxml", "Done", "output", "CopyFlags", "ResultsX"]
        for SubDirName in SubDirNames:
            SubDir = os.path.join(self.ScratchDir, SubDirName)
            self.MakeDirectory(SubDir)

class JCVIGridEnv(ClusterEnv):
    def __init__(self, projectName=None):
        ClusterEnv.__init__(self)
        self.ScratchDir = "/usr/local/projects/PGP/run/%s/" % self.USER_NAME
#        self.ScratchDir = "/tmp/t/%s/" % self.USER_NAME
        if projectName:
            self.ScratchDir = os.path.join( self.ScratchDir, os.path.basename(projectName))

        self.CopyFlagDir   = self.ScratchDir + "CopyFlags"
        self.MZXMLDir      = self.ScratchDir + "mzxml"
        self.DoneDir       = self.ScratchDir + "Done"
        self.InspectDir    = self.ScratchDir + "Inspect"
        self.GenomeZipPath = self.ScratchDir + "IPIIPI.zip"
        self.ResultsXDir   = self.ScratchDir + "ResultsX"
        self.OutputDir     = self.ScratchDir + "output"
        self.JobDir        = self.ScratchDir + "jobs"
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
