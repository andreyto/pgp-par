"""
This script contains some miscellaneous constants and utilities that are helpful for
running jobs on the grid.  It's imported ON the grid.
"""
import getpass
import os


# Temp: wait only a short while, to shake bugs out of the process:
ONE_HOUR = 60 # * 60

# Blind block size is 1000 spectra now, not 10000:  
BLIND_BLOCK_SIZE = 1000

USER_NAME = getpass.getuser()
if USER_NAME[:4] == "fwg.":
    USER_NAME = USER_NAME[4:]
    
if os.path.exists("/nas3/%s"%USER_NAME):
    System = "nbcr"
else:
    System = "fwgrid"
    
if System == "nbcr":
    HomeDir = "/home/%s"%USER_NAME
    CopyFlagDir = "/nas3/%s/CopyFlags"%USER_NAME
    MZXMLDir = "/nas3/%s/mzxml"%USER_NAME
    DoneDir = "/nas3/%s/Done"%USER_NAME
    InspectDir = "/nas3/%s/Inspect"%USER_NAME
    GenomeZipPath = "/nas3/%s/IPIIPI.zip"%USER_NAME
    ResultsXDir = "/nas3/%s/ResultsX"%USER_NAME
    OutputDir = "/nas3/%s/output"%USER_NAME
    JobDir = "/nas3/%s/jobs"%USER_NAME
    ScratchDir = "/nas3/%s"%USER_NAME
    MAX_SIMULTANEOUS_COPIES = 3
    MAXIMUM_JOBS = 63
else:
    HomeDir = "/home/%s"%USER_NAME
    CopyFlagDir = "/scratch/%s/CopyFlags"%USER_NAME
    MZXMLDir = "/scratch/%s/mzxml"%USER_NAME
    DoneDir = "/scratch/%s/Done"%USER_NAME
    InspectDir = "/scratch/%s/Inspect"%USER_NAME
    GenomeZipPath = "/scratch/%s/IPIIPI.zip"%USER_NAME
    ResultsXDir = "/scratch/%s/ResultsX"%USER_NAME
    OutputDir = "/scratch/%s/output"%USER_NAME
    JobDir = "/scratch/%s/jobs"%USER_NAME
    ScratchDir = "/scratch/%s"%USER_NAME
    MAX_SIMULTANEOUS_COPIES = 5
    MAXIMUM_JOBS = 63

class JobClass:
    """
    A Job represents a single search that should be launched.
    (NOT YET IMPLEMENTED)
    """
    def __init__(self):
        pass

def MakeDirectory(Dir):
    try:
        os.makedirs(Dir)
    except:
        pass

def MakeGridDirectories():
    """
    Create some directories for input files, results, etc.
    """
    Dir = os.path.join(ScratchDir, USER_NAME)
    MakeDirectory(Dir)
    SubDirNames = ["Inspect", "jobs", "mzxml", "Done", "output", "CopyFlags", "ResultsX"]
    for SubDirName in SubDirNames:
        SubDir = os.path.join(Dir, SubDirName)
        MakeDirectory(SubDir)

def GetRunningJobCount():
    TempFilePath = "TempJobCount.txt"
    Command = "qstat | grep -c %s > %s"%(USER_NAME, TempFilePath)
    #print Command
    os.system(Command)
    File = open(TempFilePath, "rb")
    FileLine = File.readline()
    File.close()
    return int(FileLine)
