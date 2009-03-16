"""
Simple script:
 - Read either NBCRFileNames.txt or FWGridFileNames.txt to get a list of search result files.
 - For each file that we don't have already, copy it over.
"""
import os
import sys
import getopt

UsageInfo = """
GridMirror.py:
-s [Server]: Nickname of server to copy files from
-x [String]: Required substring in filename
"""

class MirrorMaster:
    def __init__(self):
        self.RequiredString = None
    def ParseCommandLine(self, Arguments):
        (Options, Args) = getopt.getopt(Arguments, "s:x:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-s":
                self.SetServer(Value)
            elif Option == "-x":
                self.RequiredString = Value
        if not OptionsSeen.has_key("-s"):
            print UsageInfo
            sys.exit(-1)
    def SetServer(self, ServerName):
        ServerName = ServerName.lower()
        if ServerName == "nbcr" or ServerName == "kry":
            self.PSCPPath = "kry:/nas3/stanner/ResultsX/"
            self.ListingFileName = "NBCRFileNames.txt"
        elif ServerName == "fwgrid":
            self.PSCPPath = "grid:/scratch/stanner/ResultsX/"
            self.ListingFileName = "FWGridFileNames.txt"
        elif ServerName == "fwgridx":
            self.PSCPPath = "grid:/scratch/stanner/mzxml/"
            self.ListingFileName = "FWGridMZXMLNames.txt"
        else:
            raise ValueError, "*** Unknown server: %s"%ServerName
    def MirrorFiles(self):
        print "Copy listing-file %s from server..."%self.ListingFileName
        Command = "pscp %s%s ."%(self.PSCPPath, self.ListingFileName)
        print Command
        os.system(Command)
        if not os.path.exists(self.ListingFileName):
            print "(Not mirroring files from '%s')"%self.ListingFileName
            return
        # PresentFiles is a dictionary of (lower-cased) file names that we have
        # locally, and needn't copy:
        PresentFiles = {}
        PresentFiles[self.ListingFileName.lower()] = 1
        for FileName in os.listdir("."):
            Stub = os.path.splitext(FileName)[0]
            PresentFiles[FileName.lower()] = 1
        File = open(self.ListingFileName, "rb")
        for FileLine in File.xreadlines():
            FileName = FileLine.strip()
            # Ignore any files which do not contain the required string:
            if self.RequiredString != None and FileName.find(self.RequiredString) == -1:
                continue
            if PresentFiles.has_key(FileName.lower()):
                # Note that this file can be deleted from the grid:
                print "rm %s"%FileName
                continue
            Command = "pscp %s%s ."%(self.PSCPPath, FileName)
            print Command
            os.system(Command)
        File.close()

"""
pscp grid:/scratch/stanner/ResultsX/FWGridFileNames.txt .
pscp sam:/scratch/spayne/ResultsX/SamFileNames.txt .
pscp kry:/nas3/stanner/ResultsX/NBCRFileNames.txt .
pscp samk:/nas3/spayne/ResultsX/SamKFileNames.txt .

"""

if __name__ == "__main__":
    Master = MirrorMaster()
    Master.ParseCommandLine(sys.argv[1:])
    Master.MirrorFiles()
    
