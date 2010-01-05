UsageInfo="""VerifyResultsForSpectraFiles.py
This just checks that a zip file contains an (output) txt
file for every mzXML in a directory

Required Options:
 -z [zip file] Zipped Inspect Results
 -d [Directory] Directory containing the mzXML files
"""

import os
import os.path
import sys
import getopt

class WrapperClass:
    def __init__ (self):
        self.ZipPath = None
        self.SpectraDir = None
        
    def Main(self):
        UnzippedLocation = self.UnzipMyZip()
        self.CompareFiles(UnzippedLocation)
        #clean up zip
        RmCommand = "rm -r %s"%UnzippedLocation
        print "about to rm: %s"%RmCommand
        os.system(RmCommand)
        
    def CompareFiles(self, UnzippedLocation):
        ## just make sure that all the spectra have a match
        OutputFilesFullPath = os.listdir(UnzippedLocation) # a list
        OutputFileNames = {}
        for FullPath in OutputFilesFullPath:
            (Path, FileName) = os.path.split(FullPath) #/path/to/A123.txt or /path/to/A123.rescore.txt
            Bits = FileName.split(".") # KLUDGE, take first only assume no . in real filename (to get around the "rescore" moniker
            OutputFileNames[Bits[0]] = 1 # dummy value, I only care about existence
        AllSpectra = os.listdir(self.SpectraDir) #just gives me files, not full paths
        MissingCount =0 
        for FileName in AllSpectra:
            (Stub, Ext) = os.path.splitext(FileName)
            #print "%s->%s, %s"%(FileName, Stub, Ext)
            #now check my filename
            if not OutputFileNames.has_key(Stub):
                print "Missing results file for %s (%s)"%(Stub, FileName)
                MissingCount += 1
        print "There were %s missing files"%MissingCount
        
        
    def UnzipMyZip(self):
        ## unzip the file into a temp directory, in the current directory
        TempDirName = "Delme.UnzippedResults"
        CWD = os.getcwd()
        TempDirPath = os.path.join(CWD, TempDirName)
        if os.path.exists(TempDirPath):
            #crapola. we don't want to mess with stuff, because we're going to nuke this dir
            #after we finish.  ABORT. ABORT.
            print "The directory %s already exists, abort."%TempDirPath
            sys.exit(1)
        print "About to make %s"%TempDirPath
        os.mkdir(TempDirPath)
        (ZipCurrDir, ZipFileName) = os.path.split(self.ZipPath)
        CopiedZipPath = os.path.join(TempDirPath, ZipFileName)
        CopyCommand = "cp %s %s"%(self.ZipPath, CopiedZipPath)
        print "About to copy: %s"%CopyCommand
        os.system(CopyCommand)
        #now we unzip this file (should be in the directory)
        UnzipCommand = "unzip %s -d %s"%(CopiedZipPath, TempDirPath)
        print "About to unzip: %s"%UnzipCommand
        os.system(UnzipCommand)
        return TempDirPath

    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "z:d:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-z":
                # -r results file(s)
                if not os.path.exists(Value):
                    print "** Error: couldn't find zip file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.ZipPath = os.path.abspath(Value)
            if Option == "-d":
                # -r results file(s)
                if not os.path.exists(Value):
                    print "** Error: couldn't find spectra directory '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.SpectraDir = os.path.abspath(Value)
        if not OptionsSeen.has_key("-d")  or not OptionsSeen.has_key("-z"):
            print UsageInfo
            sys.exit(1)


if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(psyco not found - running in non-optimized mode)"
    DoStuff = WrapperClass()
    DoStuff.ParseCommandLine(sys.argv[1:])
    DoStuff.Main()        