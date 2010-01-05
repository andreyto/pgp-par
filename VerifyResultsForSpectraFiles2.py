UsageInfo="""VerifyResultsForSpectraFiles2.py
Given a guide file for gels we care about, find all raw files in the given 
folder and then checks that a zip file contains an (output) txt 
file for every raw

Required Options:
 -z [zip file directory] Directory containing zipped Inspect results
 -d [Directory] Directory containing the folders with raw files
 -g [file] Guide file with gel ids that we care about
 -t [temp dir] temporary directory to unzip stuff
"""

import os
import os.path
import sys
import getopt

class WrapperClass:
    def __init__ (self):
        self.ZipDir = None
        self.SpectraSuperDir = None
        self.GuideFile = None
        self.TempDir = None
        self.OutFile ="RenameYourOutput.txt"
        
    def Main(self):
        RelevantGelId = self.GetGelsFromGuide(self.GuideFile)
        #now I'm going to get all the gel cuts from this dir, so I only
        #have to mess with the directory and iterate through it once
        AllGelCutsInSuperDir = os.listdir(self.SpectraSuperDir)
        Handle = open (self.OutFile, "wb")
        for Gel in RelevantGelId:
            #we have zips in A712.rescore.zip which contain A712A and A712B
            # so let's unzip all these now
            ZipFileName = "%s.rescore.zip"%Gel
            ZipPath = os.path.join(self.ZipDir, ZipFileName)
            print "Looking at Gel %s, and unzipping %s"%(Gel, ZipPath)
            UnzippedResultsFiles = self.UnzipMyZip(ZipPath) #unzips into the temp dir
            CutsOfThisGel = self.MatchGelAndCuts(Gel, AllGelCutsInSuperDir)
            for Cut in CutsOfThisGel:
                FullPathToCut = os.path.join(self.SpectraSuperDir, Cut)
                RawFiles = os.listdir(FullPathToCut)
                self.FindMissingResultsFiles(Cut, RawFiles, UnzippedResultsFiles, Handle)
                
        
    def FindMissingResultsFiles(self, Cut, RawFiles, ResultsFiles, Handle):
        """Given a cut 'A712A' and all the raw files in it's directory
        and a bunch of results files, I compare them 
        compare files.  quicky check to see if num files matches
        """
        Handle.write("now comparing cut %s\n"%Cut)
        ResultsStubs = {} #without the .rescore.txt, really 
        MissingCount =0 
        #make a list of files I have, then iterate through the raws noting absences
        for File in ResultsFiles:
            if File.find(Cut) == -1: # this is to weed out all the A712B from our A712A 
                continue
            ReplaceText = ".rescore.txt" #change this to .txt if need arises
            ResultsStub = File.replace(ReplaceText, "")
            ResultsStubs[ResultsStub] = 1
        for File in RawFiles:
            if File.find(".RAW") == -1:
                continue #for .sld or other crap
            
            RawStub = File.replace( ".RAW", "") 
            if not ResultsStubs.has_key(RawStub):
                Handle.write("##Missing results file for %s\n"%(File))
                MissingCount += 1
        Handle.write("##There were %s missing files for %s\n"%(MissingCount,Cut))
        


        
    def MatchGelAndCuts(self, Gel, Cuts):
        """tying to match a Gel 'A712' to a cut 'A712A'
        """
        ReturnList = []
        for Cut in Cuts:
            if Cut.find(Gel) == 0:
                ReturnList.append(Cut)
        return ReturnList
        
    def GetGelsFromGuide(self, FilePath):
        """This gets lines like
        A712     A,B,C,D  othercrap
        All I care about is the A712 identifier.  After that it's your job
        to find all the relevant gel cuts.  this is just the *gel*.
        """
        Handle = open (FilePath, "rb")
        List = []
        for Line in Handle.xreadlines():
            if not Line.strip():
                continue #SNAFU
            Bits = Line.strip().split("\t")
            AGel = Bits[0]
            List.append(AGel)
        return List
        
    def UnzipMyZip(self, ZipPath):
        #now we unzip this file into a directory
        self.NukeAllTextFiles(self.TempDir)
        UnzipCommand = "unzip %s -d %s"%(ZipPath, self.TempDir)
        print "About to unzip: %s"%UnzipCommand
        os.system(UnzipCommand)
        UnzippedFiles = os.listdir(self.TempDir)
        return UnzippedFiles

    def NukeAllTextFiles(self, Directory):
        #1. nuke whatever is in this directory
        AllFiles = os.listdir(Directory)
        for File in AllFiles:
            if not File.find(".txt") == -1:
                FullPath = os.path.join(Directory, File)
                os.remove(FullPath)
        

    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "z:d:g:t:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-z":
                # -r results file(s)
                if not os.path.exists(Value):
                    print "** Error: couldn't find zip directory '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.ZipDir = os.path.abspath(Value)
            if Option == "-d":
                # -r results file(s)
                if not os.path.exists(Value):
                    print "** Error: couldn't find spectra directory '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.SpectraSuperDir = os.path.abspath(Value)
            if Option == "-t":
                # -r results file(s)
                if not os.path.exists(Value):
                    print "** Error: couldn't find temp directory '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.TempDir = os.path.abspath(Value)
            if Option == "-g":
                # -r results file(s)
                if not os.path.exists(Value):
                    print "** Error: couldn't find guide file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.GuideFile = os.path.abspath(Value)
        if not OptionsSeen.has_key("-d") or not OptionsSeen.has_key("-z") or not OptionsSeen.has_key("-t") or not OptionsSeen.has_key("-g"):
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