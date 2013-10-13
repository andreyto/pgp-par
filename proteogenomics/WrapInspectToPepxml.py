UsageInfo="""WrapInspectToPepxml.py
This script is a wrapper for the Pepxml conversion.  I can
write this faster than alter their script.

Required Options:
 -s [Directory] the top directory for a set of mzxml spectra 
 -r [Directory] Directory of Inspect results
 -w [Directory] Directory to dump all the pepxml files
 -m Flag for asking to make a mapping of the files
"""

import os
import sys
import getopt

TemplateString ="""spectra,%s
instrument,ESI-ION-TRAP
protease,Trypsin
mod,+57,C,fix
db,/data/BulkProteomics/YersiniaPestis/Databases/ForPepXMLConversion/NC_004088.6frame.RS.trie
db,/data/BulkProteomics/YersiniaPestis/Databases/ForPepXMLConversion/NC_004836.6frame.RS.trie
db,/data/BulkProteomics/YersiniaPestis/Databases/ForPepXMLConversion/NC_004837.6frame.RS.trie
db,/data/BulkProteomics/YersiniaPestis/Databases/ForPepXMLConversion/NC_004838.6frame.RS.trie
db,/data/BulkProteomics/YersiniaPestis/Databases/ForPepXMLConversion/Common.RS.trie
TagCount,25
PMTolerance,3.0
"""



class WrapperClass:
    def __init__ (self):
        self.ResultsDir = None
        self.SpectraDir = None
        self.OutDir = None
        self.MakeMappingFlag =0 #this is a flag for mapping between filenames like 123.msgf -> spectrum1.msgf
        self.FileNameMaps = {} #
        
    def Main(self):
        if self.MakeMappingFlag:
            self.MakeMapping()
        AllSpectraDirs = self.FindSubDirs(self.SpectraDir)
        self.WrapIt(AllSpectraDirs)
        
    def MakeMapping(self):
        """Parameters: non
        return: none
        Description: 0ur pipeline renames everything to be 1.txt and 2.txt for ease of
        working with the grid.  That means that the inspect (and other) results files no
        longer have the same name as their corresponding mzXML file.  So let's make that
        mapping
        """
        for ResultFile in os.listdir(self.ResultsDir):
            ResultFilePath = os.path.join(self.ResultsDir, ResultFile)
            Handle = open(ResultFilePath, "r") #like the msgf result file /path/to/1.msgf
            Data = Handle.readline() # the inspect header perhaps
            if Data[0] == "#":
                Data = Handle.readline()
            #some check here for empty results
            Bits = Data.strip().split("\t")
            FullPath = Bits[0]
            (Path, FullFileName) = os.path.split(FullPath)
            (FileName, Ext) = os.path.splitext(FullFileName)
            self.FileNameMaps[FileName] = ResultFilePath #want the full path
        
    def WrapIt(self, Directories):
        """Parameters: a list of directories
        Return: none
        Description: go through the directories, find mzxml files, or their gz. confirm
        the results file, make the parameters from the template and call the converter.        
        """
        for Directory in Directories:
            for SpectraFile in os.listdir(Directory):
                GZip = 0
                SpectraPath = os.path.join(Directory, SpectraFile)
                if SpectraPath[-3:] == ".gz":
                    GZip = 1 #so we can rezip it at the end
                    #command to gunzip it
                    
                    Command = "gunzip %s"%SpectraPath
                    print Command
                    os.system(Command)
                    SpectraPath = SpectraPath[:-3] #strip off the .gz now
                #now check to see if it's an mzxml or not
                (Path, Ext) = os.path.splitext(SpectraPath)
                if not Ext == ".mzXML":
                    continue
                ResultsPath = self.VerifyResultsFile(SpectraPath)
                if not ResultsPath:
                    print "WARNING: No results found for %s"%SpectraPath
                    continue
                ParamPath = self.GenerateParameterFile(SpectraPath)
                OutPath = self.MakeOutPath(SpectraPath)
                ## system call
                ConvertCommand = "python InspectToPepXML.py -i %s -o %s -p %s -m %s"%(ResultsPath, OutPath, ParamPath, Directory)
                print ConvertCommand
                os.system(ConvertCommand)
                if GZip:
                    RegzipCommand = "gzip %s"%SpectraPath
                    print RegzipCommand
                    os.system(RegzipCommand)
                print

    def MakeOutPath(self, ResultsPath):
        (Path, ResultsFile) = os.path.split(ResultsPath)
        (File, Ext) = os.path.splitext(ResultsFile)
        NewFile =  "%s.%s"%(File, "pepXML")
        NewPath = os.path.join(self.OutDir, NewFile)
        return NewPath

    def GenerateParameterFile(self, SpectraPath):
        String = TemplateString%SpectraPath
        Path = "./Temp.Inspect.in"
        Handle = open(Path, "wb")
        Handle.write(String)
        Handle.close()
        return Path

    def VerifyResultsFile(self, SpectraPath):
        """Parameters: path to a spectrum
        Return: Path to the corresponding results file
        Description: root around in the results dir to find the file
        """
        (Path, SpectrumFile) = os.path.split(SpectraPath)
        (File, Ext) = os.path.splitext(SpectrumFile)
        if self.MakeMappingFlag:
            #now we just look for this spectrumfile name in our mapping
            if not self.FileNameMaps.has_key(File):
                return None
            return self.FileNameMaps[File] #the full path where this thing lives
        
        Extension = "rescore.txt"
        
        ResultsFile = "%s.%s"%(File, Extension)
        ResultsPath = os.path.join(self.ResultsDir, ResultsFile)
        if not os.path.exists(ResultsPath):
            return None
        return ResultsPath

    def FindSubDirs(self, Directory):
        "cycle through and if it's a directory add to the list."
        ReturnDirs = []
        ReturnDirs.append(Directory) #temp hack.  We should figure this out
        print "trying abspath on %s"%Directory
        FullPath = os.path.abspath(Directory)
        for SubDir in os.listdir(FullPath):
            FullPathOfSubDir = os.path.join(FullPath, SubDir)
            if os.path.isdir(FullPathOfSubDir):
                ReturnDirs.append(FullPathOfSubDir)
        return ReturnDirs
    

    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:s:w:m")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-r":
                if not os.path.exists(Value):
                    print "** Error: couldn't find Inspect results directory '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.ResultsDir = Value
            if Option == "-s":
                # -r results file(s)
                if not os.path.exists(Value):
                    print "** Error: couldn't find spectra dir '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.SpectraDir = Value
            if Option == "-w":
                if not os.path.exists(Value):
                    print "** Error: couldn't find output dir '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.OutDir = Value
            if Option == "-m":
                self.MakeMappingFlag = 1
                
        if not OptionsSeen.has_key("-r") or not OptionsSeen.has_key("-s") or not OptionsSeen.has_key("-w")  :
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