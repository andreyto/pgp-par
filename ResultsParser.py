"""
Constants and methods for parsing (Inspect) search results
"""
import os
import random
class Columns:
    "Constants for which columns contain which data"
    SpectrumFile = 0
    ScanNumber = 1
    Annotation = 2
    ProteinName = 3
    Charge = 4
    MQScore = 5
    Length = 6
    NTT = 12
    PValue = 13
    FScore = 14
    DeltaScoreAny = 15
    DeltaScore = 16
    ProteinID = 17
    DBPos = 18
    FileOffset = 19
    LFDR = 20
    
class SpectrumOracleMixin:
    def __init__(self):
        self.SpectrumOracle = {}
    def FixSpectrumPath(self, Path):
        FileName = os.path.split(Path)[-1]
        Stub = os.path.splitext(FileName)[0]
        return self.SpectrumOracle.get(Stub, Path)
    def PopulateSpectrumOracle(self, RootDirectory):
        """
        Used when mzxml files are spread over multiple subdirectories.
        MZXMLOracle[Stub] = full path to the corresponding MZXML file
        Used with -M option (not with -s option)
        """
        if not RootDirectory or not os.path.exists(RootDirectory):
            return
        print "Populate oracle from %s..."%RootDirectory
        for SubFileName in os.listdir(RootDirectory):
            # Avoid expensive iteration through results directories:
            if SubFileName[:7] == "Results":
                continue
            SubFilePath = os.path.join(RootDirectory, SubFileName)
            if os.path.isdir(SubFilePath):
                self.PopulateSpectrumOracle(SubFilePath)
                continue
            (Stub, Extension) = os.path.splitext(SubFileName)
            Extension = Extension.lower()
            if Extension == ".mzxml":
                self.SpectrumOracle[Stub] = os.path.join(RootDirectory, SubFileName)
            elif Extension == ".mgf":
                self.SpectrumOracle[Stub] = os.path.join(RootDirectory, SubFileName)
            elif Extension == ".ms2":
                self.SpectrumOracle[Stub] = os.path.join(RootDirectory, SubFileName)
                
class ResultsParser:
    def __init__(self, *args, **kw):
        self.Columns = Columns
    def ProcessResultsFiles(self, FilePath, Callback, MaxFilesToParse = None, QuietFlag = 0):
        """
        Function for applying a Callback function to one search-reuslts file, or to every
        search-results file in a directory.
        """
        print "ResultsParser:%s"%FilePath
        FileCount = 0
        if os.path.isdir(FilePath):
            FileNames = os.listdir(FilePath)
            random.shuffle(FileNames)
            for FileNameIndex in range(len(FileNames)):
                FileName = FileNames[FileNameIndex]
                if not QuietFlag:
                    print "(%s/%s) %s"%(FileNameIndex, len(FileNames), FileName)
                (Stub, Extension) = os.path.splitext(FileName)
                if Extension.lower() not in (".txt", ".filtered", ".res", ".csv", ".out"):
                    continue
                FileCount += 1
                SubFilePath = os.path.join(FilePath, FileName)
                apply(Callback, (SubFilePath,))
                # Don't parse every single file, that will take too long!
                if MaxFilesToParse != None and FileCount > MaxFilesToParse:
                    break 
        else:
            apply(Callback, (FilePath,))
    
