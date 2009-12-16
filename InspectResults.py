"""
Constants and methods for parsing (Inspect) search results
"""
import os
import bz2
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
    TotalPRM = 7
    MedianPRM = 8
    FractionY = 9
    FraxtionB = 10
    Intensity = 11
    NTT = 12
    PValue = 13
    FScore = 14
    DeltaScoreAny = 15
    DeltaScore = 16
    ProteinID = 17
    DBPos = 18
    FileOffset = 19
    LFDR = 20
    
class Parser():
    def __init__(self, FilePath, wantedCols=[getattr(Columns,x) for x in dir(Columns) if not x[0] == '_'],
            MaxFilesToParse = None, QuietFlag = 0):
        self.Columns = Columns
        self.extensions = (".txt", ".filtered", ".res", ".csv", ".out")
        wantedCols.sort()
        self.wantedColumns = wantedCols
        self.maxFiles = MaxFilesToParse
        self.quiet = QuietFlag
        self.filePath = FilePath
        self.FileCount = 0

    def __iter__(self):
        """
        Function for applying a Callback function to one search-results file, or to every
        search-results file in a directory.
        """
        print "ResultsParser:%s" % self.filePath
        if os.path.isdir(self.filePath):
            FileNames = os.listdir(self.filePath)
            random.shuffle(FileNames)
        else:
            FileNames = [self.filePath]

        for FileNameIndex in range(len(FileNames)):
            FileName = FileNames[FileNameIndex]
            if not self.quiet:
                print "(%s/%s) %s"%(FileNameIndex, len(FileNames), FileName)

            fileHandle = None
            (Stub, Extension) = os.path.splitext(FileName)
            if Extension.lower() in self.extensions:
                fileHandle = open( FileName )
            elif Extension.lower() == '.bz2':
                (prefix, Extension) = os.path.splitext(Stub)
                if Extension.lower() in self.extensions:
                    fileHandle = bz2.BZ2File( FileName )
            else:
                continue

            self.FileCount += 1

            for line in fileHandle:
                if line[0] == '#':
                    continue
                cols = line.strip().split("\t")
                if len(cols) != Columns.LFDR:
                    raise Exception("Error in number of columns, got %d:\n%s\n" % (len(cols),cols))
                wanted = []
                for i in range(len(cols)):
                    if i in self.wantedColumns:
                        wanted.append( cols[i] )
                yield wanted

            fileHandle.close()

            # Don't parse every single file, that will take too long!
            if self.maxFiles != None and self.FileCount > self.maxFiles:
                raise StopIteration 
