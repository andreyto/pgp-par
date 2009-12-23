"""
Constants and methods for parsing (Inspect) search results
"""
import os
import bz2
import random

class Row:
    "Object representing one row of inspect output."
    def __init__(self):
        self.SpectrumFile = None
        self.ScanNumber = None
        self.Annotation = None
        self.ProteinName = None
        self.Charge = None
        self.MQScore = None
        self.Length = None
        self.TotalPRM = None
        self.MedianPRM = None
        self.FractionY = None
        self.FraxtionB = None
        self.Intensity = None
        self.NTT = None
        self.PValue = None
        self.FScore = None
        self.DeltaScoreAny = None
        self.DeltaScore = None
        self.ProteinID = None
        self.DBPos = None
        self.FileOffset = None
        self.LFDR = None

    def __str__(self):
        outData = [ self.SpectrumFile,
        str(self.ScanNumber),
        self.Annotation,
        self.ProteinName,
        str(self.Charge),
        str(self.MQScore),
        str(self.Length),
        str(self.TotalPRM),
        str(self.MedianPRM),
        str(self.FractionY),
        str(self.FraxtionB),
        str(self.Intensity),
        str(self.NTT),
        str(self.PValue),
        str(self.FScore),
        str(self.DeltaScoreAny),
        str(self.DeltaScore),
        str(self.ProteinID),
        str(self.DBPos),
        str(self.FileOffset)]
        if self.LFDR != None:
            outData.append( str(self.LFDR) )
        return "\t".join(outData) + "\n"

    def populateFromString(self, inspectLine):
        cols = inspectLine.strip().split("\t")
        numCols = len(cols)
        if numCols not in [20,21]:
            raise Exception("Error in number of columns, got %d:\n%s\n" % (len(cols),cols))
        self.SpectrumFile = cols[0]
        self.ScanNumber   = int(cols[1])
        self.Annotation   = cols[2]
        self.ProteinName  = cols[3]
        self.Charge       = int(cols[4])
        self.MQScore      = float(cols[5])
        self.Length       = float(cols[6])
        self.TotalPRM     = float(cols[7])
        self.MedianPRM    = float(cols[8])
        self.FractionY    = float(cols[9])
        self.FraxtionB    = float(cols[10])
        self.Intensity    = float(cols[11])
        self.NTT          = float(cols[12])
        self.PValue       = float(cols[13])
        self.FScore       = float(cols[14])
        self.DeltaScoreAny= float(cols[15])
        self.DeltaScore   = float(cols[16])
        self.ProteinID    = int(cols[17])
        self.DBPos        = int(cols[18])
        self.FileOffset   = int(cols[19])
        if numCols == 21:
            self.LFDR     = float(cols[20])
    
class Parser():
    """Class to iterate through an inspect results file, or every
        search-results file in a directory. Supports bz2 compressed files,
        as well as a method to mirror the input directory structure into
        a new location on output.
    """
    def __init__(self, FilePath, MaxFilesToParse = None, QuietFlag = 0, inputMirrorTo = None):
        """Constructor arguments: FilePath is input file or dir of inspect results.
        MaxFilesToParse: sets a hard limit on how many files from a dir to read.
        QuietFlag: True for quiet False for verbose
        inputMirrorTo: an output directory to mirror the input dir structure to.
            obj.mirrorOutHandle will be the handle to write out to the current file.
        """
        self.extensions = (".txt", ".filtered", ".res", ".csv", ".out")
        self.maxFiles = MaxFilesToParse
        self.quiet = QuietFlag
        self.filePath = FilePath
        self.FileCount = 0
        self.header = ''
        self.currentFileName = None
        self.mirrorOutHandle = None
        self.outputMirrorInput = inputMirrorTo
        if inputMirrorTo and not os.path.exists( inputMirrorTo ):
            os.makedirs( inputMirrorTo )

    def __createHandles__( self, FileName ):
        """Internal method that creates the input and optional
        output file handles for the iterator.
        """
        fileHandle = None
        (dirName,baseName) = os.path.split( FileName )
        (Stub, Extension) = os.path.splitext(FileName)
        if Extension.lower() in self.extensions:
            fileHandle = open( FileName )
            if self.outputMirrorInput:
                self.mirrorOutHandle = open( os.path.join(
                    self.outputMirrorInput, baseName ), "w")

        elif Extension.lower() == '.bz2':
            Extension = os.path.splitext(Stub)[1]
            if Extension.lower() in self.extensions:
                fileHandle = bz2.BZ2File( FileName )
                if self.outputMirrorInput:
                    self.mirrorOutHandle = bz2.BZ2File( os.path.join(
                        self.outputMirrorInput, baseName ),"w")

        return fileHandle

    def __iter__(self):
        print "ResultsParser:%s" % self.filePath
        if os.path.isdir(self.filePath):
            FileNames = [os.path.join(self.filePath,x) for x in os.listdir(self.filePath)]
            random.shuffle(FileNames)
        else:
            FileNames = [self.filePath]

        for FileNameIndex in range(len(FileNames)):
            FileName = FileNames[FileNameIndex]
            if not self.quiet:
                print "(%s/%s) %s"%(FileNameIndex, len(FileNames), FileName)

            fileHandle = self.__createHandles__( FileName )
            if not fileHandle:
                continue

            self.FileCount += 1
            self.currentFileName = FileName

            for line in fileHandle:
                if line[0] == '#':
                    self.header = line
                    continue
                row = Row()
                row.populateFromString(line)
                yield row

            fileHandle.close()
            if self.outputMirrorInput:
                self.mirrorOutHandle.close()

            # Don't parse every single file, that will take too long!
            if self.maxFiles != None and self.FileCount > self.maxFiles:
                raise StopIteration 
