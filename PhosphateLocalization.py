"""PhosphateLocalization.py

This script is a glorified wrapper for Label.py. It calls label and
calculates the PLS score for each spectral annotation in the input set.

1. read in input data.  If it is not in native inspect format, then we
send it to the GetByteOffset part, so that we can use it like Inspect.

2. Label all possible annotations of the string and get their peptide score (think binomial)

3. Find the difference between the top two scores, report. print.
The results are reported by appending two extra columns to the data from the input
file.  These correspond to the top annotation, and it's PLS.
"""

UsageInfo = """PhosphateLocalization.py
Calculates the Phosphate Localization Score (PLS) for each spectral
annotation in the input file.  Make sure to read the tutorial so
that you understand how to use it correctly.

Required Options:
 -r [FileName] File of formatted annotations.
 -m [Directory] Directory containing spectra files (not filename)
 -w [FileName] Output of this program

Additional Options:
 -d [Directory] Directory for the images and annotated peak lists
      created during the label process.  Default "LabelSpewage"
 -R Replace input annotation with better annotation.  By default the
     better annotation is simply put in a new column

"""

import os
import sys
import getopt
import ResultsParser
import GetByteOffset
import string
import Label

class DetectiveClass(ResultsParser.ResultsParser):
    def __init__(self):
        self.InputFilePath = None
        self.OutputFilePath = None
        self.LabeledAnnotationsDir = "LabelSpewage" # for Labeled output
        self.MZXMLDir = None
        self.InspectFormat = 0
        self.ParsedAllSpectra = 0
        self.ReplacementFlag = 0 #whether we replace the annotation with better ones, or simply print a new column
        self.ScanOffset = {} # potentially large dictionary for storing the byte offset of each spectrum
        self.OldInspectResults = {} #(file, scan) => (MQScore, Annotation)  #file name only, not! path
        self.PLSDict = {} # self.PLSDict[(SpectrumFile, Scan)] = (PLS, NewPeptide)
        ResultsParser.ResultsParser.__init__(self)
        
    def Main(self):
        "this can now work with single file inputs, or directories"
        ##1. If a directory is passed in, then create an output directory
        if os.path.isdir(self.InputFilePath):
            MakeDirectory(self.OutputFilePath)
        self.ProcessResultsFiles(self.InputFilePath, self.DoPLSForFile)
            
    def DoPLSForFile(self, InputFilePath):
        """This takes an input file, and does the whole round of computation
        for that file.  read input, label it, make output
        """
        #print "DoPLSForFile"
        AcceptableFormat = self.CheckInputFormat(InputFilePath)
        if not AcceptableFormat:
            # get out of here, and quick.
            return
        if not self.InspectFormat:
            #hacky, but this gets called on an entire directory of spectra, so no use doing it over again.
            if not self.ParsedAllSpectra:
                self.GetByteOffsetsForSpectra()
                self.ParsedAllSpectra = 1
        if os.path.isdir(self.InputFilePath):
            # we should make separate labeled spewage directories for each file
            (Path, FileWithExt) = os.path.split(InputFilePath)
            FileName = os.path.splitext(FileWithExt)[0]
            LabelSpewageDir = "%sFor%s"%(self.LabeledAnnotationsDir, FileName)
            print "Making LabelSpewageDir %s"%LabelSpewageDir
            #LabelSpewageDir = os.path.join(self.LabeledAnnotationsDir, FileName)
            OutputPath = os.path.join(self.OutputFilePath, FileWithExt)
        else:
            LabelSpewageDir = self.LabeledAnnotationsDir
            OutputPath = self.OutputFilePath
        MakeDirectory(LabelSpewageDir)
        self.LabelMe(InputFilePath, LabelSpewageDir)
        self.MakeOutput(InputFilePath, LabelSpewageDir, OutputPath)


    def MakeOutput(self, InputFilePath, LabelSpewageDir, OutputFilePath):
        """The results of Label.py have been put into a folder, and we now have to parse
        those and put them back into the file that people gave us.
        """
        
        ## get all the stuff from Label
        self.ProcessResultsFiles(LabelSpewageDir, self.ParseLabelSpewage)
        # start putting it into the output
        Handle = open(InputFilePath, "rb")
        OutHandle = open(OutputFilePath, "wb")
        for Line in Handle.xreadlines():
            if not Line.strip():
                continue
            if Line[0] == "#":
                #header
                if self.ReplacementFlag:
                    OutHandle.write("%s\tPLS\n"%Line.strip())
                else:
                    OutHandle.write("%s\tPLS\tBetterAnnotation\n"%Line.strip())
                continue
            Bits = Line.strip().split("\t")
            SpectrumFullPath = Bits[self.Columns.SpectrumFile]
            SpectrumFile = os.path.split(SpectrumFullPath)[1]
            Scan = Bits[self.Columns.ScanNumber]
            Annotation = Bits[self.Columns.Annotation]
            Tuple = (SpectrumFile, Scan)
            if not self.PLSDict.has_key(Tuple):
                print "NO KEY, %s, %s"%(SpectrumFile, Scan)
                continue
            (PLS, NewPeptideAnnotation, BetterMQScore) = self.PLSDict[(SpectrumFile, Scan)]
            #now write stuff out
            Bits.append("%s"%PLS)  # this actually happens regardless of the replacement flag
            if self.ReplacementFlag:
                # here we want to replace the old Annotation, and the MQScore if the
                # input format was Inspect
                if NewPeptideAnnotation:
                    Bits[self.Columns.Annotation] = NewPeptideAnnotation
                    if self.InspectFormat:
                        Bits[self.Columns.MQScore] = BetterMQScore
            else:
                Bits.append("%s"%NewPeptideAnnotation)
            String = "\t".join(Bits)
            OutHandle.write("%s\n"%String)
        OutHandle.close()

    def ParseLabelSpewage(self, FilePath):
        """In each file I am going to grep out
        filename, scan number, PLS, better peptide if such exists
        """
        
        ##in the filename are the scan number and mzxml filename
        if not FilePath[-3:] == "txt":
            return #skip png images
        (Path, FileName) = os.path.split(FilePath)
        Pos = FileName.find("mzXML") + 5
        SpectrumFile = FileName[:Pos]
        Dot = FileName.find(".", Pos+1)
        Scan = FileName[Pos+1:Dot] # string value, not int
        NewPeptide = None
        BetterMQScore = None
        Handle= open(FilePath, "rb")
        PLS = "N/A" #default, shoudl get overridden for every file
        for Line in Handle.xreadlines():
            Line = Line.strip()
            #hard coded magic
            if Line[:10] == "Phosphate ":
                #Phosphate Localization Score: 52.2
                Colon = Line.find(":")
                PLS = Line[Colon + 1:]
                #print Line
                #print "I parsed out %s"%PLS
            if Line[:7] == "WARNING":
                #parse out new peptide
                ToSplit = Line.replace("WARNING: Better annotation than input.", "")
                (BetterMQScore, NewPeptide) = ToSplit.split(",")
                NewPeptide = NewPeptide.strip()
            if Line[:2] == "b2":
                #this means we've started to get into the rest of the verbose output
                # and past what we care about
                break
        Handle.close()
        Tuple = (SpectrumFile, Scan)
        self.PLSDict[Tuple] = (PLS, NewPeptide, BetterMQScore)


    def GetByteOffsetsForSpectra(self):
        "Read mzXML from either a single file, or directory, creating the self.ScanOffset dictionary"
        Abacus = GetByteOffset.Abacus()
        if os.path.isdir(self.MZXMLDir):
            for FileName in os.listdir(self.MZXMLDir):
                (Stub, Extension) = os.path.splitext(FileName)
                if Extension.lower() == ".mzxml":
                    Path = os.path.join(self.MZXMLDir, FileName)
                    ScanOffsetSingleFile = Abacus.GetByteOffset(Path)
                    for (ScanNumber, ScanOffset) in ScanOffsetSingleFile.items():
                        self.ScanOffset[(FileName, ScanNumber)] = (Path, ScanOffset)
        else:
            ScanOffsetSingleFile = Abacus.GetByteOffset(self.MZXMLDir)
            FileName = os.path.split(self.MZXMLDir)[1]
            for (ScanNumber, ScanOffset) in ScanOffsetSingleFile.items():
                self.ScanOffset[(FileName, ScanNumber)] = (self.MZXMLDir, ScanOffset)
                #print "Storing value (%s,%s) with key (%s, %s)"%(self.MZXMLDir, ScanOffset, FileName, ScanNumber)

    def LabelMe(self, InputFilePath, LabelSpewageDir):
        Handle = open(InputFilePath, "rb")
        Dymo = Label.LabelClass()
        Count = 0
        GoodScoreCount = 0
        WrongChargeCount = 0
        ScoredWorseCount = 0
        for Line in Handle.xreadlines():
            if Line[0] == "#":
                continue
            if not Line.strip():
                continue
            Bits = list(Line.strip().split("\t"))
            #Charge = int (Bits[self.Columns.Charge])  I don't thin I need this anymore
            Count +=1
            Annotation = Bits[self.Columns.Annotation]
            #print "Annotation :%s:"%Annotation
            FileName = Bits[self.Columns.SpectrumFile]
            Scan = int(Bits[self.Columns.ScanNumber])
            if not self.InspectFormat:
                FileNameMinusPath = os.path.split(FileName)[1]
                (FullPathDummy, ByteOffset) = self.ScanOffset[(FileNameMinusPath, Scan)]
                #print (FullPathDummy, ByteOffset)
                #continue
            else:
                ByteOffset = int(Bits[self.Columns.FileOffset])
            (Path,File) = os.path.split(FileName)
            FileName = os.path.join(self.MZXMLDir, File)
            VerboseFileName = "%s.%s.%s.verbose.txt"%(File, Scan, Annotation[2:-2])
            ImageFileName = "%s.%s.%s.png"%(File, Scan, Annotation[2:-2])
            VerboseFilePath = os.path.join(LabelSpewageDir, VerboseFileName)
            ImageFilePath = os.path.join(LabelSpewageDir, ImageFileName)
            ## as we've got a single Dymo object, we must be passing in full args list
            ## -p to suppress the image popup, and -P for the PLS score
            Args = " -r %s -b %d -a %s -v %s -w %s -p -P"%(FileName, ByteOffset, Annotation, VerboseFilePath, ImageFilePath)
            ArgsList = Args.split()
            print "LabelMe for %s, scan %s"%(FileName, Scan)
            print "Args: %s"%Args
            Dymo.ParseCommandLine(ArgsList)
            Dymo.Main()

        Handle.close()
        
    def CheckInputFormat(self, FilePath):
        """This method serves to catch input files that are not in the
        proper Inspect format.  If this is the case, then we must convert the
        files to Inspect format.  This basically means that we put a byte offset at the
        end.
        Expected format. (tab delimited, 3 columns)
        Spectrum File          Spectrum Number (int)            Annotation (string, no! numbers!)
        """
        Handle = open (FilePath, "rb")
        ## 1. get the first line and see if it's already in Inspect Format
        Line = Handle.readline()
        try:
            Bits = Line.strip().split("\t")
        except:
            print "####################################################"
            print "%s in improper format. Please read tutorial."%FileName
            return 0
        if not len(Bits) < self.Columns.FileOffset:
            self.InspectFormat = 1
            return 1# in inspect format.  it's okay
        ## 2. Check to see if each line of the input file has the proper format
        Reject = 0
        for Line in Handle.xreadlines():
            if Line[0] == "#":
                continue
            try:
                Bits = Line.strip().split("\t")
            except:
                print "####################################################"
                print "%s in improper format. Please read tutorial."%FileName
                return 0
            #now check to see if column 1 is a number, and 2 is a string (with no brackets)
            try:
                SpectrumNumber = int(Bits[1])
            except:
                Reject = 1
                print "Second column must be a integer representing the spectrum number"
            Annotation = Bits[2]
            AcceptArray = string.ascii_letters
            AcceptArray += "."  #for delimiting the prefix/suffix
            AcceptArray += "*"  # for the beginning/end of a protein. should only be in prefix/suffix
            AcceptArray += string.digits
            for Index in range(len(Annotation)):
                if not Annotation[Index] in AcceptArray:
                    print "This annotation is in an improper format %s"%Annotation
                    Reject = 1
                    break
            if Reject:
                print "####################################################"
                print "There were formatting problems with %s"%FileName
                print "We cannot proceed.  Please read the tutorial."
                return 0
        print "Input file %s received in the correct format"%FileName
    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:w:m:d:R")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-r":
                # -r results file(s)
                if not os.path.exists(Value):
                    print "** Error: couldn't find results file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.InputFilePath = Value
            if Option == "-d":
                self.LabeledAnnotationsDir = Value
            if Option == "-w":
                self.OutputFilePath = Value
            if Option == "-m":
                self.MZXMLDir = Value
            if Option == "-R":
                #flag to replace the current annotation.  
                self.ReplacementFlag = 1
            
        if not OptionsSeen.has_key("-r") or not OptionsSeen.has_key("-m") or not OptionsSeen.has_key("-w"):
            print UsageInfo
            sys.exit(1)


def MakeDirectory(Dir):
    if os.path.exists(Dir):
        return 
    try:
        os.makedirs(Dir)
    except:
        raise
    


if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(psyco not found - running in non-optimized mode)"
    MacGyver = DetectiveClass()
    MacGyver.ParseCommandLine(sys.argv[1:])
    MacGyver.Main()