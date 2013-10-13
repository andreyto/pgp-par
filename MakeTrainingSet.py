"""MakeTrainingSet.py
This script takes in a set of annotation for some spectra.  It is assumed
that several programs have all annotated these spectra, and that their
annotations have been converted to Inspect format, and that the annotations
listed are quality annotations.  There will be no filtering in this program.
It is further assumed that the first part of the filename corresponds
to the program eg. Inspect.results.txt or Sequest.results.txt

"""
UsageInfo = """ MakeTrainingSet.py
Takes annotations from a variety of programs and compiles
them to make a training set.

Required Options
 -r [Directory] Directory where the results files are kept.
     Assumed that all files in this directory contain results.
 -w [OutputFile] Name of the output file
 -m [Directory] Directory for the mzxml files
Additional Options
 -P [number] Phosphorylation flag.  default 0
 -n [number] The number of programs that must annotate a scan in
         order to be included in the output. Default = 1, union
 -v [Directory] Directory for results of each Venn segment

"""

import sys
import os
import getopt
import traceback
import ResultsParser
from Utils import *
Initialize()

class CompileClass(ResultsParser.ResultsParser, ResultsParser.SpectrumOracleMixin):
    def __init__(self):
        self.ResultsDir = None
        self.OutputPath = None # the training set
        self.VenFilesDir = None #user requested
        ResultsParser.ResultsParser.__init__(self)
        ResultsParser.SpectrumOracleMixin.__init__(self)
        #hopefully all programs will annotate scans with the same peptide, but it's no guarentee
        self.AllAnnotations = {} # (filename, scan) => [(Program1, Annotation1), (program2,Annotation2), ..]
        self.ProgramNames = []
        self.PhosFlag = 0
        self.InspectMissedPath = "InspectMissed.txt"
        self.SpectraDir = None
        self.WitnessMinimum = 1
        
    def Main(self):
        if self.VenFilesDir:
            MakeDirectory(self.VenFilesDir)
        self.ProcessResultsFiles(self.ResultsDir, self.ParseFileCallback)
        self.WriteResults()

    def WriteResults(self):
        """I have decided that I am only going to deal with results
        that are confirmed by two programs.  that gives me a reasonable
        number.  All one hit wonders are routed to another function, and
        honestly forgotten.
        """
        
        "print out and also make Ven counts"
        Ven = {} #key = (Program1,Program2,) value = Count
        VenFileHandles = {} # key = (Program1,Program2,), value = file handle
        VenKeys = self.MakeVenKeys()
        for Key in VenKeys:
            Ven[tuple(Key)] = 0
            if self.VenFilesDir:
                KeyName = ".".join(Key) + ".txt"    
                VenFilePath = os.path.join(self.VenFilesDir, KeyName)
                VenFileHandles[tuple(Key)] = open(VenFilePath, "wb")
        Handle = open(self.OutputPath, "wb")
        InspectMissed = open(self.InspectMissedPath, "wb")
        for Key in self.AllAnnotations.keys():
            ValueList = self.AllAnnotations[Key]
            ## keep track of which programs contain this scan, for the ven diagram
            ScanPresence = [0]*len(self.ProgramNames)
            (FileName,Scan) = Key
            PrevAminos = None
            ProgramList = []
            Conflict = 0
            for (Program, Line) in ValueList:
                ProgramIndex = self.ProgramNames.index(Program)
                ScanPresence[ProgramIndex] = 1
                ProgramList.append(Program)
                ## keep track of all the annotations for this scan
                Bits = list(Line.split("\t"))
                try:
                    Annotation = Bits[self.Columns.Annotation]
                except:
                    continue # SNAFU
                ## check if it is a conflicting annotation
                if PrevAminos:
                    Peptide = GetPeptideFromModdedName(Annotation)
                    if not Peptide.Aminos == PrevAminos:
                        Conflict = 1
                        #print "conflict between %s and %s"%(PrevAminos, Peptide.Aminos)
                        #print Line
                else:
                    Peptide = GetPeptideFromModdedName(Annotation)
                    PrevAminos = Peptide.Aminos
            ## 1. exclude annotations that disagree between programs
            if Conflict: 
                self.ProcessConflict(ValueList)
            ## 2. Make Ven
            VenKey = []
            for Index in range(len(self.ProgramNames)):
                if ScanPresence[Index]:
                    VenKey.append(self.ProgramNames[Index])
            Ven[tuple(VenKey)] += 1
            ## 3. Write to appropiate Places
            if self.VenFilesDir:
                VenFileHandles[tuple(VenKey)].write(Line)
            if len(ValueList) >= self.WitnessMinimum:
                Handle.write(Line)
            if ProgramList.count("Inspect") == 0:
                #no inspect, put spectrum and annotation in separate file
                InspectMissed.write(Line)
            
        Handle.close()
        InspectMissed.close()
        #done with all the annotations.  Print out the ven stuff
        for Key in VenKeys:
            print "%s\t%d"%(Key,Ven[tuple(Key)])

    def ProcessConflict (self, ValueList):
        """I'm not sure at the moment how to process conflicting annotations.
        I do know that I don't want them included in the numbers or data
        for the training set.
        """
        #currently doing nothing.
        pass
              
    def MakeVenKeys(self):
        """makes all possible pairings for the program names
        Uses a bit string to calculate all possible combinations.
        for 3 programs, it calculates all numbers up to 2^3 as bit strings
        000  001  010  011  100  101  110  111
        The bit string is interpret to mean, 1 = include item in list at this index
        """
        Len = len(self.ProgramNames)
        N = pow(2,Len)
        AllComboList = []
        BitNumbers = []
        for I in range(Len):
            BitNumbers.append(pow(2,I))
        for Number in range(N):
            ## this number is perhaps 001100, so let's find which programs to include
            SingleComboList = []
            for Index in range(Len):
                BitNumber = BitNumbers[Index]
                if Number & BitNumber:
                    SingleComboList.append(self.ProgramNames[Index])
            #done creating the singleComboList for this number.
            #add it to the all combo list
            AllComboList.append(SingleComboList)
        #done with all the numbers.
        #for Item in AllComboList:
        #    print Item
        return AllComboList
    
    def ParseFileCallback(self, FilePath):
        """Called by the ResultsParser.ProcessResultsFiles
        Here I parse out lines and get annotations, Putting them in a list
        """
        FileName = os.path.split(FilePath)[-1]
        Program = FileName.split(".")[0]
        self.ProgramNames.append(Program)
        Handle = open(FilePath, "rb")
        Count = 0
        for Line in Handle.xreadlines():
            #print "Parsing line %s"%Line
            if Line[0] == "#":
                continue # comment line
            if not Line.strip():
                continue
            Bits = list(Line.split("\t"))
            try:
                Annotation = Bits[self.Columns.Annotation]
                SpectrumPath = self.FixSpectrumPath(Bits[self.Columns.SpectrumFile])
                SpectrumFile = os.path.split(SpectrumPath)[1]
                FixedPath = os.path.join(self.SpectraDir, SpectrumFile)
                #print "SpectrumFile, hopefully just the file %s"%SpectrumPath
                Spectrum = (SpectrumFile, Bits[self.Columns.ScanNumber])
            except:
                traceback.print_exc()
                continue # SNAFU
            #print Spectrum
            NumPhos = Annotation.count("phos")
            if NumPhos == 0:
                continue
            if NumPhos > self.PhosFlag:
                continue #exclude annotations with too many phosphorylation.
            Count +=1
            if not self.AllAnnotations.has_key(Spectrum):
                self.AllAnnotations[Spectrum] = []
            Bits[0] = FixedPath
            FixedLine = "\t".join(Bits) 
            #print Line
            #print FixedLine
            self.AllAnnotations[Spectrum].append((Program, FixedLine)) #keep the whole line so I can write it back out again
        print "%d lines parsed from %s"%(Count, FilePath)
        Handle.close()

    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:w:P:m:n:v:l:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-r":
                # -r results file(s)
                if not os.path.exists(Value):
                    print "** Error: couldn't find results file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.ResultsDir = Value
            if Option == "-w":
                self.OutputPath = Value
            if Option == "-P":
                self.PhosFlag = int (Value)
            if Option == "-m":
                self.SpectraDir = Value
            if Option == "-n":
                self.WitnessMinimum = int (Value)
            if Option == "-v":
                self.VenFilesDir = Value
        if not OptionsSeen.has_key("-r") or not OptionsSeen.has_key("-w") or not OptionsSeen.has_key("-m"):
            print UsageInfo
            sys.exit(1)


if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(psyco not found - running in non-optimized mode)"
    Gather = CompileClass()
    Gather.ParseCommandLine(sys.argv[1:])
    Gather.Main()