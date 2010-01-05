UsageInfo = """ PeptidomeProteinReport.py
This program is designed to produce a table for peptidome
submissions, at least the protein portion of the submission
Required Options
 -r [FileName] Filename or directory of annotations
 -d [TrieFile] Database for the Search
 -w [FileName] Output filename
 -f            2 peptide per protein flag (default off)
 
"""


import os
import getopt
import sys
import ResultsParser
import SelectProteins
import ProteinStatistics
from Utils import *
Initialize()


class QuickyProteinObject:
    def __init__(self, Name):
        self.Name = Name
        self.Peptides = {} #key = aminos value = [filename:spectrum, filename:spectrum, ..]

    def AddPeptide (self, Aminos, SpectraAndFiles):
        self.Peptides[Aminos] = SpectraAndFiles
        #since we have already filtered by unique peptides (see the parse file method),
        # then we will never need to do any checking


class AbacusClass(ResultsParser.ResultsParser):
    def __init__(self):
        self.InputFile = None
        self.OutputFile = "RenameYourOutput.txt"
        # i keep a list of the unique peptides, mostly because there is bound
        #to be redundancy, and if I only have to look for the location once, then
        # I save perhaps a lot of time.
        self.UniquePeptides = {} # key = aminos value = [filename:spectrum, filename:spectrum, ..]
        self.DBPath = []  # array for possible multiple files
        self.ProteinObjects = {} # key = proteinID (like returned from find peptide location), value = object
        self.FilterTwoPeptideFlag = 0
        ResultsParser.ResultsParser.__init__(self)
        
    def Main(self):
        self.ProteinPicker = SelectProteins.ProteinSelector()
        self.ProteinPicker.LoadMultipleDB(self.DBPath)
        self.ProcessResultsFiles(self.InputFile, self.ParseFile)
        self.MakeProteins()
        self.PrintProteins()              
                    
        
    def ParseFile(self, FileName):
        """This adds peptides to the list.
        """
        Handle = open(FileName, "rb")
        LineCount = 0
        for Line in Handle.xreadlines():
            LineCount += 1
            if LineCount % 1000 == 0:
                print "Parsing Line %s"%LineCount
            if not Line.strip():
                continue
            Line = Line.strip()
            if Line[0] == "#":
                continue
            Bits = Line.split("\t")
            SpectrumFile = Bits[self.Columns.SpectrumFile]
            ScanNumber = Bits[self.Columns.ScanNumber]
            Annotation = Bits[self.Columns.Annotation]
            File = self.CleanFileName(SpectrumFile)
            DictValue = "%s:%s"%(File, ScanNumber)
            Peptide = GetPeptideFromModdedName(Annotation)
            if not self.UniquePeptides.has_key(Peptide.Aminos):
                self.UniquePeptides[Peptide.Aminos] = []
            self.UniquePeptides[Peptide.Aminos].append(DictValue)
        Handle.close()

    def CleanFileName(self, ProbablyAPathAndName):
        (Path, Name) = os.path.split(ProbablyAPathAndName)
        return Name

    def MakeProteins(self):
        """After reading in all of the files in our set, we fill in proteins
        and create the quicky objects that it needs.
        """
        for Aminos in self.UniquePeptides.keys():
            Locations = self.ProteinPicker.FindPeptideLocations(Aminos)
            #for all of the locations, put data into the quicky protein object
            for (ProteinID, StartResidue) in Locations:
                # StartResidue is ZERO based, most people think 1 based
                # just remember that for printing out and stuff
                if not self.ProteinObjects.has_key(ProteinID):
                    ProteinName = self.ProteinPicker.ProteinNames[ProteinID]
                    self.ProteinObjects[ProteinID] = QuickyProteinObject(ProteinName)
                # now put in the peptides and spectra
                self.ProteinObjects[ProteinID].AddPeptide(Aminos, self.UniquePeptides[Aminos])

    def PrintProteins(self):
        """put stuff in the format that peptidome wants
        """
        Handle = open(self.OutputFile, "wb")
        PrintHeader = "Protein\tPeptide\tSpectrum Files\n"
        Handle.write(PrintHeader)
        for Protein in self.ProteinObjects.values():
            ### first line
            ## protein      peptide           file:spectrum, file:spectrum
            NumPeptides = len(Protein.Peptides)
            if (NumPeptides == 1) and self.FilterTwoPeptideFlag:
                continue
            First =1
            for (Aminos, Spectra) in Protein.Peptides.items():
                if First:
                    Line = "%s\t%s\t%s\n"%(Protein.Name, Aminos, Spectra)
                    Handle.write(Line)
                    First =0
                else:
                    Line = " \t%s\t%s\n"%(Aminos, Spectra)
                    Handle.write(Line)

        Handle.close()
                    

    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:d:w:f")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-r":
                # -r results file(s)
                if not os.path.exists(Value):
                    print "** Error: couldn't find results file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.InputFile = Value
            elif Option == "-d":
                self.DBPath.append(Value)
            elif Option == "-w":
                self.OutputFile = Value
            elif Option == "-f":
                self.FilterTwoPeptideFlag = 1
        if not OptionsSeen.has_key("-r") or not OptionsSeen.has_key("-d"):
            print UsageInfo
            sys.exit(1)
            

if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(psyco not found - running in non-optimized mode)"
    BeanCounter = AbacusClass()
    BeanCounter.ParseCommandLine(sys.argv[1:])
    BeanCounter.Main()                
    
