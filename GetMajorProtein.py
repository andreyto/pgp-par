UsageInfo = """ GetMajorProtein.py
This program is designed for processing the results of 2D gels, in that
it looks for the main protein found in each mzXML file, assuming that 
each mzXML file corresponds to a single protein.
Required Options
 -r [FileName] Filename or directory of annotations
 -d [TrieFile] Database for the Search
 -w [FileName] Output filename
"""


import os
import getopt
import sys
import ResultsParser
import SelectProteins
from Utils import *
Initialize()

class AbacusClass(ResultsParser.ResultsParser):
    def __init__(self):
        self.InputFile = None
        self.OutputFile = None
        self.Header = None
        self.UniquePeptides = {} # Peptide => (MQScore, Bits)
        self.UniqueSites = {} #site in DB => (PLS, Bits)
        self.DBPath = []  # array for possible multiple files
        self.TrypsinCount = 0
        self.KeratinCount = 0

        ResultsParser.ResultsParser.__init__(self)
        
    def Main(self):
        self.ProteinPicker = SelectProteins.ProteinSelector()
        self.ProteinPicker.LoadMultipleDB(self.DBPath)
        self.OutputHandle = open(self.OutputFile, "wb")
        self.OutputHandle.write("File\tDominantProtein\tSpectra\tPercent\tContaminatingSpectra\tPercent\tTotalSpectra\n")
        self.ProcessResultsFiles(self.InputFile, self.ParseFile)
        #self.ProcessResultsFiles(self.InputFile, self.CountContaminants)
        self.OutputHandle.close()
        print "trypsin %s, keratin %s"%(self.TrypsinCount, self.KeratinCount)
        
    def CountContaminants(self, FileName):    
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
            ProteinName = Bits[self.Columns.ProteinName]
            if ProteinName.find("trypsin") > -1:
                self.TrypsinCount += 1
            if ProteinName.find("Keratin") > -1:
                self.KeratinCount += 1
    
    def ParseFile(self, FileName):
        """This is a wrapper for the Assess ProteinComponents, and goes in chunks of mzxml files
        This is meant to fix for the case that a single result file contains the results of several
        mzxml files
        """
        CurrentFile = None
        PeptidesSeen = {}
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
            if not CurrentFile == SpectrumFile:
                #new file found, so do something
                if not CurrentFile:
                    #this happens if CurrentFile == None
                    CurrentFile = SpectrumFile
                    #basically we do nothing here
                    continue
                else:
                    self.AssessProteinComponents(CurrentFile, PeptidesSeen)
                    #now start vars over
                    CurrentFile = SpectrumFile
                    PeptidesSeen.clear()
            Annotation = Bits[self.Columns.Annotation]
            Peptide = GetPeptideFromModdedName(Annotation)
            if not PeptidesSeen.has_key(Peptide.Aminos):
                PeptidesSeen[Peptide.Aminos] =0
            PeptidesSeen[Peptide.Aminos] +=1

    def AssessProteinComponents(self, FileName, PeptidesDict):
        """Given the results of a single mzXML file, we try and find the most common protein
        represented in those results.  We will also keep track of minor components.  I think that
        we will keep track of coverage of the major component, perhaps minor components also.
        I think that I will also keep track of hits to common contaminants separately.  perhaps
        reporting those as a QC for the experiment.
        """
        ProteinCoverage = {}  # cleaned out each function call.  key = proteinID, value = array for length of protein all zero, then + if covered
        ProteinSpectrumCount = {}
        LineCount = 0
        for (Aminos, Spectra) in PeptidesDict.items():
            #print "AssessProteinComponents, %s"%Aminos
            Locations = self.ProteinPicker.FindPeptideLocations(Aminos)
            #for all of the locations, put in the protein coverage, and spectrum count
            for (ProteinID, StartResidue) in Locations:
                # StartResidue is ZERO based, most people think 1 based
                # just remember that for printing out and stuff
                ## check for presence
                if not ProteinSpectrumCount.has_key(ProteinID):
                    ProteinSpectrumCount[ProteinID] = 0
                    Length = len(self.ProteinPicker.ProteinSequences[ProteinID])
                    ProteinCoverage[ProteinID] = [0,] * Length
                ## now put in data
                ProteinSpectrumCount[ProteinID] += Spectra
                #for Index in range (StartResidue, StartResidue + len(Aminos)):
                #    ProteinCoverage[ProteinID][Index] += Spectra
        # Now we're done with the lines from the file.  Let's assess the coverage
        DominantProteinBySpectrum = None
        MaxSpectra = 0
        TotalSpectra = 0
        CommonContaminantSpectra = 0
        for (ProteinID, SpectraCount) in ProteinSpectrumCount.items():
            TotalSpectra += SpectraCount
            ProteinName = self.ProteinPicker.ProteinNames[ProteinID]
            IsContaminant = self.CheckContaminants(ProteinName)
            if IsContaminant:
                CommonContaminantSpectra += SpectraCount
                continue
            #now for the rest of the stuff
            if SpectraCount > MaxSpectra:
                MaxSpectra = SpectraCount
                DominantProteinBySpectrum = ProteinID
        if DominantProteinBySpectrum:
            DominantProteinName = self.ProteinPicker.ProteinNames[DominantProteinBySpectrum]
        else:
            DominantProteinName = "None" #it found nothing!?!>!
        PercentInDominant = MaxSpectra/float(TotalSpectra)
        PercentInContaminants = CommonContaminantSpectra/float(TotalSpectra)
        OutLine = "%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(FileName, DominantProteinName, MaxSpectra, PercentInDominant, CommonContaminantSpectra, PercentInContaminants, TotalSpectra)
        self.OutputHandle.write(OutLine)
            

    def CheckContaminants(self, ProteinName):
        """For a protein name, check to see if its a common contaminant
        You can add things to the list here for customization
        """
        CommonContaminants = ["trypsin", "Keratin"]
        IsContaminant = 0
        for Contaminant in CommonContaminants:
            if ProteinName.find(Contaminant) > -1:
                IsContaminant = 1
                break
        return IsContaminant

    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:d:w:")
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