UsageInfo = """ PeptidomeProteinReport.py
This program is designed to produce a table for peptidome
submissions, at least the protein portion of the submission
Required Options
 -r [FileName] Filename or directory of annotations
 -d [TrieFile] Database for the Search
 -p [float] PValue cutoff
 -w [FileName] Output filename
 -f            2 peptide per protein flag (default off)
 -t            fully tryptic flag (default off)
 
"""


import os
import getopt
import sys
import traceback
import ResultsParser
import SelectProteins
import ProteinStatistics
from Utils import *
Initialize()
import InspectResults

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
        self.AllPeptides = {}
        self.PeptideSources = {}
        self.DBPath = []  # array for possible multiple files
        self.ProteinObjects = {} # key = proteinID (like returned from find peptide location), value = object
        self.FilterTwoPeptideFlag = 0
        self.TrypticFlag = 0
        self.PValueCutoff = 0.005 #something relevant for inspect
        ResultsParser.ResultsParser.__init__(self)
        
        
    def Main(self):
        self.ProteinPicker = SelectProteins.ProteinSelector()
        self.ProteinPicker.LoadMultipleDB(self.DBPath)
        #self.ProcessResultsFiles(self.InputFile, self.ParseFile)
        self.ParseInspect(self.InputFile)
        #return
        self.MakeProteins()
        self.PrintProteins()              

    def ParseInspect(self, FilePath):
        """Here I parse out Inspect Results to get peptide annotations, 
        Putting them in a hash for safe keeping
        """
        FalseAminos = []
        SpectrumCount = 0
        inspectParser = InspectResults.Parser( FilePath )
        for result in inspectParser:
            try:
                Annotation = result.Annotation
                Peptide = GetPeptideFromModdedName(Annotation)
                Aminos = Peptide.Aminos
                PValue = result.PValue
                InspectMappedProtein = result.ProteinName
                FilePath = result.SpectrumFile
                Spectrum = result.ScanNumber
                (Path, File) = os.path.split(FilePath)
                DictValue = "%s:%s"%(File, Spectrum)

            except:
                traceback.print_exc()
                continue # SNAFU
            #here's some hacking that needs to be fixed.  I currently want to filter stuff by lfdr, but
            #that may not always be present
            if result.LFDR == None: #meaning that I don't have the column in question
                if PValue > self.PValueCutoff: #substitute pvalue for LFDR
                    continue
            else:
                LFDR = result.LFDR
                PValue = LFDR
                if LFDR > self.PValueCutoff:
                    continue
            #now we check for tryptic (if you want)
            if self.TrypticFlag:
                if not self.FullyTryptic(Annotation):
                    continue
            #everybody passed this line gets a cookie (you passed pvalue cutoff)
            SpectrumCount += 1
            #just a little damage control here.  We want to count the number of false positive peptides
            if InspectMappedProtein[:3] == "XXX":
                #this is a true negative.  let's count them
                if not Aminos in FalseAminos:
                    FalseAminos.append(Aminos)
                continue

            if not self.AllPeptides.has_key(Aminos):
                self.AllPeptides[Aminos] = PValue
                self.PeptideSources[Aminos] = []
                self.PeptideSources[Aminos].append(DictValue)
            else:
                self.PeptideSources[Aminos].append(DictValue)
                if PValue <  self.AllPeptides[Aminos]:
                    self.AllPeptides[Aminos] = PValue

        print "I got %s truedb peptides, and %s decoy peptides (%s spectra)"%(len(self.AllPeptides), len(FalseAminos), SpectrumCount)


    def FullyTryptic(self, Annotation):
        "expecting a full annotation, not just aminos.  so R.ABSCDR.F"
        if not Annotation[0] in ["R", "K"]:
            return 0
        if not Annotation[-3] in ["R", "K"]:
            return 0
        return 1
        
    def MakeProteins(self):
        """After reading in all of the files in our set, we fill in proteins
        and create the quicky objects that it needs.
        """
        for Aminos in self.PeptideSources.keys():
            Locations = self.ProteinPicker.FindPeptideLocations(Aminos)
            #for all of the locations, put data into the quicky protein object
            for (ProteinID, StartResidue) in Locations:
                # StartResidue is ZERO based, most people think 1 based
                # just remember that for printing out and stuff
                if not self.ProteinObjects.has_key(ProteinID):
                    ProteinName = self.ProteinPicker.ProteinNames[ProteinID]
                    self.ProteinObjects[ProteinID] = QuickyProteinObject(ProteinName)
                # now put in the peptides and spectra
                self.ProteinObjects[ProteinID].AddPeptide(Aminos, self.PeptideSources[Aminos])

    def PrintProteins(self):
        """put stuff in the format that peptidome wants
        """
        Handle = open(self.OutputFile, "wb")
        PrintHeader = "Protein\tPeptide\tSpectrum Files\n"
        Handle.write(PrintHeader)
        print "Processing results for %s proteins"%len(self.ProteinObjects)
        Count = 0
        for Protein in self.ProteinObjects.values():
            ### first line
            ## protein      peptide           file:spectrum, file:spectrum
            NumPeptides = len(Protein.Peptides)
            if (NumPeptides == 1) and self.FilterTwoPeptideFlag:
                continue
            First =1
            Count += 1
            for (Aminos, Spectra) in Protein.Peptides.items():
                if First:
                    Line = "%s\t%s\t%s\n"%(Protein.Name, Aminos, Spectra)
                    Handle.write(Line)
                    First =0
                else:
                    Line = " \t%s\t%s\n"%(Aminos, Spectra)
                    Handle.write(Line)
        print "%s proteins passed the input filters"%Count
        Handle.close()
                    

    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:d:w:ftp:")
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
            elif Option == "-t":
                self.TrypticFlag = 1
            elif Option == "-p":
                self.PValueCutoff = float(Value)
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
    
