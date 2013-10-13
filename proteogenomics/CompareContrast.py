"""
Compares multiple groups of InspectResults.  It looks at both similarities and differences.
Caveats -----
1. Internally *ALL* invitro modifications are removed.  They are evaporated!  InVitro
modifications are defined in the file InVitroModifications.txt, included in the
Inspect package.  These common modifications are not of interest to me, and so
they are removed.  That means that the annotation "S+43AMMY" is the same as
"SAMMY" from this progam's point of view.
2. It is our experience that looking for two modifications per peptide yields
junky results.  Almost all annotations are unreliable, and their number
is so few that it is not worth the trouble.  Therefore, this program assumes
a max of ONE modification per peptide.  It will not deal with two, in fact
that would probably break it.
"""
import os
import sys
import getopt
import string
import struct
from Utils import *
Initialize()



UsageInfo = """CompareContrast.py
    Compares and contrasts groups of Inspect results.  Two reports are
    saved to the "Reports" folder in the current working directory.  First,
    <userStem>.all.txt is a report that compares the basic protein statictics
    between the groups (e.g. coverage, number of peptides, number of spectra)
    <userStem>.mods.txt is a report that compares sites of modification between
    the groups of results (e.g. the modification type, number of spectra with
    the modification, number of spectra total over the site)

Required options:
 -r [FileName] - The file explaining the results groupings - see code for example
 -d [TrieFileName] - The name of the database file searched. 
 -w [FileNameStem] - Generates the two results pages with the FileNameStem
 
Additional options:
 -p [Value] - Cutoff p-value.  Annotations with inferior p-value are ignored.
    Defaults to 0.05.
 -m [Value] Minimum number of unique peptides per protein.  Proteins with fewer
    assigned peptides will be removed. Default vallue is 2
    
"""

GroupingExample = """
This code is to explain what a file for the -r option needs to be like. It should have
lines with a *FULL PATH* to the results folders, one per line.  If I had two groups,
and wanted to compare them, a possible file could look like follows 
---------------
C:\CancerResults\Case
C:\CancerResults\Control

----------------
Caveats:
1. All Results in a single group should be in a single folder. The program
will only search that folder (no subfolders) for results files.

"""


class SummaryClass:
    def __init__(self):
        "constructor"
        self.PValueCutoff = 0.05
        self.SpectrumCount = 0
        self.MinimumPeptideLength = 6
        self.GroupingsFile = None
        self.NumGroupings = None
        self.ResultsDirs = []
        self.DatabasePath = None
        self.OutputFileStem = None
        self.MinimumProteinHits = 2 #need at least two peptides per protein to pass
        self.ReportsDir = "Reports"
        self.MakeReportFlag =0
        self.AllProteins = {} #key = Name, Value = [Protein Object0,1,2,3...]
            #ProteinClass Object 0 (above) is instantiated for all proteins in the database.
            #It has NO peptides.  It is merely a holder so that I can do stuff
        
    def ReadSearchResults(self):
        "Read search results from each of the specified directories"
        DirCounter=0
        for Dir in self.ResultsDirs:
            DirCounter +=1
            for FileName in os.listdir(Dir):
                (Stub, Extension) = os.path.splitext(FileName)
                if Extension.lower() == ".txt":
                    Path = os.path.join(Dir, FileName)
                    self.ReadSearchResultsFromFile(Path,DirCounter)
            
        print "File input complete.  Using %d spectra"%self.SpectrumCount
        #now we cut out the proteins who have too few peptides
        print "Paring the Protein database.  Keeping those with %d or more peptides"%self.MinimumProteinHits
        AllProteinNames = self.AllProteins.keys()
        for ProteinName in AllProteinNames:
            #1. Delete proteins which have a len =1.  This means that no others were added
            if(len(self.AllProteins[ProteinName] == 1):
               del self.AllProteins[ProteinName]
            else:   
                #Keep a protein if at least one of the cell types has enough
                KeepMe =0
                for ProteinObject in self.AllProteins[ProteinName]
                    if len (ProteinObject.Peptides) >= self.MinimumProteinHits:
                        KeepMe =1
                        break
                if not KeepMe:
                    del self.AllProteins[ProteinName]
        NumProteinsKept = len(self.AllProteins.keys())
        print "Keeping %d proteins"%NumProteinsKept

    def ReadSearchResultsFromFile(self, FileName, CellType):
        "Populate and self.ProteinNames from an Inspect output file"
        File = open(FileName, "r")
        
        OldSpectrum = None
        LineNumber = 0
        print "Read search results from %s..."%FileName
        for FileLine in File.xreadlines():
            Bits = FileLine.split("\t")
            LineNumber += 1
            if LineNumber%10000 == 0:
                print "Read %s line %s..."%(FileName, LineNumber)
            try:
                PValue = float(Bits[10])
                ProteinName = Bits[3]
                Spectrum = (Bits[0], Bits[1])
            except:
                continue
            if Spectrum == OldSpectrum:
                continue #only take a single annotation per spectrum
            if PValue > self.PValueCutoff:
                continue
            if ProteinName[:3] == "XXX":
                continue
            OldSpectrum = Spectrum
            Annotation = Bits[2]
            if Annotation[1] == ".":
                Annotation = Annotation[2:-2]
            Sequence = ""  # the pure AminoAcid sequence, without any annotations
            for I in range(len(Annotation)):
                if Annotation[I] in string.uppercase:
                    Sequence += Annotation[I]
            if len(Sequence) < self.MinimumPeptideLength:
                continue
            self.SpectrumCount += 1

            ## 1. See if the reported ProteinName has any instances yet, if not then make them
            if len(self.AllProteins[ProteinName]) ==1:
                for Group in self.ResultsDirs:
                    Empty = self.AllProteins[ProteinName][0]
                    Protein = ProteinClass(Empty.Sequence,Empty.DBStart,Empty.DBEnd,Group)
                    self.AllProteins[ProteinName].append(Protein)
            ## 2. Test to see if the peptide is already included
            ## 3. If not, then add it to the Inspect annotated "ProteinName"
            ###     as well as any other proteins that it has an exact match to
            FoundPeptide =0
            for Peptide in self.AllProteins[ProteinName][CellType].Peptides:
                if Peptide.IsMe(Sequence):
                    #the query sequence matches the peptide objec.
                    #add annotation and break out!
                    Peptide.AddAnnotation(Annotation)
                    FoundPeptide=1
                    break
            if not FoundPeptide:
                NumLocations =1
                Peptide = UniversalPeptide(Annotation, NumLocations)
                self.AllProteins[ProteinName][CellType].Peptides.append(Peptide)

    def FinalizeProteinStructures(self):
        """
        This method does some finishing of the protein objects, 
        """
        for ProteinName in self.AllProteins.keys():
            Protein = self.AllProteins[ProteinName][0]
            Protein.GenerateModList()
            Protein.GenerateSequenceCoverage()
            Protein = self.AllProteins[ProteinName][1]
            Protein.GenerateModList()
            Protein.GenerateSequenceCoverage()
    
    def ReadProteinDatabase(self):
        """
        This reads in the trie and index database files.  The trie file is a long string
        and the index file is a packed struct of (SourceFilePos, SquishedFilePos, ID)
        where source file is the fasta file, and squished file is the trie file.  This
        is the start position. 
        """
        print "Reading Database"
        self.IndexPath = os.path.splitext(self.DatabasePath)[0] + ".index"
        if not os.path.exists(self.IndexPath):
            print "* Error: Database index file '%s' not found!"%self.IndexPath
        # Read protein sequence, with I->L substitution:
        File = open(self.DatabasePath, "rb")
        self.DB = File.read()
        File.close()
        # Read protein names from the index file:
        File = open(self.IndexPath, "rb")
        self.ProteinNames = []
        self.ProteinStarts = []
        while (1):
            Data = File.read(92)
            if not Data:
                break
            Tuple = struct.unpack("<qi80s", Data)
            Name = Tuple[2]
            NullPos = Name.find(chr(0))
            if NullPos != -1:
                Name = Name[:NullPos]
            self.ProteinNames.append(Name)
            self.ProteinStarts.append(Tuple[1])
        NumProteins = len(self.ProteinNames)
        self.ProteinIndiciesDict = {}
        for I in range(NumProteins):
            Name = self.ProteinNames[I]
            Start = self.ProteinStarts[I]
            if Name[:3] == "XXX":
                continue #not making objects for fake entries
            if I == NumProteins-1: #last entry
                DBLen = len(self.DB)
                End = DBLen-1
            else:
                End = self.ProteinStarts[I+1]-1
            Sequence = self.DB[Start:End]
            EmptyProtein = ProteinClass(Sequence,Start,End,"Empty")
            self.AllProteins[Name] = []
            self.AllProteins[Name].append(EmptyProtein)
    
    
    def OpenDir(self):
        """
        This method creates the Reports dir.  It is called because the user wants
        some reports made.  We won't delete the directory if it exists, but that
        means that old files will be mixed with new.  It is your responsbility
        as the user to know which ones you care about
        """
        try:
            os.mkdir(self.ReportsDir) #make a directory for these files in the current directory
            print "Made directory \"Reports\" for all output"
        except:
            print "Reports directory already exists.  Using it for all reports"
            
        
    def GenerateModResultsReport(self):
        """
        Generate the most basic report comparing the different cell types
        
        GeneName    Position    Modification    ModSpectra(Group1)   TotalSpectra(Group1)   (Group2) ...
        ----------------------------------------------------------------------------------
        rpl5        62,K        Acetylation            15                25               
        """
        ModResultsPath = os.path.join(self.ReportsDir, "%s.mod.txt"%self.OutputFileStem)
        ModResultsFileHandle = open(ModResultsPath, "wb")
        print "Writing File %s"%ModResultsPath
        PrintString = "ProteinName\tPosition\tModification"
        for Group in self.ResultsDirs:
            PrintString += "\tModSpectra(%s)\tTotalSpectra(%s)"%(Group,Group)
        PrintString += "\n"
        ModResultsFileHandle.write(PrintString)
        for ProteinName in self.AllProteins.keys():
            PrintString = "%s"%ProteinName
            ### 1. get a huge key dictionary - a list of all keys for the
            ### positionModSpectraDict for this protein (in all groups)
            AllKeys = {}
            for Index in range (1, self.NumGroupings +1):
                ModDictionary = self.AllProteins[ProteinName][Index].PositionModSpectraDict
                for Key in ModDictionary.keys():
                    AllKeys[Key] = AllKeys.get(Key,0) +1
            ### 2. Now go through AllKeys and get the needed information
            for Key in AllKeys.keys():
                (ModName,Pos)= Key
                Residue = self.AllProteins[ProteinName][0].Sequence[Pos]
                PrintString += "\t%s,%d\t%s"%(Residue,Pos,ModName)
                for Index in range (1, self.NumGroupings +1):
                    ModSpectraCount = self.AllProteins[ProteinName][Index].PositionModSpectraDict[Key]
                    TotalSpectraCount = self.AllProteins[ProteinName][Index].SequenceCoverage[Pos]
                    PrintString += "\t%d\t%d"%(ModSpectraCount,TotalSpectraCount)
                PrintString += "\n"
                ModResultsFileHandle.write(PrintString)
        ModResultsFileHandle.close()

    def GenerateAllResultsReport(self):
        """
        Makes a report showing proteins found in the two samples
        
        ProteinName  MW      Coverage(Group1)  Peptides(Group1)   Spectra(Group1)       Coverage(Group2) ...   
        ----------------------------------------------------------------------------------------------------
          rps15     102434          .75             42                  256                 .45    
        """
        AllResultsPath = os.path.join(self.ReportsDir, "%s.all.txt"%self.OutputFileStem)
        AllResultsFileHandle = open(AllResultsPath, "wb")
        print "Writing File %s"%AllResultsPath
        PrintString = "ProteinName\tMW"
        for Group in self.ResultsDirs:
            PrintString += "\tCoverage(%s)\tPeptides(%s)\tSpectra(%s)"%(Group,Group,Group)
        PrintString += "\n"
        AllResultsFileHandle.write(PrintString)
        for ProteinName in self.AllProteins.keys():
            Sequence = self.AllProteins[ProteinName][0].Sequence
            try:
                Mass = GetMass(Sequence)
            except:
                Mass = 0
                print "Your sequence has a unknown amino acid character, Protein %s"%ProteinName
            PrintString = "%s\t%d"%(ProteinName,Mass)
            for Index in range (1,self.NumGroupings +1):
                Coverage = self.AllProteins[ProteinName][Index].Coverage
                PeptideList = self.AllProteins[ProteinName][Index].Peptides
                NumPeps = len(PeptideList)
                NumSpectra = 0
                for Peptide in PeptideList:
                    NumSpectra += Peptide.TotalSpectraCount
                PrintString += "\t%f\t%d\t%d"%(Coverage,NumPeps,NumSpectra)
            PrintString += "\n"
            AllResultsFileHandle.write(PrintString)
        AllResultsFileHandle.close()

    def ParseResultsGroupings(self):
        """
        Given a results groupings file (as explained above) this
        parses out the folders with contain groups.
        """
        FileHandle = open(self.GroupingsFile, "rb")
        for Line in FileHandle.xreadlines():
            Line = Line.rstrip()
            if os.path.isdir(Line):
                self.ResultsDirs.append(Line)
        self.NumGroupings = len(self.ResultsDirs)

if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(psyco not found - running in non-optimized mode)"
    Summarizer = SummaryClass()
    (Options, Args) = getopt.getopt(sys.argv[1:], "r:d:p:m:w:")
    OptionsSeen = {}
    for (Option, Value) in Options:
        OptionsSeen[Option] = 1
        if Option == "-r":
            # -r results file(s)
            if not os.path.exists(Value):
                print "** Error: couldn't find results file '%s'\n\n"%Value
                print UsageInfo
                sys.exit(1)
            Summarizer.GroupingsFile = Value
        elif Option == "-d":
            #-d the Database file
            if not os.path.exists(Value):
                print "** Error: couldn't find Database file '%s'\n\n"%Value
                print UsageInfo
                sys.exit(1)
            Summarizer.DatabasePath = Value            
        elif Option == "-m":
            # -m Minimum number of spectra for a new protein
            Summarizer.MinimumProteinHits = int(Value)
        elif Option == "-p":
            # -p p-value cutoff
            Summarizer.PValueCutoff = float(Value)
        elif Option == "-w":
            #-t tab-delimited outputfile
            Summarizer.OutputFileStem = Value
            Summarizer.MakeReportFlag=1
    # Error out, if we didn't see required options:
    if not OptionsSeen.has_key("-r") or not OptionsSeen.has_key("-w") or not OptionsSeen.has_key("-d"):
        print "*************************************************************"
        print "** Error: Please specify the required options\n"
        print UsageInfo
        sys.exit(1)
    Summarizer.ReadProteinDatabase()
    Summarizer.ParseResultsGroupings()
    Summarizer.ReadSearchResults()
    Summarizer.FinalizeProteinStructures()
    if Summarizer.MakeReportFlag:
        Summarizer.OpenDir()
    if Summarizer.OutputFileStem:
        Summarizer.GenerateAllResultsReport()
        Summarizer.GenerateModResultsReport()
        
    