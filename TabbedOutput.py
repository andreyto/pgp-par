"""
This file is used to generate some tabbed delimited outputs from Inspect output.
Caveats -----
1. Internally *ALL* invitro modifications are removed prior to any parsing.
InVitro modifications are defined in the file InVitroModifications.txt.
Thus the annotation "S+43AMMY" is the same as "SAMMY" .
2. It is our experience that looking for two modifications per peptide yields
poor results.  When sorted by pvalue, almost all annotations are unreliable.
Therefore, this program assumes a max of ONE modification per peptide.
It will not deal with two, in fact that would probably break it.
"""
import os
import sys
import getopt
import string
import struct
import ResultsParser
import SelectProteins
from Utils import *
Initialize()



UsageInfo = """TabbedOutput.py - An alternative results representation. This
    script generates tab delimited output, saved to the "Reports" folder in
    the current working directory. Assumes a maximum of one PTM per peptide
    Example reports are:
    TotalResultsPage: Generates a basic report on every gene.
    GeneResultsPage: Generates a report on the a specific gene.
    ModificationTypePage: Generates a report of a specific type of
        modification, eg. phosphorylation

Required options:
 -r [FileName] - The name of the results file to parse.  If a directory is
    specified, then all .txt files within the directory will be combined into
    one report
 -d [TrieFileName] - The name of the database file searched.
 
Additional options:
 -t [FileName] - Generates the TotalResultsPage
 -g [FileName] - Generates a GeneResultsPage for every gene in the input
    file.  The input file lists genes, one per line, as a unique string from
    the database header.
 -o [FileName] - Generates a PTM Results Page, listing all the proteins with
    the specified PTM.  The input file has a list of PTMs, which are copied
    from the InVivoModifications.txt file.
 -p [Value] - Cutoff p-value.  Annotations with inferior p-value are ignored.
    Defaults to 0.05.
 -w Write out the annotations used to a file with the file extension .kept.txt
 -m [Value] Minimum number of unique peptides per protein.  Proteins with fewer
    assigned peptides will be removed. Default vallue is 2
    
Examples:
    TabbedOutput.py -r ResultsDir -d Database%sDicty.trie -t AllResults.txt
    TabbedOutput.py -r ResultsDir -d Database%sDicty.trie -g GeneList.txt
    
"""%(os.sep, os.sep)

#Global variables, because I want to be able to access them without having to pass them
InVivoMods = []
InVitroMods = []


class SummaryClass(ResultsParser.ResultsParser):
    def __init__(self):
        "constructor"
        self.PValueCutoff = 0.05
        self.SpectrumCount = 0
        self.MinimumPeptideLength = 6
        self.ResultsFileName = None
        self.DatabasePath = None
        self.TotalResultsFileName = None
        self.GeneInputFileName = None
        self.PTMInputFileName = None
        self.WriteKeptOutputFlag = 0 
        self.GenesForReport = [] # a list of geneNames for which to make the gene report
        self.PTMsForReport = [] # a list of PTMs for reporting
        self.MinimumProteinHits = 2 #need at least two peptides per protein to pass
        self.ReportsDir = "Reports"
        self.MakeReportFlag =0
        self.AllProteins = {} #key = Name, Value = Protein Object
        ResultsParser.ResultsParser.__init__(self)

    def Main(self):
        "after reading in the arguments, do stuff!"
        self.ReadModifications()
        self.ProteinSelector = SelectProteins.ProteinSelector()
        self.ProteinSelector.PValueCutoff = self.PValueCutoff
        print "LoadDB:"
        self.ProteinSelector.LoadDB(self.DatabasePath)
        self.ProteinSelector.RetainRepresentativeCount = 1 #so that we keep track
        self.ProcessResultsFiles(self.ResultsFileName, self.ProteinSelector.ParseAnnotations)
        self.ProteinSelector.ChooseProteins()
        self.CreateProteins()
        ## Reports
        if self.MakeReportFlag:
            MakeDirectory(self.ReportsDir) # in Utils.py
        if self.TotalResultsFileName:
            self.GenerateTotalResultsReport()
        if self.GeneInputFileName:
            self.ReadGeneFile()
            self.GenerateGeneReports()
        """
        if self.PTMInputFileName:
            self.ReadPTMInputFile()
            self.GeneratePTMReports()
        """
    def ReadGeneFile(self):
        """
        Read in a file denoting genes to be used in the gene report.  Gene names
        are one per line, and must be a unique string within the fasta lines of the
        database.
        """
        File = open (self.GeneInputFileName, "rb")
        GenesInput = []
        for Line in File.xreadlines():
            Line = Line.rstrip()
            GenesInput.append(Line)
        File.close()
        #check to see if the requested GeneIDs are unique
        AllProteinNames = []
        for ProteinID in self.AllProteins.keys():
            AllProteinNames.append(self.ProteinSelector.ProteinNames[ProteinID])
        ##
        for Gene in GenesInput:
            HitProtein = None
            for ProteinID in self.AllProteins.keys():
                ProteinName = self.ProteinSelector.ProteinNames[ProteinID]
                if ProteinName[:3] == "XXX":
                    continue #skiping all fake genes.  Override if you don't like it
                if ProteinName.find(Gene) >=0:
                    #found the desired gene name
                    if HitProtein: #if I've already made a match. Non-unique ID -> skip it!
                        HitProtein = "MultipleHits"    
                        break
                    HitProtein = (ProteinID, Gene) #Gene is the shortened unique version of proteinName
            #done with the for AllProteinnames loop
            if not HitProtein:
                print "Error: The gene ID %s was not found.  No output generated."%Gene
                continue
            if HitProtein == "MultipleHits":
                print "Error: The gene ID %s is not unique.  No output generated."%Gene
                continue
            self.GenesForReport.append(HitProtein)
        

    def CreateProteins(self):
        """
        This method creates protein objects, with the list of protein assignments from ChooseProteins
        """
        ## 1. create objects
        for (Annotation, RepList) in self.ProteinSelector.BestRepresentatives.items():
            if len(RepList) < 1:
                continue
            Peptide = RepList[0][1]
            ProteinID = self.ProteinSelector.PeptideProteins.get(Peptide.Aminos, None)
            ProteinName = self.ProteinSelector.ProteinNames[ProteinID]
            #print "Peptide %s belongs to protein %s\n\n"%(Peptide.GetFullModdedName(),ProteinName)
            if not self.AllProteins.has_key(ProteinID):
                "create it, need sequence"
                ProteinSequence = self.ProteinSelector.ProteinSequences[ProteinID]
                NewProtein = ProteinClass(ProteinSequence)
                self.AllProteins[ProteinID] = NewProtein
            SpectrumCounts = self.ProteinSelector.AnnotationSpectrumCounts[Peptide.GetFullModdedName()]
            self.AllProteins[ProteinID].AddAnnotation(Peptide,SpectrumCounts)
        ## 2. Process some characteristics
        for (ProteinID,Protein) in self.AllProteins.items():
            ProteinName = self.ProteinSelector.ProteinNames[ProteinID]
            Protein.GenerateSequenceCoverage()
            #print "sequence coverage for %s is %f"%(ProteinName,Protein.Coverage)
            Protein.GenerateModList()
            
    
    def GenerateTotalResultsReport(self):
        """
        Print off all results at the gene level.  This
        report lists all proteins with the following information: GeneName, Mass,
        Number of unique peptides, Total number of spectra, Number of modified
        peptides, number of spectra for the modified peptides.
        """
        TotalResultsPath = os.path.join(self.ReportsDir, self.TotalResultsFileName)
        TotalResultsFileHandle = open(TotalResultsPath, "wb")
        print "Writing File %s"%TotalResultsPath
        PrintString = "ProteinName\tMW\tCoverage\tUniquePeptides\tSpectraCount\tMod Peptides\tMod SpectraCount\n"
        TotalResultsFileHandle.write(PrintString)
        for Protein in self.AllProteins.keys():
            ProteinName = self.ProteinSelector.ProteinNames[Protein]
            Sequence = self.AllProteins[Protein].Sequence
            try:
                Mass = GetMass(Sequence)
            except:
                Mass = 0
                print "Your sequence has a unknown amino acid character, Protein %s"%Protein
            Coverage = self.AllProteins[Protein].Coverage
            PeptideList = self.AllProteins[Protein].Peptides # a list of UnmodifiedPeptides associated with this protein
            NumPeps = len(PeptideList)
            TotalSpectra =0
            ModifiedPeptideCount =0
            UnmodifiedSpectra =0
            for Peptide in PeptideList: #all UnmodifiedPeptide objects in the protein
                TotalSpectra += Peptide.TotalSpectraCount
                UnmodifiedSpectra += Peptide.UnmodifiedSpectraCount
                for Peptide in Peptide.Peptides: #all Peptide Objects in the UnmodifiedPeptide object
                    if len(Peptide.Modifications) > 0: #this peptide has a modification
                        ModifiedPeptideCount += 1
            ModifiedSpectraCount = TotalSpectra - UnmodifiedSpectra
            #finished tabulating.  Now lets print it to a file
            PrintString = "%s\t%d\t%f\t%d\t%d\t%d\t%d\n"%(ProteinName,Mass,Coverage,NumPeps,TotalSpectra,ModifiedPeptideCount,ModifiedSpectraCount)
            TotalResultsFileHandle.write(PrintString)
        TotalResultsFileHandle.close()
            
    def GenerateGeneReports(self):
        """
        This method goes through all the IDs in the GeneReportList, and finds
        the correct Protein, and then outputs the following information
        GeneName
        Peptide1Sequence:TotalSpectra:UnmodedSpectra:[modSequence:SpectraCount:Position]
        Peptide2 ...
        """
        for (ProteinID, ProteinName) in self.GenesForReport:
            FullSequence = self.AllProteins[ProteinID].Sequence
            GeneFileName = ProteinName + ".txt"
            GeneFilePath = os.path.join(self.ReportsDir, GeneFileName)
            print "Writing file %s for %s"%(GeneFilePath, ProteinName)
            GeneFileHandle = open(GeneFilePath, "wb")
            PrintString = "%s\nCoverage %f\n"%(ProteinName,self.AllProteins[ProteinID].Coverage)
            GeneFileHandle.write(PrintString)
            GeneFileHandle.write("Unmodified Sequence\tPeptidePosition\tTotalSpectra\tUnmodifiedSpectra\tMod Sequence\tResiduePositions\tSpectra\n")
            for UPeptide in self.AllProteins[ProteinID].Peptides: #peptide is a UnmodifiedPeptide objects
                ##Gather the Information about the Peptide
                PrintString = "%s\t"%UPeptide.Aminos
                PeptidePositionInProtein = FullSequence.find(UPeptide.Aminos)
                PrintString += "%d\t"%PeptidePositionInProtein
                PrintString += "%d\t%d\t"%(UPeptide.TotalSpectraCount,UPeptide.UnmodifiedSpectraCount)
                
                for Peptide in UPeptide.Peptides: #Peptide objects
                    if len(Peptide.Modifications) == 0:
                        continue #don't reprint the unmodified version
                    ## Determine the position of the modification in the sequence.
                    ModifiedSequence = Peptide.GetModdedName() #not always modified, but nonetheless
                    SpectraCount = UPeptide.SpectraCount[ModifiedSequence]
                    ModificationPosition = [] #list because there may be multiple
                    for (AminoIndex, ModificationList) in Peptide.Modifications.items():
                        ModificationPosition.append(PeptidePositionInProtein + AminoIndex)
                    PrintString += "%s\t%s\t%s\t"%(ModifiedSequence, SpectraCount, ModificationPosition)
                PrintString += "\n"
                GeneFileHandle.write(PrintString)
            GeneFileHandle.close()

    def ReadModifications(self):
        """
        This method reads in two files: InVivoModifications.txt and InVitroModifications.txt
        It makes a ModificationTypeObject out of each mod listed in the files
        (except fixed mods).  These input files are expected to be of the format
        mod,14,KR    TAB       #methylation.
        mod,DELTA_MASS,AMINOS,POSITION__TAB__#Modification name
        
        """
        InVivoHandle = open("InVivoModifications.txt","rb")
        for Line in InVivoHandle.xreadlines():
            Line = Line.rstrip()
            Data = Line.split("\t")
            Name = Data[1][1:] #should get rid of the '#'
            Latin = "InVivo"
            InspectInput = Data[0].rstrip() #get rid of any right side junk
            Data = InspectInput.split(",")
            DeltaMass = int (Data[1])
            Residues = Data[2]
            if len(Data) > 3:
                Position = Data[3]
            else:
                Position = None
            Mod = ModificationTypeObject(Latin,Name,DeltaMass,Residues,Position)
            InVivoMods.append(Mod)
        InVivoHandle.close()
        InVitroHandle = open("InVitroModifications.txt", "rb")
        for Line in InVitroHandle.xreadlines():
            if Line.find("fixed") >= 0:
                #this is a fixed modification, skip it
                continue
            Line = Line.rstrip()
            Data = Line.split("\t")
            Name = Data[1][1:] #should get rid of the '#'
            Latin = "InVivo"
            InspectInput = Data[0].rstrip() #get rid of any right side junk
            Data = InspectInput.split(",")
            DeltaMass = int (Data[1])
            Residues = Data[2]
            if len(Data) > 3:
                Position = Data[3]
            else:
                Position = None
            Mod = ModificationTypeObject(Latin,Name,DeltaMass,Residues,Position)
            InVitroMods.append(Mod)
        InVitroHandle.close()
        
    def ReadPTMInputFile(self):
        """
        Read in the PTMs that the user wants genes for
        """
        File = open (self.PTMInputFileName, "rb")
        for Line in File.xreadlines():
            Line = Line.rstrip()
            self.PTMsForReport.append(Line)
        File.close()
        
    def GeneratePTMReports(self):
        """
        Go through all the PTMs in the list, and print off genes for the PTM
        """
        AllProteins = self.AllProteins.keys()
        PrintIt=0
        for Mod in self.PTMsForReport:
            print "Scanning Annotaions for PTMs of type %s" %Mod
            PTMFileName = Mod + ".txt"
            PTMFilePath = os.path.join(self.ReportsDir, PTMFileName)
            PTMFileHandle = open(PTMFilePath, "wb")
            print "Results will be put into %s"%PTMFilePath
            PTMFileHandle.write("Genes and Annotations for %s\n"%Mod)
            PTMFileHandle.write("Gene\tAnnotation\tSpectrumCount\n")
            #mod is a string of a modification
            for Protein in AllProteins:
                #get the modification list, and see if any of them are it
                PeptideList = self.AllProteins[Protein].Peptides
                for Peptide in PeptideList:
                    #get the list of all modified versions of this peptide
                    ModList = Peptide.ModsList # list of ModifiedPeptide Objects
                    if PrintIt:
                        print Peptide.UnmodifiedSequence
                    for PepModObject in ModList:
                        (ModName,Position) = PepModObject.SpecificModification
                        if ModName == Mod:
                            #Print this Gene, and Peptide!!!
                            GeneName = Protein
                            ModifiedSequence = PepModObject.MSequence
                            PTMFileHandle.write("%s\t%s\t%d\n"%(GeneName,ModifiedSequence,PepModObject.SpectraCount))

    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:d:p:m:t:g:o:w")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-r":
                # -r results file(s)
                if not os.path.exists(Value):
                    print "** Error: couldn't find results file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                Summarizer.ResultsFileName = Value
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
            elif Option == "-w":
                # -w Write out annotations to fiel
                Summarizer.WriteKeptOutputFlag =1
                Summarizer.MakeReportFlag=1
            elif Option == "-p":
                # -p p-value cutoff
                Summarizer.PValueCutoff = float(Value)
            elif Option == "-t":
                #-t tab-delimited outputfile
                Summarizer.TotalResultsFileName = Value
                Summarizer.MakeReportFlag=1
            elif Option == "-g":
                #-g tab delmited gene files
                if not os.path.exists(Value):
                    print "** Error: couldn't find Gene Input file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                Summarizer.GeneInputFileName = Value
                Summarizer.MakeReportFlag=1
            elif Option == "-o":
                # -o PTM specific file
                if not os.path.exists(Value):
                    print "** Error: couldn't find PTM Input file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                Summarizer.PTMInputFileName = Value
                Summarizer.MakeReportFlag=1
        # Error out, if we didn't see required options:
        if not OptionsSeen.has_key("-r") or not OptionsSeen.has_key("-d"):
            print "*******************************************************************"
            print "** Error: Please specify results file (-r) and a database file (-d)\n"
            print UsageInfo
            sys.exit(1)

                        
if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(psyco not found - running in non-optimized mode)"
    Summarizer = SummaryClass()
    Summarizer.ParseCommandLine(sys.argv[1:])
    Summarizer.Main()
    sys.exit(1)
        
    