#!/usr/bin/env python

UsageInfo = """ProteogenomicsPostProcessing.py
Given a set of annotation for some spectra, a genome, and a proteome
we perform a variety of proteogenomic analyses. 
1. map peptides to the genome
2. cluster and filter
3. Find mis-predicted proteins
4. Cleavage Analysis

Required Parameters
 -r [FileName] File or directory of Peptide/Spectrum Match results.
 -g [FileName] GFF file of mapped peptides (alternate with -r) (single file)
 -o [Trie file] ORF Database used for search
 -d [Trie file] Proteome Database used for search
 -w [FileName]  Output file
 
Additional Parameters
 -i [int] Interpeptide distance cutoff for clustering (default 1000 bp)
 -u   Flag for requiring peptide uniqueness within the database.  Default
     behavior is to allow any peptide.
 -v   Flag to print a more verbose trace of the program progression
     through various methods and whatnot
 -p [float] Pvalue cutoff (Default 0.05)
 -M    Flag to skip looking for mis-predicted proteins
 -C    Flag to skip cleavage analysis
 -G [FileName] write out the peptide mappings in GFF format
 -W    Flag to print out warnings (as opposed to the verbose program progress)

"""

import sys
import os
import getopt
import traceback
import ResultsParser
import PeptideMapper
import GenomicLocations
import PGORFFilters
import PGCleavageAnalysis
import BasicStats
import GFFIO
from Utils import *
Initialize()


class FinderClass(ResultsParser.ResultsParser):
    def __init__(self):
        self.ReferenceResults = None
        self.MappedGFFResults = None #alternate input form, mostly used to shunt mapping, because that takes so long
        self.OutputPath = "RenameYourOutput.txt"
        self.ProteomeDatabasePaths = [] #possibly multiple
        self.ORFDatabasePaths = [] #possibly multiple
        self.AllPeptides = {} # AminoSequence -> best pvalue, nulled out in MapAllPeptides
        self.AllLocations = [] #list of PeptideMapper.GenomicLocationForPeptide
        self.AllPredictedProteins = {} #predictedProteinName ->GenomicLocationForPeptide Object, used in CreateORFs
        self.AllORFs = {} #ORF name -> GenomeLocationForORF object
        self.ProteomicallyObservedORFs = [] # this is populated when peptides are mapped, and deleted after ORF Objects are created
        self.UniquenessFlag = 0
        self.SpectrumCount = 0
        self.ORFPeptideMapper = PeptideMapper.PeptideMappingClass()
        self.InterPeptideDistanceMax = 1000 # sensible default
        self.PValueLimit = 0.05 #pvalues are 0.00 (good) to 1.0 (bad) 
        self.Verbose = 0
        self.VerboseWarnings = 0
        self.SearchForMispredictions = 1
        self.SearchForCleavage = 1
        self.OutputPeptidesToGFF = 0
        self.GFFOutputPath = "RenameYourOutput.gff"
        ResultsParser.ResultsParser.__init__(self)
        
    def Main(self):
        self.ORFPeptideMapper.LoadDatabases(self.ORFDatabasePaths)
        if self.ReferenceResults:
            self.ProcessResultsFiles(self.ReferenceResults, self.ParseInspectCallback)
            print "I found %s peptides from %s spectra"%(len(self.AllPeptides), self.SpectrumCount)
            self.MapAllPeptides()
            
        else :
            self.LoadResultsFromGFF(self.MappedGFFResults)
        ##now map all the peptides
        self.MapPredictedProteins()
        self.CreateORFs()
    
        ## Now start the analyses
        #self.FindOverlappingDubiousGeneCalls()
        self.FilterORFs()
        if self.OutputPeptidesToGFF:
            self.WritePeptideGFFFile()
        if self.SearchForMispredictions:
            self.FindMiscalls()
        if self.SearchForCleavage:
           self.AnalyzeCleavage()

    
    def FindOverlappingDubiousGeneCalls(self):
        """Parameters: none
        Return: none
        Description: There are genomic regions for which two gene calls overlap. 
        For some badly predicted genomes, the overlap is substantial (like >50 bp).
        I believe that most of these are bad gene calls, and I want to filter them 
        out.  This method calls the PRORFFilters method FindOverlappingDubiousGeneCalls
        which does that.  First we have to build the right dictionary,  so we do that here
        """
        ProteinDictionary = {} #proteinname->(start, stop)
        MaxOverlap = 50 #the base pairs
        for ORF in self.AllORFs.values():
            if not ORF.ProteinPrediction:
                continue #don't work with those lacking any protein.
            #print "%s"%ORF.ProteinPredictionName
            #ORF.ProteinPrediction.PrintMe()
            #print "\n"
            Name = ORF.ProteinPredictionName
            Start = ORF.ProteinPrediction.StartNucleotide 
            #'Start' is the small number. ALWAYS.  5' refers to something genic, but start is always just the small number
            Stop = ORF.ProteinPrediction.StopNucleotide
            #these coords don't yet take the stop codon into account, which NCBI does
            if ORF.Strand == "+":
                Stop += 3
            else:
                Start -= 3
            ProteinDictionary[Name] = (Start, Stop)
        OverlappingList = PGORFFilters.FindOverlappingDubiousGeneCalls(ProteinDictionary, MaxOverlap)
        #now go through and see whether we have peptide evidence for some overlappers
        for ORF in self.AllORFs.values():
            if not ORF.ProteinPrediction:
                continue #don't work with those lacking any protein.
            if len(ORF.PeptideLocationList) < 1:
                continue
            #ORF has both a protein and peptides, good.
            Name = ORF.ProteinPredictionName
            if Name in OverlappingList:
                print "%s is on the overlapping list, and has peptide representation"%Name


    def FilterORFs(self):
        """
        Parameters: None
        Return: None
        Description: glorified wrapper for calling the filter function. Afterwards
        we do some cleanup of ORFs that lack any peptides
        """
        if self.Verbose:
            print "ProteogenomicsPostProcessing.py:FilterORFs"
        for ORF in self.AllORFs.values():
            PGORFFilters.FilterThisORF(ORF)
        DeleteMeList = []
        for (ORFName, ORF) in self.AllORFs.items():
            if len(ORF.PeptideLocationList) == 0:
                DeleteMeList.append(ORFName)
        #can't delete from a list during the iteration.  so now we clean up
        for Name in DeleteMeList:
            del self.AllORFs[Name]
        print "\t%s ORFs left after filtering"%len(self.AllORFs)
        if self.Verbose:
            for (ORFName, ORF) in self.AllORFs.items():
                ORF.PrintMe(0,1)

    def AnalyzeCleavage(self):
        """
        Parameters: None
        Return: None
        Description: glorified wrapper for calling the filter function. 
        """
        Enzymes = ["Trypsin", ]
        if self.Verbose:
            print "ProteogenomicsPostProcessing.py:AnalyzeCleavage"
        for ORF in self.AllORFs.values():
            PGCleavageAnalysis.Analysis(ORF, Enzymes)



    def LoadResultsFromGFF(self, GFFFile):
        """Parameters: GFF file path
        Return: None
        Description: Take a gff file of peptides mapped to the genome, and use
        that to populate our results, as opposed to parsing Inspect results
        and mapping them (which takes time)
        """
        Handle = open (GFFFile, "rb")
        for Line in Handle.xreadlines():
            Dictionary = GFFIO.ParseGFFLine(Line)
            if not Dictionary:
                continue #in case they didn't return anything
            ##now here is some hacking, which should be improved.  Currently I'm only expecting
            ## to parse peptides, as opposed to ORFs or anything like that.  So I'm just
            ## blindly populating it to Peptides, but in the future, we should take care of this
            Location = GenomicLocations.GenomicLocationForPeptide()
            Location.FillFromGFF(Dictionary)
            MappedORF = Location.ProteinName
            if not MappedORF in self.ProteomicallyObservedORFs:
                self.ProteomicallyObservedORFs.append(MappedORF)
            self.AllLocations.append(Location)

    def WritePeptideGFFFile(self):
        """
        Parameters: None
        Return: None
        Description: This takes the peptide objects from self.AllLocations and
        makes a GFF File containing all of them.  
        """
        Handle = open(self.GFFOutputPath, "wb")
        GlobalLocationCount = 0
        for Location in self.AllLocations:
            Line = Location.GetGFF3LineNew(GlobalLocationCount)
            Handle.write(Line)#remember the Line has its own newline
            GlobalLocationCount += 1
        Handle.close()
        
    def FindMiscalls(self):
        """This function goes through the ORFs and calls their method for
        determining if they were mispredicted. just a wrapper really
        """
        if self.Verbose:
            print "ProteogenomicsPostProcessing.py:FindMiscalls"
        Count = 0
        TotalORFCount = len(self.AllORFs)
        for ORF in self.AllORFs.values():
            #if not ORF.ORFName == "Protein35275":
            #    continue
            ORF.AmIMispredicted()
            Count +=1
            if (Count %1000) == 0 and self.Verbose:
                print "Processed %s / %s ORF Objects"%(Count, TotalORFCount)
            

    def CreateORFs(self):
        """This function is to make ORF objects (GenomeLocationForORFs)
        We should put peptides and predictedProteins in here. 
        """
        DebugString = "Protein9325"
        if self.Verbose:
            print "ProteogenomicsPostProcessing.py:CreateORFs"
        Count = 0
        TotalORFCount = len(self.ORFPeptideMapper.ProteinPicker.ProteinNames.keys())
        for ID in self.ORFPeptideMapper.ProteinPicker.ProteinNames.keys():
            FastaLine = self.ORFPeptideMapper.ProteinPicker.ProteinNames[ID]
            #if FastaLine.find(DebugString) == -1:
            #    continue
            if not self.IsThisORFObserved(FastaLine):
                continue
            ORFSequence = self.ORFPeptideMapper.ProteinPicker.ProteinSequences[ID]
            ORFObject = GenomicLocations.GenomicLocationForORF(FastaLine, len(ORFSequence))
            self.AllORFs[ORFObject.ORFName] = ORFObject
            #now we go through and put all the protiens on the ORFS, then the peptides on the ORFs
            for (ProteinName, ProteinObject) in self.AllPredictedProteins.items():
                if ProteinObject.ProteinName == ORFObject.ORFName:
                    ## here it might be easy to confuse the ProteinName and ProteinObject.ProteinName
                    ## remember that PRoteinObject.ProteinName was set when a protein sequence was mapped onto an ORF
                    ## so really it should be called ORF Name, but that's they way things are
                    ## ProteinName, a key from self.AllPredictedProteins is actually the name of the predicted protein from the genome annotation
                    ORFObject.ProteinPrediction = ProteinObject
                    ORFObject.ProteinPredictionName = ProteinName
                    break #we're done here
            #now that we've finished running through the protein list, we run through the peptide list
            for Location in self.AllLocations:
                #print "trying to get locations with name %s"%Location.ProteinName
                if Location.ProteinName == ORFObject.ORFName:
                    ORFObject.PeptideLocationList.append(Location)
            #now done assigning both proteins and peptides to this ORF
            Count += 1
            if (Count %1000) == 0 and self.Verbose:
                print "\tCreated %s / %s ORF Objects"%(Count, TotalORFCount)
        print "\tCreated %s ORF Objects"%Count

    def IsThisORFObserved(self, FastaLine):
        """Parameters: fasta line of an ORF from the six frame translation
        Return: 0/1 
        Description: We want to answer whether this ORF has some peptides associated with
        it in this dataset
        """
        if FastaLine.find("XXX.") == 0:
            FastaLine = FastaLine.replace("XXX.", "XXX")
        InfoBits = FastaLine.split(".")
        #chromosome, strand
        ORFName = InfoBits[0]
        if ORFName in self.ProteomicallyObservedORFs:
            return 1
        return 0
                
    def MapAllPeptides(self):
        """Go through the list in self.AllPeptides, and find the location for each
        peptide. It is true that we could have parsed this out of the Inspect file, but
        some peptides have multiple locations, so we have to look it up. This sets up
        the self.AllLocations variable
        """
        if self.Verbose:
            print "ProteogenomicsPostProcessing.py:MapAllPeptides"
        Count = 0
        LocationCount=0
        for (Aminos, PValue) in self.AllPeptides.items():
            Count += 1
            if (Count %1000) == 0 and self.Verbose:
                print "Mapped %s / %s peptides"%(Count, len(self.AllPeptides))
            #print Aminos
            GenomicLocations = self.ORFPeptideMapper.MapMe(Aminos, PValue)
            if self.UniquenessFlag and (len(GenomicLocations) > 1):
                continue #skip out on adding it to the list
            self.AllLocations.extend(GenomicLocations)
            LocationCount += len(GenomicLocations)
            #now put into the observedORF stuff
            for Location in GenomicLocations:
                MappedORF = Location.ProteinName
                if not MappedORF in self.ProteomicallyObservedORFs:
                    self.ProteomicallyObservedORFs.append(MappedORF)

        self.AllPeptides = {}  # set to null just for the memory savings
        if self.Verbose:
            print "ProteogenomicsPostProcessing:MapAllPeptides - mapped %s peptides to %s genomic locations"%(Count, LocationCount)

    def MapPredictedProteins(self):
        """This will look very much like MapAllPeptides, because really a predicted protein can fit into an 
        ORF just like a MS/MS peptide. So we do the same kind of thing. 
        NOTE: because of the alternate start sites, we have a bit of a problem, in that most proteins
        don't map to the six frame translation, specifically the initial M.  See the following example
            RDVLNRVMYYIILARFINYRLISLSCRSKRMRIFQGVVCGMALFLA  (six frame translation)
                  MMYYIILARFINYRLISLSCRSKRMRIFQGVVCGMALFLA  (protein sequence)
       to remedy this, we are going to axe off the initial M, and then kludge back
       the start site of the mapping
        
        """
        if self.Verbose:
            print "ProteogenomicsPostProcessing::MapPredictedProteins"
        PredictedProteinReader = PeptideMapper.PeptideMappingClass() #just to use the LoadDB parts.  
        # I could do that here, but then I would be importing SelectProteins in this also, and why
        # should I do that when i've already written the code elsewhere
        PredictedProteinReader.LoadDatabases(self.ProteomeDatabasePaths)
        for ProteinID in PredictedProteinReader.ProteinPicker.ProteinNames.keys():
            ProteinName = PredictedProteinReader.ProteinPicker.ProteinNames[ProteinID]
            ProteinSequence = PredictedProteinReader.ProteinPicker.ProteinSequences[ProteinID]
            #This is where it gets a bit tricky. see note in method comment for explanation of the weirdness
            SearchableSequence = ProteinSequence[1:]
            ProteinLocation = self.ORFPeptideMapper.MapMe(SearchableSequence, 1)
            if len(ProteinLocation) == 0:
                if self.VerboseWarnings:
                    print "WARNING: protein %s does not map to any ORF"%ProteinName
                continue #always
            if len(ProteinLocation) > 1 and self.VerboseWarnings:
                print "WARNING: Protein %s maps to multiple ORFs"%ProteinName
            for Location in ProteinLocation: #potentially many
                Location.AddOneAminoAcidFivePrime()
                #Location.Aminos = ProteinSequence# replace the full sequence,because it got snipped in the beginning
            
            self.AllPredictedProteins[ProteinName] = ProteinLocation[0] #explicit assumption that there is only ONE
        
    def ParseInspectCallback(self, FilePath):
        """Called by the ResultsParser.ProcessResultsFiles
        Here I parse out Inspect Results to get peptide annotations, 
        Putting them in a hash for safe keeping
        """
        Handle = open(FilePath, "rb")
        for Line in Handle.xreadlines():
            #print "Parsing line %s"%Line
            if Line[0] == "#":
                continue # comment line
            if not Line.strip():
                continue
            Bits = list(Line.split("\t"))
            try:
                Annotation = Bits[self.Columns.Annotation]
                Peptide = GetPeptideFromModdedName(Annotation)
                Aminos = Peptide.Aminos
                PValue = float(Bits[self.Columns.PValue])
            except:
                traceback.print_exc()
                continue # SNAFU
            #here's some hacking that needs to be fixed.  I currently want to filter stuff by lfdr, but
            #that may not always be present
            if len(Bits) < self.Columns.LFDR: #meaning that I don't have the column in question
                if PValue > self.PValueLimit: #substitute pvalue for LFDR
                    continue
            else:
                LFDR = float(Bits[self.Columns.LFDR])
                PValue = LFDR
                if LFDR > self.PValueLimit:
                    continue
            if not self.AllPeptides.has_key(Aminos):
                self.AllPeptides[Aminos] = PValue
            else:
                if PValue <  self.AllPeptides[Aminos]:
                    self.AllPeptides[Aminos] = PValue
            self.SpectrumCount += 1
        Handle.close()

    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:g:d:w:uvi:o:p:CMG:W")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-r":
                # -r results file(s)
                if not os.path.exists(Value):
                    print "** Error: couldn't find results file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.ReferenceResults = Value
            if Option == "-g":
                if not os.path.exists(Value):
                    print "** Error: couldn't find results file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.MappedGFFResults = Value
            if Option == "-d":
                if not os.path.exists(Value):
                    print "** Error: couldn't find database file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.ProteomeDatabasePaths.append( Value)
            if Option == "-o":
                if not os.path.exists(Value):
                    print "** Error: couldn't find database file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.ORFDatabasePaths.append( Value)
            if Option == "-w":
                self.OutputPath = Value
            if Option == "-u":
                self.UniquenessFlag = 1
            if Option == "-i":
                self.InterPeptideDistanceMax = int(Value)
            if Option == "-p":
                self.PValueLimit = float(Value)
            if Option == "-v":
                self.Verbose = 1
            if Option == "-M":
                self.SearchForMispredictions = 0
            if Option == "-C":
                self.SearchForCleavage = 0
            if Option == "-G":
                self.OutputPeptidesToGFF = 1
                self.GFFOutputPath = Value
            if Option == "-W":
                self.VerboseWarnings =1
        if not OptionsSeen.has_key("-w") or not OptionsSeen.has_key("-d") or not OptionsSeen.has_key("-o"):
            print UsageInfo
            sys.exit(1)
        if not OptionsSeen.has_key("-r")  and not OptionsSeen.has_key("-g"):
            print UsageInfo
            sys.exit(1)


if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(psyco not found - running in non-optimized mode)"
    Gumshoe = FinderClass()
    Gumshoe.ParseCommandLine(sys.argv[1:])
    Gumshoe.Main()
