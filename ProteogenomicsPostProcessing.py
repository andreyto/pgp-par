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
 -b [FileName] The genbank file for the genome includes proteins and nucleotides
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


import getopt
import traceback
import InspectResults
import GenomicLocations
import PeptideMapper
import PGPeptide
import PGORFFilters
import PGPrimaryStructure
import SelectProteins
import PGCleavageAnalysis
import BasicStats
import GFFIO
from Utils import *
Initialize()


class FinderClass():
    def __init__(self):
        self.InspectResultsPath = None #Peptide/Spectrum Matches from Inspect
        self.GFFInputPath = None #alternate input form, mostly used to shunt mapping, because that takes so long
        self.GenbankPath  = None
        self.OutputPath = "RenameYourOutput.txt" #really a stub for output, as it gets modified to *.novel.faa, *.novel.info, etc.
        self.ProteomeDatabasePaths = [] #possibly multiple
        self.ORFDatabasePaths = [] #possibly multiple
        self.AllPeptides = {} # AminoSequence -> best pvalue, nulled out in MapAllPeptides
        self.AllLocatedPeptides = [] #list of PeptideMapper.GenomicLocationForPeptide
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

    def Main(self):
        self.ORFPeptideMapper.LoadDatabases(self.ORFDatabasePaths) #a handle for the 6frame translations db (called ORF)
        chromReader = PGPeptide.GenbankGenomeReader(self.GenbankPath, self.ORFDatabasePaths)
        genome = chromReader.makeGenomeWithProteinORFs()
        #1. we map peptides, either from Inspect, or from pre-mapped GFFs
        if self.InspectResultsPath:
            self.ParseInspect( self.InspectResultsPath )
            print "I found %s peptides from %s spectra"%(len(self.AllPeptides), self.SpectrumCount)
            self.MapAllPeptides(genome) # modify to add too genome/chromsome

        else :
            gffReader = PGPeptide.GFFPeptide( self.GFFInputPath )
            gffReader.generateORFs( self.ORFDatabasePaths, genome )

        self.FilterORFs(genome)

        if self.OutputPeptidesToGFF:
            self.WritePeptideGFFFile()
        if self.SearchForMispredictions:
            self.FindMiscalls(genome)
        if self.SearchForCleavage:
            self.AnalyzeCleavage()

    def CheckComplexity(self):
        if self.Verbose:
            print "ProteogenomicsPostProcessing.py:CheckComplexity"
        LowMW = ["G", "A"]
        List = []
        for ORF in self.AllORFs.values():
            PeptideString = ""
            for PeptideObject in ORF.PeptideLocationList:
                PeptideString += PeptideObject.Aminos
            Count = 0
            for Letter in PeptideString:
                if Letter in LowMW:
                    Count +=1
            Normalized = Count / float (len(PeptideString))
            List.append(Normalized)
        Histogram = BasicStats.CreateHistogram(List, 0, 0.05)
        BasicStats.PrettyPrintHistogram(Histogram, None)

    def FindOverlappingDubiousGeneCalls(self):
        """Parameters: none
        Return: none
        Description: There are genomic regions for which two gene calls overlap. 
        For some badly predicted genomes, the overlap is substantial (like >50 bp).
        I believe that most of these are bad gene calls, and I want to filter them 
        out.  This method calls the PROilters methodFindOverlapnGeneCalls
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
            #'Start' is the small number. ALWAYS.  5' refers t but start is always just the small number
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


    def FilterORFs(self,genome):
        """
        Parameters: None
        Return: None
        Description: glorified wrapper for calling the filter function. Afterwards
        we do some cleanup of ORFs that lack any peptides
        """
        if self.Verbose:
            print "ProteogenomicsPostProcessing.py:FilterORFs"
            print "Genome has %d chroms totaling %d simple orfs, %d pepOnly." % (
                genome.numChromosomes(), genome.numOrfs('Simple'), genome.numOrfs('PepOnly') )

        #this is the way we do filters.  We create all the filters that
        #we want to use, and then put them into th efilter list. It does the 
        #magic for us.
        SequenceComplexity = PGORFFilters.SequenceComplexityFilter()
        MinPeptide = PGORFFilters.MinPeptideFilter(2)
        Uniqueness = PGORFFilters.UniquenessFilter()
        FilterList = PGORFFilters.FilterList([SequenceComplexity, Uniqueness, MinPeptide]) # an anonymous list
        genome.filterORFs( FilterList )

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
        Return: None, but we populate ORF objects
        Description: Take a gff file of peptides mapped to the genome, and use
        that to populate our results, as opposed to parsing Inspect results
        and mapping them (which takes time)
        """
        #GFFReader = PGPeptide.GFFPeptide( GFFFile )
        #ORFDict   = GFFReader.generateORFs( Test.SixFrame )
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
            self.AllLocatedPeptides.append(Location)

    def WritePeptideGFFFile(self):
        """
        Parameters: None
        Return: None
        Description: This takes the peptide objects from self.AllLocatedPeptides and
        makes a GFF File containing all of them.  
        """
        GFFOut = PGPeptide.GFFPeptide(self.GFFOutputPath, "w") #creates my GFF converter. 'w' parameter is for writing

        for ORF in self.AllORFs.values(): # we cycle through the ORFs, writing all their located peptides
            GFFOut.writeORFPeptides(ORF)
            #print "%s"%ORF
        GFFOut.close()

    def FindMiscalls(self, genome):
        """This function goes through the ORFs and calls their method for
        determining if they were mispredicted. just a wrapper really
        """
        if self.Verbose:
            print "ProteogenomicsPostProcessing.py:FindMiscalls"
            print "Genome has %d chroms totaling %d simple orfs, %d pepOnly." % (
                genome.numChromosomes(), genome.numOrfs('Simple'), genome.numOrfs('PepOnly') )
        Count = 0
        NovelCount = 0
        UnderPredictedCount = 0
        for (chromName, chrom) in genome.chromosomes.items():
            BagChecker = PGPrimaryStructure.PrimaryStructure(
                self.OutputPath, str(chrom.sequence) )
            if self.Verbose:
                print "ORFs for chromosome %s of len %d" % (chromName, len(chrom.sequence))

            for (orfName,ORF) in chrom.simpleOrfs.items() + chrom.pepOnlyOrfs.items():

                Status = BagChecker.CheckStructure(ORF)
                if Status == "NOVEL":
                    NovelCount += 1
                if Status == "UNDERPREDICTED":
                    UnderPredictedCount += 1
                Count +=1
                if (Count %1000) == 0 and self.Verbose:
                    print "Processed %s ORF Objects" % (Count)
        #done with loop
        print "Processed %s. %s novel and %s underpredicted"%(Count, NovelCount, UnderPredictedCount)
#        print "Named %s, hypothetical %s"%(BagChecker.NamedCount, BagChecker.HypotheticalCount)
#        BagChecker.OutputGCFiles()
#        BagChecker.OutputLenFiles()


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

    def MapAllPeptides(self,genome):
        """
        Parameters: Genome object to add peptides too
        Return: None
        Description: Go through self.AllPeptides, and find their location 
        within ORFs. We make LocatedProtein objects out of these. This sets 
        up the self.AllLocatedPeptides variable. We also make a list of 
        observed ORFs to that when we create ORFs, we only do so for those
        that will matter.
        """
        if self.Verbose:
            print "ProteogenomicsPostProcessing.py:MapAllPeptides"
        Count = 0
        LocationCount=0
        inGenomeCount = 0
        for (Aminos, PValue) in self.AllPeptides.items():
            Count += 1
            if (Count %1000) == 0 and self.Verbose:
                print "Mapped %s / %s peptides"%(Count, len(self.AllPeptides))
            #print Aminos
            LocatedPeptides = self.ORFPeptideMapper.MapPeptide(Aminos, PValue)
            if self.UniquenessFlag and (len(LocatedPeptides) > 1):
                continue #skip out on adding it to the list because it's not unique.  And that's what you asked for
            self.AllLocatedPeptides.extend(LocatedPeptides)
            LocationCount += len(LocatedPeptides)

            #now put into the observedORF stuff
            for Location in LocatedPeptides:

                orfInGenome = genome.getOrf( Location.ORFName )
                if orfInGenome:
                    inGenomeCount += 1
                    orfInGenome.addLocatedPeptide( Location )
                else:
                    # ORF not in genome, create a Peptide only orf
                    orf = PGPeptide.OpenReadingFrame(name=Location.ORFName)
                    orf.location = Location.location
                    genome.addOrf( orf, 'PepOnly' )
                    orf.addLocatedPeptide( Location )

                if not Location.ORFName in self.ProteomicallyObservedORFs:
                    self.ProteomicallyObservedORFs.append(Location.ORFName)

        genome.addSeqToPepOnlyOrfs( self.ORFDatabasePaths )

        self.AllPeptides = {}  # set to null just for the memory savings
        if self.Verbose:
            print "ProteogenomicsPostProcessing:MapAllPeptides - mapped %s peptides to %s genomic locations %s found in genome"%(Count, LocationCount,inGenomeCount)

    def ParseInspect(self, FilePath):
        """Here I parse out Inspect Results to get peptide annotations, 
        Putting them in a hash for safe keeping
        """
        FalseAminos = []
        inspectParser = InspectResults.Parser( FilePath )
        for result in inspectParser:
            try:
                Annotation = result.Annotation
                Peptide = GetPeptideFromModdedName(Annotation)
                Aminos = Peptide.Aminos
                PValue = result.PValue
                InspectMappedProtein = result.ProteinName
            except:
                traceback.print_exc()
                continue # SNAFU
            #here's some hacking that needs to be fixed.  I currently want to filter stuff by lfdr, but
            #that may not always be present
            if result.LFDR == None: #meaning that I don't have the column in question
                if PValue > self.PValueLimit: #substitute pvalue for LFDR
                    continue
            else:
                LFDR = result.LFDR
                PValue = LFDR
                if LFDR > self.PValueLimit:
                    continue

            if InspectMappedProtein[:3] == "XXX":
                #this is a true negative, let's count them
                if not Aminos in FalseAminos:
                    FalseAminos.append(Aminos)
                continue
            if not self.AllPeptides.has_key(Aminos):
                self.AllPeptides[Aminos] = PValue
            else:
                if PValue <  self.AllPeptides[Aminos]:
                    self.AllPeptides[Aminos] = PValue
            self.SpectrumCount += 1
        print "I got %s truedb peptides, and %s decoy peptides"%(len(self.AllPeptides), len(FalseAminos))

    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "b:r:g:d:w:uvi:o:p:CMG:Wn:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-r":
                # -r results file(s)
                if not os.path.exists(Value):
                    print "** Error: couldn't find results file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.InspectResultsPath = Value
            if Option == "-b":
                if not os.path.exists(Value):
                    print "** Error: couldn't find input file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.GenbankPath = Value
            if Option == "-g":
                if not os.path.exists(Value):
                    print "** Error: couldn't find results file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.GFFInputPath = Value
            if Option == "-d":
                if not os.path.exists(Value):
                    print "** Error: couldn't find database file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.ProteomeDatabasePaths.append( Value)
            if Option == "-o":
                if not os.path.exists(Value):
                    print "** Error: couldn't find openreadingframe database file '%s'\n\n"%Value
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
        #if not OptionsSeen.has_key("-w") or not OptionsSeen.has_key("-d") or not OptionsSeen.has_key("-o"):
        if not OptionsSeen.has_key("-w") or not OptionsSeen.has_key("-o"):
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
