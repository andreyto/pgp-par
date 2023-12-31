#!/usr/bin/env python
###############################################################################
#                                                                             # 
#       Copyright (c) 2009 J. Craig Venter Institute.                         #     
#       All rights reserved.                                                  #
#                                                                             #
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.    #
#                                                                             #
###############################################################################





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
import PGCleavageAnalysis
import GFFIO
from Utils import *
Initialize()
from itertools import combinations # for the combinations of double for loops

class FinderClass():
    def __init__(self):
        self.InspectResultsPath = None #Peptide/Spectrum Matches from Inspect
        self.GFFInputPath = None #alternate input form, mostly used to shunt mapping, because that takes so long
        self.GenbankPath  = None
        self.OutputPath = "RenameYourOutput.txt" #really a stub for output, as it gets modified to *.novel.faa, *.novel.info, etc.
        self.ProteomeDatabasePaths = [] #possibly multiple
        self.ORFDatabasePaths = [] #possibly multiple
        self.AllPeptides = {} # AminoSequence -> best pvalue, nulled out in MapAllPeptides
        self.PeptideSources = {} #aminosequence->[(file,spectrum), (file,spectrum), ...]
        self.AllLocatedPeptides = [] #list of PeptideMapper.GenomicLocationForPeptide
        self.AllPredictedProteins = {} #predictedProteinName ->GenomicLocationForPeptide Object, used in CreateORFs
        self.ProteomicallyObservedORFs = [] # this is populated when peptides are mapped, and deleted after ORF Objects are created
        self.UniquenessFlag = 0
        self.InterPeptideDistanceMax = 1000 # sensible default
        self.PValueLimit = 0.05 #pvalues are 0.00 (good) to 1.0 (bad) 
        self.Verbose = 0
        self.VerboseWarnings = 0
        self.SearchForMispredictions = 1
        self.SearchForCleavage = 1
        self.OutputPeptidesToGFF = 0
        self.GFFOutputPath = "RenameYourOutput.gff"
        self.Report = PGPReport()

    def Main(self):
        chromReader = PGPeptide.GenbankGenomeReader(self.GenbankPath, self.ORFDatabasePaths)
        genome = chromReader.makeGenomeWithProteinORFs()
        self.Report.SetValue("MappedProteins", genome.numOrfs('Simple'))
        self.Report.SetValue("UnmappedProteinsComplex", genome.numOrfs('Complex'))
        self.Report.SetValue("UnmappedProteinsSNAFU", genome.numOrfs('Other'))
        #1. we map peptides, either from Inspect, or from pre-mapped GFFs
        if self.InspectResultsPath:
            self.ParseInspect( self.InspectResultsPath )
            self.MapAllPeptides(genome) # modify to add too genome/chromsome

        else :
            gffReader = PGPeptide.GFFPeptide( self.GFFInputPath )
            gffReader.generateORFs( self.ORFDatabasePaths, genome )

        self.FilterORFs(genome)

        #now look for overlaps.  This is an effort at quality of assignments
        # or confusion in the genome
        self.CheckForOverlaps(genome, self.OutputPath)

        if self.OutputPeptidesToGFF:
            self.WritePeptideGFFFile(genome)
        if self.SearchForMispredictions:
            self.FindMiscalls(genome)
        if self.SearchForCleavage:
            self.AnalyzeCleavage(genome)
        #write out the simple protein inference
        self.WriteProteinInference(genome)
        #write out the report
        (Path, Ext) = os.path.splitext(self.OutputPath)
        ReportPath = "%s.report.txt"%Path
        self.Report.SetFileName(ReportPath)        
        self.Report.WriteReport()
        
    def WriteProteinInference(self, genome):
        """Parameters: NOne
        Return: NOne
        Description: write a simple protein inference file, just proteins 
        and their mapped peptides, separated out into unique and shared
        """
        Header = "Protein Name\tNum Peptides\tUnique Count\tUniquePeptides\tShared Peptides\n"
        (Path, Ext) = os.path.splitext(self.OutputPath)
        InferencePath = "%s.proteininference.txt"%Path
        Handle = open(InferencePath, "w")
        Handle.write(Header)
        for (chromName, chrom) in genome.chromosomes.items():
            for (orfName,ORF) in chrom.simpleOrfs.items() + chrom.pepOnlyOrfs.items():
                if ORF.numPeptides() == 0: #remember, we keep these around for the overlap comparison
                    continue
                ORF.WriteSimpleProteinInference(Handle)
        Handle.close()
        
        #now we write a peptidome protein submission file
        PeptidomePath = "%s.peptidome.txt"%Path
        Handle = open(PeptidomePath, "w")
        Header = "Protein\tPeptide\tSpectrum\n"
        Handle.write(Header)
        for (chromName, chrom) in genome.chromosomes.items():
            for (orfName,ORF) in chrom.simpleOrfs.items() + chrom.pepOnlyOrfs.items():
                if ORF.numPeptides() == 0: #remember, we keep these around for the overlap comparison
                    continue
                ORF.WritePeptidome(Handle)
        Handle.close()
        


    def CheckForOverlaps(self, genome, OutputPath):
        """Parameters: none
        Return: none
        Description: There are genomic regions for which two gene calls overlap. 
        For some badly predicted genomes, the overlap is substantial (like >50 bp).
        I believe that most of these are bad gene calls, and I want to filter them 
        out.  This method calls the PROilters methodFindOverlapnGeneCalls
        which does that.  First we have to build the right dictionary,  so we do that here
        """
        #first let's get our output path fixed up
        (Path, Ext) = os.path.splitext(OutputPath)
        NewPath = "%s.%s"%(Path, "conflictreport.txt")
        OutHandle = open(NewPath, "w")
        ConflictLevelCount = [0,0,0,0,0,0,0]
        if self.Verbose:
            print "ProteogenomicsPostProcessing.py:CheckForOverlaps"
        for (chromName, chrom) in genome.chromosomes.items():
            if self.Verbose:
                print "ORFs for chromosome %s of len %d" % (chromName, len(chrom.sequence))
            AllORFs = chrom.simpleOrfs.values() + chrom.pepOnlyOrfs.values()
            for (ORF1, ORF2) in combinations(AllORFs, 2): #this is from itertools and gives all combinations (symetry cancelled)
                #first we check that at least one of these dudes has some proteomic coverage
                #because we really don't care if stuff is not supported 
                if ORF1.numPeptides() == 0 and ORF2.numPeptides() == 0:
                    continue
                #now we have to get the locations right for these things.  Let's do the simple stuff
                #first.  We just get the overlap of gene coords (or coverage coords for pepOnlyOrfs)
                #this can help us with level 0,1,2 then we do some more difficult stuff
                Location1 = self.GetLocationForOverlapComparison(ORF1)
                Location2 = self.GetLocationForOverlapComparison(ORF2)
                Result = Location1.overlap(Location2)
                ##Result could be 'None' meaning no overlap
                if not Result:
                    ConflictLevelCount[0] += 1
                    continue
                (StartOverlap, StopOverlap) = Result
                Len = StopOverlap - StartOverlap
                if Len < 10:
                    ConflictLevelCount[1] += 1
                elif Len < 40: 
                    ConflictLevelCount[2] += 1
                else:
                    State = self.ComplexOverlapAnalysis(ORF1, ORF2, Len, OutHandle)
                    ConflictLevelCount[State] += 1
        String = "Conflict States\n"
        String += "0: no conflict, %s\n"%ConflictLevelCount[0]
        String += "1: overlap < 10bp, %s\n"%ConflictLevelCount[1]
        String += "2: overlap < 40bp, %s\n"%ConflictLevelCount[2]
        String += "3: overlap > 40bp, overlap an unsupported hypothetical, %s\n"%ConflictLevelCount[3]
        String += "4: overlap > 40bp, overlap an unsupported named protein, %s\n"%ConflictLevelCount[4]
        String += "5: overlap > 40bp, no peptides in overlap region, %s\n"%ConflictLevelCount[5]
        String += "6: overlap > 40bp, peptides in overlap region, %s\n"%ConflictLevelCount[6]
        OutHandle.write(String)
        print String
                
    def ComplexOverlapAnalysis(self, ORF1, ORF2, Len, Handle):
        """Parameters: two OpenReadingFrame objects, an open file handle
        Return: none
        Description: for ORFs with more than 40bp of overlap, we do some more 
        indepth analysis.  We see if they both have proteomic support, and
        whether the supporting regions overlap
        """
        Peps1 = ORF1.numPeptides()
        Peps2 = ORF2.numPeptides()
        if Peps1 == 0:
            #no peptides for ORF1.  Is it named "hypothetical"
            ProteinName1 = ORF1.GetProteinName()
            if ProteinName1.find("hypothetical") > -1:
                Handle.write("Level 3 conflict: overlap (%s bp) an unsupported hypothetical\n"%Len )
                Handle.write("%s\n%s\n\n"%(ORF1, ORF2))# uses PGPeptide.OpenReadingFrame.__string__
                ORF1.SetConflictString("Level 3 conflict: overlap an unsupported hypothetical\n%s\n"%ORF2)
                ORF2.SetConflictString("Level 3 conflict: overlap an unsupported hypothetical\n%s\n"%ORF1)
                return 3
            Handle.write("Level 4 conflict: overlap (%s bp) an unsupported named protein\n"%Len )
            Handle.write("%s\n%s\n\n"%(ORF1, ORF2))# uses PGPeptide.OpenReadingFrame.__string__
            ORF1.SetConflictString("Level 4 conflict: overlap an unsupported named protein\n%s\n"%ORF2 )
            ORF2.SetConflictString("Level 4 conflict: overlap an unsupported named protein\n%s\n"%ORF1 )
            return 4
        #now same for ORF2
        if Peps2 == 0:
            #print "%s"%ORF2
            ProteinName2 = ORF2.GetProteinName()
            if ProteinName2.find("hypothetical") > -1:
                Handle.write("Level 3 conflict: overlap (%s bp) an unsupported hypothetical\n"%Len  )
                Handle.write("%s\n%s\n\n"%(ORF1, ORF2))# uses PGPeptide.OpenReadingFrame.__string__
                ORF1.SetConflictString("Level 3 conflict: overlap an unsupported hypothetical\n%s\n"%ORF2)
                ORF2.SetConflictString("Level 3 conflict: overlap an unsupported hypothetical\n%s\n"%ORF1)
                return 3
            Handle.write("Level 4 conflict: overlap (%s bp) an unsupported named protein\n"%Len  )
            Handle.write("%s\n%s\n\n"%(ORF1, ORF2))# uses PGPeptide.OpenReadingFrame.__string__
            ORF1.SetConflictString("Level 4 conflict: overlap an unsupported named protein\n%s\n"%ORF2 )
            ORF2.SetConflictString("Level 4 conflict: overlap an unsupported named protein\n%s\n"%ORF1 )
            return 4
        #now we're both represented. that's not good news
        #we need to make new locations, which are simply the proteomics-stop.  
        #then we compare if those areas have a conflict.
        (Start1, Stop1) = ORF1.GetObservedDNACoords()
        (Start2, Stop2) = ORF2.GetObservedDNACoords()
        Strand1 = ORF1.GetStrand()
        Strand2 = ORF2.GetStrand()
        Location1 = PGPeptide.GenomicLocation(Start1, Stop1, Strand1)
        Location2 = PGPeptide.GenomicLocation(Start2, Stop2, Strand2)
        Result = Location1.overlap(Location2)
        if not Result:
            #no overlap
            Handle.write("Level 5 conflict: overlap of two supported proteins (%s bp), but overlapping region not supported\n"%Len  )
            Handle.write("%s\n%s\n\n"%(ORF1, ORF2))# uses PGPeptide.OpenReadingFrame.__string__
            ORF1.SetConflictString("Level 5 conflict: overlap of two supported proteins, but overlapping region not supported\n%s\n"%ORF2 )
            ORF2.SetConflictString("Level 5 conflict: overlap of two supported proteins, but overlapping region not supported\n%s\n"%ORF1 )
            return 5
        Handle.write("Level 6 conflict: overlap of two supported proteins (%s bp) , overlapping region has peptide support\n"%Len  )
        Handle.write("%s\n%s\n\n"%(ORF1, ORF2))# uses PGPeptide.OpenReadingFrame.__string__
        ORF1.SetConflictString("Level 6 conflict: overlap of two supported proteins, overlapping region has peptide support\n%s\n"%ORF2 )
        ORF2.SetConflictString("Level 6 conflict: overlap of two supported proteins, overlapping region has peptide support\n%s\n"%ORF1 )
        return 6   
        
                    
    def GetLocationForOverlapComparison(self, ORF):
        """Parameters: an OpenReadingFrame object
        Return: a GenomicLocation object
        Description: we get the location of a protein or peptideset for comparinson
        in the overlap contest.  When an ORF has a locatedprotein, then we use that
        if not, then we give the bounding box of the peptides-stop.
        """
        Protein = ORF.GetLocatedProtein()
        if Protein:
            return Protein.GetLocation()
        #now we make this up.  It has no protein, so it must have peptides
        (Start, Stop) = ORF.GetObservedDNACoords()
        Strand = ORF.GetStrand()
        Location = PGPeptide.GenomicLocation(Start, Stop, Strand)
        return Location

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
        PeptideDistance = PGORFFilters.PeptideDistance()
        MinPeptide = PGORFFilters.MinPeptideFilter(2)
        Uniqueness = PGORFFilters.UniquenessFilter()
        Tryptic = PGORFFilters.TrypticFilter()

        MinPepONE = PGORFFilters.MinPeptideFilter(1)
        #FilterList = PGORFFilters.FilterList( [MinPeptide, ], self.OutputPath) # an anonymous list
        FilterList = PGORFFilters.FilterList([ # an anonymous list
            SequenceComplexity, PeptideDistance, Uniqueness, Tryptic, MinPeptide
                                             ], self.OutputPath)
        genome.filterORFs( FilterList )
        #last thing is to set the report for filters
        ListString = FilterList.GetListString()
        self.Report.SetValue("FiltersUsed", ListString)


    def AnalyzeCleavage(self, genome):
        """
        Parameters: None
        Return: None
        Description: glorified wrapper for calling the filter function.
        I get returned the following
        0 = tryptic nterminal first peptide
        1 + non-tryptic nterminal first peptide
        1 = too short
        2 = too long
        3 = no patch
        4 = patch 
        5 = all three signals
        """

        #make ourselves a file handle here
        (PathStub, Ext) = os.path.splitext(self.OutputPath)
        SignalPeptidePath = "%s.%s"%(PathStub, "signalpeptide.info")
        Handle = open (SignalPeptidePath, "w")
        Header = PGCleavageAnalysis.GetSignalPeptideHeader()
        Handle.write(Header)
        
        ProteinsWithNonTrypticNTermTotal = 0
        ProteinsWithHydrophobicPatch = 0
        ProteinLackHydrophobicPatch = 0
        ProteinsWithItAll = 0
        WrongLength = 0
        if self.Verbose:
            print "ProteogenomicsPostProcessing.py:AnalyzeCleavage"
        for (chromName, chrom) in genome.chromosomes.items():
            if self.Verbose:
                print "ORFs for chromosome %s of len %d" % (chromName, len(chrom.sequence))

            for (orfName,ORF) in chrom.simpleOrfs.items() + chrom.pepOnlyOrfs.items():
                #first we do some checking.  
                if ORF.numPeptides() == 0:
                    # we dont want to deal with things that don't have peptides,
                    continue

                if not ORF.CDS: ## temporary maybe.  I should not be looking
                    ## at things without a protein prediction, or I should just get the protein
                    ## sequence of the whole ORF.  Another decision postponed.
                    continue
                
                (Decision, String) = PGCleavageAnalysis.Analysis(ORF, "Trypsin")
                if Decision: #all things not zero
                    ProteinsWithNonTrypticNTermTotal += 1
                if Decision == 1 or Decision == 2:
                    WrongLength +=  1
                if Decision == 3:
                    ProteinLackHydrophobicPatch += 1
                if Decision >= 4:
                    ProteinsWithHydrophobicPatch += 1
                if Decision == 5:
                    ProteinsWithItAll += 1
                if Decision >= 3:
                    #these are thing I want to print out, conveniently they also have a non NULL String
                    Handle.write(String)
        print "These are the counts"
        print "all those with a non-tryptic: %s"%ProteinsWithNonTrypticNTermTotal
        print "wrong length: %s"%WrongLength
        print "no patch: %s"%ProteinLackHydrophobicPatch
        print "Have Patch: %s"%ProteinsWithHydrophobicPatch
        print "Have All required motifs: %s"%ProteinsWithItAll
        self.Report.SetValue("SignalPeptides", ProteinsWithItAll)
        Handle.close()

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

    def WritePeptideGFFFile(self, genome):
        """
        Parameters: None
        Return: None
        Description: This takes the peptide objects from self.AllLocatedPeptides and
        makes a GFF File containing all of them.  
        """
        if self.Verbose:
            print "ProteogenomicPostProcessing.py::WritePeptideGFFFile"

        GFFOut = PGPeptide.GFFPeptide(self.GFFOutputPath, "w") #creates my GFF converter. 'w' parameter is for writing

        for (ChromName, Chrom) in genome.chromosomes.items():
            for (ORFName, ORF) in Chrom.simpleOrfs.items() + Chrom.pepOnlyOrfs.items():
                # we cycle through the ORFs, writing all their located peptides
                GFFOut.writeORFPeptides(ORF)

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
        print "Processed %s for potential miscalls. %s novel and %s underpredicted"%(Count, NovelCount, UnderPredictedCount)
        self.Report.SetValue("NovelORFCount", NovelCount)
        self.Report.SetValue("ORFWrongStartCount", UnderPredictedCount)


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
        Count = 0 #just a raw total.  no higher meaning
        LocationCount=0
        inGenomeCount = 0
        NonMappingCount = 0
        MappingCount = 0 #those that do map to our databases
        ORFPeptideMapper = PeptideMapper.PeptideMappingClass()
        ORFPeptideMapper.LoadDatabases(self.ORFDatabasePaths) #a handle for the 6frame translations db (called ORF)
        for (Aminos, PValue) in self.AllPeptides.items():
            Count += 1
            if (Count %1000) == 0 and self.Verbose:
                print "Mapped %s / %s peptides"%(Count, len(self.AllPeptides))
            #print Aminos
            MSMSSources = self.PeptideSources[Aminos] #the are all the spectra that are for a amino acid string
            LocatedPeptides = ORFPeptideMapper.MapPeptide(Aminos, PValue, MSMSSources)
            
            if len(LocatedPeptides) == 0:
                #this peptide sequence did not map to the database.  That might mean that
                #it is a common contaminant, like trypsin and I'm not caring about those. 
                #but in any case it is something that we want to keep track of
                NonMappingCount += 1
            else:
                MappingCount += 1

            if self.UniquenessFlag and (len(LocatedPeptides) > 1):
                continue #skip out on adding it to the list because it's not unique.  And that's what you asked for
            self.AllLocatedPeptides.extend(LocatedPeptides)
            LocationCount += len(LocatedPeptides)

            #now put into the observedORF stuff
            for Location in LocatedPeptides:

                orfInGenome = genome.getOrf( Location.ORFName, Location.chromosome )
                if orfInGenome:
                    inGenomeCount += 1
                    orfInGenome.addLocatedPeptide( Location )
                else:
                    # ORF not in genome annotation, create a Peptide only orf
                    orf = PGPeptide.OpenReadingFrame(name=Location.ORFName)
                    orf.location = Location.location
                    genome.addOrf( orf, 'PepOnly' )
                    orf.addLocatedPeptide( Location )

                if not Location.ORFName in self.ProteomicallyObservedORFs:
                    self.ProteomicallyObservedORFs.append(Location.ORFName)

        genome.addSeqToPepOnlyOrfs( self.ORFDatabasePaths )

        self.AllPeptides = {}  # set to null just for the memory savings
        #now we set variables in the report
        PeptidesMappingOutsideAnnotation = LocationCount - inGenomeCount
        self.Report.SetValue("UnmappedPeptides", NonMappingCount)
        self.Report.SetValue("MappedPeptideLocationCount", LocationCount)
        self.Report.SetValue("MappedPeptideCount", MappingCount)
        self.Report.SetValue("PeptidesMappingOutsideAnnotation", PeptidesMappingOutsideAnnotation)
        
        #we're now done with the self.PeptideSources array, let's clean it up just to be nice for memory footprint
        del self.PeptideSources
        if self.Verbose:
            print "ProteogenomicsPostProcessing:MapAllPeptides - mapped %s peptides to %s genomic locations %s found in current annotation"%(MappingCount, LocationCount,inGenomeCount)


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
                SourceTuple = (File,Spectrum)
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
                self.PeptideSources[Aminos].append(SourceTuple)
            else:
                self.PeptideSources[Aminos].append(SourceTuple)
                if PValue <  self.AllPeptides[Aminos]:
                    self.AllPeptides[Aminos] = PValue

        print "I got %s truedb peptides, and %s decoy peptides (%s spectra)"%(len(self.AllPeptides), len(FalseAminos), SpectrumCount)
        self.Report.SetValue("TruePeptides", len(self.AllPeptides))
        self.Report.SetValue("DecoyPeptides", len(FalseAminos))
        self.Report.SetValue("SpectraProcessed", SpectrumCount)

    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "b:r:g:d:w:uvi:o:p:CMG:W")
        OptionsSeen = {}
        #set our report
        self.Report.SetValue("CommandLine", Arguments)
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


class PGPReport():
    """Takes information about the data run and formats it so that you can
    have a pretty report for later.
    """
    def __init__(self):
        self.Info = {}
        self.Info["CommandLine"] = None
        self.Info["SpectraProcessed"] = None
        self.Info["TruePeptides"] = None
        self.Info["DecoyPeptides"] = None
        self.Info["MappedProteins"] = None
        self.Info["UnmappedProteinsSNAFU"] = None
        self.Info["UnmappedProteinsComplex"] = None
        self.Info["MappedPeptideCount"] = None
        self.Info["MappedPeptideLocationCount"] = None
        self.Info["UnmappedPeptides"] = None
        self.Info["ORFCountPreFilter"] =None
        self.Info["ORFCountPostFilter"] = None
        self.Info["NovelORFCount"] = None
        self.Info["ORFWrongStartCount"] = None
        self.Info["PGPVersion"] = None
        self.Info["FiltersUsed"] = None
        self.Info["PeptidesMappingOutsideAnnotation"]= None
        self.Info["SignalPeptides"] = None
        self.FileName = "renameyourreport.txt"

    def WriteReport(self):
        #create the long string.  Sooo boring
        String = ""
        String += "Proteogenomics version %s\n"%self.Info["PGPVersion"]
        String += "Command Line: %s\n"%self.Info["CommandLine"]
        String += "Spectra Passing PValue: %s\n"%self.Info["SpectraProcessed"]
        String += "True Peptides %s, False Peptides %s\n"%(self.Info["TruePeptides"], self.Info["DecoyPeptides"])
        String += "Mapped %s proteins (%s failed complex, %s failed other)\n"%(self.Info["MappedProteins"],
                                                                               self.Info["UnmappedProteinsComplex"],
                                                                               self.Info["UnmappedProteinsSNAFU"])
        String += "Mapped %s peptides to %s locations\n"%(self.Info["MappedPeptideCount"], 
                                                          self.Info["MappedPeptideLocationCount"])
        String += "\t%s Unmappable peptides\n"%self.Info["UnmappedPeptides"]
        String += "\t%s Peptide locations outside of current annotation\n"%self.Info["PeptidesMappingOutsideAnnotation"]
        String += "ORFs analyzed: %s prefilter, %s postfilter\n"%(self.Info["ORFCountPreFilter"],
                                                                  self.Info["ORFCountPostFilter"])
        String += "Filters employed %s\n"%self.Info["FiltersUsed"]
        String += "Novel proteins: %s\n"%self.Info["NovelORFCount"]
        String += "Underpredicted proteins: %s\n"%self.Info["ORFWrongStartCount"]
        String += "Likely Signal Peptide cleavages: %s\n"%self.Info["SignalPeptides"]

        Handle = open(self.FileName, "w")
        Handle.write(String)
        Handle.close()

    def SetValue(self, Key, Value):
        if not self.Info.has_key(Key):
            print "WARNING: report called with bad key %s"%Key
            return
        self.Info[Key] = Value
    def SetFileName(self, FileName):
        self.FileName = FileName


if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(psyco not found - running in non-optimized mode)"
    Gumshoe = FinderClass()
    Gumshoe.ParseCommandLine(sys.argv[1:])
    Gumshoe.Main()
