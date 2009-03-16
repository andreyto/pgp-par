UsageInfo = """ProfilePeptideClusters.py
Given a set of annotation for some spectra, we measure the interpeptide 
distance for peptides which belong to the same ORF.
WARNING: We are trying to get the genomic distance, therefore
the database must have genome location in the ORF name.
WARNING: note that peptides must belong to the same ORF. 
WARNING: I have not dealt with spliced location crap yet.

Required Options
 -r [FileName] File or directory of results.
 -d [Trie file] Database used for search
 -w [FileName]  Output file
 
Required Options
 -u   Flag for requiring peptide uniqueness within the database.  Default
     behavior is to allow any peptide.


"""

import sys
import os
import getopt
import traceback
import ResultsParser
#import SelectProteins
import PeptideMapper
import BasicStats
from Utils import *
Initialize()


class ProfileClass(ResultsParser.ResultsParser):
    def __init__(self):
        self.ReferenceResults = None
        self.OutputPath = "RenameYourOutput.txt"
        self.DatabasePaths = [] #possibly multiple
        self.AllPeptides = {} # AminoSequence -> spectrumCount
        self.AllLocations = [] #list of PeptideMapper.GenomicLocationForPeptide
        self.UniquenessFlag = 0
        self.SpectrumCount = 0
        self.PeptideMapper = PeptideMapper.PeptideMappingClass()
        self.InterPeptideLengths = [] #just append to the list.
        self.PValueLimit = 1.0 #essentially lets anything through, because pvalues are 0.00 (good) to 1.0 (bad)
        self.LargestGapPerProtein = {} # protein->interpeptideDistance 
        ResultsParser.ResultsParser.__init__(self)
        
    def Main(self):
        UseTargetDB = 1 
        UseDecoyDB = 0
        self.PeptideMapper.LoadDatabases(self.DatabasePaths)
        self.ProcessResultsFiles(self.ReferenceResults, self.ParseInspectCallback)
        print "I found %s peptides from %s spectra"%(len(self.AllPeptides), self.SpectrumCount)
        ##now map all the peptides
        self.MapAllPeptides()
        self.OutHandle = open(self.OutputPath, "wb")
        self.CharacterizeProteins(UseTargetDB,UseDecoyDB) 
        self.OutHandle.close()
        
    def CharacterizeProteins(self, UseTargetDB = 1, UseDecoyDB = 0):
        """now that peptides have been mapped, we start looking at how much space
        exists inbetween them. If peptides overlap, it is possible for this space
        to be very small.  WE have two default parameters about whether we will
        use the target, decoy or both.  by default, we only work for the target db
        although you could be fancy and do something different
        """
        print "ProfilePeptideClusters.py:CharacterizeProteins"
        AllProteinNames = self. GetAllProteinNamesFromLocations(self.AllLocations)
        #go through all locations, get peptides in the same protein. then send
        #send those out to get measured
        Count =0 
        for Protein in AllProteinNames:
            TargetDB = 1
            if Protein[:3] == "XXX":
                TargetDB = 0
            #now that we determined where this came from, decide
            if TargetDB and (not UseTargetDB):
                continue
            if (not TargetDB) and (not UseDecoyDB):
                continue
            #print Protein
            LocationsForSingleProtein = []
            for Location in self.AllLocations:
                if Location.ProteinName == Protein:
                    LocationsForSingleProtein.append(Location)
                    
            #done with the loop
            #print "Protein %s has %s peptides"%(Protein, len(LocationsForSingleProtein))
            self.FindInterPeptideLengths(LocationsForSingleProtein, Protein)
            Count+= 1
            if (Count%100) == 0:
                print "characterized %s of %s proteins"%(Count, len(AllProteinNames))
        Histogram = BasicStats.CreateHistogram(self.InterPeptideLengths, 0, 1)
        BasicStats.PrettyPrintHistogram(Histogram, self.OutHandle)
            #break
            
    def FindInterPeptideLengths(self, ListOfLocations, ProteinName):
        """Given a list of peptides, we sort them and then find the distance in between neighboring pairs
        If they overlap in sequence, then there is no interpeptide length
        """
        ListOfLocations.sort(self.SortLocations)
        PrevLocation = None
        BiggestGap = -1 #used for the larget interpeptide distance per protein
        for Location in ListOfLocations:
            #Location.PrintMe()
            if PrevLocation:
                # do the comparison, first look for overlap
                Overlap = self.DoPeptidesOverlap(Location, PrevLocation)
                if Overlap:
                    #here we set the interpeptide length to zero.
                    self.InterPeptideLengths.append(0)
                    if BiggestGap < 0:
                        BiggestGap = 0
                else:
                    ##defined as B.start-A.stop.  we can do it this way, because 
                    ###start and stop are numerically ordered, not 5' associated
                    Distance = int(Location.StartNucleotide - PrevLocation.StopNucleotide)
                    self.InterPeptideLengths.append(Distance)
                    if Distance > BiggestGap:
                        BiggestGap = Distance
            PrevLocation = Location
        self.LargestGapPerProtein[ProteinName] = BiggestGap
        
    def DoPeptidesOverlap(self, A, B):
        """test whether Peptides overlap in their nucleotide sequence
        I know there is probably a better way to do this, but simply
        checking for set intersection is pretty straightforward
        """
        ANucs = range(A.StartNucleotide, A.StopNucleotide+1)
        BNucs = range(B.StartNucleotide, B.StopNucleotide+1)
        for ANuc in ANucs:
            if ANuc in BNucs:
                return 1
            
        return 0
        
        
    def SortLocations(self, A, B):
        """Used to sort. smallest start, then smallest stop
        """
        if A.StartNucleotide < B.StartNucleotide:
            return -1
        elif A.StartNucleotide > B.StartNucleotide:
            return 1
        else: # same start, sort by stop
            if A.StopNucleotide < B.StopNucleotide:
                return -1
            else:
                return 1
        
    def GetAllProteinNamesFromLocations(self, ListOfLocations):
        """Given a list, we cycle through.  simple"""
        ProteinNames = []
        Counter = 0
        for Location in ListOfLocations:
            SingleProtein = Location.ProteinName
            Counter += 1
            if not SingleProtein in ProteinNames:
                ProteinNames.append(SingleProtein)
        print "I got %s protein names from %s locations"%(len(ProteinNames),Counter)
        return ProteinNames

    def MapAllPeptides(self):
        """Go through the list in self.AllPeptides, and find the location for each
        peptide. It is true that we could have parsed this out of the Inspect file, but
        some peptides have multiple locations, so we have to look it up. This sets up
        the self.AllLocations variable
        """
        print "ProfilePeptideClusters.py:MapAllPeptides"
        Count = 0
        LocationCount=0
        count2 =0 
        for Aminos in self.AllPeptides.keys():
            #print Aminos
            GenomicLocations = self.PeptideMapper.MapMe(Aminos, 1)
            count2 += len (GenomicLocations)
            if self.UniquenessFlag and (len(GenomicLocations) > 1):
                continue #skip out on adding it to the list
            self.AllLocations.extend(GenomicLocations)
            LocationCount += len(GenomicLocations)
            Count += 1
            if (Count %500) == 0:
                print "Mapped %s peptides"%Count
        self.AllPeptides = {}  # set to null just for the memory savings
        print "I mapped %s peptides to %s (%s) genomic locations"%(Count, LocationCount, count2)
        
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
            if PValue > self.PValueLimit:
                continue
            ##Hack to get the fake Distribution, cut out all targetDB
            #Protein = Bits[self.Columns.ProteinName]
            #if not Protein[:3] == "XXX":
            #    continue
            if not self.AllPeptides.has_key(Aminos):
                self.AllPeptides[Aminos] = 0 
            self.AllPeptides[Aminos] += 1
            self.SpectrumCount += 1
        Handle.close()


    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:d:w:u")
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
            if Option == "-d":
                if not os.path.exists(Value):
                    print "** Error: couldn't find database file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.DatabasePaths.append( Value)
            if Option == "-w":
                self.OutputPath = Value
            if Option == "-u":
                self.UniquenessFlag = 1
            if Option == "-p":
                #undocumented pvalue filter
                self.PValueLimit = float(Value)
        if not OptionsSeen.has_key("-r")  or not OptionsSeen.has_key("-w") or not OptionsSeen.has_key("-d"):
            print UsageInfo
            sys.exit(1)


if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(psyco not found - running in non-optimized mode)"
    FBI = ProfileClass()
    FBI.ParseCommandLine(sys.argv[1:])
    FBI.Main()