"""PGORFFilters.py
This auxiliary set of functions is to be used to filter ORFs (not proteins.  
open reading frames).  The filters are 

1. 2 peptides/ORF
2. Clustering
3. Sequence content
4. Tryptic content of cluster

NOTE: this is a utility, and not executable from the command line

"""

class Filter:
    """Class Filter: this is a generic filter, meant to be inherited to 
    the actual specific filter classes, e.g. SequenceComplexityFilter
    Classes that inherit from this must define - 
    Variables: 
        self.name
    Functions:
        apply(ORF)
    """
    def __init__(self):
        """Parameters: none
        Return: none
        Description: trivially empty constructor for the parent class
        """        
        self.name = "Filter Parent Class"
        
    def apply(self, ORF):
        """
        Parameters: an ORF object that is filled with peptides 
        Return: don't know
        Description: this is an empty function for the parent class
        """
        return


class UniquenessFilter(Filter):
    """Class UniquenessFilter: this is an ORF level filter for proteogenomcis,
    and works to get rid of ORFs that do not have any uniquely mapping peptides
    The problem with ts is not tat we doubt the reality of their
    peptide matches, but that without any uniquely mapping peptides, we cannot
    say for certain whether it was this locus or another that needs the 
    attention.
    Rules for Filter:
    Rule 1. If there are no unique peptides, then we signal for deleting 
    the ORF.
    
    """
    def __init__(self):
        """
        Parameters: none
        Return: none
        Description: a rather small constructor
        """               
        self.name = "UniquenessFilter"
        

    def apply(self, ORF):
        """
        Parameters: an ORF object that is filled with peptides 
        Return: 0/1 keep/destroy
        Description: Apply the filter. Check to see if there are any
        unique peptides mapping within the ORF. If not, then delete
        the ORF
        """
        Save = 0
        for Peptide in ORF.PeptideLocationList:
            if Peptide.Unique: # == 1 
                Save = 1
                break
        if Save:
            return 0 #keep me around
        return 1 # delete me NOW

class TrypticFilter(Filter):
    """Class TrypticFilter: this is an ORF level filter for proteogenomcis,
    and works to get rid of ORFs that do not have any tgryptic peptides
    If trypsin is used to generate the peptides in the experiment, then
    we expect that there are SOME tryptic peptides.  We know that not all
    of them will be tryptic, and for some very abundant proteins, the 
    tryptic peptides may be in the minority, but they will be present.
    Rules for Filter:
    Rule 1. If there are no tryptic peptides, then we signal for deleting 
    the ORF.
    
    NOTE: THIS WILL NOT WORK YET, as the peptide object does not have
    the tryptic variable. But it will soon
    
    """
    def __init__(self):
        """
        Parameters: none
        Return: none
        Description: a rather small constructor
        """               
        self.name = "TrypticFilter"
        

    def apply(self, ORF):
        """
        Parameters: an ORF object that is filled with peptides 
        Return: 0/1 keep/destroy
        Description: Apply the filter. Check to see if there are any
        tryptic peptides mapping within the ORF. If not, then delete
        the ORF
        """
        Save = 0
        for Peptide in ORF.PeptideLocationList:
            if Peptide.Tryptic: # == 1 
                Save = 1
                break
        if Save:
            return 0 #keep me around
        return 1 # delete me NOW


    
class SequenceComplexityFilter(Filter):
    """Class SequenceComplexityFilter: this is an ORF level filter for 
    proteogenomics, and works to get rid of ORFs who are represented by
    peptides with an uncharacteristically low sequence complexity.  I've
    found a lot of ORFs with peptides like these: GGGGAAA, GGGAGGR, GGAGAGA,
    DPAAAAR, AAGGAPP, etc. Lots of small mol weight amino acids.  I believe
    that these are random annotations made possible by the small MW making
    it easy to match a lot of peaks.  
    Rules for Filter: 
    Rule 1. I have taken a survey of proteomics data on hand, and
    found that the average normalized MW per peptide is XXX, with Stdev XXX
    If I find an ORF that hasmostly G,A,S then the it will be much lower than 
    the average.  I will do something like - if you have 90% G or A, then delete
    Maybe also look for those way outside the stdev. we'll see about that later
    Rule 2: if the sequence is entirely G,A, then I am going to just delete the
    peptide, and let it float through the rest of the filters
    
    Variables:
        self.name
    Functions:
        apply (ORF)
    """
    def __init__(self):
        """
        Parameters: none
        Return: none
        Description: a rather small constructor
        """               
        self.name = "SequenceComplexityFilter"
        self.MaximumLowMWContent = 0.9
        self.LowMW = ["G", "A"]
        self.GA = ["G", "A"]
    def apply(self, ORF):
        """
        Parameters: an ORF object that is filled with peptides 
        Return: don't know
        Description: Apply the filter. Look for sequence complexity
        that is outside of the acceptable range.
        """
        self.RemoveExclusivelyLowMWPeptides(ORF)
        return
        #1. Get a big string of all the peptides in the ORF
        PeptideString = ""
        for PeptideObject in ORF.PeptideLocationList:
            PeptideString += PeptideObject.Aminos
        #2. Count the small aminos
        Count = 0
        for Letter in PeptideString:
            if Letter in self.LowMW:
                Count +=1
        Normalized = Count / float (len(PeptideString))
        if Normalized > self.MaximumLowMWContent:
            return 1 # delete ME
        return 0 #keep ME
            
    def RemoveExclusivelyLowMWPeptides(self, ORF):
        """
        Parameters: an ORF object that is filled with peptides
        Return: None
        Description: This looks for peptides that are just
        LowMW residues and removes them entirely.  They are crappy
        annotations and I don't want to see them!
        """
        Count = len(ORF.PeptideLocationList)
        DeleteList = []
        for Index in range (Count):
            PeptideObject = ORF.PeptideLocationList[Index]
            Saved = 0
            #now cycle through and see if it's only G and A
            for Letter in PeptideObject.Aminos:
                if not Letter in self.GA:
                    Saved = 1
                    break
            if not Saved:
                DeleteList.append(Index)
        #now we have to go through the peptides in reverse order when deleting, so that
        #we dont' screw up our numbering, if we delete index=1, what was at index2 will now be at index1
        DeleteList.reverse()
        for Index in DeleteList:
            del ORF.PeptideLocationList[Index] 
        


MinimumPeptidesPerORF = 2

def FilterThisORF(ORF):
    """Parameters: ORF object from GenomicLocations.py
    Return: none
    Description: this is the main entry for filtering.  it does all the 
    rerouting and figuring things out. If it finds peptides to filter
    it will delete them, then you check afterwards if the ORF now contains
    any peptides
    """
    ## 1. Meet the minimum peptide requirement
    if len(ORF.PeptideLocationList) < MinimumPeptidesPerORF:
        del ORF.PeptideLocationList[:] #remove all peptides
        return 
    #now we start the more complicated filters, which I will do later.

def FindOverlappingDubiousGeneCalls(Dictionary, MaxOverlap):
        """Parameters: none
        Return: list of offending genes
        Description: There are genomic regions for which two gene calls overlap. 
        For some badly predicted genomes, the overlap is substantial (like >50 bp).
        I believe that most of these are bad gene calls, and I want to filter them 
        out.  This method finds such cases and calls them to attention.
        """
        #1. first sort the proteins by their start and stop
        TwoLayerTuples = Dictionary.items()
        DictLen = len(Dictionary)
        print "I foiund %s items in the dictionary"%DictLen
        TwoLayerTuples.sort(lambda (k1, (b1,e1)), (k2, (b2,e2)): SortStartStop(b1,e1,b2,e2))
        Overlappers = []
        # I now have a sorted list that looks like this
        # (protein1 (start, stop)), (protein2 (start, stop))...
        #so now I look for overlaps, is start2< stop1.  If so, by like 50 bases?
        Len = len(TwoLayerTuples)
        for Index in range(Len):
            P1 = TwoLayerTuples[Index]
            (P1Name, (P1Start, P1Stop)) = P1
            #print "P1 %s, %s"%(P1Start, P1Name)
            for Jndex in range(Index+1, Len):
                P2 = TwoLayerTuples[Jndex]
                #now do the compares.  Do P1 and P2 overlap?
                (P2Name, (P2Start, P2Stop)) = P2
                #print "Comparing \n\t%s,%s \n\t%s,%s"%(P1Start, P1Name, P2Start, P2Name)
                if P1Stop > P2Start:
                    #this is overlap, now get the amount of overlap
                    Overlap = P1Stop - P2Start
                    if Overlap > MaxOverlap:
                        print "%s\t%s\t%s"%(Overlap, P1Name, P2Name)
                        if not P1Name in Overlappers:
                            Overlappers.append(P1Name)
                        if not P2Name in Overlappers:
                            Overlappers.append(P2Name)
                    #if they do overlap, then it's possible (although crazy) that another P2 will also overlap
                    #so don't break out 
                else:
                    #there is no overlap with this, and it's a sorted list, so time to advance p1
                    break
                
        return Overlappers
                
                    
                    
                 
                
                
            
        
def SortStartStop(Begin1, End1, Begin2, End2):
    if Begin1 > Begin2:
        return 1
    elif Begin1 == Begin2:
        #would normally return 0 here, but we then sort by the end coord
        if End1 > End2:
            return 1
        else: #we tacitly assume that there are no things with the same start and stop
            return -1
    else:
        #begin1 < Begin2
        return -1
