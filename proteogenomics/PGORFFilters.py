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



"""PGORFFilters.py
THis is a set of classes designed to help filter out spurious 
ORFs from the proteogenomics pipeline.  Each new filter should
be added as a new class, inheriting from 'Filter' and implementing
a method called 'filterORF' where the magic gets done.
 
NOTE: this is a utility, and not executable from the command line

"""
#import math #used to calculate logs
import os

class FilterList:
    """Class FilterList: This is a container object for a list
    of filters.  It wraps the sequential calling of each filter
    (called in order that they were added to the list).
    Variables:
        self.List
    Functions: ApplyAllFilters(ORFList)
    """
    def __init__(self, ListOfFilterObjects, OutputPath):
        """Parameters: List of filter objects, or classes that inherit from 'Filter' 
        Return: none
        Description: trivially constructor 
        """
        self.List = ListOfFilterObjects
        (Path, Ext) = os.path.splitext(OutputPath)
        FilterPath = "%s.filterreport.txt"%Path
        self.Handle = open(FilterPath, "w")
        
    def ApplyAllFilters(self, DictionaryOfORFs):
        """
        Parameters: A dict of ORF objects, where the key = name, value = ORF object
        Return: a filtered dict of ORF objects
        Description: run the filters one at a time, deleting
        the ORF objects that the filters tell me to delete
        """
        KilledPeptideCount =0
        KilledProteinCount = 0
        KillList = [] # a list of ORF.name to kill
        FiltersKillingORFs = {} #filter.name -> [list of orfs killed]
        ORFStrings = {} #orf.name -> orf.__string__
        for (Name, ORF) in DictionaryOfORFs.items():
            #now cycle through all our filters.  We use the ordering
            #of the list, meaning that you should have throught about
            #what order you wanted to filter them in when you created 
            #me as an object
            DeleteMe = 0 #assume innnocence
            CurrFilter = None
            for aFilter in self.List:
                CurrFilter = aFilter.name
                DeleteMe = aFilter.filterORF(ORF)
                if DeleteMe:
                    break #quit cycling through all the filters already, we know it sucks
            # so here's the deal.  We have created ORFs for two reasons
            #1. It contains peptides
            #2. It contains a predicted protein, which is important to keep track
            #    of when we are resolving conflicts.  So we don't really want to 
            #    kill ORFs for no reason.  In fact it would actually be better for 
            #    if we just deleted all the peptides of stuff on the list.
            if DeleteMe:
                KilledPeptideCount += ORF.numPeptides()
                KilledProteinCount += 1
                #before we nuke this thing, let's get the string version for our report
                ORFStrings[ORF.GetName()] = "%s"%ORF # uses the __string__ method we cleverly built in
                ORF.DeleteAllPeptides()
                #now after killing all peptides, we keep track of which filter killed which
                #protein
                if not FiltersKillingORFs.has_key(CurrFilter):
                    FiltersKillingORFs[CurrFilter] = [] #empty list
                FiltersKillingORFs[CurrFilter].append(ORF.GetName())
                #we only truly delete something if it has no peptides AND no protein
                if ORF.GetLocatedProtein() == None:
                    KillList.append(Name)

        for Name in KillList:
            del DictionaryOfORFs[Name]
        #now we write out to the report file
        #first the skinny
        for (FilterName, ProteinList) in FiltersKillingORFs.items():
            self.Handle.write("Filter %s eliminated %s\n"%(FilterName, len(ProteinList)))
        
        #now the long
        for (FilterName, ProteinList) in FiltersKillingORFs.items():
            self.Handle.write("Filter %s results\n"%FilterName)
            #now we want to list out all the things that got killed
            for ORFName in ProteinList:
                self.Handle.write("%s\n"%ORFStrings[ORFName])
            #one final return just to help out
            self.Handle.write("\n")
        return DictionaryOfORFs

    def GetListString(self):
        """return the string names of all the filters in this list"""
        String = ""
        for Item in self.List:
            String += " %s,"%Item.name
        
        return String


class Filter:
    """Class Filter: this is a generic filter, meant to be inherited to 
    the actual specific filter classes, e.g. SequenceComplexityFilter
    Classes that inherit from this must define - 
    Variables:
        self.name
    Functions:
        filterORF(ORF)
    """
    def __init__(self):
        """Parameters: none
        Return: none
        Description: trivially empty constructor for the parent class
        """
        self.name = "Filter Parent Class"

    def filterORF(self, ORF):
        """
        Parameters: an ORF object that is filled with peptides 
        Return: don't know
        Description: this is an empty function for the parent class
        """
        return


class UniquenessFilter(Filter):
    """Class UniquenessFilter: this is an ORF level filter for proteogenomcis,
    and works to get rid of ORFs that do not have any uniquely mapping peptides
    The problem with ts is not that we doubt the reality of their
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
        Filter.__init__(self)
        self.name = "UniquenessFilter"


    def filterORF(self, ORF):
        """
        Parameters: an ORF object that is filled with peptides 
        Return: 0/1 keep/destroy
        Description: Apply the filter. Check to see if there are any
        unique peptides mapping within the ORF. If not, then delete
        the ORF
        """
        Save = 0
        for Peptide in ORF.peptideIter():
            if Peptide.isUnique: # == 1 
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

    Peptides Trypticness should be set when peptides are mapped

    """
    def __init__(self):
        """
        Parameters: none
        Return: none
        Description: a rather small constructor
        """
        Filter.__init__(self)
        self.name = "TrypticFilter"


    def filterORF(self, ORF):
        """
        Parameters: an ORF object that is filled with peptides 
        Return: 0/1 keep/destroy
        Description: Apply the filter. Check to see if there are any
        tryptic peptides mapping within the ORF. If not, then delete
        the ORF
        """
        Save = 0
        for Peptide in ORF.peptideIter():
            if Peptide.isFullyTryptic(): # == 1 
                Save = 1
                break
        if Save:
            return 0 #keep me around
        return 1 # delete me NOW

class MinPeptideFilter(Filter):
    """Class MinPeptideFilter: this is an ORF level filter for proteogenomcis,
    and works to get rid of ORFs thathave too few peptides.  Although Pavel
    may disagree, we sometimes use this heuristic, and it's pretty effective
    We actually should want to look closely at the true significance of matches
    with only a single peptide, but that is not the purpose of this filter
    Rules for Filter:
    Rule 1. If there are fewer peptides than our cutoff, then we signal for 
    deleting the ORF.


    """
    def __init__(self, MinimumPeptides = 2):
        """
        Parameters: The minimum number of peptides that you require. defalut = 2 
        Return: none
        Description: a rather small constructor
        """
        Filter.__init__(self)
        self.name = "MinPeptideFilter"
        self.MinimumPeptides = MinimumPeptides

    def filterORF(self, ORF):
        """
        Parameters: an ORF object that is filled with peptides 
        Return: 0/1 keep/destroy
        Description: Apply the filter. Check to see if there are any
        tryptic peptides mapping within the ORF. If not, then delete
        the ORF
        """
        if ORF.numPeptides() >= self.MinimumPeptides:
            return 0 #keep me around
        return 1 # delete me NOW


class PeptideDistance(Filter):
    """Class PeptideDistance: an ORF level filter that removes peptides
    that are greater then some specified distance from another peptide
    in the ORF.
    """
    def __init__(self,maxDistance=750):
        "Parameter: the maximum distance a peptide can be from another peptide."
        Filter.__init__(self)
        self.name = "PeptideDistance"
        self.maxDistance = maxDistance

    def filterORF(self, ORF):
        """
        Parameters: an ORF object that is filled with peptides
        Return: 0 keep the ORF, this filter only deletes peptides.
        Description: Apply the filter. Check to see if there are any
        peptides greater then maxDistance away from each other, and if so
        delete the isolated peptide.
        """
        numPeptides = ORF.numPeptides()
        # Make sure the peptides are sorted by start
        prevPep = ORF.GetFivePrimePeptide()
        # The first peptide is always to far on it's left from the previous
        left2far = True
        # these next 2 are set up here so they are still in scope after the loop
        lower = 0
        prevStop = 0
        for peptide in ORF.peptideIter():
            # Skip the 1st peptide, since we're need to look at both sides
            if peptide == prevPep:
                continue
            lower = peptide.GetStart()
            prevStop = prevPep.GetStop()
            right2far = abs(lower - prevStop) > self.maxDistance
#            print "left,right2far %s,%s lower %s prevStop %s %d" % (
#                left2far,right2far,lower,prevStop,lower-prevStop)

            # if the previous peptide is 2 far from others on both sides, delete
            if left2far and right2far:
                ORF.deletePeptide( prevPep )
                numPeptides -= 1
            prevPep = peptide
            left2far = right2far

        if numPeptides > 0:
#            print "left2far %s lower %s prevStop %s %d\n%s" % (
#                left2far,lower,prevStop,lower-prevStop,ORF)

            # For the last peptide the right is always too far, so just check left
            if left2far:
                ORF.deletePeptide( prevPep )
                numPeptides -= 1

        # We want to keep the ORF, we only delete peptides so return 0
        return 0

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
        filterORF (ORF)
    """
    def __init__(self):
        """
        Parameters: none
        Return: none
        Description: a rather small constructor
        """
        Filter.__init__(self)
        self.name = "SequenceComplexityFilter"
        self.MaximumLowMWContent = 0.9
        self.LowMW = ["G", "A"]

    # declaring this as a static method so that it doesn't need an object instance
    # to work on, ie self in the arugment list
    @staticmethod
    def lowComplexFilter(PeptideObject):
        """
        Parameters: an Peptide object 
        Return: True to include peptide, False to exclude peptide from list
        Description: This looks for peptides that are just
        LowMW residues and removes them entirely.  They are crappy
        annotations and I don't want to see them!
        
        """
        Saved = 0
        GA = ["G", "A"]
        Len = len(PeptideObject.aminos)
        LenAboveReproach = 10 #total magic number.  But the short peptides are the ones giving me fits
        if Len >= LenAboveReproach:
            return True #this passes filter, because it's above reproach
        MinNonGA = Len / 3 # was /2 but I'm testing to get bad results
        Count =0
        #now cycle through and see if it's only G and A
        for Letter in PeptideObject.aminos:
            if not Letter in GA:
                Count += 1
        #now that Iv'e tallied my non-ga, see if we meet the criteria
        if Count >= MinNonGA:
            Saved =1
        if Saved:
            
            return True
        else:
            return False

    def filterORF(self, ORF):
        """
        Parameters: an ORF object that is filled with peptides 
        Return: None
        Description: Apply the filter. Look for sequence complexity
        that is outside of the acceptable range.
        """
        
        #end debugging/research
        ORF.filterPeptides( self.lowComplexFilter )
        return #this is an intentionally blank return.  We are filtering individual peptides here
        #so we don't return a call of this particular ORF or not.  
        #that is as long as we're doing the cheap hack and not the full entropy based model
        
        ###### OLD CODE for ENTROPY ########
        #some upfront debugging
        ### this should be reworked to call ProteinStatistics code which calculated entropy
        Sequence =  ORF.GetProteinSequence()
        if Sequence:
            Entropy = self.SequenceEntropy(Sequence)

        #1. Get a big string of all the peptides in the ORF
        PeptideString = ""
        for PeptideObject in ORF.peptideIter():
            PeptideString += PeptideObject.aminos
        #2. Count the small aminos
        Count = 0
        for Letter in PeptideString:
            if Letter in self.LowMW:
                Count +=1
        Normalized = Count / float (len(PeptideString))
        if Normalized > self.MaximumLowMWContent:
            return 1 # delete ME
        return 0 #keep ME



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
        print "I found %s items in the dictionary"%DictLen
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
