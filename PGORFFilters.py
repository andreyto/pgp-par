"""PGORFFilters.py
This auxiliary set of functions is to be used to filter ORFs (not proteins.  
open reading frames).  The filters are 

1. 2 peptides/ORF
2. Clustering
3. Sequence content
4. Tryptic content of cluster

NOTE: this is a utility, and not executable from the command line

"""


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
