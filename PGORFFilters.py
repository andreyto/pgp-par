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

