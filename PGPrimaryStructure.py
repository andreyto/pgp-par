"""PGPrimaryStructure.py
THis is a set of classes designed to help analyze the primary
structure of proteins, and whether they agree with the proteomic
evidences
 
NOTE: this is a utility, and not executable from the command line

"""

class PrimaryStructure:
    """Class PrimaryStructure: This is an analysis object to help
    us determine whether the predicted protein is in harmony with
    the peptides observed in proteomics experiments
    Variables:
        -
    Functions: CheckStructure(PGPeptide.OpenReadingFrame)
    """
    def __init__(self):
        """Parameters: none 
        Return: none
        Description: trivial constructor 
        """        
        pass
        
    def CheckStructure(self, ORF):
        """
        Parameters: a PGPeptide.OpenReadingFrame object
        Return: ?
        Description: various QA/QC checks
        """
        Novel = self.IsItNovel(ORF)
        if Novel:
            return "NOVEL"
        UnderPredicted = self.IsItUnderPredicted(ORF)
        if UnderPredicted:
            return "UNDERPREDICTED"
        
    def IsItNovel(self, ORF):
        """
        Parameters: a PGPeptide.OpenReadingFrame object
        Return: true/false
        Description: check to see if this ORF has peptides but not a 
        predicted protein
        
        NOTE: in the old version of this, I had a minimum ORF length.  Should
        I put that in here? I'm not inclined yet, ,mostly because I don't 
        remember why I put that in there, and don't have any really solid
        research compelling me to do so.
        """
        PredictedProtein = ORF.GetLocatedProtein()
        if PredictedProtein:
            return False #there is something here, so no it's not novel
        #well, if we get here, then we have something interesting
        #let's try and get the observed sequence.  That's firstpeptide->stop
        ObservedSequence = ORF.GetObservedSequence()
        print "I got this observed sequence from %s\n%s\n\n"%(ORF, ObservedSequence)
        return True
        
    def IsItUnderPredicted(self, ORF):
        """
        Parameters: a PGPeptide.OpenReadingFrame object
        Return: true/false
        Description: Check to see if this ORF has peptides upstream 
        of the currently predicted start (at least two amino acids
        upstream). We should also include whether there is a start 
        site upstream to use
        """
        #1. we get the nuc of the 5' peptide
        FirstObservedPeptide = ORF.GetFivePrimePeptide()
        FirstObservedNucleotide = FirstObservedPeptide.GetFivePrimeNucleotide()
        StartCodon = ORF.GetNucleotideStartOfTranslation()
        Strand = ORF.GetStrand()
        
        #2. Strand switch that we use all over the place
        if Strand == "+":
            if FirstObservedNucleotide + 3 < StartCodon: # do the plus three
                #because we want more than a single amino acid upstream.
                UpstreamExtent = StartCodon - FirstObservedNucleotide
                print "Peptide %s is %s bases upstream of protein in %s\n\n"%(FirstObservedPeptide, UpstreamExtent, ORF)
                return True
        else:
            if FirstObservedNucleotide > StartCodon + 3:
                UpstreamExtent = FirstObservedNucleotide - StartCodon
                print "Peptide %s is %s bases upstream of protein in %s\n\n"%(FirstObservedPeptide, UpstreamExtent, ORF)
                return True
                
        #3. Is the actual start codon observed as a L or V, instead of M
        if FirstObservedNucleotide == StartCodon:
            FirstObservedAminoAcid = FirstObservedPeptide.aminos[0]
            if not FirstObservedAminoAcid == "M":
                print "Peptide %s is at the start codon of protein %s\n\n"%(FirstObservedPeptide, ORF)
                return True
        