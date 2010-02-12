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
        self.IsItNovel(ORF)
        
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
        NumPeptides = ORF.numPeptides() #not sure really why I care about this, because it's already passed the filters
        #let's try and get the observed sequence.  That's firstpeptide->stop
        ObservedSequence = ORF.GetObservedSequence()
        print "I got this observed sequence from %s\n%s\n\n"%(ORF, ObservedSequence)
        
