"""PGCleavageAnalysis.py
This auxiliary set of functions is to be used to analyze peptide
cleavage.  We take as input a GenomicLocationForORF object that 
is fully populated (with peptides and a amino acid sequence).  So far we are characterizing

1. Covered peptides
2. Covered cleavage sites
3. Abutting peptides

Convention: a cleavage break point (cut point) is listed by the 
preceding amino acid index.  e.g. for the peptide "OTEINS" on the 
protein sequence 

123456789012345678
MAPROTEINSEQUENCE*

the peptide maps to (5-10) inclusive and the two recorded
cleavage points are index 4 and 10.  index 4 is considered
to follow tryptic proteolysis rules, but 10 is not.  The hash
we use for a container of peptide locations looks something like
(5,10) -> (Y, N)
where key is the peptide location, and the value is whether that
break follows proteolysis rules.

NOTE: this is a utility, and not executable from the command line

"""

## these cuts may be a bit confusing.  but here's the way we define them
## for the residues on either side of the cut, we state if they are on the
##n-term side of the cut, or c-term.  e.g. trypsin cuts after lysine, leaving
##lysine of the n-term half of the cut.  so we define it as n-term specificity
## warning, these are very simplistic enzymes, no motif, just single letter
## enzymes, like trypsin

NTermCuts = {}
NTermCuts["Trypsin"] = ["K", "R"]  

CTermCuts = {}
CTermCuts["Trypsin"] = [] #non-specific to the c-term residue


def Analysis(ORFObject, ProteolysisEnzymes):
    """
    Parameters: GenomicLocationForORF object
    Return: ???
    Description: do some analysis on the peptide cleavage, see note at
    top of the file.
    """
    if not ORFObject.ProteinPrediction: ## temporary maybe.  I should not be looking
        ## at things without a protein prediction, or I should just get the protein
        ## sequence of the whole ORF.  Another decision postponed.
        print "the below ORF has no protein"
        ORFObject.PrintMe(0,1)
    else:
        MapPeptidesIntoProteinSpace(ORFObject, ProteolysisEnzymes)



def MapPeptidesIntoProteinSpace(ORFObject, ProteolysisEnzymes):
    """Parameters: GenomicLocationForORF object, List of proteolysis enzymes
    Return: Mapping Dictionary
    Description: the ORF object, as typically populated will map peptides
    onto their nucleotide position on the chromosome (or plasmid). But for
    us it is very, very useful to get these mapped onto the amino acid 
    sequence space.
    """
    ProteinSequence = ORFObject.ProteinPrediction.Aminos
    MappingDictionary = {} #key = (startamino, stopamino) value = (follows proteolytic rule n, and c)
    for PeptideObject in ORFObject.PeptideLocationList:
        Aminos = PeptideObject.Aminos
        StartResidue = ProteinSequence.find(Aminos)
        StopResidue = StartResidue + len(Aminos)
        ## now we look to understand the specificity of the cut points
        ## first we look at the front of the peptide
        StartCutFollowsRules = "N" #for no
        if StartResidue == 0: # does it start with the protein start
            StartCutFollowsRules = "Y"
        else:
            StartNResidue = ProteinSequence[StartResidue -1] #-1 to get the n side of the cut
            AcceptableNResidues = [] 
            for Enzyme in ProteolysisEnzymes:
                AcceptableNResidues.extend(NTermCuts[Enzyme])
            #now check
            if StartNResidue in AcceptableNResidues:
                StartCutFollowsRules = "Y"
            #here's space for doing the c-term residue, but I don't want to do that now.  rush, rush ASMS
        StopCutFollowsRules = "N"
        if StopResidue == (len(ProteinSequence) -1): #does it end at the protein end
            StopCutFollowsRules = "Y"
        else:
            StopNResidue = ProteinSequence[StopResidue -1] #-1 to get the n side of the cut
            AcceptableNResidues = [] 
            for Enzyme in ProteolysisEnzymes:
                AcceptableNResidues.extend(NTermCuts[Enzyme])
            #now check
            if StopNResidue in AcceptableNResidues:
                StopCutFollowsRules = "Y"
            #here's space for doing the c-term residue, but I don't want to do that now.  rush, rush ASMS
        #now put in the hash
        Key = (StartResidue, StopResidue)
        Value = (StartCutFollowsRules, StopCutFollowsRules)
        MappingDictionary[Key] = Value
            
    
    
    
    
    
