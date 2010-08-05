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

import ProteinStatistics

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
    return EvaluateSignalPeptide(ORFObject)
        



def MapPeptidesIntoProteinSpace(ORFObject, ProteolysisEnzymes):
    """Parameters: GenomicLocationForORF object, List of proteolysis enzymes
    Return: Mapping Dictionary
    Description: the ORF object, as typically populated will map peptides
    onto their nucleotide position on the chromosome (or plasmid). But for
    us it is very, very useful to get these mapped onto the amino acid 
    sequence space.
    """
    
    """ #############DEPRECATED>  DO I USE THIS ??????????? #####################"""
    
    ProteinSequence = ORFObject.annotatedProtein.Aminos
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
            
def GetSignalPeptideHeader():
    """ Return a string header corresponding to the output of EvaulateSignalPeptide.  
    If you change one but not the other, your name is MUD
    The header contains column names for the information columns, but not for the
    columns representing solely the hyrdophobicity score of a given index
    
    
    """
    Header = "ProteinName\tFirstPeptideIndex\tPrefixSequence\tSuffixSequence\tHasHydrophobicPatch\tHasAxBMotif\tHasBasicEarly\n"
    return Header

    
def EvaluateSignalPeptide(ORF):
    """Parameters: an openreadingframe object
    Return: Decision (Int), String of info
    Description: we look at whether a protein has evidence of a signal
    peptide.  We assume that everything passed to us HAS PEPTIDES.
    We return a decision on whether there is information supporting
    a signal peptide cleavage, and then we also return a string
    which is full of information.  Return Values for Decision:
        0 = tryptic nterminal first peptide - no evidence
        1 + non-tryptic nterminal first peptide - perhaps evidence
        1 = too short - evidence against
        2 = too long - no evidence
        3 = no patch - evidence against
        4 = patch  - evidence supporting signal peptide cleavage
    NOTE: if decision <3, then the string will be None, the NULL variable, 
    not the string 'None'
    """

    HPlot = ProteinStatistics.HyrdopathyPlot()
    FirstObservedPeptide = ORF.GetFivePrimePeptide(Unique=1)
    FirstResidue = ORF.aaseq.find(FirstObservedPeptide.aminos) #the index into the OpenReadingFrame of the first observed peptide
    AnnotatedProteinStart = ORF.ProteinAminoAcidOffset #index within the OpenReadingFrame
    
    # level 1 - is it non-tryptic on the N-terminus
    if FirstObservedPeptide.IsTrypticNterm():
        return (0, None) # this means N-term tryptic.  the none is because there is no string to return
    #Level 2 - the length filter, which is somewhat stupid, because a 
    ##protein could be mispredicted, but it's what we have to go with 
    ##and we expect signal peptides to fall within a certain length range
    PeptideOffsetIntoProtein = FirstResidue - AnnotatedProteinStart
    MaxLen = 50
    MinLen = 15
    if PeptideOffsetIntoProtein < MinLen:
        return (1, None) # this means too short
    if PeptideOffsetIntoProtein > MaxLen:
        return (2, None) #this means too long
    #meet's level 2
    PrefixSequence = ORF.GetProteinSequence(0, PeptideOffsetIntoProtein) # zero is the start.
    AfterCut = ORF.GetProteinSequence(PeptideOffsetIntoProtein, PeptideOffsetIntoProtein + 2)
    ProteinName = ORF.GetProteinName()
    #now some trickery
    Plot = HPlot.MakePlot(PrefixSequence)
    LongerSignal = HPlot.IsConsistentSignal(Plot) #default to 10
    #### WARNING  !!!!!!!!!!!!!! MAGIC NUMBER
    #in the line below there is a '3 +' which is a magic number
    #it is because the Plot variable is not the hydrophobic index of a single
    #amino acid, but a 5 residue floating window, centered over the middle.  Thus the 
    #array 'Plot' is indexed +3 into the array "PrefixSequence"
    ShorterSignalStartIndex = 3 + HPlot.IsConsistentSignal(Plot, 0.5, 8)
    ShorterSignal = 0 #assume that there is no shorter signal
    if ShorterSignalStartIndex > -1:
        ShorterSignal = 1 #reset if you get one
    AxBMotif = ProteinStatistics.HasSignalPeptidaseMotif(PrefixSequence)
    #we look for the basic residue
    HasBasicEarlyResidue = 0 # assume guilt
    if ShorterSignal:
        HasBasicEarlyResidue = ProteinStatistics.HasBasicResidue(PrefixSequence, 0, ShorterSignalStartIndex)
    #I decided to go with the shorter hydrophobic patch length, but I kept the other var in there just incase
    String = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t"%(ProteinName, PeptideOffsetIntoProtein, PrefixSequence, AfterCut, ShorterSignal, AxBMotif, HasBasicEarlyResidue)
    PlotString =""
    # need to front pad these with zeros
    PadLen = MaxLen - len(Plot)
    for i in range(PadLen):
        PlotString += "\t0"
    
    for Item in Plot:
        PlotString += "\t%s"%Item
    ReturnString = "%s%s\n"%(String, PlotString) #plotString has an explicit \t at start. .Put \n to make it completely self enclosed
    if not ShorterSignal:
        return (3, ReturnString) #lacks the hydrophobic patch
    return  (4, ReturnString) #it has a hyrophobic patch 
            


    
    
