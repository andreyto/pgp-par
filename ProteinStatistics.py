"""These are for functions on protein sequences.  Right now it's only hydropathy, but maybe more will come up
"""
import BasicStats
import Global
import string
import math

class HydrophobicIndex:
    """scores Kyte-Doolittle, grouping a la Lehninger"""
    HI = {}
    #NonPolar 
    HI["G"] = -0.4 #Glycine
    HI["A"] = 1.8 #Alanine
    HI["P"] = -1.6 #Proline  #### Typo in Lenhinger!!!
    HI["V"] = 4.2 #Valine
    HI["L"] = 3.8 #Leucine
    HI["I"] = 4.5 #Isoleucine
    HI["M"] = 1.9 #Methionine
    #Aromatic 
    HI["F"] = 2.8 #Phenylalanine
    HI["Y"] = -1.3 #Tyrosine
    HI["W"] = -0.9 #Tryptophan
    #polar, uncharged 
    HI["S"] = -0.8 #Serine
    HI["T"] = -0.7 #Threonine
    HI["C"] = 2.5 #Cysteine
    HI["N"] = -3.5 #Asparagine
    HI["Q"] = -3.5 #Glutamine
    #polar, charged 
    HI["K"] = -3.9 #Lysine
    HI["H"] = -3.2 #Histidine
    HI["R"] = -4.5 #Arginine
    HI["D"] = -3.5 #Aspartate
    HI["E"] = -3.5 #Glutamate
    
class HyrdopathyPlot:
    def __init__(self):
        self.Metric = HydrophobicIndex
    
    def MakePlot(self, Sequence, BinLen = 5):
        """given a sequence, do the plot.  you can supply a bin len, but it's your
        job to figure out what a good length is. I read somewhere 5-7.
        """
        EffectiveLength = len(Sequence) - BinLen + 1 #plus one to get the last bin
        Plot = []
        for Index in range(EffectiveLength):
            #now look at the kmer starting at Index
            Values = []
            for Jndex in range(Index, Index+BinLen):
                Letter = Sequence[Jndex]
                if not Letter in (string.uppercase):
                    print "I'm barfing on %s"%Sequence
                    continue
                Values.append( self.Metric.HI[Letter])
            Mean = BasicStats.GetMean(Values)
            Plot.append(Mean)
        return Plot

    def IsConsistentSignal(self, List, Min = 0.5, Len = 10):
        """I am trying to see if here is a consistent hydrophobic patch in this
        plot.  So for standards, we count positive values, at least 10 in a row
        """
        Counter  = 0
        Index = -1 #start at -1 so I can increment right off the bat
        StartIndex = -1 #start of the hydrophobic patch
        for Item in List:
            Index += 1 #index into the array
            
            if Item > Min:
                if not Counter:
                    #this is the first one in a row
                    StartIndex = Index
                Counter += 1
                if Counter >= Len:
                    return StartIndex
            else:
                Counter =0
                StartIndex =-1 
        return -1 # not found




def GetMW(Aminos):
    Aminos = Aminos.upper()
    MW = 0
    for Amino in Aminos:
        MW += Global.AminoMass.get(Amino, 0)
    MW += 19
    return MW

def HasSignalPeptidaseMotif(Aminos):
    """Looking for the end of the peptide (the prefix to the observed protein) to have
    AxB.. where A in [ILVAGS] and B in [AGS]
    """
    APosition = Aminos[-3]
    BPosition = Aminos[-1]
    AcceptableA = ["I","L","V","A","G","S"]
    AcceptableB = ["A","G","S"]
    if not APosition in AcceptableA:
        return 0
    if not BPosition in AcceptableB:
        return 0
    return 1

def HasBasicResidue(Sequence, Start = 0, End = None):
    """
    return 0/1 if there is a basic residue in the sequence given,
    and the bracketed subsequence
    """
    SubSequence = Sequence
    if End:
        SubSequence = Sequence[Start:End]
    else:
        SubSequence = Sequence[Start:]
    Basic = ["R", "K"]
    for Letter in SubSequence:
        if Letter in Basic:
            return 1
        
    return 0

class SequenceComplexity:
    """this class is meant to measure a growing number of sequence complexity traits
    in protein sequences.
    """
    
    def __init__(self):
        pass

    
    def CheckComplexity(self, genome):
        if self.Verbose:
            print "ProteogenomicsPostProcessing.py:CheckComplexity"
        LowMW = ["G", "A"]
        List = []
        for (chromName, chrom) in genome.chromosomes.items():
            for (orfName,ORF) in chrom.simpleOrfs.items() + chrom.pepOnlyOrfs.items():
        
                PeptideString = ""
                for PeptideObject in ORF.PeptideLocationList:
                    PeptideString += PeptideObject.Aminos
                Count = 0
                for Letter in PeptideString:
                    if Letter in LowMW:
                        Count +=1
                Normalized = Count / float (len(PeptideString))
                List.append(Normalized)
        Histogram = BasicStats.CreateHistogram(List, 0, 0.05)
        BasicStats.PrettyPrintHistogram(Histogram, None)


    def SequenceComplexityHack(self, genome):
        """hack"""
        
        ORFList = []
        PeptideList = []
        DiffList = []
        RealProteinList = []
        for (chromName, chrom) in genome.chromosomes.items():
            for (orfName,ORF) in chrom.simpleOrfs.items() + chrom.pepOnlyOrfs.items():
        
                #we try for the predicted protein first
                Sequence =  ORF.GetProteinSequence()
                if Sequence:
                    Entropy = self.SequenceEntropy(Sequence)
                    RealProteinList.append(Entropy)
                if not Sequence:
                    Sequence = ORF.GetObservedSequence()
                ProteinEntropy = self.SequenceEntropy(Sequence)
                #now do for the peptides
                PeptideCat = ""
                for Peptide in ORF.peptideIter():
                    PeptideCat += Peptide.GetAminos()
                PeptideEntropy = self.SequenceEntropy(PeptideCat)
                
                #put stuff in lists
                ORFList.append(ProteinEntropy)
                PeptideList.append(PeptideEntropy)
                Diff = ProteinEntropy - PeptideEntropy
                DiffList.append(Diff)
            
        ORFHandle = open("ORFEntropy.txt", "wb")
        ORFLine = "\t".join(map(str, ORFList))
        ORFHandle.write(ORFLine)
        ORFHandle.close()
    
        PeptideHandle = open("PeptideEntropy.txt", "wb")
        Line = "\t".join(map(str, PeptideList))
        PeptideHandle.write(Line)
        PeptideHandle.close()
        
        DiffHandle = open("DiffEntropy.txt", "wb")
        Line = "\t".join(map(str, DiffList))
        DiffHandle.write(Line)
        DiffHandle.close()
        
        RealHandle = open("RealProteinEntropy.txt", "wb")
        Line = "\t".join(map(str, RealProteinList))
        RealHandle.write(Line)
        RealHandle.close()
        
   
    def SequenceEntropy(self, Sequence):
        """Parameters: An amino acid sequence
        Return: the H(x) entropy
        Description: Use the classic information entropy equation to calculate
        the entropy of the input sequence.
        H(x) = SUM p(xi) * log(1/ p(xi))
        xi = letter of the sequence
        e.g. GGGAS
        x1 = G, p(G) = 3/5
        x2 = A, p(A) = 1/5
        X3 = S, p(S) = 1/5
        """
        ProbTable = self.GetProbabilityTable(Sequence)
        Sum =0 
        for (Letter, Probability) in ProbTable.items():
            LogValue = math.log(1 / Probability) #currently the natural log.  not sure the base of the log matters
            Sum += (Probability * LogValue)
        return Sum

    def GetProbabilityTable(self, Sequence):
        """Parameters: an amino acid sequence
        Return: a dictionary of probability (frequence/n) for each letter
        Description: just convert counts to probability.  easy.
        """
        CountDict = {}
        ProbabilityDict = {}
        for Letter in Sequence:
            if not CountDict.has_key(Letter):
                CountDict[Letter] = 0 #initialize
            CountDict[Letter] += 1
        Len = float(len(Sequence)) #cast to float so we can do real division
        for (Key, Value) in CountDict.items():
            Probability = Value / Len
            ProbabilityDict[Key] = Probability
        return ProbabilityDict
            

