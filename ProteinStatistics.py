"""These are for functions on protein sequences.  Right now it's only hydropathy, but maybe more will come up
"""
import BasicStats
import Global
import string

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
        Found = 0
        for Item in List:
            if Item > Min:
                Counter += 1
                if Counter >= Len:
                    Found =1
            else:
                Counter =0 
        return Found


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
