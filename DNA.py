"""
Various functions to handle DNA sequences.
"""

GeneticCode = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
               "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
               "TAT":"Y", "TAC":"Y", "TAA":"X", "TAG":"X",
               "TGT":"C", "TGC":"C", "TGA":"X", "TGG":"W",
               "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
               "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
               "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
               "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
               "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
               "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
               "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
               "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
               "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
               "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
               "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
               "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
               }
RCDict = {"A":"T", "G":"C", "T":"A", "C":"G",
          "a":"T", "g":"C", "t":"A", "c":"G"}

#A    |    x
#C A G|G U x A G U . . . . . . . . C A G|G 
#3 4 5 6 7 8 9 1011                15161718
SpliceBoundaryProfile = [
      {"A":0.3501  ,"C":0.3420  ,"G":0.1873  ,"T":0.1207},
      {"A":0.6085  ,"C":0.1094  ,"G":0.1220  ,"T":0.1601},
      {"A":0.0956  ,"C":0.0381  ,"G":0.7933  ,"T":0.0730},
      {"A":0.0028  ,"C":0.0010  ,"G":0.9949  ,"T":0.0014},
      {"A":0.0028  ,"C":0.0124  ,"G":0.0012  ,"T":0.9836},
      {"A":0.5512  ,"C":0.0304  ,"G":0.3846  ,"T":0.0338},
      {"A":0.7011  ,"C":0.0764  ,"G":0.1142  ,"T":0.1084},
      {"A":0.0787  ,"C":0.0545  ,"G":0.8026  ,"T":0.0642},
      {"A":0.1570  ,"C":0.1569  ,"G":0.1926  ,"T":0.4935},
      {"A":0.0560  ,"C":0.6400  ,"G":0.0094  ,"T":0.2947},
      {"A":0.9937  ,"C":0.0011  ,"G":0.0013  ,"T":0.0039},
      {"A":0.0011  ,"C":0.0014  ,"G":0.9943  ,"T":0.0032},
      {"A":0.2453  ,"C":0.1438  ,"G":0.4948  ,"T":0.1160},
    ]

def GetSpliceSignalScore(SpliceSignal):
    Score = 0
    for Pos in range(len(SpliceSignal)):
        Frequency = SpliceBoundaryProfile[Pos].get(SpliceSignal[Pos], 0.25)
        Score += math.log(Frequency)
    return Score
def ReverseComplement(DNA):
    "Returns the reverse complement of a DNA sequence."
    Str = ""
    for Index in range(len(DNA) - 1, -1, -1):
        Str += RCDict.get(DNA[Index], DNA[Index])
    return Str

def Translate(DNA):
    "Returns the peptide translation of a sequence."
    Peptide = ""
    for Index in range(0, len(DNA) - 2, 3):
        Codon = DNA[Index:Index+3].upper()
        AA = GeneticCode.get(Codon, "")
        Peptide += AA
        # Include a character for stop codons, and then end it:
        #if AA == "X":
        #    break
    return Peptide

def RCTranslate(DNA):
    "Take the reverse complement of a DNA sequence, then translate it"
    RCDNA = ReverseComplement(DNA)
    return Translate(RCDNA)
