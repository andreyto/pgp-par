"""These are for functions on nucleotide sequences.  
"""
#from DNA import Frequency
class ProteinTranslationClass:
    def __init__(self):
        """Simply the translation table, so it can be accessed
        """
        self.Table = {}
        self.Table["ATT"] = "I" #Isoleucine      I    ATT, ATC, ATA
        self.Table["ATC"] = "I"
        self.Table["ATA"] = "I"
        self.Table["CTT"] = "L" #Leucine         L    CTT, CTC, CTA, CTG, TTA, TTG
        self.Table["CTC"] = "L"
        self.Table["CTA"] = "L"
        self.Table["CTG"] = "L"
        self.Table["TTA"] = "L"
        self.Table["TTG"] = "L"
        self.Table["GTT"] = "V" #Valine          V    GTT, GTC, GTA, GTG
        self.Table["GTC"] = "V"
        self.Table["GTA"] = "V"
        self.Table["GTG"] = "V"
        self.Table["TTT"] = "F" #Phenylalanine   F    TTT, TTC
        self.Table["TTC"] = "F"
        self.Table["ATG"] = "M" #Methionine      M    ATG
        self.Table["TGT"] = "C" #Cysteine        C    TGT, TGC
        self.Table["TGC"] = "C"
        self.Table["GCT"] = "A" #Alanine         A    GCT, GCC, GCA, GCG
        self.Table["GCC"] = "A"
        self.Table["GCA"] = "A"
        self.Table["GCG"] = "A"
        self.Table["GGT"] = "G" #Glycine         G    GGT, GGC, GGA, GGG
        self.Table["GGC"] = "G"
        self.Table["GGA"] = "G"
        self.Table["GGG"] = "G"
        self.Table["CCT"] = "P" #Proline         P    CCT, CCC, CCA, CCG
        self.Table["CCC"] = "P"
        self.Table["CCA"] = "P"
        self.Table["CCG"] = "P"
        self.Table["ACT"] = "T" #Threonine       T    ACT, ACC, ACA, ACG
        self.Table["ACC"] = "T"
        self.Table["ACA"] = "T"
        self.Table["ACG"] = "T"
        self.Table["TCT"] = "S" #Serine          S    TCT, TCC, TCA, TCG, AGT, AGC
        self.Table["TCC"] = "S"
        self.Table["TCA"] = "S"
        self.Table["TCG"] = "S"
        self.Table["AGT"] = "S"
        self.Table["AGC"] = "S"
        self.Table["TAT"] = "Y" #Tyrosine        Y    TAT, TAC
        self.Table["TAC"] = "Y"
        self.Table["TGG"] = "W" #Tryptophan      W    TGG
        self.Table["CAA"] = "Q" #Glutamine       Q    CAA, CAG
        self.Table["CAG"] = "Q"
        self.Table["AAT"] = "N" #Asparagine      N    AAT, AAC
        self.Table["AAC"] = "N"
        self.Table["CAT"] = "H" #Histidine       H    CAT, CAC
        self.Table["CAC"] = "H"
        self.Table["GAA"] = "E" #Glutamic acid   E    GAA, GAG
        self.Table["GAG"] = "E"
        self.Table["GAT"] = "D" #Aspartic acid   D    GAT, GAC
        self.Table["GAC"] = "D"
        self.Table["AAA"] = "K" #Lysine          K    AAA, AAG
        self.Table["AAG"] = "K"
        self.Table["CGT"] = "R" #Arginine        R    CGT, CGC, CGA, CGG, AGA, AGG
        self.Table["CGC"] = "R"
        self.Table["CGA"] = "R"
        self.Table["CGG"] = "R"
        self.Table["AGA"] = "R"
        self.Table["AGG"] = "R"
        self.Table["TAA"] = "*" #Stop codons     *    TAA, TAG, TGA
        self.Table["TAG"] = "*" 
        self.Table["TGA"] = "*" 

class CodonCount:
    def __init__(self):
        self.Table = {}
        self.Table["ATT"] = 0 #Isoleucine      I    ATT, ATC, ATA
        self.Table["ATC"] = 0
        self.Table["ATA"] = 0
        self.Table["CTT"] = 0 #Leucine         L    CTT, CTC, CTA, CTG, TTA, TTG
        self.Table["CTC"] = 0
        self.Table["CTA"] = 0
        self.Table["CTG"] = 0
        self.Table["TTA"] = 0
        self.Table["TTG"] = 0
        self.Table["GTT"] = 0 #Valine          V    GTT, GTC, GTA, GTG
        self.Table["GTC"] = 0
        self.Table["GTA"] = 0
        self.Table["GTG"] = 0
        self.Table["TTT"] = 0 #Phenylalanine   F    TTT, TTC
        self.Table["TTC"] = 0
        self.Table["ATG"] = 0 #Methionine      M    ATG
        self.Table["TGT"] = 0 #Cysteine        C    TGT, TGC
        self.Table["TGC"] = 0
        self.Table["GCT"] = 0 #Alanine         A    GCT, GCC, GCA, GCG
        self.Table["GCC"] = 0
        self.Table["GCA"] = 0
        self.Table["GCG"] = 0
        self.Table["GGT"] = 0 #Glycine         G    GGT, GGC, GGA, GGG
        self.Table["GGC"] = 0
        self.Table["GGA"] = 0
        self.Table["GGG"] = 0
        self.Table["CCT"] = 0 #Proline         P    CCT, CCC, CCA, CCG
        self.Table["CCC"] = 0
        self.Table["CCA"] = 0
        self.Table["CCG"] = 0
        self.Table["ACT"] = 0 #Threonine       T    ACT, ACC, ACA, ACG
        self.Table["ACC"] = 0
        self.Table["ACA"] = 0
        self.Table["ACG"] = 0
        self.Table["TCT"] = 0 #Serine          S    TCT, TCC, TCA, TCG, AGT, AGC
        self.Table["TCC"] = 0
        self.Table["TCA"] = 0
        self.Table["TCG"] = 0
        self.Table["AGT"] = 0
        self.Table["AGC"] = 0
        self.Table["TAT"] = 0 #Tyrosine        Y    TAT, TAC
        self.Table["TAC"] = 0
        self.Table["TGG"] = 0 #Tryptophan      W    TGG
        self.Table["CAA"] = 0 #Glutamine       Q    CAA, CAG
        self.Table["CAG"] = 0
        self.Table["AAT"] = 0 #Asparagine      N    AAT, AAC
        self.Table["AAC"] = 0
        self.Table["CAT"] = 0 #Histidine       H    CAT, CAC
        self.Table["CAC"] = 0
        self.Table["GAA"] = 0 #Glutamic acid   E    GAA, GAG
        self.Table["GAG"] = 0
        self.Table["GAT"] = 0 #Aspartic acid   D    GAT, GAC
        self.Table["GAC"] = 0
        self.Table["AAA"] = 0 #Lysine          K    AAA, AAG
        self.Table["AAG"] = 0
        self.Table["CGT"] = 0 #Arginine        R    CGT, CGC, CGA, CGG, AGA, AGG
        self.Table["CGC"] = 0
        self.Table["CGA"] = 0
        self.Table["CGG"] = 0
        self.Table["AGA"] = 0
        self.Table["AGG"] = 0
        self.Table["TAA"] = 0 #Stop codons     *    TAA, TAG, TGA
        self.Table["TAG"] = 0
        self.Table["TGA"] = 0

def GetLength(Sequence):
    return len(Sequence)

def GetGC(Sequence):
    """Parameters: DNA sequence
    Return: integer percent GC, eg 45 for 45%
    Description: count g and c
    """
    GCCount =0
    TotalCount = len(Sequence)
    for Letter in Sequence:
        if Letter in ["G", "C", "g", "c"]:
            GCCount += 1
    Fraction = GCCount / float (TotalCount)
    Percent = int (Fraction *100)
    return Percent

def CodonUsageFractions(Sequence):
    """Parameters: Sequence of nucleotides
    Return: a 64 space vector (the alhpabetical list of codons, e.g. AAA, AAC, AAG ...)
    Description: return the fractional usage of each of the 64 codons.
    """
    #first I have to make sure that this is a simple string, otherwise things fail
    if type(Sequence) == "str":
        #this is a naked, plain string, 
        pass
    else:
        try :
            Sequence = str(Sequence)
        except:
            print "You gave me something which does not transform into a string. ABORT"
    Table = CodonCount()
    TotalCount = 0.0  #with a decimal, it's now a float and will divide correctly with out casting
    while 1:
        #print "My sequence is %s"%Sequence
        
        Codon = Sequence[:3]
        Sequence = Sequence[3:]
        if not Table.Table.has_key(Codon):
            print "SNAFU in codon counting.  No codon exists for :%s:"%Codon
        else:
            Table.Table[Codon] += 1

        TotalCount += 1
        if len(Sequence) < 3:
            break
    #now I'm done with the loop I need to normalize everything
    TempHash = {}
    for Key in Table.Table.keys():
        Count = Table.Table[Key]
        Frequency = Count / TotalCount
        TempHash[Key] = Frequency
    ToReturn = [] # a list of frequency values sorted by their codon
    Keys = Table.Table.keys()
    Keys.sort()
    for Key in Keys:
        ToReturn.append(TempHash[Key])
    return ToReturn

 
 
 
if __name__ == "__main__":
    String = "AAABBBCCCDDDEEEFFFGGGHHHIIIJJJKKK"
    CodonUsageFractions(String)
    
