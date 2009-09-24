"""
Given the output of a splice-tolerant search, find all peptides that aren't present in a
reference "non-spliced" database.
"""

import os
import sys
import string
from Utils import *
Global.FixedMods = {"C":57.0518}
Initialize()
MQCutoff = 0.5

MatchScores = {}
for AA in "ACDEFGHIKLMNPQRSTVWY":
    for AA2 in "ACDEFGHIKLMNPQRSTVWY":
        MatchScores[(AA, AA2)] = -0.5
    MatchScores[(AA, AA)] = 1
##    MatchScores[(AA, "*")] = -99
##    MatchScores[("*", AA)] = -99
##    MatchScores[(AA, "X")] = -99
##    MatchScores[("X", AA)] = -99
##    MatchScores[(AA, "O")] = -99
##    MatchScores[("O", AA)] = -99
##    
MatchScores[("I", "L")] = 0.95
MatchScores[("L", "I")] = 0.95
MatchScores[("Q", "K")] = 0.80
MatchScores[("K", "Q")] = 0.80

def PrintClosestHomolog(Database, Peptide):
    ScoreTable = {}
    DirTable = {}
    BestScore = -99
    BestPos = None
    for Y in range(0, len(Database)):
        DBAA = Database[Y]
        for X in range(0, len(Peptide)):
            PepAA = Peptide[X]
            AAScore = MatchScores.get((DBAA, PepAA), -99)
            BackScore = ScoreTable.get((X-1, Y-1), None)
            if BackScore > 0:
                DirTable[(X, Y)] = 1
                ScoreTable[(X, Y)] = AAScore + BackScore
            else:
                DirTable[(X, Y)] = 0
                ScoreTable[(X, Y)] = AAScore
            if ScoreTable[(X, Y)] > BestScore:
                BestScore = ScoreTable[(X, Y)]
                BestPos = (X, Y)
    (BestX, BestY) = BestPos
    StrA = ""
    StrB = ""
    StrC = ""
    MassDelta = 0
    for Y in range(BestY-X-10, BestY+(len(Peptide)-X)+10):
        if (Y>=0 and Y < len(Database)):
            DBAA = Database[Y]
        else:
            DBAA = " "
        StrA += DBAA
        X = Y + (BestX-BestY)
        if (X>=0 and X<len(Peptide)):
            PepAA = Peptide[X]
        else:
            PepAA = " "
        StrB += PepAA
        if (DBAA!=" " and DBAA == PepAA):
            StrC += "*"
        elif (DBAA!=" " and PepAA!=" "):
            StrC += "-"
            MassDelta += Global.AminoMass.get(PepAA,0) - Global.AminoMass.get(DBAA,0)
        else:
            StrC += " "
    #print "%s matched with score %s"%(Peptide, BestScore)
    print StrA
    if MassDelta:
        print StrC, "Mass delta %s"%MassDelta
    else:
        print StrC
    print StrB , "Match score %s"%BestScore

def GetOldAnnotations():
    Dict = {}
    #Dir = r"E:\ms\PeptideAtlas\A8_IP_searched"
    #Dir = r"E:\ms\PeptideAtlas\lgmem2_searched"
    Dir = r"E:\ms\PeptideAtlas\ICAT-3B-S1_searched"
    for FileName in os.listdir(Dir):
        if os.path.splitext(FileName)[1].lower() != ".html":
            continue
        print "Read legacy annotations from %s..."%FileName
        Path = os.path.join(Dir, FileName)
        File = open(Path, "r")
        for FileLine in File.xreadlines():
            Pos = FileLine.find("HREF=")
            if Pos == -1:
                continue
            Pos2 = FileLine.find("\">", Pos+1)
            Slashy = FileLine[Pos:Pos2].split("/")
            Dotty = Slashy[-1].split(".")
            Key = (Dotty[0]+".ms2", int(Dotty[1]))
            Pos2 = FileLine.rfind("</A>")
            Pos = FileLine.rfind(">", 0, Pos2 - 1)
            Peptide = FileLine[Pos+1:Pos2]
            #print Pos, Pos2, Peptide
            if Dict.has_key(Key):
                Dict[Key].append(Peptide)
            else:
                Dict[Key] = [Peptide]
    return Dict
    #      1  <A TARGET="Win1" HREF="/cgi-bin/sequest-tgz-out.cgi?OutFile=/regis/sbeams/archive/nking/Colorado/A8_IP/HsIPI_v3.01/A8_013_07_X08_J5/./A8_013_07_X08_J5.0012.0012.1.out">./A8_013_07_X08_J5.0012.0012.1</A>  1320.4 (-0.7)  0.6628  0.004     68.6    2  <A TARGET="Win1" HREF="/cgi-bin/sequest-tgz-plot.cgi?Dta=/regis/sbeams/archive/nking/Colorado/A8_IP/HsIPI_v3.01/A8_013_07_X08_J5/A8_013_07_X08_J5.0012.0012.1.dta&amp;MassType=0&amp;NumAxis=1&amp;DMass2=16.000000&amp;MassC=160.1390&amp;Pep=VLHAAASAGEHEK">  5/ 24</A>  <A TARGET="Win1" HREF="/cgi-bin/comet-fastadb.cgi?Ref=IPI00178203&amp;Db=/dbase/IPI/ipi.HUMAN.fasta.v3.01&amp;NucDb=0&amp;Pep=VLHAAASAGEHEK&amp;MassType=0">IPI00178203</A>  <A TARGET="Win1" HREF="/cgi-bin/comet-fastadb.cgi?Db=/dbase/IPI/ipi.HUMAN.fasta.v3.01&amp;NucDb=0&amp;Pep=VLHAAASAGEHEK&amp;MassType=0">+1</A>  R.<A TARGET="Win1" HREF="http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Web&amp;LAYOUT=TwoWindows&amp;AUTO_FORMAT=Semiauto&amp;ALIGNMENTS=50&amp;ALIGNMENT_VIEW=Pairwise&amp;CDD_SEARCH=on&amp;CLIENT=web&amp;COMPOSITION_BASED_STATISTICS=on&amp;DATABASE=nr&amp;DESCRIPTIONS=100&amp;ENTREZ_QUERY=(none)&amp;EXPECT=1000&amp;FILTER=L&amp;FORMAT_OBJECT=Alignment&amp;FORMAT_TYPE=HTML&amp;I_THRESH=0.005&amp;MATRIX_NAME=BLOSUM62&amp;NCBI_GI=on&amp;PAGE=Proteins&amp;PROGRAM=blastp&amp;SERVICE=plain&amp;SET_DEFAULTS.x=41&amp;SET_DEFAULTS.y=5&amp;SHOW_OVERVIEW=on&amp;END_OF_HTTPGET=Yes&amp;SHOW_LINKOUT=yes&amp;QUERY=VLHAAASAGEHEK">VLHAAASAGEHEK</A>.V                                    


def Main():
    OldAnnotations = GetOldAnnotations()
    DatabaseFile = open(sys.argv[2], "r")
    RawDatabase = DatabaseFile.read()
    DatabaseFile.close()
    Database = RawDatabase.replace("I","L").replace("Q","K")

    Peptides = {} # from sequence to a list of spectra
    OutputFile = open(sys.argv[1], "r")
    OldSpectrum = None
    PendingPeptide = None
    DBExactMatchFound = 0
    LineNumber = 0
    for FileLine in OutputFile.xreadlines():
        LineNumber += 1
        if LineNumber % 1000 == 0:
            print LineNumber
        Bits = FileLine.strip().split("\t")
        Spectrum = (Bits[0], Bits[1])
        if Spectrum != OldSpectrum:
            # Finish old spectrum, maybe:
            if PendingPeptide and not DBExactMatchFound:
                if not Peptides.has_key(PendingPeptide):
                    Peptides[PendingPeptide] = []
                Peptides[PendingPeptide].append(PendingPeptideLine)
            Spectrum = OldSpectrum
            PendingPeptide = None
            DBExactMatchFound = 0
            
        Bits = FileLine.split("\t")
        try:
            MQScore = float(Bits[5])
        except:
            MQScore = -9999
        if MQScore < MQCutoff:
            continue
        Aminos = ""
        for Char in Bits[2][2:-2]:
            if Char in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
                Aminos += Char
        if RawDatabase.find(Aminos)!=-1:
            DBExactMatchFound = 1
            continue
        if not PendingPeptide:
            PendingPeptide = Aminos
            PendingPeptideLine = FileLine.strip()
        
    # Finish old spectrum, maybe:
    if PendingPeptide and not DBExactMatchFound:
        if not Peptides.has_key(PendingPeptide):
            Peptides[PendingPeptide] = []
        Peptides[PendingPeptide].append(PendingPeptideLine)

    SortedList = []
    for Peptide in Peptides.keys():
        SortedList.append((len(Peptides[Peptide]), Peptide))
    SortedList.sort()
    SortedList.reverse()
    for (SpectrumCount, Peptide) in SortedList:
        # Check to see that there's at least one line with a *good* MQScore
        # Also, check to see if the peptide's already been matched:
        MQScore = -99
        AlreadyMatched = 0
        for Line in Peptides[Peptide]:
            Bits = Line.split("\t")
            MQScore = max(MQScore, float(Bits[5]))
            Key = (Bits[0].split("\\")[-1], int(Bits[1]))
            Legacy = OldAnnotations.get(Key, [])
            if Bits[2][2:-2] in Legacy:
                AlreadyMatched = 1
                break
        if MQScore < 1.8:
            continue  
        if AlreadyMatched:
            continue
        print "Peptide %s matched for %s spectra."%(Peptide, SpectrumCount)
        PrintClosestHomolog(RawDatabase, Peptide)
        for Line in Peptides[Peptide]:
            Bits = Line.split("\t")
            Key = (Bits[0].split("\\")[-1], int(Bits[1]))
            Legacy = OldAnnotations.get(Key, [])
            NewBits = [Bits[0], Bits[20], Bits[1], Bits[2], str(Legacy), Bits[3], Bits[4], Bits[5]]
            print string.join(NewBits, "\t")
            #Bits.extend(OldAnnotations.get(Key, []))
            #print Line


if __name__ == "__main__":
    import psyco
    psyco.full()
    Main()

##    File = open("Database\\NewISB.trie", "r")
##    DB = File.read()
##    File.close()
##    PrintClosestHomolog(DB, "VKEAMAPK")
##    print 
##    PrintClosestHomolog(DB, "VKEAMGPK")
    #Main()