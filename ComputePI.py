"""
Generate a tab-delimited file specifying the PI and hydrophobicity and molecular weight
for a protein sequence.  Indicate whether the pattern of signal peptide cleavage is
known for the protein.
"""
import sys
import os
import traceback
from Utils import *
Initialize()

BriefFlag = 0 # for debugging

# No cysteine in these dictionaries, becuase cysteines are modified
# with +57 protecting group, abolishing acidic activity:
AAPKA = {"Nterm":8, "Cterm":3.1, "D":4.4, "R":12, "E":4.4, "H":6.5, "Y":10, "K":9.8, "C":8.5}
AAPKASign = {"Nterm":1, "Cterm":-1, "D":-1, "R":1, "E":-1, "H":1, "Y":-1, "K":1, "C":-1}

def GetPI(Aminos, MaskCysteines = 0):
    Counts = {"Nterm":1, "Cterm":1}

    # If the n-terminus is acetylated, then it's not going to get charged as usual:
    if Aminos[0] != Aminos[0].upper():
        Counts["Nterm"] = 0
        
    # If MaskCysteines is true, then cysteines have a protecting group that abolishes
    # their basic activity.  To capture that effect on PI, we skip all cysteines:
    if MaskCysteines:
        Aminos = Aminos.replace("C","X")
    for Amino in Aminos:
        if Amino in ("D", "R", "C", "E", "H", "Y", "K"):
            Counts[Amino] = Counts.get(Amino, 0) + 1
    PrevCharge = 99        
    for PH in range(1,13):
        Charge = 0
        for (Key, PKA) in AAPKA.items():
            Diff = AAPKASign[Key] * (PKA - PH)
            ChargeContribution = AAPKASign[Key] * Counts.get(Key, 0) * 10**Diff / float(10**Diff + 1)
##            if Counts.get(Key):
##                print AAPKASign[Key], Counts[Key], Diff, 10**Diff, 10**Diff + 1, 10**Diff / float(10**Diff + 1)
##                print "Charge contribution from %s: %s"%(Key, ChargeContribution)
            Charge += ChargeContribution
        if Charge <= 0:
            DY = PrevCharge - Charge
            return (PH - 1) + (PrevCharge / DY)
        PrevCharge = Charge
    return 13

AAHydropathy = {"I":4.5, "V":4.2, "L":3.8, "F":2.8, "C":2.5,
                "M":1.9, "A":1.8, "G":-.4, "T":-0.7, "W":-0.9,
                "S":-0.8, "Y":-1.3, "P":-1.6, "H":-3.2, "E":-3.5,
                "Q":-3.5, "D":-3.5, "N":-3.5, "K":-3.9, "R":-4.5
                }

def GetHydropathy(Aminos):
    AvgHydro = 0
    Count = 0
    for Amino in Aminos:
        AvgHydro += AAHydropathy.get(Amino, 0)
        Count += 1
    return AvgHydro / float(Count)

def GetMW(Aminos):
    Aminos = Aminos.upper()
    MW = 0
    for Amino in Aminos:
        MW += Global.AminoMass.get(Amino, 0)
    MW += 19
    return MW

class ProteinClass:
    def __init__(self):
        self.IPIID = ""
        self.Sprot = ""
        self.SprotNice = ""
        self.Maturation = "Unknown"
    def Copy(self):
        Protein = ProteinClass()
        Protein.IPIID = self.IPIID
        Protein.Sprot = self.Sprot
        Protein.SprotNice = self.SprotNice
        Protein.Maturation = self.Maturation
        Protein.PI = self.PI
        Protein.Hydropathy = self.Hydropathy
        Protein.MW = self.MW
        Protein.MaturePI = self.MaturePI
        Protein.MatureHydropathy = self.MatureHydropathy
        Protein.MatureMW = self.MatureMW
        return Protein
    def ComputeParams(self, Sequence, SetMature = 1):
        self.PI = GetPI(Sequence)
        self.Hydropathy = GetHydropathy(Sequence)
        self.MW = GetMW(Sequence)
        if SetMature:
            self.MaturePI = self.PI
            self.MatureHydropathy = self.Hydropathy
            self.MatureMW = self.MW
    def ComputeMatureParams(self, MatureSequence):
        self.MaturePI = GetPI(MatureSequence)
        self.MatureHydropathy = GetHydropathy(MatureSequence)
        self.MatureMW = GetMW(MatureSequence)
    def SetMatureForm(self, Start, End):
        if End == None:
            self.MatureSequence = self.MatureSequence[Start-1:]
        elif Start == None:
            self.MatureSequence = self.MatureSequence[:End-1]
        else:
            self.MatureSequence = self.MatureSequence[Start-1:End-1]
    def Describe(self):
        Str = "%s\t%s\t%s\t"%(self.IPIID, self.Sprot, self.SprotNice)
        Str += "%.2f\t%.2f\t%s\t"%(self.PI, self.Hydropathy, self.MW)
        Str += "%.2f\t%.2f\t%s\t"%(self.MaturePI, self.MatureHydropathy, self.MatureMW)
        Str += "%s\t"%self.Maturation
        print Str

class PIReporterClass:
    def __init__(self):
        self.Proteins = {}
        self.ProteinList = []
    def FinishIPIProtein(self, CurrentID, CurrentSprot, CurrentSprotNice, CurrentSequence):
        if not CurrentSequence:
            return
        #print CurrentID, CurrentSprot, CurrentSprotNice
        if self.Proteins.has_key(CurrentSprot):
            Protein = self.Proteins[CurrentSprot]
            if Protein.IPIID:
                Protein = Protein.Copy()
                self.ProteinList.append(Protein)
        elif self.Proteins.has_key(CurrentSprotNice):
            Protein = self.Proteins[CurrentSprotNice]
            if Protein.IPIID:
                Protein = Protein.Copy()
                self.ProteinList.append(Protein)
        else:
            Protein = ProteinClass()
            Protein.ComputeParams(CurrentSequence)
            Protein.Sprot = CurrentSprot
            Protein.SprotNice = CurrentSprotNice
            self.ProteinList.append(Protein)
        Protein.IPIID = CurrentID
        self.Proteins[CurrentID] = Protein
    def FinishSprotProtein(self, CurrentSprot, CurrentSprotNice, Sequence, CurrentFeatures):
        if not Sequence:
            return
        Protein = ProteinClass()
        self.Proteins[CurrentSprot] = Protein
        if CurrentSprotNice:
            self.Proteins[CurrentSprotNice] = Protein
        Protein.Sprot = CurrentSprot
        Protein.SprotNice = CurrentSprotNice
        Protein.ComputeParams(Sequence)
        Maturation = "Unknown"
        for Feature in CurrentFeatures:
            Bits = Feature.split()
            FeatureType = Bits[1]
            try:
                FeatureStart = int(Bits[2])
                FeatureEnd = int(Bits[3])
            except:
                continue # multiline feature, probably.
            if FeatureType == "SIGNAL":
                Sequence = ("X" * FeatureEnd) + Sequence[FeatureEnd:]
                if len(Bits)>4:
                    Maturation = Bits[4]
                    if Maturation == "By":
                        Maturation = "By Similarity"
                else:
                    Maturation = ""
            elif FeatureType == "DISULFID":
                Sequence = Sequence[:FeatureStart-1] + "c" + Sequence[FeatureStart:FeatureEnd-1] + "c" + Sequence[FeatureEnd:]
                ##disulfid 10 15
                ##01234567890123456789
                ##1234567890123456789
                ##         x    x
            elif FeatureType == "MOD_RES":
                # Most common: N-terminal acetylation
                if Bits[4][:8] == "N-acetyl":
                    Sequence = Sequence[:FeatureStart - 1] + Sequence[FeatureStart - 1].lower() + Sequence[FeatureStart:]
        Protein.ComputeMatureParams(Sequence)
        Protein.Maturation = Maturation
        self.Proteins[CurrentSprot] = Protein
        self.Proteins[CurrentSprotNice] = Protein
        self.ProteinList.append(Protein)
    def ReadSprot(self, FileName):
        """
        Read proteins from the swiss-prot database
        """
        File = open(FileName, "r")
        CurrentSequence = ""
        CurrentSprot = ""
        CurrentSprotNice = ""
        CurrentFeatures = []
        LineNumber = 0
        for FileLine in File.xreadlines():
            LineNumber += 1
            if BriefFlag and LineNumber > 10000:
                break
            if FileLine[:3] == "ID ":
                self.FinishSprotProtein(CurrentSprot, CurrentSprotNice, CurrentSequence, CurrentFeatures)
                CurrentSequence = ""
                CurrentFeatures = []
                CurrentSprotNice = FileLine.split()[1]
            elif FileLine[:3] == "AC ":
                CurrentSprot = FileLine.strip().split()[1]
                if CurrentSprot[-1] == ";":
                    CurrentSprot = CurrentSprot[:-1]
            elif FileLine[:3] == "FT ":
                # Features that may be of interest:
                # FT   SIGNAL        1     30       Potential.
                # FT   DISULFID     41     77       By similarity.
                # FT   MOD_RES      22     22       Pyrrolidone carboxylic acid.
                CurrentFeatures.append(FileLine)
            elif FileLine[0] == " ":
                CurrentSequence += FileLine.strip().replace(" ","")
        self.FinishSprotProtein(CurrentSprot, CurrentSprotNice, CurrentSprot, CurrentFeatures)
    def ReadIPI(self, FileName):
        "Read proteins from the IPI database"
        # Called *after* we have read in swiss-prot.
        File = open(FileName, "r")
        CurrentSequence = ""
        CurrentID = ""
        CurrentSprot = ""
        CurrentSprotNice = ""
        LineNumber = 0
        for FileLine in File.xreadlines():
            LineNumber += 1
            if BriefFlag and LineNumber > 10000:
                break
            # ID   IPI00000001.1         IPI;      PRT;   577 AA.
            if FileLine[:3] == "ID ":
                self.FinishIPIProtein(CurrentID, CurrentSprot, CurrentSprotNice, CurrentSequence)
                CurrentSprot = ""
                CurrentSprotNice = ""
                CurrentSequence = ""
                CurrentID = FileLine.split()[1]
                if CurrentID[-2] == ".":
                    CurrentID = CurrentID[:-2]
            # DR   UniProtKB/Swiss-Prot; O95793-1; STAU_HUMAN; M.
            if FileLine[:3] == "DR " and FileLine.find("Swiss-Prot")!=-1:
                CurrentSprot = FileLine.split()[2]
                CurrentSprotNice = FileLine.split()[3]
                if CurrentSprot[-1] == ";":
                    CurrentSprot = CurrentSprot[:-1]
                if CurrentSprot[-2] == "-":
                    CurrentSprot = CurrentSprot[:-2]
                if CurrentSprotNice[-1] == ";":
                    CurrentSprotNice = CurrentSprotNice[:-1]
            # Only sequence lines start with a blank:
            if FileLine[0] == ' ':
                CurrentSequence += FileLine.strip().replace(" ","")
        self.FinishIPIProtein(CurrentID, CurrentSprot, CurrentSprotNice, CurrentSequence)
    def Report(self):
        for Protein in self.ProteinList:
            Protein.Describe()

def LoadPITable(FileName = "IPIPI.txt"):
    "Read the PI table that was written out by PIReporterClass::Report()"
    File = open(FileName, "r")
    IPIPI = {} # cute variable name...ipipipipipi longstocking
    for FileLine in File.xreadlines():
        Bits = FileLine.split("\t")
        if len(Bits) < 10:
            continue
        PI = float(Bits[6])
        Hydropathy = float(Bits[7])
        MW = float(Bits[8])
        IPIPI[Bits[0]] = (PI, Hydropathy, MW, Bits[9])
    return IPIPI
            
if __name__ == "__main__":
    PIReporter = PIReporterClass()
    PIReporter.ReadSprot("Database\\Sprot47.2.dat")
    PIReporter.ReadIPI("Database\\ipi.HUMAN.v3.11.dat")
    PIReporter.Report()