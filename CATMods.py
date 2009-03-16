"""
Given the CAT (consensus annotation table), report on modification sites and their presence
in different samples.  For Larry David lab lens spectra.
"""
import os
import sys
import traceback
import math
import string
from Utils import *
Initialize()
from LDUtil import *

def PrintConsensusHiscore(ModCoverage):
    OutFile = open("TrueConsensusAnnotationTable.txt", "w")
    for ModInfoList in ModCoverage:
        SortedList = [] # sort by score!
        for ModInfo in ModInfoList:
            if (ModInfo[9] and ModInfo[10]) and ModInfo[11] and ModInfo[12]:
                SortedList.append((ModInfo[13], ModInfo))
        if not SortedList:
            continue
        SortedList.sort()
        SortedList.reverse()
        (Score, ModInfo) = SortedList[0]
        if Score < 1.0:
            continue
        Line = ModInfo[-1]
        Bits = Line.split("\t")
        FileName = Bits[3]
        (Stub, Dummy) = os.path.splitext(FileName)
        Path = LDGetFilePath(FileName)
        Str = "%s\t%s\t%s\t"%(FileName, ModInfo[14], Score)
        OutFile.write(Str + "\n")
        Command = "Label.py %s %s LDConsensus\%s.png"%(Path, ModInfo[14], Stub)
        print Command
        os.system(Command)
        Command = "copy %s LDConsensus"%Path
        print Command
        os.system(Command)
        
        
def GetGroups(A, B, C):
    """
    Return a list of all groups a spectrum is contained in.  All spectra are in group "all".
    Example: For (93cat, Insoluble, QTOF) a spectrum is in groups all, cat, 93cat, Insoluble, QTOF.
    """
    Groups = ["all"]
    if C == "QTOF":
        Groups.append("qtof")
        if B == "Soluble":
            Groups.append("QSoluble")
        else:
            Groups.append("QInsoluble")
        return Groups # That's *all*, for QTOF
    elif C == "IonOld":
        Groups.append("ltq")
        if B == "Soluble":
            Groups.append("OldIonSoluble")
        else:
            Groups.append("OldIonInsoluble")
        return Groups # That's *all*, for QTOF
    else:
        Groups.append("ltq")
        Groups.append(B)
    Groups.append(A)
    if A != "3day":
        Groups.append("aged")
        if B == "Soluble":
            Groups.append("AgedSol")
        else:
            Groups.append("AgedInsol")
    if A in ("3day", "70cont"):
        Groups.append("cont")
        if B == "Soluble":
            Groups.append("ContSol")
        else:
            Groups.append("ContInsol")
    else:
        Groups.append("cat")
        if B == "Soluble":
            Groups.append("CatSol")
        else:
            Groups.append("CatInsol")
        
    if A == "93cat" and B == "Insoluble":
        Groups.append("93CatInsol")
    if A == "93cat" and B == "Soluble":
        Groups.append("93CatSol")
    if A == "70cat" and B == "Insoluble":
        Groups.append("70CatInsol")
    if A == "70cat" and B == "Soluble":
        Groups.append("70CatSol")
    if A == "70cont" and B == "Soluble":
        Groups.append("70ContSol")
    if A == "70cont" and B == "Insoluble":
        Groups.append("70ContInsol")
    if A == "3day" and B == "Soluble":
        Groups.append("3daySol")
    #print A, B, C, Groups
    return Groups

def PlotProteinCoverage(GroupCoverage):
    ReportGroups = ("93CatInsol", "93CatSol", "70CatInsol", "70CatSol", "70ContInsol", "70ContSol","3daySol")
    Str = "DBPos\tProtein\tResidue\tAA\t"
    for Group in ReportGroups:
        Str += "%s\t"%Group
    for Group in ReportGroups:
        Str += "%sNorm\t"%Group
    print Str
    for ProteinIndex in range(min(25, len(DBProteinNames))):
        GroupTotal = {}
        for Group in ReportGroups:
            GroupTotal[Group] = 1
        ProteinName = DBProteinNames[ProteinIndex]
        ProteinStart = DBProteinPositions[ProteinIndex]
        if ProteinIndex < len(DBProteinNames) - 1:
            ProteinEnd = DBProteinPositions[ProteinIndex + 1]
        else:
            ProteinEnd = len(DB) - 1
        for Pos in range(ProteinStart, ProteinEnd):
            Count = GroupCoverage[Group][Pos]
            for Group in ReportGroups:
                GroupTotal[Group] = max(GroupTotal[Group], Count)
        for Pos in range(ProteinStart, ProteinEnd):
            Str = "%s\t%s\t%s\t%s\t"%(Pos, ProteinName, Pos - ProteinStart + 1, DB[Pos])
            for Group in ReportGroups:
                Str += "%s\t"%GroupCoverage[Group][Pos]
            for Group in ReportGroups:
                Str += "%s\t"%(GroupCoverage[Group][Pos] / float(GroupTotal[Group]))
            print Str
        print

def PrintCoverageTable(GroupCoverage):
    NextPos = DBProteinPositions[1]
    Name = DBProteinNames[0]
    ProteinIndex = 0
    ProteinPos = 1
    File = open("LDCoverage.txt", "w")
    for Pos in range(len(DB)):
        if Pos >= NextPos:
            ProteinIndex += 1
            if ProteinIndex < len(DBProteinNames):
                Name = DBProteinNames[ProteinIndex]
                if ProteinIndex < len(DBProteinPositions) - 1:
                    NextPos = DBProteinPositions[ProteinIndex + 1]
                ProteinPos = 1
        else:
            ProteinPos += 1
        Amino = DB[Pos]
        if Amino != "*":
            File.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n"%(Pos, Name, ProteinPos, DB[Pos],
                                                 GroupCoverage["all"][Pos],
                                                 GroupCoverage["3daySol"][Pos],
                                                 GroupCoverage["cat"][Pos],
                                                 GroupCoverage["CatInsol"][Pos],
                                                 ))
    File.close()            

def PrintModsTable(HighlightUnexplained = 0):
    #Coverage = [0] * (len(DB) + 1)
    GroupCoverage = {}
    GroupNTermini = {}
    GroupCTermini = {}
    for GroupName in ("all", "Soluble", "Insoluble", "QSoluble", "QInsoluble",
                      "OldIonSoluble", "OldIonInsoluble", 
                      "93cat", "70cat", "70cont", "3day", "cat", "cont",
                      "qtof", "ltq", "93CatInsol", "70CatInsol", "3daySol",
                      "93CatSol", "70CatSol", "70ContSol", "70ContInsol",
                      "CatInsol", "CatSol", "ContInsol", "ContSol", "aged",
                      "AgedSol", "AgedInsol"):
        GroupCoverage[GroupName] = [0]*(len(DB) + 1)
        GroupNTermini[GroupName] = [0]*(len(DB) + 1)
        GroupCTermini[GroupName] = [0]*(len(DB) + 1)
    ModCoverage = []
    for X in range(len(DB) + 1):
        ModCoverage.append([])
    MinSpectrumCount = 1 #  2
    MinMQScore = -1.0 # 0
    MinBestScore = 0.5 # 1
    File = open("e:\\ms\\LarryDavid\\CAT3.txt")
    LineNumber = 0
    for FileLine in File.xreadlines():
        LineNumber += 1
        #if LineNumber > 100000:
        #    break 
        Bits = FileLine.split("\t")
        if len(Bits)<15:
            continue
        try:
            Score = float(Bits[6])
        except:
            continue
        if Score < MinMQScore:
            continue
        Peptide = GetPeptideFromModdedName(Bits[4][2:-2])
        Pos = DB.find(Peptide.Aminos)
        if Pos == -1:
            continue 
        # Take note of which groups this scan is in:
        if Bits[0] == "?":
            continue
        Groups = GetGroups(Bits[0], Bits[1], Bits[2])
        (ProteinName, Residue) = GetLensProtein(Peptide.Aminos)
        for Group in Groups:
            GroupNTermini[Group][Pos] += 1
            GroupCTermini[Group][Pos + len(Peptide.Aminos) - 1] += 1        
        for Index in range(len(Peptide.Aminos)):
            for Group in Groups:
                GroupCoverage[Group][Pos + Index] += 1
            if Peptide.Modifications.has_key(Index):
                for Mod in Peptide.Modifications[Index]:
                    ModInfo = (Peptide.Aminos, int(round(Mod.Mass)), Bits[0], Bits[1], Bits[2],
                               int(Bits[11]), int(Bits[12]), int(Bits[13]), int(Bits[14]),
                               int(Bits[15]), int(Bits[16]), int(Bits[17]), int(Bits[18]),
                               Score, Peptide.GetModdedName(), ProteinName, Residue + Index, FileLine)
                    ModCoverage[Pos + Index].append(ModInfo)
    return PrintConsensusHiscore(ModCoverage) 
    # Summarize coverage:
    #PrintCoverageTable(GroupCoverage)
    #PlotProteinCoverage(GroupCoverage)
    #PlotProteinCoverage(GroupNTermini)
    #PlotProteinCoverage(GroupCTermini)
    # Now summarize the sites:
    # Header:
    GroupPairings = [("3day", "aged"),
                     ("AgedInsol", "AgedSol"),
                     ("70CatInsol", "3daySol"),
                     ("93CatInsol", "3daySol"),
                     ("93CatInsol", "3daySol"),
                     ("70cat", "70cont"),
                     ("70CatInsol", "70CatSol"),
                     ("93CatInsol", "70ContSol"),
                     ("93CatInsol", "93CatSol"),
                     ("cat", "cont"),
                     ("Soluble", "Insoluble"),
                     ]
    Str = "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t"
    for (A, B) in GroupPairings:
        Str += "%s vs %s\t\t\t\t"%(A, B)
    print Str
    Str = "Protein\tResidue\tMod\tPeptides\tBestScore\tMatch\tConsensus\tConsensusSameSpectrum\tConsensusSite\tConsensusSiteSameSpectrum\tRate\tModAll\tAll\tModCat\tCat\tModCont\tCont\tModSol\tSol\tModInsol\tInsol\tMod93cat\t93cat\tMod70cat\t70cat\tMod70cont\t70cont\tMod3day\t3day\tModQTOF\tQTOF\tModLTQ\tLTQ\t"
    for X in range(len(GroupPairings)):
        Str += "RateDiff\tRateChange\tChiSquared\tDetails\t"
    print Str
    # Print a grop of rows for each locus, a row for each modification on the site.
    for Pos in range(len(ModCoverage)):
        Mods = ModCoverage[Pos]
        if len(Mods) < MinSpectrumCount:
            continue
        # Report the modifications from commonest to rarest:
        MassCounts = {}
        for Mod in Mods:
            MassCounts[Mod[1]] = MassCounts.get(Mod[1], 0) + 1
        SortedList = []
        for (Mass, Count) in MassCounts.items():
            SortedList.append((Count, Mass))
        SortedList.sort()
        for (ThisMassCount, Mass) in SortedList:
            if ThisMassCount < MinSpectrumCount:
                continue
            # Loop over modifications:
            # Count the number of peptides.
            # Count the number of hits for all groups.
            # Determine whether the mod is triply-confirmed overall,
            # or triply confirmed in one spectrum
            GroupCounts = {}
            Peptides = {}
            PepNames = {}
            Consensus = [0, 0, 0]
            SSConsensus = [0, 0, 0]
            ConsensusTrue = 0
            SSConsensusTrue = 0
            SpectrumCount = 0
            BestScore = -999
            for Mod in Mods:
                if Mod[1]!=Mass:
                    continue
                Peptides[Mod[0]] = Peptides.get(Mod[0], 0) + 1
                GroupList = GetGroups(Mod[2], Mod[3], Mod[4])
                for Group in GroupList:
                    GroupCounts[Group] = GroupCounts.get(Group, 0) + 1
                Consensus[0] += Mod[5] + Mod[6]
                Consensus[1] += Mod[7]
                Consensus[2] += Mod[8]
                SSConsensus[0] += Mod[9] + Mod[10]
                SSConsensus[1] += Mod[11]
                SSConsensus[2] += Mod[12]
                Flags = (Mod[5] or Mod[6]) + Mod[7] + Mod[8]
                if Flags == 3:
                    ConsensusTrue = 1
                Flags = (Mod[9] or Mod[10]) + Mod[11] + Mod[12]
                if Flags == 3:
                    SSConsensusTrue = 1
                SpectrumCount += 1
                if Mod[13] > BestScore:
                    BestScore = Mod[13]
                PepNames[Mod[14]] = PepNames.get(Mod[14], 0) + 1
            PepSorted = []
            for (Name, Count) in PepNames.items():
                PepSorted.append((Count, Name))
            PepSorted.sort()
            BestName = PepSorted[-1][1]
            if Consensus[0] and Consensus[1] and Consensus[2]:
                Consensus = 1
            else:
                Consensus = 0
            if not Consensus and BestScore < MinBestScore:
                continue # Not interesting!
            if SSConsensus[0] and SSConsensus[1] and SSConsensus[2]:
                SSConsensus = 1
            else:
                SSConsensus = 0
            if HighlightUnexplained:
                if BestScore > 2 and SpectrumCount>1 and not Consensus:
                    for Mod in Mods:
                        print Mod[-1].strip()
                continue
            Str = "%s\t%s\t"%(Mod[15], Mod[16])
            Str += "%s\t%s\t%.2f\t%s\t"%(Mass, len(Peptides.keys()), BestScore, BestName)
            Str += "%s\t%s\t%s\t%s\t"%(Consensus, ConsensusTrue, SSConsensus, SSConsensusTrue)
            OverallRate = GroupCounts.get("all", 0) / float(GroupCoverage["all"][Pos])
            Str += "%.1f\t"%(100*OverallRate)
            for GroupName in ("all", "3day", "aged", "AgedInsol", "AgedSol", "70cat", "70cont", "cat", "cont"):
                Str += "%s\t%s\t"%(GroupCounts.get(GroupName, 0), GroupCoverage[GroupName][Pos])
            # Perform a simple statistical test to see whether the PTM seems to have different
            # rates between samples:
            DiffString = ""
            for (GroupA, GroupB) in GroupPairings:
                #("cat", "cont"), ("93cat", "70cat"), ("70cont", "3day"), ("Soluble", "Insoluble"),
                #("QSoluble", "QInsoluble"), ("OldIonSoluble", "OldIonInsoluble"),
                #("70cat","70cont")]:
                ModA = GroupCounts.get(GroupA, 0)
                CountA = GroupCoverage[GroupA][Pos]
                ModB = GroupCounts.get(GroupB, 0)
                CountB = GroupCoverage[GroupB][Pos]
                # Require that (a) both samples have spectral coverage >10, and (b) one sample has
                # at least 10 modified spectra
    ##            if max(ModA, ModB) < 5:
    ##                Str += "\t\t\t\t"
    ##                continue
    ##            if min(CountA, CountB) < 10:
    ##                Str += "\t\t\t\t"
    ##                continue
                AllMod = ModA + ModB
                AllModless = (CountA + CountB) - AllMod
                if AllModless == 0:
                    Str += "\t\t\t\t"
                    continue
                RateA = ModA / max(1, float(CountA))
                RateB = ModB / max(1, float(CountB))
                RateOverall = (ModA + ModB) / float(CountA + CountB)
    ##            #print ModA, CountA, RateA, ModB, CountB, RateB
    ##            if RateA > 0.99 or RateA < 0.01:
    ##                DifferentLP = 0
    ##            else:
    ##                DifferentLP = ModA * math.log(RateA) + (CountA - ModA) * math.log(1 - RateA)
    ##            if RateB > 0.99 or RateB < 0.01:
    ##                DifferentLP += 0
    ##            else:
    ##                DifferentLP += ModB * math.log(RateB) + (CountB - ModB) * math.log(1 - RateB)
                ChiSquared = 0
                Expect1 = RateOverall * CountA
                if Expect1:
                    ChiSquared += (ModA - Expect1)**2 / float(Expect1)
                Expect2 = (1.0 - RateOverall) * CountA
                if Expect2:
                    ChiSquared += (CountA - ModA - Expect2)**2 / float(Expect2)
                Expect3 = RateOverall * CountB
                if Expect3:
                    ChiSquared += (ModB - Expect3)**2 / float(Expect3)
                Expect4 = (1.0 - RateOverall) * CountB
                if Expect4:
                    ChiSquared += (CountB - ModB - Expect4)**2 / float(Expect4)
                if min(Expect1, Expect2, Expect3, Expect4)<=5:
                    ChiSquared = ""
                #Ratio = DifferentLP - SameLP
                RateDiff = 100*(RateA - RateB)
                if max(RateA, RateB):
                    RateChange = 100*(RateA - RateB)/float(max(RateA, RateB))
                else:
                    RateChange = 0
                
                Str += "%.2f\t%.2f\t%s\t%s/%s (%.1f) vs %s/%s (%.1f)\t"%(RateDiff, RateChange, ChiSquared, ModA, CountA, 100*RateA, ModB, CountB, 100*RateB)
                #Str += "%s vs %s: %s/%s vs %s/%s (%.2f)\t"%(GroupA, GroupB,
                #    ModA, CountA, ModB, CountB, Ratio)
                #Str += DiffString
            print Str

def PrintDiffs():
    File = open("Mods2.txt", "rb")
    SortedLists = [[],[],[],[]]
    for FileLine in File.xreadlines():
        Bits = FileLine.split("\t")
        try:
            Score = float(Bits[4])
        except:
            continue
        Consensus = int(Bits[6])
        if not Consensus:
            continue
        # Decide whether the line is "interesting":
        for X in range(4):
            DiffCol = 33 + X*4
            RateChangeCol = 34 + X*4
            RatioTestCol = 35 + X*4
            if RatioTestCol >= len(Bits):
                continue
            try:
                Diff = float(Bits[DiffCol])
            except:
                continue
            if abs(Diff) < 3:
                continue
            Change = float(Bits[RateChangeCol])
            if abs(Change) < 25:
                continue
            RatioTest = float(Bits[RatioTestCol])
            if RatioTest < 3.84:
                continue
            SortedLists[X].append((Bits[2],FileLine.strip()))
    for X in range(4):
        print "\n\n"
        SortedLists[X].sort()
        for (Dummy, Line) in SortedLists[X]:
            print Line

if __name__ == "__main__":
    PrintModsTable(0)
    #PrintDiffs()
