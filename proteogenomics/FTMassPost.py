"""
Post-process the output of FTMass.py.  Compute the average *monoisotopic* mass
for each (common) modification type.  For example, for oxidized methionine, the
average should be ~15.9994 (according to UNIMOD)!
"""
import os
import sys
import BasicStats

class PostProcessor:
    def __init__(self):
        self.MassesByType = {}
        self.MassesByType = {}
        self.TypesToGroups = {}
        self.TypesToGroups["oxidation M"] = "Oxidation"
        self.TypesToGroups["oxidation W"] = "Oxidation"
        self.TypesToGroups["CAM H"] = "CAM HK"
        self.TypesToGroups["CAM K"] = "CAM HK"
        self.TypesToGroups["dehydration D"] = "dehydration"
        self.TypesToGroups["dehydration E"] = "dehydration"
        self.TypesToGroups["methyl ester Y"] = "methyl ester"
        self.TypesToGroups["methyl ester E"] = "methyl ester"
        self.TypesToGroups["double oxidation W"] = "double oxidation"
        self.TypesToGroups["double oxidation M"] = "double oxidation"
        self.Groups = ["Oxidation", "CAM HK", "dehydration",
                       "double oxidation", "methyl ester"]
        # Fix masses - even if the isotopic distribution confuses us!
        self.OverrideMasses = {}
        self.OverrideMasses["Intramolecular disulfide"] = -116
        self.OverrideMasses["methyl ester Y"] = 14
        self.OverrideMasses["methyl ester E"] = 14
        self.OverrideMasses["hydroxyproline"] = 16
        self.OverrideMasses["oxidation W"] = 16
        self.OverrideMasses["double oxidation W"] = 32
        self.OverrideMasses["disulfide/sulfone"] = -84
        self.OverrideMasses["carbamylation"] = 43
        self.OverrideMasses["methylglutamine"] = 14
        self.OverrideMasses["formylation"] = 28
        self.OverrideMasses["C-57"] = -57
        self.OverrideMasses["Hydroxylation"] = 16
        self.OverrideMasses["CAM+carb"] = 100
        self.OverrideMasses["form+CAM"] = 84
        self.OverrideMasses["CAM+CAM"] = 114
        self.OverrideMasses["dehydration D"] = -18
        self.OverrideMasses["dehydration E"] = -18
        self.OverrideMasses["dehydration"] = -18
        self.OverrideMasses["CAM H"] = 57
        self.OverrideMasses["CAM K"] = 57
        self.OverrideMasses["CAM HK"] = 57
        self.OverrideMasses["NT+CAM"] = 57
        self.OverrideMasses["form+carb"] = 71
        self.OverrideMasses["M+16-64"] = -48
        self.OverrideMasses["persulfide"] = 32
    def LoadMasses(self, FileName):
        File = open(FileName, "rb")
        LineNumber = 0
        for FileLine in File.xreadlines():
            LineNumber += 1
            Bits = FileLine.split("\t")
            try:
                PTMType = Bits[0]
                AA = Bits[1]
                MassType = int(Bits[2])
                DBPos = int(Bits[3])
                Mass = float(Bits[4]) - 1.0078
            except:
                continue # header or crap
            if not self.MassesByType.has_key(PTMType):
                self.MassesByType[PTMType] = []
            self.MassesByType[PTMType].append(Mass)
##            # Remember *+Mass
##            Key = (AA, MassType)
##            if not self.MassesByType.has_key(Key):
##                self.MassesByType[Key] = []
##            self.MassesByType[Key].append(Mass)
##            # Remember AA+Mass
##            Key = (None, MassType)
##            if not self.MassesByType.has_key(Key):
##                self.MassesByType[Key] = []
##            self.MassesByType[Key].append(Mass)
        File.close()
    def ReportHelper(self, AllMasses, Name):
        OverrideMass = self.OverrideMasses.get(Name, None)
        IsotopeCounts = {}
        for Mass in AllMasses:
            IntMass = int(round(Mass))
            IsotopeCounts[IntMass] = IsotopeCounts.get(IntMass, 0) + 1
        BestIsotope = None
        BestIsotopeCount = None
        for (IntMass, Count) in IsotopeCounts.items():
            if Count > BestIsotopeCount:
                BestIsotopeCount = Count
                BestIsotope = IntMass
        # OVERRIDE:
        if OverrideMass:
            BestIsotope = OverrideMass 
        SelectedMassList = []
        for Mass in AllMasses:
            Diff = abs(Mass - BestIsotope)
            if Diff > 0.5:
                continue
            SelectedMassList.append(Mass)
        (Mean, StdDev) = BasicStats.GetMeanStdDev(SelectedMassList)
        Str = "%s\t"%(Name)
        Str += "%s\t%s\t"%(Mean, StdDev)
        Str += "%s\t%s\t"%(len(AllMasses), len(SelectedMassList))
        Str += "%s\t"%str(IsotopeCounts)
        print Str        
    def ReportImportantMasses(self):
        AllMasses = []
        for PTMType in self.MassesByType.keys():
            AllMasses = self.MassesByType[PTMType]
            self.ReportHelper(AllMasses, PTMType)
        # Report grouped PTM types:
        for PTMGroup in self.Groups:
            AllMasses = []
            for (PTMType, MappedGroup) in self.TypesToGroups.items():
                if MappedGroup == PTMGroup:
                    print "Group %s gets %s"%(PTMGroup, PTMType)
                    AllMasses.extend(self.MassesByType.get(PTMType, []))
            self.ReportHelper(AllMasses, PTMGroup)                
    def xReportImportantMasses(self):
        ImportantMasses = [
                           [("M",  16),("W",  16)],
                           [("W",  16),],
                           [("M",  16),],
                           [("M",  43),],
                           [(None, 43),],
                           [("Q",  -17),],
                           [(None, 28),],
                           [("H",  57),],
                           [(None, 57),],
                           [("K",  57),],
                           [("K",  57),("H",  57)],
                           [("E",  14),],
                           [("Y",  14),],
                           [("E",  14),("Y",  14)],
                           [("M",  32),("W",  32)],
                           [("N",  -17),],
                           [(None, -116),],
                           [(None, 85),],
                           [(None, 53),],
                           [("C",  77),],
                           [("D", -18),("E", -18),],
                           [("C",  109),],
                           [(None, 114),],
                           [(None, 100),],
                           ]
        for KeyList in ImportantMasses:
            AllMasses = []
            for Key in KeyList:
                MassList = self.MassesByType.get(Key, None)
                if MassList:
                    AllMasses.extend(MassList)
                else:
                    print "%s\t%s\tNo hits!"%(Key[0], Key[1])
            IsotopeCounts = {}
            for Mass in AllMasses:
                IntMass = int(round(Mass))
                IsotopeCounts[IntMass] = IsotopeCounts.get(IntMass, 0) + 1
            BestIsotope = None
            BestIsotopeCount = None
            for (IntMass, Count) in IsotopeCounts.items():
                if Count > BestIsotopeCount:
                    BestIsotopeCount = Count
                    BestIsotope = IntMass
            # OVERRIDE:
            BestIsotope = Key[1]
            SelectedMassList = []
            for Mass in AllMasses:
                Diff = abs(Mass - BestIsotope)
                if Diff > 0.5:
                    continue
                SelectedMassList.append(Mass)
            (Mean, StdDev) = BasicStats.GetMeanStdDev(SelectedMassList)
            Str = "%s\t"%(KeyList)
            Str += "%s\t%s\t"%(Mean, StdDev)
            Str += "%s\t%s\t"%(len(AllMasses), len(SelectedMassList))
            Str += "%s\t"%str(IsotopeCounts)
            print Str
                
if __name__ == "__main__":
    Bob = PostProcessor()
    Bob.LoadMasses(sys.argv[1])
    Bob.ReportImportantMasses()