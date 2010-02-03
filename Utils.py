"""
Various support functions and constants to support ms/ms algorithms.
Amino acid mass importer, ion type list generation.
"""
import Global
import os
import sys
import string

if hasattr(os, "sysconf"):
    IS_WINDOWS = 0
else:
    IS_WINDOWS = 1


# IDs for all chromosome files.  Usually ID == chromosome number, but we also
# give numbers to X, Y, and all of the unlocalized ("random" as in "random access") sequences
ChromosomeMap = {"chr1":1, "chr2":2, "chr3":3, "chr4":4,
                 "chr5":5, "chr6":6, "chr7":7, "chr8":8,
                 "chr9":9, "chr10":10, "chr11":11, "chr12":12,
                 "chr13":13, "chr14":14, "chr15":15, "chr16":16,
                 "chr17":17, "chr18":18, "chr19":19, "chr20":20,
                 "chr21":21, "chr22":22, "chrX":23, "chrY":24,
                 "chrM":25, "chr1_random":26, "chr2_random":27, "chr3_random":28,
                 "chr4_random":29, "chr5_random":30, "chr6_random":31, "chr7_random":32,
                 "chr8_random":33, "chr9_random":34, "chr10_random":35, "chr11_random":36,
                 "chr12_random":37, "chr13_random":38, "chr14_random":39, "chr15_random":40,
                 "chr16_random":41, "chr17_random":42, "chr18_random":43, "chr19_random":44,
                 "chr20_random":45, "chr21_random":46, "chr22_random":47, "chrX_random":48,
                 "chrx":23, "chry":24, "chrm":25, "chrx_random":48
                 }

ReverseChromosomeMap = ["", "chr1", "chr2", "chr3", "chr4", "chr5",
    "chr6", "chr7", "chr8", "chr9", "chr10",
    "chr11", "chr12", "chr13", "chr14", "chr15",
    "chr16", "chr17", "chr18", "chr19", "chr20",
    "chr21", "chr22", "chrX", "chrY", "chrM",
    "chr1_random", "chr2_random", "chr3_random", "chr4_random", "chr5_random",
    "chr6_random", "chr7_random", "chr8_random", "chr9_random", "chr10_random",
    "chr10_random", "chr12_random", "chr13_random", "chr14_random", "chr15_random",
    "chr16_random", "chr17_random", "chr18_random", "chr19_random", "chr20_random",
    "chr21_random", "chr22_random", "chrX_random",]

MassPTMods = {}

class Bag:
    "Generic argument-container class"
    pass

def FixupPath(Path):
    if IS_WINDOWS:
        return Path.replace("/", "\\")
    else:
        return Path.replace("\\", "/")

class PTModClass:
    InstanceCount = 0
    "A type of post-translational modification, such as phosphorylation."
    def __init__(self, Name):
        self.Name = Name
        # Peptides that this modification can affect:
        self.Bases = {}
        self.BaseString = ""
        self.Mass = 0.0
        self.Score = 0.0
        PTModClass.InstanceCount += 1
        try:
            Mass = int(Name)
            self.Mass = Mass
        except:
            pass
    def __del__(self):
        if PTModClass:
            PTModClass.InstanceCount -= 1
    def __str__(self):
        return "<PTMod '%s'>"%self.Name

def LoadPTMods():
    """
    Read the definitions of post-translational modifications from PTMods.txt.
    Line format is tab-delimited, like this: "Alkylation	14.01564	CKRHDENQ"
    (This is rarely used in practice, but is useful for search results that are annotated with names instead of masses)
    """
    File = open(checkConfigPath("PTMods.txt"),"r")
    for FileLine in File.xreadlines():
        FileLine = FileLine.strip()
        if (not FileLine) or FileLine[0]=="#":
            continue
        Bits = FileLine.split("\t")
        PTMod = PTModClass(Bits[0])
        PTMod.Mass = float(Bits[1])
        PTMod.BaseString = Bits[2]
        for Char in Bits[2]:
            PTMod.Bases[Char] = 1
        Global.PTMods[PTMod.Name.lower()] = PTMod
        Global.PTModByShortName[PTMod.Name[:3].lower()] = PTMod
        Global.PTModByShortName[PTMod.Name[:4].lower()] = PTMod
        Global.PTModList.append(PTMod)
    File.close()

class ProteinClass:
    """
    Class representing a single protein: a collection of peptides.
    Can compute sequence coverage, as well as the list of modifications (PTM)
    """
    def __init__(self,Sequence, Type= None):
        self.CellType = Type #use by sam for comparing cell types
        self.Sequence = Sequence
        self.SequenceCoverage = [0]*len(self.Sequence) #counts of spectra at each residue
        self.Peptides = [] #list of UnmodifiedPeptide objects
        self.PositionModSpectraDict = {} #key = (Name,Position) value = ModdedSpectra
        self.Coverage = 0.0
    def GenerateSequenceCoverage(self):
        """
        Goes through all the peptides and increments TotalSpectraCount for the
        Residues that it covers.
        """
        for UPeptide in self.Peptides: #UnmodifiedPeptide
            Start = self.Sequence.find(UPeptide.Aminos)
            Len = len(UPeptide.Aminos)
            for I in range(Start,Start+Len,1):
                self.SequenceCoverage[I] += UPeptide.TotalSpectraCount
        Covered =0
        for I in range(len(self.Sequence)):
            if self.SequenceCoverage[I] > 0:
                Covered +=1
        Coverage = Covered/ float(len(self.Sequence))
        self.Coverage = Coverage
    def GenerateModList(self):
        """
        This method runs through all the peptides in self.Peptides
        and generates a list of modifications on residues
        """
        for UPep in self.Peptides:
            PeptidePosition = self.Sequence.find(UPep.Aminos)
            TotalSpectra = UPep.TotalSpectraCount
            for Peptide in UPep.Peptides: # PeptideClass Object
                for (AminoIndex, ModList) in Peptide.Modifications.items():
                    SpectraThisPeptide = UPep.SpectraCount[Peptide.GetModdedName()]
                    ModificationPosition = AminoIndex + PeptidePosition
                    for Mod in ModList: #PTModClass objects
                        Name = Mod.Name
                        Key = (Name,ModificationPosition)
                        if self.PositionModSpectraDict.has_key(Key):
                            self.PositionModSpectraDict[Key] += SpectraThisPeptide
                        else:
                            self.PositionModSpectraDict[Key] = SpectraThisPeptide

    def AddAnnotation(self,Peptide,SpectrumCounts):
        "Add a Peptide to the protein"
        Found = 0
        for UPeptide in self.Peptides:
            if UPeptide.Aminos == Peptide.Aminos: # same aminos already a UPep object,
                    #perhaps this is a modified version of the same
                UPeptide.AddAnnotation(Peptide,SpectrumCounts)
                Found =1
                break
        if not Found:
            ToAdd = UnmodifiedPeptide(Peptide,SpectrumCounts)
            self.Peptides.append(ToAdd)

class UnmodifiedPeptide:
    """A wrapper for the PeptideClass, it contains all modified states of the
    same amino acid sequence. Useful sometimes.
    """
    def __init__(self,Peptide,SpectrumCounts): #PeptideClass below
        "constructor"
        self.Aminos = Peptide.Aminos
        self.UnmodifiedSpectraCount =0
        self.TotalSpectraCount =0
        self.Peptides = [] #an array of PeptideClass objects
        self.SpectraCount ={} #key = fullname of peptide (no prefix/suffix), value = spectracount
        self.AddAnnotation(Peptide,SpectrumCounts)
        #maybe some sort of modification list of distinct modifications

    def IsMe(self,Peptide):
        if self.Aminos == Peptide.Aminos:
            return 1
        return 0

    def PrintMe(self):    
        print "UnmodifiedPeptide Object: %s"%self.UnmodifiedSequence
        print "Total Spectra %d, UnmodifiedSpectra %d"%(self.TotalSpectraCount,self.UnmodifiedSpectraCount)
        
    def AddAnnotation(self,NewPeptide,SpectrumCounts):
        "adds an annotation to my list, and updates tallies"
        self.TotalSpectraCount += SpectrumCounts
        if len (NewPeptide.Modifications) == 0:
            if self.UnmodifiedSpectraCount == 0:
                self.Peptides.append(NewPeptide)
            self.UnmodifiedSpectraCount += SpectrumCounts
        else:
            #determine if we've already got this modification
            Found =0
            for MyPeptide in self.Peptides:
                if MyPeptide.GetModdedName() == NewPeptide.GetModdedName():
                    Found = 1
                    break
            if not Found:
                self.Peptides.append(NewPeptide)
        self.SpectraCount[NewPeptide.GetModdedName()] = SpectrumCounts
    

class ModificationTypeObject:
    """
    This holds information about a specific type of Modification.  Remember the format
    mod,14,KR    TAB       #methylation.
    mod,DELTA,AMINOS,POSITION,whatever TAB #Name
    """
    def __init__(self,Latin,Name,DMass,Residues,Position):
        self.inLatin = Latin #inVivo, inVitro
        self.Name = Name
        self.DeltaMass = DMass
        self.Residues = Residues
        if self.Residues == "*":
            self.Residues = "ACDEFGHIKLMNPQRSTVWY"
        self.Position = Position
        self.InspectID = ""
        ##Do a littl processing.  the way a modification shows up in the Inspect
        ## output is with a +43, or a -17, or possibly phos.  So I need to generate
        ## the InspectIdentifier that can be compared later in the RemoveSelf function
        if self.Name == "phosphorylation":
            self.InspectID = "phos"
        elif self.DeltaMass < 0:
            self.InspectID = "%s"%self.DeltaMass
        elif self.DeltaMass > 0:
            "this is for positive values of DeltaMass which present a problem"
            self.InspectID = "+" + "%s"%self.DeltaMass
            #try:
            #    if self.DeltaMass[0] == "+":
            #       self.InspectID = "%s"%self.DeltaMass
            #       self.DeltaMass = int(self.DeltaMass[1:]) 
            #except:
            #   self.InspectID = "+" + "%s"%self.DeltaMass
        
    def PrintMe(self):
        "Simple Debugging printer"
        print "I am a ModificationTypeObject for %s"%self.Name
        print "InspectID %s"%self.InspectID
        print "AcceptableResidues %s"%self.Residues
        
    def RemoveSelf(self,Annotation):
        """
        This method takes an input string and looks for modifications which correspond
        to its identity.  If any are found, it removes them from the string and returns it.
        It will remove all copies of itself from the String
        """
        InMod = 0
        StartIndex = -1
        ModString = ""
        I = 0 #loop iterater
        while I < len(Annotation):
            Letter = Annotation[I]
            if not Letter in string.uppercase:
                if not InMod:
                    StartIndex=I
                InMod = 1
                ModString += Letter
            elif InMod:
                #this is the first upper case letter after a modification.
                if ModString == self.InspectID:
                    #this is my Identifier, Check position and residue
                    PositionCheck = 0 #false
                    ResidueCheck = 0
                    if self.Position == "nterminal" and StartIndex == 1:
                        PositionCheck = 1
                    elif not self.Position:
                        PositionCheck = 1 #no position specified (it should be None)
                    ###   Add other position things in here as you get them ###
                    ModifiedResidue = Annotation[StartIndex-1]
                    if self.Residues.find(ModifiedResidue) >= 0:
                        #found the modified residue in self.residue string
                        ResidueCheck = 1
                    if PositionCheck and ResidueCheck:
                        Front = Annotation[:StartIndex]
                        EndIndex = StartIndex + len(self.InspectID)
                        Back = Annotation[EndIndex:]
                        Annotation = Front + Back
                        I = StartIndex-1 ##### VERY IMPORTANT to go back once the Annotation has been reset.
                #regardless of whether this was actually me or not, still reset the vars below
                InMod = 0
                ModString = ""
            I += 1 #must increment for the while loop 
        return Annotation

    def IsMe(self, Identifier,Residue, Position):
        """
        Check to see if all the criteria match
        """
        if not Identifier == self.InspectID:
            return 0
        if self.Position == "nterminal":
            if Position > 0: #zero indexed string
                return 0
        ### add other position identifiers if you have them
        if self.Residues.find(Residue) < 0:
            return 0 # returned a -1 for "not found"
        return 1

def LoadModifications():
    """
    This method reads in two files: InVivoModifications.txt and InVitroModifications.txt
    It makes a ModificationTypeObject out of each mod listed in the files
    (except fixed mods).  These input files are expected to be of the format
    mod,14,KR    TAB       #methylation.
    mod,DELTA_MASS,AMINOS,POSITION__TAB__#Modification name
    """
    LoadModificationsFromFile("InVivoModifications.txt", Global.InVivoMods, "InVivo")
    LoadModificationsFromFile("InVitroModifications.txt", Global.InVitroMods, "InVitro")
    
def LoadModificationsFromFile(FileName, ModificationList, ChemistryType):
    try:
        File = open(checkConfigPath(FileName),"rb")
    except:
        #print "File '%s' not found - not loading mods"%FileName
        return
    for Line in File.xreadlines():
        Line = Line.rstrip()
        Data = Line.split("\t")
        Name = Data[1][1:] #should get rid of the '#'
        Latin = "InVivo"
        InspectInput = Data[0].rstrip() #get rid of any right side junk
        Data = InspectInput.split(",")
        DeltaMass = int (Data[1])
        Residues = Data[2]
        if len(Data) > 3:
            Position = Data[3]
        else:
            Position = None
        Mod = ModificationTypeObject(ChemistryType, Name, DeltaMass, Residues, Position)
        ModificationList.append(Mod)
    File.close()

class PeptideClass:
    """
    Class representing one peptide, possibly with modifications.  We get one PeptideClass instance
    for every match from the trie-based search.  A PeptideClass instance can also (if its PrefixMass
    and SuffixMass members are set) represent a tag.
    """
    # Track number of live instances:
    InstanceCount = 0 
    def __init__(self, Aminos = ""):
        "Constructor - if we have amino acids, get our masses now."
        self.Aminos = Aminos
        self.Masses = []
        # Modifications map amino acid indices to a list of PTModClass instances
        self.Modifications = {}
        self.Score = None
        self.ID = None
        self.RecordNumber = None
        self.PValue = 0
        self.DeltaCN = 0
        self.DeltaCNOther = 0
        if Aminos:
            self.ComputeMasses()
        PeptideClass.InstanceCount += 1
    def GetPTMBeforeAfter(self, Mass):
        PTMBefore = {}
        PTMAfter = {}
        for (AminoIndex, List) in self.Modifications.items():
            for Entry in List:
                if Entry.Mass == Mass:
                    for OtherIndex in range(AminoIndex, len(self.Aminos)+1):
                        PTMBefore[OtherIndex] = 1
                    for OtherIndex in range(0, AminoIndex+1):
                        PTMAfter[OtherIndex] = 1
        return (PTMBefore, PTMAfter)
        
    def GetPhosphoBeforeAfter(self):
        PhosBefore = {}
        PhosAfter = {}
        for (AminoIndex, List) in self.Modifications.items():
            for Entry in List:
                if Entry.Name == "Phosphorylation":
                    for OtherIndex in range(AminoIndex, len(self.Aminos)+1):
                        PhosBefore[OtherIndex] = 1
                    for OtherIndex in range(0, AminoIndex+1):
                        PhosAfter[OtherIndex] = 1
        return (PhosBefore, PhosAfter)
    def __del__(self):
        if PeptideClass:
            PeptideClass.InstanceCount -= 1
    def GetParentMass(self):
        if not self.Masses:
            self.ComputeMasses()
        return 19 + self.Masses[-1]
    def ComputeMasses(self):
        """
        Populate our Masses list, based upon Aminos and Modifications.  Must be called,
        if self.Modifications is edited!
        """
        self.Masses = [0]
        Mass = 0
        for Index in range(len(self.Aminos)):
            Amino = self.Aminos[Index]
            AminoMass = Global.AminoMass.get(Amino, None)
            if AminoMass == None:
                if Amino == "X":
                    print "** Warning: Peptide '%s' contains wild-card amino X, mass is probably wrong."%(self.Aminos)
                    AminoMass = 0
                else:
                    raise ValueError, "Bad amino '%s' in peptide '%s'"%(Amino, self.Aminos)
            Mass += AminoMass
            Mass += Global.FixedMods.get(Amino, 0)
            for Mod in self.Modifications.get(Index, []):
                Mass += Mod.Mass
                # Warn, but don't fail here. (The trick case: We generate tag GVQ instead of GVK,
                # and biotin can't attach to Q.  Bah!)
                #if not Mod.Bases.has_key(Amino):
                #    print "Warning: Amino '%s' in peptide '%s' has illegal modification %s at %s"%(Amino, self.Aminos, Mod.Name, Index)
            self.Masses.append(Mass)
    def GetPTMCount(self):
        Total = 0
        for Key in self.Modifications.keys():
            Total += len(self.Modifications[Key])
        return Total
    def GetFullModdedName(self):
        return "%s.%s.%s"%(self.Prefix, self.GetModdedName(), self.Suffix)
    def GetModdedName(self):
        "Returns the amino sequence with modifications included, like this: EAM+16APK"
        Str = ""
        for Index in range(len(self.Aminos)):
            Amino = self.Aminos[Index]
            Str += "%s"%(Amino)
            for Mod in self.Modifications.get(Index, []):
                Str += "%s"%(Mod.Name[:4].lower())
        return Str
    def __str__(self):
        return "<Peptide '%s'>"%self.Aminos
    def IsValidTag(self, TagPeptide, Epsilon = 2.0):
        """
        Returns true if TagPeptide is a valid tag for this (full-length) peptide
        """
        TotalResidueMass = self.Masses[-1]
        TagLength = len(TagPeptide.Aminos)
        TagAminos = TagPeptide.Aminos.replace("I", "L").replace("Q", "K")
        Aminos = self.Aminos.replace("I", "L").replace("Q", "K")
        for Pos in range(len(self.Masses)):
            PrefixMass = self.Masses[Pos]
            # Check flanking mass:
            if abs(PrefixMass - TagPeptide.PrefixMass) > Epsilon:
                #print "Pos %s: Invalid (prefix %s vs %s)"%(Pos, PrefixMass, TagPeptide.PrefixMass)
                continue
            # Check amino acids:
            if Aminos[Pos:Pos + TagLength] != TagAminos:
                #print "Pos %s: Invalid (aminos %s vs %s)"%(Pos, Aminos[Pos:Pos + TagLength], TagAminos)
                continue
            # Check suffix mass:
            SuffixMass = TotalResidueMass - self.Masses[Pos + TagLength]
            if abs(SuffixMass - TagPeptide.SuffixMass) > Epsilon:
                #print "Pos %s: Invalid (suffix %s vs %s)"%(Pos, SuffixMass, TagPeptide.SuffixMass)
                continue
            return 1
        #Mass = TagPeptide.PrefixMass + TagPeptide.SuffixMass + GetMass(TagPeptide.Aminos)
        
    def IsSame(self, OtherPeptide):
        SubstDict = {"Q": "K", "I": "L"}
        if len(self.Aminos) != len(OtherPeptide.Aminos):
            return 0
        for AminoIndex in range(len(self.Aminos)):
            OurAmino = self.Aminos[AminoIndex]
            TheirAmino = OtherPeptide.Aminos[AminoIndex]
            OurMods = []
            TheirMods = []
            for Mod in self.Modifications.get(AminoIndex, []):
                if Mod.Name[1:3] == "->":
                    OurAmino = Mod.Name[-1].upper()
                else:
                    OurMods.append(Mod.Mass)
            for Mod in OtherPeptide.Modifications.get(AminoIndex, []):
                if Mod.Name[1:3] == "->":
                    TheirAmino = Mod.Name[-1].upper()
                else:
                    TheirMods.append(Mod.Mass)
            OurAmino = SubstDict.get(OurAmino, OurAmino)
            TheirAmino = SubstDict.get(TheirAmino, TheirAmino)
            if OurAmino != TheirAmino:
                return 0
            OurMods.sort()
            TheirMods.sort()
            if OurMods != TheirMods:
                return 0
        return 1
    def __cmp__(self, OtherPeptide):
        if (not isinstance(OtherPeptide, PeptideClass)):
            return 1
        # Sort by score, best to worst:
        if self.Score > OtherPeptide.Score:
            return -1
        if self.Score < OtherPeptide.Score:
            return 1
        return 0
    def GetNTT(self):
        """
        Returns the number of tryptic termini.  (assumes self.prefix and self.suffix
        are set)
        """
        NTT = 0
        if self.Prefix in ("-*X"):
            NTT += 1
        elif (self.Prefix in ("KR")) and (self.Aminos[0] !="P"):
            NTT += 1
        if self.Suffix in ("-*X"):
            NTT += 1
        elif (self.Aminos[-1] in "KR") and (self.Suffix != "P"):
            NTT += 1
        return NTT
    def IsFullyTryptic(self):
        if self.Prefix in ("-", "*"):
            pass
        elif (self.Prefix in ("K", "R")) and self.Aminos[0] != "P":
            pass
        else:
            return 0
        if self.Suffix in ("-", "*"):
            pass
        elif self.Aminos[-1] in ("K", "R") and self.Suffix != "P":
            pass
        else:
            return 0
        return 1
    def GetNiceAnnnotation(self):
        """
        Return an annotation suitable for a filename.  *.ABC.D turns into -.ABC.D
        """
        Str = "%s.%s.%s"%(self.Prefix, self.Aminos, self.Suffix)
        return Str.replace("*", "-")
    
def GetPeptideFromModdedName(TagName):
    """
    Parse a tag with form like "ATphosQ", adding PTMs at the correct spots.
    """
    StringPos = 0
    Peptide = PeptideClass()
    # If the name has the form K.ABCDER.G, then strip off the prefix and suffix:
    if len(TagName) > 4 and TagName[1] == "." and TagName[-2] == ".":
        Peptide.Prefix = TagName[0]
        Peptide.Suffix = TagName[-1]
        TagName = TagName[2:-2]
    try:
        while (1):
            if StringPos >= len(TagName):
                break
            if TagName[StringPos] in string.uppercase:
                Peptide.Aminos += TagName[StringPos]
                StringPos += 1
            else:
                # It's a modification:
                ModName = ""
                while (StringPos<len(TagName) and TagName[StringPos] not in string.uppercase) and len(ModName)<4:
                    if ModName and ModName[0] in ("-","+") and TagName[StringPos] not in "0123456789k":
                        break
                    ModName += TagName[StringPos]
                    StringPos += 1
                Mod = Global.PTModByShortName.get(ModName)
                if len(ModName)<2:
                    print "!???", TagName, ModName
                if not Mod and ModName[-2]==">": #Mutation is annotated as "a->g", etc.
                    Mod = PTModClass(ModName)
                    Mod.Mass = Global.AminoMass[ModName[-1].upper()] - Global.AminoMass[ModName[0].upper()]
                if not Mod and ModName[0] in ("-","+"):
                    ModName = ModName.replace("(","")
                    # Keep a cache of "mass mods":
                    ModMass = int(ModName)
                    Mod = MassPTMods.get(ModMass, None)
                    if not Mod:
                        Mod = PTModClass(ModName)
                        Mod.Mass = ModMass
                        MassPTMods[ModMass] = Mod
                if Mod:
                    Pos = len(Peptide.Aminos) - 1
                    if not Peptide.Modifications.has_key(Pos):
                        Peptide.Modifications[Pos] = []
                    Peptide.Modifications[Pos].append(Mod)
                else:
                    print "** Warning: Unknown mod '%s' in '%s'"%(ModName, TagName)
    except:
        print TagName
        raise
    Peptide.ComputeMasses()
    return Peptide

class AminoClass:
    def __init__(self, Name, ShortName, Abbreviation, LeftMass, RightMass):
        self.Name = Name # "Histidine"
        self.ShortName = ShortName # "His"
        self.Abbreviation = Abbreviation # "H"
        self.LeftMass = LeftMass
        self.RightMass = RightMass
        self.RequiredModification = None 
        
def LoadAminoAcids():
    """
    Read in the masses of all amino acids.
    Populate dictionaries AminoMass, AminoMassRight and list AminoMasses
    """
    File = open(checkConfigPath("AminoAcidMasses.txt"),"r")
    for FileLine in File.xreadlines():
        # Line is whitespace-delimited.  Pieces are:
        # Long, short, abbrev, left-mass, right-mass
        # Example: "Glycine Gly G 57.02146 57.0520"
        FileLine = FileLine.strip()
        if FileLine[0] == "#":
            continue
        Bits = FileLine.split(" ")
        if len(Bits)<5:
            continue
        LeftMass = float(Bits[3])
        RightMass = float(Bits[4])
        Global.AminoMass[Bits[2]] = LeftMass
        Global.AminoMassRight[Bits[2]] = RightMass
        Global.AminoMasses.append(LeftMass)
        # Put the Amino object into Global.AminoAcids:
        Amino = AminoClass(Bits[0], Bits[1], Bits[2], LeftMass, RightMass)
        Global.AminoAcids[Amino.Abbreviation] = Amino
    File.close()
    Global.AminoMasses.sort()
    

def DebugPrintPTMods():
    Keys = Global.PTMods.keys()
    Keys.sort()
    print "--PTMods--"
    for Key in Keys:
        PTMod = Global.PTMods[Key]
        BaseString = ""
        for Base in PTMod.Bases.keys():
            BaseString += Base
        print "  %s mass %s bases '%s'"%(PTMod.Name, PTMod.Mass, BaseString)
    print "-----"
    
class IonClass:
    """
    Each IonClass corresponds to an ion type, such as b or y-nh3.
    Each spectral peak gives rise to one PRM peak for each ion type;
    these PRM peaks remember their associated ion class
    """
    def __init__(self, Name):
        self.Name = Name
        self.Opposite = None
        self.Charge = 1
        self.Score = 1.0
    def __str__(self):
        return "<ion '%s'>"%self.Name
    def GetPRMMass(self, Mass, ParentMass):
        """
        Returns the prm peak for a spectrum peak of the given mass.  For instance,
        for b ions, GetPRMMass() returns the peak mass minus 1.  (Because the spectral peak
        appears 1amu to the right of the actual prefix mass)
        """
        return None
    def GetPeakMass(self, Mass, ParentMass):
        """
        Returns the peak for a PRM of the given mass.  Inverse of GetPRMMass.
        For instance, for b ions, GetPeakMass() returns the PRM plus 1.
        """
        return None



AllIons = []
Global.AllIonDict = {}
def DefineIons():
    """
    Define all the ion types we care about.  
    (This function is repetitive, but easy enough to maintain since the zoo of ion types
    is pretty small...the scores should be in a datafile, though!)
    """
    IonB = IonClass("b")
    IonB.GetPeakMass = lambda L, P:L+1
    IonB.GetPRMMass = lambda M, P:M-1
    AllIons.append(IonB)
    #
    IonBH = IonClass("b-h2o")
    IonBH.GetPeakMass = lambda L, P:L-17
    IonBH.GetPRMMass = lambda M, P:M+17
    AllIons.append(IonBH)
    #
    IonBN = IonClass("b-nh3")
    IonBN.GetPeakMass = lambda L, P:L-16
    IonBN.GetPRMMass = lambda M, P:M+16
    AllIons.append(IonBN)
    #
    Ion = IonClass("b-h2o-h2o")
    Ion.GetPeakMass = lambda L, P:L-17-18
    Ion.GetPRMMass = lambda M, P:M+17+18
    AllIons.append(Ion)
    #
    Ion = IonClass("b-h2o-nh3")
    Ion.GetPeakMass = lambda L, P:L-16-18
    Ion.GetPRMMass = lambda M, P:M+16+18
    AllIons.append(Ion)
    #
    Ion = IonClass("b-p'")
    Ion.GetPeakMass = lambda L, P:L-79
    Ion.GetPRMMass = lambda M, P:M+79
    AllIons.append(Ion)
    #
    Ion = IonClass("b-p")
    Ion.GetPeakMass = lambda L, P:L-97
    Ion.GetPRMMass = lambda M, P:M+97
    AllIons.append(Ion)
    #
    Ion = IonClass("b-p-h2o")
    Ion.GetPeakMass = lambda L, P:L-97-18
    Ion.GetPRMMass = lambda M, P:M+97+18
    AllIons.append(Ion)
    #
    Ion = IonClass("b-p-nh3")
    Ion.GetPeakMass = lambda L, P:L-97-17
    Ion.GetPRMMass = lambda M, P: M+97+17
    AllIons.append(Ion)
    # for oxidized methionine:
    Ion = IonClass("b-*")
    Ion.GetPeakMass = lambda L, P:L-63
    Ion.GetPRMMass = lambda M, P:M+63
    AllIons.append(Ion)
    #
    IonY = IonClass("y")
    IonY.GetPeakMass = lambda L, P:P-L
    IonY.GetPRMMass = lambda M, P:P-M
    AllIons.append(IonY)
    #
    IonYH = IonClass("y-h2o")
    IonYH.GetPeakMass = lambda L, P:P-(L+18)
    IonYH.GetPRMMass = lambda M, P:(P-M)-18
    AllIons.append(IonYH)
    #
    IonYN = IonClass("y-nh3")
    IonYN.GetPeakMass = lambda L, P:P-(L+17)
    IonYN.GetPRMMass = lambda M, P:(P-M)-17
    AllIons.append(IonYN)
    #
    Ion = IonClass("y-h2o-nh3")
    Ion.GetPeakMass = lambda L, P:P-(L+17+18)
    Ion.GetPRMMass = lambda M, P:(P-M)-17-18
    AllIons.append(Ion)
    #
    Ion = IonClass("y-h2o-h2o")
    Ion.GetPeakMass = lambda L, P:P-(L+18+18)
    Ion.GetPRMMass = lambda M, P:(P-M)-18-18
    AllIons.append(Ion)
    #
    Ion = IonClass("y-p'")
    Ion.GetPeakMass = lambda L, P:(P-L)-80
    Ion.GetPRMMass = lambda M, P:P-(M+80)
    AllIons.append(Ion)
    #
    Ion = IonClass("y-p")
    Ion.GetPeakMass = lambda L, P:(P-L)-98
    Ion.GetPRMMass = lambda M, P:P-(M+98)
    AllIons.append(Ion)
    # For oxidized methionine:
    Ion = IonClass("y-*")
    Ion.GetPeakMass = lambda L, P:(P-L)-64
    Ion.GetPRMMass = lambda M, P:P-(M+64)
    AllIons.append(Ion)
    #
    IonA = IonClass("a")
    IonA.GetPeakMass = lambda L, P:L-27
    IonA.GetPRMMass = lambda M,P:M+27
    AllIons.append(IonA)
    #
    IonAN = IonClass("a-nh3")
    IonAN.GetPeakMass = lambda L, P:L-27-17
    IonAN.GetPRMMass = lambda M,P:M+27+17
    AllIons.append(IonAN)
    #
    IonAH = IonClass("a-h2o")
    IonAH.GetPeakMass = lambda L, P:L-27-18
    IonAH.GetPRMMass = lambda M,P:M+27+18
    AllIons.append(IonAH)
    #
    Ion = IonClass("b2")
    Ion.GetPeakMass = lambda L,P:(L/2)+1
    Ion.GetPRMMass = lambda M,P:(M-1)*2
    Ion.Charge = 2
    AllIons.append(Ion)
    #
    Ion = IonClass("b2-h2o")
    Ion.GetPeakMass = lambda L,P:(L/2)+1 - 9
    Ion.GetPRMMass = lambda M,P:(M-1)*2 + 18
    Ion.Charge = 2
    AllIons.append(Ion)
    #
    Ion = IonClass("b2-nh3")
    Ion.GetPeakMass = lambda L,P:(L/2)+1 - 8.5
    Ion.GetPRMMass = lambda M,P:(M-1)*2 + 17 
    Ion.Charge = 2
    AllIons.append(Ion)
    #
    Ion = IonClass("b2-nh3-h2o")
    Ion.GetPeakMass = lambda L,P:(L/2)+1 - 17.5
    Ion.GetPRMMass = lambda M,P:(M-1)*2 + 35 
    Ion.Charge = 2
    AllIons.append(Ion)
    #
    Ion = IonClass("b2-p")
    Ion.GetPeakMass = lambda L,P:(L/2)+1 - 49
    Ion.GetPRMMass = lambda M,P:(M-1)*2 + 98 
    Ion.Charge = 2
    AllIons.append(Ion)
    #
    Ion = IonClass("y2")
    Ion.GetPeakMass = lambda L,P:(P-L+1.0078)/2
    Ion.GetPRMMass = lambda M,P:P - (M*2 - 1.0078)
    Ion.Charge = 2
    AllIons.append(Ion)
    #
    Ion = IonClass("y2-h2o")
    Ion.GetPeakMass = lambda L,P:(P-L+1.0078 - 18)/2
    Ion.GetPRMMass = lambda M,P:P - (M*2 - 1.0078 + 18)
    Ion.Charge = 2
    AllIons.append(Ion)
    #
    Ion = IonClass("y2-nh3")
    Ion.GetPeakMass = lambda L,P:(P-L+1.0078 - 17)/2
    Ion.GetPRMMass = lambda M,P:P - (M*2 - 1.0078 + 17)
    Ion.Charge = 2
    AllIons.append(Ion)
    #
    Ion = IonClass("y2-nh3-h2o")
    Ion.GetPeakMass = lambda L,P:(P-L+1.0078 - 17 - 18)/2
    Ion.GetPRMMass = lambda M,P:P - (M*2 - 1.0078 + 17 + 18)
    Ion.Charge = 2
    AllIons.append(Ion)
    #
    Ion = IonClass("y2-p")
    Ion.GetPeakMass = lambda L,P:(P-L+1.0078 - 98)/2
    Ion.GetPRMMass= lambda M,P:P - (M*2 - 1.0078 + 98)
    Ion.Charge =2
    AllIons.append(Ion)
    # For oxidized methionine:
    Ion = IonClass("y2-*")
    Ion.GetPeakMass = lambda L,P:(P-L+1.0078 - 64)/2
    Ion.GetPRMMass = lambda M,P:P - (M*2 - 1.0078 + 64)
    Ion.Charge = 2
    AllIons.append(Ion)
    #
    
    for Ion in AllIons:
        Global.AllIonDict[Ion.Name] = Ion

def GetMass(Str):
    "Return the mass of a string of amino acids.  Useful in interactive mode."
    Mass = 0
    for Char in Str:
        Mass += Global.AminoMass[Char]
        Mass += Global.FixedMods.get(Char, 0)
    return Mass

        

def GetIsotopePatterns():
    Global.IsotopeWeights = {}
    File = open(checkConfigPath("IsotopePatterns.txt"), "r")
    for FileLine in File.xreadlines():
        Bits = FileLine.split("\t")
        if len(Bits) < 2:
            continue
        Global.IsotopeWeights[int(Bits[0])] = float(Bits[1])
        

INITIALIZED = 0
DummyIon = None
def Initialize():
    global INITIALIZED
    global DummyIon
    if INITIALIZED:
        return 
    DefineIons()

    # dummy ion type, for the spectral edge peaks we put at mass 0 and at parent-mass:
    DummyIon = IonClass("")
    DummyIon.GetPeakMass = lambda L,P:L
    DummyIon.GetPRMMass = lambda M,P:M

    # Do this initialization once, up front:
    LoadAminoAcids()
    LoadPTMods()
    LoadModifications()
    GetIsotopePatterns()
    INITIALIZED = 1

def MakeDirectory(Dir):
    if os.path.exists(Dir):
        return 
    try:
        os.makedirs(Dir)
    except:
        raise
    
def checkConfigPath(path):
    if os.path.exists(path):
        return path

    alt = os.path.join(sys.path[0],path)
    if os.path.exists(alt):
        return alt
    alt = os.path.join( "..",path)
    if os.path.exists(alt):
        return alt
   
    raise IOError("Can't open %s or %s" % (path,alt))
