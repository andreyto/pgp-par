""" ScorpionUtils.py
This holds a bunch of data that is swapped inbetween all the scorpion scripts
and includes a class that can be subclassed so that everyone has access to the
same easy information.
"""

class PhosphorylationIons:
    #columns
    Y = 0
    B = 1
    Y2 = 2
    B2 = 3
    YMinusP = 4
    BMinusP = 5
    YMinusH2O = 6
    BMinusH2O = 7
    YMinusNH3 = 8
    BMinusNH3 = 9
    BMinusPH2O = 10
    BMinusPNH3 = 11
    Y2MinusH2O = 12
    B2MinusH2O = 13
    Y2MinusNH3 = 14
    B2MinusNH3 = 15
    Y2MinusP = 16
    B2MinusP = 17
    NumIons = 18
    IonTypeNames = ["y","b","y2","b2","y-p","b-p","y-h2o","b-h2o","y-nh3","b-nh3", "b-p-h2o","b-p-nh3", 
                        "y2-h2o", "b2-h2o", "y2-nh3", "b2-nh3", "y2-p","b2-p"]
    #LevelingFunction = [IntensityLeveling for ALL]
    #three columns for each: Mass, Ratio, Skew
    #note starting with a TAB
    String = "\tyMass\tyRatio\tySkew\tbMass\tbRatio\tbSkew"
    String += "\ty2Mass\ty2Ratio\ty2Skew\tb2\tb2\tb2"
    String += "\ty-p\ty-p\ty-p\tb-p\tb-p\tb-p"
    String += "\ty-h2o\ty-h2o\ty-h2o\tb-h2o\tb-h2o\tb-h2o"
    String += "\ty-nh3\ty-nh3\ty-nh3\tb-nh3\tb-nh3\tb-nh3"
    String += "\tb-p-h2o\tb-p-h2o\tb-p-h2o\tb-p-nh3\tb-p-nh3\tb-p-nh3"
    String += "\ty2-h2o\ty2-h2o\ty2-h2o\tb2-h2o\tb2-h2o\tb2-h2o"
    String += "\ty2-nh3\ty2-nh3\ty2-nh3\tb2-nh3\tb2-nh3\tb2-nh3"
    String += "\ty2-p\ty2-p\ty2-p\tb2-p\tb2-p\tb2-p"

class RegularIons:
    #columns
    Y = 0
    B = 1
    Y2 = 2
    B2 = 3
    YMinusH2O = 4
    BMinusH2O = 5
    YMinusNH3 = 6
    BMinusNH3 = 7
    Y2MinusH2O = 8
    B2MinusH2O = 9
    Y2MinusNH3 = 10
    B2MinusNH3 = 11
    A = 12
    NumIons = 13
    IonTypeNames = ["y","b","y2","b2","y-h2o","b-h2o","y-nh3","b-nh3",  
                        "y2-h2o", "b2-h2o", "y2-nh3", "b2-nh3", "a"]
    #LevelingFunction = [IntensityLeveling for ALL]
    #three columns for each: Mass, Ratio, Skew
    #note starting with a TAB
    String = "\tyMass\tyRatio\tySkew\tbMass\tbRatio\tbSkew"
    String += "\ty2Mass\ty2Ratio\ty2Skew\tb2\tb2\tb2"
    String += "\ty-h2o\ty-h2o\ty-h2o\tb-h2o\tb-h2o\tb-h2o"
    String += "\ty-nh3\ty-nh3\ty-nh3\tb-nh3\tb-nh3\tb-nh3"
    String += "\ty2-h2o\ty2-h2o\ty2-h2o\tb2-h2o\tb2-h2o\tb2-h2o"
    String += "\ty2-nh3\ty2-nh3\ty2-nh3\tb2-nh3\tb2-nh3\tb2-nh3"
    String += "\ta\ta\ta"

class BreakInfoColumns:
    Peptide = 0
    File = 1
    BreakIndex = 2
    PrefixAA = 3
    SuffixAA = 4
    NBase = 5
    NAcid = 6
    CBase = 7
    CAcid = 8
    Sector = 9
    Charge = 10
    PhosOnN = 11 #whether the phosphate resides on the nterminus of the break (b ion) or not
                #set to zero for non-phos searches. thus containing no information
    String = "Peptide\tFile\tBreakIndex\tPrefixAA\tSuffixAA\tNBase\tNAcid\tCBase\tCAcid\tSector\tCharge\tPhosOnN"
    LevelingScheme = ["None","None","None","AminoAcid","AminoAcid","AcidBase","AcidBase","AcidBase","AcidBase","Sector","None","None"]
    ColumnNames = String.split("\t")
    NumCols = 12
    
class ScorpionParser:
    def __init__(self, *args, **kw):
        self.SpecialCase = None
        self.FeaturedIons = RegularIons
        self.BreakInfoColumns = BreakInfoColumns
        #self.FeaturedIons = PhosphorylationIons

    def SetSpecialCase(self,SpecialCase):
        #print "I'm setting a special case, %s"%SpecialCase
        if SpecialCase == "Phosphorylation":
            self.FeaturedIons = PhosphorylationIons
            self.SpecialCase = SpecialCase
           