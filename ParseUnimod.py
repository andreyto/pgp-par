"""
Parse the Unimod database dump from XML into something useful
"""
from xml.dom.minidom import parse, parseString

Classifications = ["", "", "Post-translational", "Co-translational", "Pre-translational",
                   "Chemical derivative", "Artefact", "N-linked glycosylation", "O-linked glycosylation",
                   "Other glyco", "synthetic protecting group", "Isotopic label", "Non-standard residue",
                   "Multiple", "Other"]

def GetText(NodeList):
    String = ""
    for Node in NodeList:
        if Node.nodeType == Node.TEXT_NODE:
            String += Node.data
    return String

class ModClass:
    def __init__(self, Mass, Name):
        self.Mass = Mass
        self.Name = Name
        self.Specificity = [] # tuples of form (AA, Terminus)
        #self.AA = None
        #self.Terminus = None

# Print header:
print "Source\tID\tMassDelta\tName\tCodeName\tAA\tTerminus\tClass\tNotes\t"

File = open("Unimod.xml", "rb")
Text = File.read()
File.close()
DOM = parseString(Text)
Modifications = {}
Notes = ""
# Parse modifications:
for Row in DOM.getElementsByTagName("modifications_row"):
    Mass = float(Row.getAttribute("avge_mass"))
    Name = str(Row.getAttribute("full_name"))
    CodeName = str(Row.getAttribute("code_name"))
    #GetText(Row.getElementsByTagName("full_name")[0])
    RecordID = int(Row.getAttribute("record_id"))
    #int(GetText(Row.getElementsByTagName("record_id")[0]))
    Modifications[RecordID] = ModClass(Mass, Name)
    NotesNode = Row.getElementsByTagName("misc_notes")[0]
    Notes = GetText(NotesNode.childNodes)
    Notes = Notes.replace("\r"," ").replace("\n", " ")
    Modifications[RecordID].Notes = Notes
    Modifications[RecordID].CodeName = CodeName
    #print Mass, Name
    

# Parse specificity:    
for Row in DOM.getElementsByTagName("specificity_row"):
    #RecordID = int(Row.getAttribute("record_id")) #int(GetText(Row.getElementsByTagName("record_id")[0]))
    RecordID = int(Row.getAttribute("mod_key"))
    if not Modifications.has_key(RecordID):
        print "** Warning: Specificity found for nonexistent modification %s"%RecordID
        continue
    Mod = Modifications[RecordID]
    Terminus = ""
    AA = Row.getAttribute("one_letter") #GetText(Row.getElementsByTagName("one_letter")[0])
    if AA in ("N-term", "C-term"):
        Terminus = AA
        AA = ""
    ClassNumber = int(Row.getAttribute("classifications_key"))
    Classification = Classifications[ClassNumber]
    Mod.Specificity.append((AA, Terminus, Classification))
##    if Mod.AA or Mod.Terminus:
##        print "ADDITIONAL specificity for %s %s"%(Mod.Name, Mod.Mass)
##    else:
##        Mod.AA = AA
##        Mod.Terminus = Terminus
##        print "%s%s specific to %s%s"%(Mod.Name, Mod.Mass, AA, Terminus)
for (ModID, Mod) in Modifications.items():
    for Specificity in Mod.Specificity:
        Str = "UniMOD\t%s\t%s\t%s\t%s\t"%(ModID, Mod.Mass, Mod.Name, Mod.CodeName)
        Str += "%s\t%s\t%s\t%s\t"%(Specificity[0], Specificity[1], Specificity[2], Mod.Notes)
        print Str
    
    