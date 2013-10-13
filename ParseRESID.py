import xml.dom.minidom
import string
import sys
import os
import traceback

def GetText(NodeList):
    TextString = ""
    for Node in NodeList:
        if Node.nodeType == Node.TEXT_NODE:
            TextString += Node.data
    return TextString

def ParseXML():
    File = open(sys.argv[1], "rb")
    Text = File.read()
    File.close()
    DOM = xml.dom.minidom.parseString(Text)
    Entries = DOM.getElementsByTagName("Entry")
    for Entry in Entries:
        EntryID = str(Entry.getAttributeNode("id").value)
        # Get the name:
        Name = Entry.getElementsByTagName("Name")[0]
        Name = GetText(Name.childNodes)
        # Get the mass:
        Mass = 0
        DeltaMass = 0
        CorrectionBlocks = Entry.getElementsByTagName("CorrectionBlock")
        if CorrectionBlocks:
            CorrectionBlock = CorrectionBlocks[0]
            WeightNodes = CorrectionBlock.getElementsByTagName("Weight")
            if WeightNodes:
                DeltaMass = GetText(WeightNodes[0].childNodes)
        FormulaBlocks = Entry.getElementsByTagName("FormulaBlock")
        if FormulaBlocks:
            FormulaBlock = FormulaBlocks[0]
            WeightNodes = FormulaBlock.getElementsByTagName("Weight")
            if WeightNodes:
                Mass = GetText(WeightNodes[0].childNodes)
        # Get the amino acid specificity:
        AANodes = Entry.getElementsByTagName("SequenceSpec")
        if AANodes:
            AADict = {}
            AAString = GetText(AANodes[0].childNodes)
##            for Char in AAString:
##                if Char in "ACDEFGHIKLMNPQRSTUVWY":
##                    AADict[Char] = 1
##            Keys = AADict.keys()
##            Keys.sort()
##            AA = string.join(Keys, "")
            AA = AAString
            
        print "%s\t%s\t%s\t%s\t%s\t"%(EntryID, Name, AA, Mass, DeltaMass)
ParseXML()