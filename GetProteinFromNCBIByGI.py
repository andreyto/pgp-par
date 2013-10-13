"""
Grab protein records from NCBI based on their GI number
"""
import urllib
import os
import sys
import string
import traceback

class GrabberClass:
    def __init__(self, OutputFileName):
        self.File = open(OutputFileName, "w")
        print "OPENED:", OutputFileName
    def ParseName(self, Text):
        NamePos = Text.find("DEFINITION") + len("DEFINITION")
        AccessionPos = Text.find("ACCESSION")
        Name = string.join(Text[NamePos:AccessionPos].split(), " ")
        return Name
    def ParseSequence(self, Text):
        OriginPos = Text.find("ORIGIN") + len("ORIGIN")
        Sequence = ""
        Chars = Text[OriginPos:].upper()
        SkippingFlag = 0 
        for Char in Chars:
            # Skip over all tags:
            if Char == ">":
                SkippingFlag = 0
                continue
            if Char == "<":
                SkippingFlag = 1
                continue
            if SkippingFlag:
                continue
            # Keep amino acid letters:
            if Char in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
                Sequence += Char
        return Sequence
    def ParseAccession(self, ID, Text):
        SourcePos = Text.find("DBSOURCE")
        if SourcePos==-1:
            return None
        KeyPos = Text.find("KEYWORDS")
        print "%s\t%s"%(ID, Text[SourcePos+len("DBSOURCE"):KeyPos].replace("\n","  "))
        StartPos = Text.find(">", SourcePos) + 1
        if StartPos >= KeyPos:
            return None
        EndPos = Text.find("<", StartPos)
        return Text[StartPos:EndPos]
    def Grab(self, ID):
        URL = "http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=protein&val=%s"%ID
        URLFile = urllib.urlopen(URL)
        Text = URLFile.read()
        URLFile.close()
        SubFile = open("GI%s.txt"%ID, "w")
        SubFile.write(Text)
        SubFile.close()
        PreBlockStart = Text.lower().find("<pre>") + len("<pre>")
        PreBlockEnd = Text.lower().find("</pre>")
        PreBlock = Text[PreBlockStart:PreBlockEnd]
        Accession = self.ParseAccession(ID, PreBlock)
        Name = "%s %s"%(ID, self.ParseName(PreBlock))
        Sequence = self.ParseSequence(PreBlock)
        if Accession:
            self.File.write(">%s %s\n"%(Accession, Name))
        else:
            self.File.write(">%s\n"%(Name))
        self.File.write("%s\n"%Sequence)
        self.File.flush()
        print "%s\t%s\t%s\t"%(ID, Accession, Name)
        #print "Wrote %s:'%s'"%(ID, Name)
         
if __name__ == "__main__":
    File = open("NCBIGrab2.txt", "r")
    Grabber = GrabberClass("NCBIGrabbed.txt")
    for FileLine in File.xreadlines():
        if FileLine[:3]!="gi|":
            continue
        Bits = FileLine.split("|")
        ID = Bits[1]
        try:
            Grabber.Grab(ID)
        except:
            traceback.print_exc()
    #BatchGrab()