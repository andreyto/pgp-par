UsageInfo = """MascotToInspect.py
Script for converting Mascot results into the Inspect format, using features as
they apply.  This is for converting the protein-based results output.

-r MascotFile/Directory
-m Spectrum Directory
"""


""" Here's the information that I can convert.  Everything else will contain the word "Dummy"
Not in column format
File - line tagged with "MS data path"

Column output starts after the "prot_hit_num" line
Information     Mascot Column        Inspect Column
----------      -------------       --------------
File            Not in columns      0  ResultsParser.Columns.SpectrumFile
Scan Number     6                   1  ResultsParser.Columns.ScanNumber
Score           13                  5  ResultsParser.Columns.MQScore
delta           NOT Included        16 ResultsParser.Columns.DeltaScore
Annotation      16,17,18            2  ResultsParser.Columns.Annotation
Protein         1 (only first)      3  ResultsParser.Columns.ProteinName
Charge          9                   4  ResultsParser.Columns.Charge
ByteOffset      NOT Included        19 ResultsParser.Columns.FileOffset


"""


import os
import sys
import getopt
import ResultsParser
import GetByteOffset

class ConvertClass(ResultsParser.ResultsParser):
    def __init__(self):
        self.MascotDir = None
        self.SpectrumFileName = None
        self.MZXMLDir = None
        #self.ScanOffset = {} scan => byteoffset
        ResultsParser.ResultsParser.__init__(self)

    def ParseMascotFile(self, MascotFilePath):
        ProteinName = None #must keep, because not printed on every line
        SpectrumFileName = None #ditto
        MascotHandle = open(MascotFilePath, "rb")
        (Path, MascotFileName) = os.path.split(MascotFilePath)
        (Stub, Ext) = os.path.splitext(MascotFileName)
        InspectFileName = "%s.inspectstyle.txt"%Stub
        InspectFilePath = os.path.join(Path, InspectFileName)
        InspectHandle = open(InspectFilePath, "wb")
        for Line in MascotHandle.xreadlines():
            Line = Line.strip()
            Bits = list(Line.split("\t"))
            try:
                int (Bits[0])
            except:
                #not an integer. must be meta data
                if Bits[0] == "MS data path":
                    SpectrumFilePath = Bits[1].replace("\"","")
                    SpectrumFileName = os.path.split(SpectrumFilePath)[1]
                    self.SpectrumFileName = os.path.join(self.MZXMLDir, SpectrumFileName)
                    ## get the scan offsets now so we can have them for the results lines
                    Abacus = GetByteOffset.Abacus()
                    self.ScanOffset = Abacus.GetByteOffset(self.SpectrumFileName)
                continue
            #all peptide hits begin their line number with an integer
            ScanNumber = int (Bits[6])
            Score = float (Bits[13])
            Charge = int (Bits[9])
            Annotation = "%s.%s.%s"%(Bits[16], Bits[17], Bits[18])
            FixedAnnotation = Annotation.replace("\"","")
            Aminos = Bits[17].replace("\"","")
            Length = len(Aminos)
            ByteOffset = self.ScanOffset[ScanNumber]
            #print FixedAnnotation
            if Bits[1]: #not printed on every line
                ProteinName = Bits[1].replace("\"","")
            
            Reformatted = "%s\t%s\t%s\t"%(SpectrumFileName, ScanNumber, FixedAnnotation)
            Reformatted += "%s\t%s\t%s\t%s\t"%(ProteinName, Charge, Score, Length)
            #print Reformatted
            Dummy = "Dummy\t"*12
            Reformatted += "%s%s\n"%(Dummy, ByteOffset)
            InspectHandle.write(Reformatted)
        InspectHandle.close()
        MascotHandle.close()

    def Main(self):
        self.ProcessResultsFiles(self.MascotDir, self.ParseMascotFile)
        
    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:m:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-r":
                # -r results file(s)
                if not os.path.exists(Value):
                    print "** Error: couldn't find results file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.MascotDir = Value
            if Option == "-m":
                self.MZXMLDir = Value
        if not OptionsSeen.has_key("-r") :
            print UsageInfo
            sys.exit(1)


if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(psyco not found - running in non-optimized mode)"
    GospelPrinciples = ConvertClass()
    GospelPrinciples.ParseCommandLine(sys.argv[1:])
    GospelPrinciples.Main()           