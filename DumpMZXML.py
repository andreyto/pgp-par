"""
Utility script:
- Parse scans from an .mzXML file
- Write scans out as .dta files
"""
import os
import sys
import getopt
import xml
import xml.sax
import ParseXML
import MSSpectrum

# Important note:
# GetByteOffset is tracked only ONCE for all parser objects.
# You can't pause parsing of one xml file, parse another, and
# then refer back to valid byte offsets!
# That's why we build a list of scans first.

UsageInfo ="""
DumpMZXML.py - Convert peaks from mzXML format to .dta or .mgf
Arguments:
-r [InputPath]
-w [InputPath]: Output .mgf file name, or .dta directory
-m: Output an .mgf file (default)
-d: Output .dta files
"""

class PeakFindParseStates:
    Skipping = 0
    
class MZXMLScanFinder(ParseXML.XMLDictionaryHandler):
    def __init__(self, MZXMLFile):
        self.State = PeakFindParseStates.Skipping
        self.StartHandlers = {"scan":self.StartScan,
                              "precursorMz":self.StartPrecursorMZ,
                              }
        self.EndHandlers = {"scan":self.EndScan}
        self.MZXMLFile = MZXMLFile
        self.MSLevel = 0
        self.PeakCount = 0
        self.AllScans = []
        ParseXML.XMLDictionaryHandler.__init__(self)
    def HandleCharacters(self, String):
        pass
    def StartScan(self, Attributes):
        self.ScanNumber = int(Attributes.get("num", 0))
        self.PeakCount = int(Attributes.get("peaksCount", 0))
        self.MSLevel = int(Attributes.get("msLevel", 0))
        self.ScanStartPos = self.Parser._parser.CurrentByteIndex
        #print "Scan %s starts at pos %s (%s:%s)"%(self.ScanNumber, self.ScanStartPos, self.Parser.getLineNumber(), self.Parser.getColumnNumber())
    def StartPrecursorMZ(self, Attributes):
        pass
    def EndScan(self):
        if self.MSLevel > 1 and self.PeakCount >= 10:
            self.AllScans.append((self.ScanStartPos, self.ScanNumber))
            #print "Scan %s complete at %s:%s"%(self.ScanNumber, self.Parser.getLineNumber(), self.Parser.getColumnNumber())
            #Spectrum = MSSpectrum.SpectrumClass()
            #self.MZXMLFile.seek(self.ScanStartPos)
            #ParseXML.GetSpectrumPeaksMZXML(Spectrum, self.MZXMLFile)
            #self.SpectrumCallback(Spectrum, self.ScanNumber)
        self.MSLevel = 0
        self.PeakCount = 0

class MZXMLDumper:
    def __init__(self):
        self.OutputFileName = "Scans.mgf"
        self.OutputType = "mgf"
        self.InputFileName = None
    def ParseCommandLine(self, Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:w:dm")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-r":
                self.InputFileName = Value
            elif Option == "-w":
                self.OutputFileName = Value
            elif Option == "-m":
                self.OutputType = "mgf"
            elif Option == "-d":
                self.OutputType = "dta"
    def Main(self):
        if not self.InputFileName:
            print UsageInfo
            sys.exit(-1)
        if self.OutputType == "mgf":
            CallbackFunction = self.DumpPeaksMGF
            self.OutputFile = open(self.OutputFileName, "wb")
        elif self.OutputType == "dta":
            try:
                os.makedirs(self.OutputFileName)
            except:
                pass
            CallbackFunction = self.DumpPeaksDTA
        else:
            raise ValueError, "Illegal OutputType"
        MZXMLFile = open(self.InputFileName, "rb")
        Finder = MZXMLScanFinder(MZXMLFile)
        SAXParser = xml.sax.make_parser()
        Finder.Parser = SAXParser
        SAXParser.setContentHandler(Finder)
        print "Parse scans from mzxml..."
        SAXParser.parse(MZXMLFile)
        print "Dump %s scans..."%len(Finder.AllScans)
        for (FilePos, ScanNumber) in Finder.AllScans:
            #print "Dump scan %s from %s..."%(ScanNumber, FilePos)
            Spectrum = MSSpectrum.SpectrumClass()
            MZXMLFile.seek(FilePos)
            ParseXML.GetSpectrumPeaksMZXML(Spectrum, MZXMLFile)
            #self.SpectrumCallback(Spectrum, ScanNumber)
            CallbackFunction(Spectrum, ScanNumber)
        MZXMLFile.close()
    def DumpPeaksMGF(self, Spectrum, ScanNumber):
        Spectrum.WriteMGFPeaks(self.OutputFile, "Scan %s"%ScanNumber, ScanNumber)
    def DumpPeaksDTA(self, Spectrum, ScanNumber):
        LatePeakIndex = int(len(Spectrum.Peaks) * 0.80)
        if Spectrum.Peaks[LatePeakIndex].Mass >= Spectrum.PrecursorMZ:
            Charge1Flag = 0
        else:
            Charge1Flag = 1
        if Charge1Flag:
            FilePath = os.path.join(self.OutputFileName, "Scan%s.1.dta"%ScanNumber)
            Spectrum.SetCharge(1)
            Spectrum.WritePeaks(FilePath)
        else:
            FilePath = os.path.join(self.OutputFileName, "Scan%s.2.dta"%ScanNumber)
            Spectrum.SetCharge(2)
            Spectrum.WritePeaks(FilePath)
            FilePath = os.path.join(self.OutputFileName, "Scan%s.3.dta"%ScanNumber)
            Spectrum.SetCharge(3)
            Spectrum.WritePeaks(FilePath)
        
        
if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(psyco not found)"
    SaddyDumpington = MZXMLDumper()
    SaddyDumpington.ParseCommandLine(sys.argv[1:])
    SaddyDumpington.Main()