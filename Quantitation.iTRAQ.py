UsageInfo="""Quantitation.iTRAQ.py
uses CID/PQD files.  Pulls the identifications (as defined by Inspect)
and then gets their quantitation from the neighboring PQD spectrum.

Required Options:
 -r [FileName] Results from Inspect, expected to be pvalued.
 -m [Directory] Directory (not file) where spectra files are stored

"""

"""Future
 -o [int]Ordering of CID/PQD.  Value of 0 =  MS1 -> CID -> PQD
                               Value of 1 =  MS1 -> PQD -> CID

"""

"""I asked Zhouxin, and here's the deal.  Each m/z was fragmented by two
different techniques: CID and PQD.  So the first spectrum of the
double is CID followed by another spectrum for the same m/z produced
by PQD.  They do this because CID does not get the low masses and PQD
does.  Apparently the low limit for CID is 0.28 * precursor

Table 1-1 iTRAQ Reagent specifications
Specification Reagent
114     115     116     117
Exact mass added to the peptide per primary amine
144.1059 144.0996 144.1021 144.1021
Monoisotopic MH+ of the Reporter Group (fragment seen by MS/MS analysis)
114.1112 115.1083 116.1116 117.1150
"""

import os
import sys
import getopt
import ResultsParser
import xml.sax.handler
import xml.sax
import MSSpectrum

## auxiliary for the mzxml files
class XMLHandler(xml.sax.handler.ContentHandler):
    def __init__(self):
        self.inOffset = 0
        self.CIDScan = 0
        self.mapping = {}
        self.CID2PQDMapping = {} # key = CID ScanNum, value = PQD Scan
    def startElement(self, name, attributes):
        if name == "offset":
            self.buffer = ""
            self.scan = attributes["id"]
            self.inOffset = 1
        if name == "scan":
            Level = int(attributes["msLevel"])
            if not Level == 2:
                self.CIDScan =0 #reset here too
                return
            #now let's get the mapping
            ScanNum = int(attributes["num"])
            if not self.CIDScan:
                # I dont' currently have any scan number stored in my CID register
                self.CIDScan = ScanNum
            else:
                # i do have a CID scan, since PQD always follows CID, then this must be the PQD
                PQD = ScanNum
                self.CID2PQDMapping[self.CIDScan] = PQD
                self.CIDScan =0
                
    def characters(self, data):
        if self.inOffset:
            self.buffer += data
    def endElement(self, name):
        if name == "offset":
            self.inOffset = 0
            self.mapping[self.scan] = self.buffer

class AbacusClass(ResultsParser.ResultsParser):
    def __init__ (self):
        self.OutputPath = None
        self.InspectResults = None
        self.SpectraDir = None
        self.InspectHits = {} #key = filename #value = scan of Inspect hit
        self.iTRAQMasses = []
        self.ScanOffset = {} #(file, scan) = offset
        self.CID2PQDMapping = {} #(file, CIDscan) = PQDScan
        # set up Dictionary of iTRAQ reporter masses
        self.ReporterMasses = {}
        self.ReporterMasses[114] = 114.1112
        self.ReporterMasses[115] = 115.1083
        self.ReporterMasses[116] = 116.1116
        self.ReporterMasses[117] = 117.1150
        ResultsParser.ResultsParser.__init__(self)
        
        
    def Main(self):
        self.GetSpectralPairs()
        self.ProcessResultsFiles(self.InspectResults, self.ParseResults)


    def GetSpectralPairs(self):
        for File in os.listdir(self.SpectraDir):
            (FileName, Extension) = os.path.splitext(File)
            Extension = Extension.lower()
            if not Extension == ".mzxml":
                continue
            print File
            FullPath = os.path.join(self.SpectraDir, File)
            self.GetSpectralPairsFile(FullPath)

    def GetSpectralPairsFile(self, FilePath):
        """Given a mzxml file, read through and get the spectrum pairs
        of neighboring scans with the same parent mass.  Put them in a
        dictionary of some kind.  Also make sure to get the file offset for such files.
        """
        Parser = xml.sax.make_parser()
        Handler = XMLHandler()
        Parser.setContentHandler(Handler)
        Parser.parse(FilePath)
        for (Scan, Offset) in Handler.mapping.items():
            ScanNumber = int(Scan)
            Offset = int(Offset)
            #print (Scan, Offset)
            Key = (FilePath, ScanNumber)
            self.ScanOffset[Key] = Offset
        for (CID, PQD) in Handler.CID2PQDMapping.items():
            Key = (FilePath, CID)
            self.CID2PQDMapping[Key] = PQD
            #print CID, PQD
            

    def ParseResults(self, FileName):
        Handle = open(FileName, "rb")
        OutHandlePath = FileName.replace(".txt", ".itraq.txt")
        OutHandle = open(OutHandlePath, "wb")
        self.InspectHits[FileName] = [] #empty list for scans
        for Line in Handle.xreadlines():
            if Line[0] == "#":
                Headers = Line.strip()
                OutHandle.write("%s\tPQDScan\t114\t115\t116\t117\n"%Headers)
                continue
            #print ":%s:"%Line
            Bits = Line.strip().split("\t")
            SpectrumPath = Bits[self.Columns.SpectrumFile]
            (Path, SpectrumFile)= os.path.split(SpectrumPath)
            Scan = int(Bits[self.Columns.ScanNumber])
            ByteOffset = int (Bits[self.Columns.FileOffset])
            SpectrumPath = os.path.join(self.SpectraDir, SpectrumFile)
            ## now I have the scan of interest.  If this scan is not in my
            ## CID2PQD mapping, then it must be a PQD itself. 
            ## once we have the PQD scan, we get the offset
            ## then we open up the spectrum and pull out the peaks
            PQDScan = self.CID2PQDMapping.get((SpectrumPath, Scan), None);
            if not PQDScan:
                PQDScanOffset = ByteOffset
                PQDScan = Scan
            else:
                PQDScanOffset = self.ScanOffset[(SpectrumPath, PQDScan)]
                
            iTRAQPeaks = self.GetiTRAQPeaks(SpectrumPath, PQDScanOffset)
            NewLine = "\t".join(Bits)
            NewLine += "\t%s"%PQDScan
            for Peak in iTRAQPeaks:
                NewLine += "\t%s"%Peak
            OutHandle.write("%s\n"%NewLine)
        Handle.close()
        OutHandle.close()

    def GetiTRAQPeaks(self, SpectrumFilePath, ByteOffset):
        """Given the spectrum and offset, open up the spectrum and then find the iTRAQ peaks.
        """
        MassTolerance = 0.3 # daltons I want peaks to be close to the right mass
        Peaks = []
        Spectrum = MSSpectrum.SpectrumClass()
        SpectrumFile = open(SpectrumFilePath, "rb")
        SpectrumFile.seek(ByteOffset)
        Spectrum.ReadPeaksFromFile(SpectrumFile, SpectrumFilePath)
        #get the peaks
        Peak114 = Spectrum.GetPeak(self.ReporterMasses[114], 0.3)
        Peak115 = Spectrum.GetPeak(self.ReporterMasses[115], 0.3)
        Peak116 = Spectrum.GetPeak(self.ReporterMasses[116], 0.3)
        Peak117 = Spectrum.GetPeak(self.ReporterMasses[117], 0.3)
        #some bulletproofing
        if not Peak114:
            Peak114Intensity = 0
        else:
            Peak114Intensity = Peak114.Intensity
        if not Peak115:
            Peak115Intensity = 0
        else:
            Peak115Intensity = Peak115.Intensity
        if not Peak116:
            Peak116Intensity = 0
        else:
            Peak116Intensity = Peak116.Intensity
        if not Peak117:
            Peak117Intensity = 0
        else:
            Peak117Intensity = Peak117.Intensity
            
        SpectrumFile.close()
        return (Peak114Intensity, Peak115Intensity, Peak116Intensity, Peak117Intensity)
        
    
    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "m:r:o:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-m":
                # -m spectra file(s)
                self.SpectraDir = Value
            if Option == "-w":
                # -r results file(s)
                self.OutputPath = Value
            if Option == "-r":
                self.InspectResults = Value
            if Option == "-o":
                self.Ordering = int (Value)
        if not OptionsSeen.has_key("-m") or not OptionsSeen.has_key("-r"):
            print UsageInfo
            sys.exit(1)


if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(psyco not found - running in non-optimized mode)"
    DoStuff = AbacusClass()
    DoStuff.ParseCommandLine(sys.argv[1:])
    DoStuff.Main()        
