"""
Build an MZXML file, given spectra in some other crazy format.
"""
import os
import sys
import getopt
import shutil
import MSSpectrum

UsageInfo = """
BuildMZXML: Merge several spectra into mzXML format.
Arguments:
-r [DIR]: Directory where input spectra are stored
-w [FILE]: Output file; defaults to Spectra.mzXML

"""

class MZXMLBuilder:
    def __init__(self):
        self.NextScanNumber = 1
        self.ScanCount = 0
    def ParseCommandLine(self, Arguments):
        self.InputDir = None
        self.OutputFileName = "Spectra.mzXML"
        (Options, Args) = getopt.getopt(Arguments, "r:w:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-r":
                self.InputDir = Value
            elif Option == "-w":
                self.OutputFileName = Value
        if not self.InputDir:
            print UsageInfo
            sys.exit(-1)
    def AddDTAFile(self, Path):
        Spectrum = MSSpectrum.SpectrumClass()
        Spectrum.ReadPeaks(Path)
        Spectrum.WriteMZXMLPeaks(self.OutputFile, self.NextScanNumber)
        self.NextScanNumber += 1
        self.ScanCount += 1
    def BuildMZXMLFile(self):
        self.TempOutputFileName = self.OutputFileName + ".tmp"
        self.OutputFile = open(self.TempOutputFileName, "wb")
        # Write body, while building up our footer index string:
        self.Footer = "\n<index name=\"scan\">"
        FileNames = os.listdir(self.InputDir)
        FileNames.sort()
        for FileName in FileNames:
            Path = os.path.join(self.InputDir, FileName)
            if os.path.isdir(Path):
                continue # don't recurse into subdirectories
            (Stub, Extension) = os.path.splitext(FileName)
            if Extension == ".dta":
                self.AddDTAFile(Path)
            else:
                print "* Sorry, extension '%s' isn't handled (yet)"%Extension
        self.OutputFile.close()
        # NOW we can write the header, body, and footer:
        Header = """<?xml version="1.0" encoding="ISO-8859-1"?>
        <mzXML
         xmlns="http://sashimi.sourceforge.net/schema_revision/mzXML_2.0"
         xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
         xsi:schemaLocation="http://sashimi.sourceforge.net/schema_revision/mzXML_2.0 http://sashimi.sourceforge.net/schema_revision/mzXML_2.0/mzXML_idx_2.0.xsd">
         <msRun scanCount="%s"
        startTime="PT0.00521S"
        endTime="PT150.057S">
        """%self.ScanCount
        self.OutputFile = open(self.OutputFileName, "wb")
        self.OutputFile.write(Header)
        TempFile = open(self.TempOutputFileName, "rb")
        while 1:
            Block = TempFile.read(1024*1024)
            if not Block:
                break
            self.OutputFile.write(Block)
        TempFile.close()
        os.remove(self.TempOutputFileName)
        # Write footer:
        self.Footer += "\n</index>\n"
        self.OutputFile.write(self.Footer)
        self.OutputFile.close()
        print "Wrote out a total of %s scans."%self.ScanCount
if __name__ == "__main__":
    Builder = MZXMLBuilder()
    Builder.ParseCommandLine(sys.argv[1:])
    Builder.BuildMZXMLFile()