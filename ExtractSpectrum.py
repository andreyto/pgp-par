"""
Extract spectra from larger files.
"""
import getopt
import os
import sys
import string
import traceback
import MSSpectrum

class Extractor:
    def __init__(self):
        self.OutputFileName = "Spectrum.dta" #default name
        self.ReadFileName = None
        self.ReadBytePosition = None
        self.ReadScanNumber = None
    def ExtractByBytePosition(self, FileName, BytePosition, OutputFileName):
        OutputFile = open(OutputFileName, "wb")
        File = open(FileName, "rb")
        File.seek(BytePosition)
        Stuff = File.read(20480)
        Pos = Stuff.find("END IONS")
        if Pos == -1:
            print "Extraction error %s:%s - didn't see END IONS"%(FileName, BytePosition)
            return
        OutputFile.write(Stuff[:Pos])
        OutputFile.write("END IONS\n")
        OutputFile.close()
    def FancyExtract(self, FileNames, ByteOffsets, OutputFileNames):
        for SpectrumIndex in range(len(FileNames)):
            File = open(FileNames[SpectrumIndex], "rb")
            File.seek(ByteOffsets[SpectrumIndex])
            Stuff = File.read(20480)
            TitlePos = Stuff.find("TITLE=")
            if TitlePos == -1:
                print "*** ERROR: No title!", SpectrumIndex
                continue
            Title = Stuff[TitlePos+6:].split()[0]
            Title = Title.replace("16O", "18O")
            File.close()
            OtherFilePath = FileNames[SpectrumIndex].replace("16O.mgf", "18O.mgf")
            File = open(OtherFilePath, "rb")
            Stuff = File.read()
            TitlePos = Stuff.find(Title)
            if TitlePos == -1:
                print "*** ERROR: didn't find title again!", SpectrumIndex, Title
                continue
            SpectrumStartPos = Stuff.rfind("BEGIN IONS", 0, TitlePos)
            SpectrumEndPos = Stuff.find("END IONS", TitlePos)
            print SpectrumStartPos, TitlePos, SpectrumEndPos
            OutputFile = open(OutputFileNames[SpectrumIndex], "wb")
            OutputFile.write(Stuff[SpectrumStartPos:SpectrumEndPos])
            OutputFile.write("END IONS\n")
            OutputFile.close()
    def FindScanNumberMZXML(self, File):
        ScanNumberText = '<scan num="%s"'%self.ReadScanNumber
        Data = ""
        FilePos = 0
        MEG = 1024*1024
        KEEP = MEG + 20000        
        while 1:
            Block = File.read(MEG)
            if not Block:
                break
            Data += Block
            Pos = Data.find(ScanNumberText)
            if Pos == -1:
                Len = len(Data)
                if Len > KEEP:
                    Data = Data[-KEEP:]
                    FilePos += (Len - KEEP)
                #print FilePos
            else:
                self.ReadBytePosition = FilePos + Pos
                print "Scan %s found at file position %s"%(self.ReadScanNumber, self.ReadBytePosition)
                break
    def FindScanNumberMS2(self, File):
        FilePos = File.tell()
        while (1):
            FileLine = File.readline()
            if not FileLine:
                break # eof
            if FileLine[0]==":":
                ScanNumber = int(FileLine.split(":")[1])
                if ScanNumber == self.ReadScanNumber:
                    self.ReadBytePosition = FilePos
                    break
            if FileLine[0]=="S":
                Bits = FileLine.split()
                ScanNumber = int(Bits[1])
                if ScanNumber == self.ReadScanNumber:
                    self.ReadBytePosition = FilePos
                    break
            FilePos += len(FileLine)
    def ExtractSpectrum(self):
        """
        Main method, after parsing command-line parameters.
        """
        Spectrum = MSSpectrum.SpectrumClass()
        File = open(self.ReadFileName, "rb")
        if self.ReadScanNumber:
            # Figure out the byte position for this scan number.
            Extension = os.path.splitext(self.ReadFileName)[1].lower()
            if Extension == ".ms2":
                self.FindScanNumberMS2(File)
            elif Extension == ".mzxml":
                self.FindScanNumberMZXML(File)
        if self.ReadBytePosition == None:
            print "** Error: Unable to find scan number %s!"%self.ReadScanNumber
            return
        if self.ReadBytePosition:
            File.seek(self.ReadBytePosition)
        # TODO: Handle scan numbers here:
        Spectrum.ReadPeaksFromFile(File, Bob.ReadFileName)
        File.close()
        Spectrum.WritePeaks(Bob.OutputFileName)        
def Test():
    FileNames = [r"E:\ms\AriO16O18\p19-O16-O18-10ug-each-10Da-1D1hr-082306-LTQ1_2_16O.mgf",
                 r"E:\ms\AriO16O18\p19-O16-O18-10ug-each-10Da-1D1hr-082306-LTQ1_1_16O.mgf",
                 r"E:\ms\AriO16O18\p19-O16-O18-10ug-each-10Da-1D1hr-082306-LTQ1_3_16O.mgf",
                 r"E:\ms\AriO16O18\p19-O16-O18-10ug-each-10Da-1D1hr-082306-LTQ1_2_16O.mgf",
                 r"E:\ms\AriO16O18\p19-O16-O18-10ug-each-10Da-1D1hr-082306-LTQ1_3_16O.mgf"
                 ]
    ByteOffsets = [4825334,
                   3712849,
                   667628,
                   3015159,
                   5229283,
                   ]
    OutputFileNames = ["GLSEDTTEETLK.mgf",
                       "DAVTYTEHAK.mgf",
                       "VLLDAPCSGTGVISK.mgf",
                       "EIAQDFK.mgf",
                       "SFLLDLLNATGK.mgf"
                       ]
    Bob = Extractor()
    Bob.FancyExtract(FileNames, ByteOffsets, OutputFileNames)    

UsageInfo = """ExtractSpectrum.py - Write out a single spectrum from a larger file.

Arguments:
-r [FileName]: Filename to parse the spectrum from.
-b [BytePosition]: Byte offset within the file.
-s [ScanNumber]: Byte offset within the file.  [NOT IMPLEMENTED]
-w [FileName]: Name to write the spectrum to; defaults to 'Spectrum.dta'

Options required: -r, and either -b or -s.

Spectrum file formats supported: .mgf, .mzxml, .ms2
"""

if __name__ == "__main__":
    Bob = Extractor()
    (Options, Args) = getopt.getopt(sys.argv[1:], "r:b:w:s:")
    for (Option, Value) in Options:
        if Option == "-w":
            Bob.OutputFileName = Value
        elif Option == "-r":
            Bob.ReadFileName = Value
        elif Option == "-b":
            Bob.ReadBytePosition = int(Value)
        elif Option == "-s":
            Bob.ReadScanNumber = int(Value)
    if not Bob.ReadFileName:
        print UsageInfo
        print "** Please specify a file name, and byte position or scan number."
        sys.exit(-1)
    # Special: Filenames of the form "Foo.mzxml:1234" include a byte offset:
    SplitBits = Bob.ReadFileName.split(":")
    try:
        Bob.ReadBytePosition = int(SplitBits[-1])
        Bob.ReadFileName = string.join(SplitBits[:-1], ":")
    except:
        pass
    if Bob.ReadBytePosition == None and Bob.ReadScanNumber == None:
        print UsageInfo
        print "** Please specify a file name, and byte position or scan number."
        sys.exit(-1)
    Bob.ExtractSpectrum()
