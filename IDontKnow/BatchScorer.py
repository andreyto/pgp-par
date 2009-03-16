"""
BatchScorer - simple PyInspect interface for Ari
"""
import getopt
import os
import sys
import PyInspect
import Utils

UsageInfo = """
BatchScorer.py - score all spectra in an .mgf file
Arguments:
-s [FILENAME]: .mgf file
-p [FILENAME]: Annotation file.  File should be tab-delimited, with
   one line per scan.  First column should be a scan number (numbering
   from zero), second column a peptide annotation, third column a charge.
   Blank lines, or lines starting with "#", are ignored.  
-w [FILENAME]: Output filename; defaults to "SpectrumScores.txt"
"""

class BatchScorer:
    def __init__(self):
        self.SpectrumFileName = None
        self.AnnotationFileName = None
        self.OutputFileName = "Scores.txt"
    def Main(self):
        # Check args:
        if not self.SpectrumFileName or not self.AnnotationFileName:
            print UsageInfo
            return
        self.ScanNumberFilePositions = self.EnumerateScans()
        File = open(self.SpectrumFileName, "rb")
        AnnotationFile = open(self.AnnotationFileName, "rb")
        LineNumber = 0
        ScoreCount = 0
        self.OutputFile = open(self.OutputFileName, "wb")
        for FileLine in AnnotationFile.xreadlines():
            LineNumber += 1
            if FileLine[0] == "#":
                continue
            Bits = FileLine.strip().split("\t")
            if len(Bits) < 3:
                continue
            try:
                ScanNumber = int(Bits[0])
                Peptide = Bits[1]
                Charge = int(Bits[2])
            except:
                print "* Syntax error on line %s of %s"%(LineNumber, self.AnnotationFileName)
                traceback.print_exc()
                continue
            ScanFilePos = self.ScanNumberFilePositions.get(ScanNumber, None)
            if ScanFilePos == None:
                print "* Error on line %s of %s: Scan number '%s' doesn't exist!"%(LineNumber, self.AnnotationFileName, ScanNumber)
                continue
            Spectrum = PyInspect.Spectrum(self.SpectrumFileName, ScanFilePos)
            Score = Spectrum.ScorePeptide(Peptide)
            Str = "%s\t%s\t%s\t%s\t"%(ScanNumber, Peptide, Charge, Score)
            self.OutputFile.write(Str + "\n")
            ScoreCount += 1
        print "Successfully scored %s annotations on %s file lines"%(ScoreCount, LineNumber)
    def EnumerateScans(self):
        File = open(self.SpectrumFileName, "rb")
        ScanText = File.read()
        File.close()
        NextScanNumber = 0
        LastPos = -1
        ScanPositions = {}
        while 1:
            Pos = ScanText.find("BEGIN IONS", LastPos + 1)
            print LastPos, Pos
            if Pos == -1:
                break
            ScanPositions[NextScanNumber] = Pos
            LastPos = Pos
            NextScanNumber += 1
        print "Found %s scans in file %s..."%(len(ScanPositions.keys()), self.SpectrumFileName)
        return ScanPositions
    def ParseCommandLine(self):
        (Options, Args) = getopt.getopt(sys.argv[1:], "s:p:w:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-s":
                if not os.path.exists(Value):
                    print "** Error: couldn't find spectrum file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.SpectrumFileName = Value
            elif Option == "-p":
                if not os.path.exists(Value):
                    print "** Error: couldn't find annotation file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.AnnotationFileName = Value
            elif Option == "-w":
                self.OutputFileName = Value
                
if __name__ == "__main__":
    Scorer = BatchScorer()
    Scorer.ParseCommandLine()
    Scorer.Main()