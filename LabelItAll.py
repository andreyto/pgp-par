"""LabelItAll.py

This script is a glorified wrapper for Label.py.  If you have a file of many
annotations, like from Inspect or MakeTrainingSet.py, this will make a pretty
picture of all of them.
"""

UsageInfo = """LabelItAll.py
This script is a wrapper for the verbose output of Label.py, making
pictures and text files about every annotation in the input file.

Required Options:
 -r [FileName] File of Inspect formatted annotations
 -w [Directory] Name of the directory where all the pictures and
     text from Label will be put.
 -m [Directory] Directory containing spectra files

Additional Options:
 -A Puts the Ascore of phosphopeptides into the verbose output.
"""

import os
import sys
import getopt
import ResultsParser
import string
import Label

class DetectiveClass(ResultsParser.ResultsParser):
    def __init__(self):
        self.InputFilePath = None
        self.InspectResultsDir = None
        self.LabeledAnnotationsDir = "LabelSpewage" # for Labeled output
        self.MZXMLDir = None
        self.OldInspectResults = {} #(file, scan) => (MQScore, Annotation)  #file name only, not! path
        self.AScore = 0
        ResultsParser.ResultsParser.__init__(self)
        
    def Main(self):
        if self.InspectResultsDir:
            self.ProcessResultsFiles(self.InspectResultsDir, self.ParseInspectResults)
        MakeDirectory(self.LabeledAnnotationsDir)
        self.LabelMe()

    def LabelMe(self):
        Handle = open(self.InputFilePath, "rb")
        Dymo = Label.LabelClass()
        Count = 0
        GoodScoreCount = 0
        WrongChargeCount = 0
        ScoredWorseCount = 0
        for Line in Handle.xreadlines():
            if Line[0] == "#":
                continue
            if not Line.strip():
                continue
            Bits = list(Line.split("\t"))
            Charge = int (Bits[self.Columns.Charge])
            Count +=1
            Annotation = Bits[self.Columns.Annotation]
            FileName = Bits[self.Columns.SpectrumFile]
            Scan = int(Bits[self.Columns.ScanNumber])
            ByteOffset = int(Bits[self.Columns.FileOffset])
            (Path,File) = os.path.split(FileName)
            FileName = os.path.join(self.MZXMLDir, File)
            VerboseFileName = "%s.%s.%s.verbose.txt"%(File, Scan, Annotation[2:-2])
            ImageFileName = "%s.%s.%s.png"%(File, Scan, Annotation[2:-2])
            VerboseFilePath = os.path.join(self.LabeledAnnotationsDir, VerboseFileName)
            ImageFilePath = os.path.join(self.LabeledAnnotationsDir, ImageFileName)
            ## as we've got a single Dymo object, we must be passing in full args list
            ## because it might have some things prestored, if we don't retell it (like charge is a killer)
            Args = " -r %s -b %d -a %s -v %s -w %s -p -c %d"%(FileName, ByteOffset, Annotation, VerboseFilePath, ImageFilePath, Charge)
            if self.AScore:
                Args += " -A"
            ArgsList = Args.split()
            #print "Parsing Results for %s, scan %s, charge %s"%(FileName, Scan, Charge)
            #print "Args: %s"%Args
            Dymo.ParseCommandLine(ArgsList)
            Dymo.Main()
            ## get Inspect top score, get this top score
            if Dymo.MQScore < 0: # don't bother with these
                continue
            GoodScoreCount += 1
            if self.OldInspectResults.has_key((File, Scan)):
                (OldScore, OldAnnotation, OldCharge) = self.OldInspectResults[(File, Scan)]
                if OldAnnotation == Annotation:
                    print "RawInspect got that answer too, but it did not pass the Pvalue filter"
                else:
                    print "RawInspect Score %f, %s, charge %d"%(OldScore, OldAnnotation, OldCharge)
                    print "The New labeling has a score difference of %f"%(Dymo.MQScore - OldScore)
                    if (Dymo.MQScore - OldScore) > 0:
                        ScoredWorseCount += 1
                        if not OldCharge == self.Charge:
                            WrongChargeCount += 1
            print #to break up the strings

        Handle.close()
        print "processed %d results"%Count
        if self.InspectResultsDir:
            print "Parsed %d, %d had a passable Inspect score."%(Count, GoodScoreCount)
            print "NewInspect had a worse score on %d, got %d charge wrong"%(ScoredWorseCount, WrongChargeCount)

    def ParseInspectResults(self, FileName):
        Handle = open (FileName, "rb")
        OldSpectrum = None
        for Line in Handle.xreadlines():
            if Line[0] == "#":
                continue
            if not Line.strip():
                continue
            Bits = list(Line.split("\t"))
            Charge = int (Bits[self.Columns.Charge])
            Annotation = Bits[self.Columns.Annotation]
            SpectrumPath = Bits[self.Columns.SpectrumFile]
            SpectrumFile = os.path.split(SpectrumPath)[1]
            Scan = int(Bits[self.Columns.ScanNumber])
            MQScore = float(Bits[self.Columns.MQScore])
            Spectrum = (SpectrumFile, Scan)
            if Spectrum == OldSpectrum:
                continue #only tabulate the top result
            OldSpectrum = Spectrum
            self.OldInspectResults[Spectrum] = (MQScore, Annotation, Charge)
        Handle.close()
        

    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:w:m:i:A")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-r":
                # -r results file(s)
                if not os.path.exists(Value):
                    print "** Error: couldn't find results file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.InputFilePath = Value
            if Option == "-w":
                self.LabeledAnnotationsDir = Value
            if Option == "-m":
                self.MZXMLDir = Value
            if Option == "-i":
                self.InspectResultsDir = Value
            if Option == "-A":
                self.AScore = 1
        if not OptionsSeen.has_key("-r") or not OptionsSeen.has_key("-m") or not OptionsSeen.has_key("-w"):
            print UsageInfo
            sys.exit(1)


def MakeDirectory(Dir):
    if os.path.exists(Dir):
        return 
    try:
        os.makedirs(Dir)
    except:
        raise
    


if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(psyco not found - running in non-optimized mode)"
    MacGyver = DetectiveClass()
    MacGyver.ParseCommandLine(sys.argv[1:])
    MacGyver.Main()