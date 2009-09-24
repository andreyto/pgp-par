"""
Look for search resuls on the exon graph which are not better explained by
a search of the IPI database.
"""
import os
import sys
import string

class ExonGraphFilter:
    def __init__(self):
        self.IPIResults = {}
        self.MSAResults = {}
    def HandleBatch(self, ResultsX, ResultsIPI, FiltrationOutput):
        for FileName in os.listdir(ResultsX):
            (Stub, Extension) = os.path.splitext(FileName)
            if Extension != ".txt":
                continue
            InputFilePath = os.path.join(ResultsX, FileName)
            OutputFilePath = os.path.join(FiltrationOutput, FileName)
            #if os.path.exists(OutputFilePath):
            #    continue # already filtered it!
            IPIResultsPath = os.path.join(ResultsIPI, "%s.filtered"%Stub)
            if not os.path.exists(IPIResultsPath):
                IPIResultsPath = os.path.join(ResultsIPI, "%s.txt"%Stub)
            try:
                self.ParseIPIResults(IPIResultsPath)
            except:
                continue
            self.FilterXGResults(InputFilePath, OutputFilePath)
    def FilterXGResults(self, FilePath, OutputFilePath):
        print "Filter from %s to %s..."%(FilePath, OutputFilePath)
        InFile = open(FilePath, "rb")
        OutFile = open(OutputFilePath, "wb")
        SkippingFlag = 0
        OldSpectrum = None
        SpectrumCount = 0
        FilterCount = 0
        OutputLineCount = 0
        LineCount = 0
        OldPeptideInfo = ""
        for FileLine in InFile.xreadlines():
            LineCount += 1
            Bits = list(FileLine.split("\t"))
            try:
                MQScore = float(Bits[5])
                PValue = float(Bits[10])
                ScanNumber = int(Bits[1])
            except:
                continue
            Spectrum = (Bits[0], Bits[1])
            if Spectrum == OldSpectrum:
                if SkippingFlag:
                    continue
            else:
                SkippingFlag = 0 # innocent by default!
                SpectrumCount += 1
                OldPeptideInfo = ""
                OldSpectrum = Spectrum
                OldTuple = self.IPIResults.get(ScanNumber, None)
                if OldTuple != None:
                    (OldPValue, OldMQScore, OldBits) = OldTuple
                    OldPeptideInfo = "%s %s"%(OldBits[2], OldBits[5])
                    if OldPValue < PValue * 0.75 and (OldMQScore > MQScore + 0.5):
                        # Throw out search hits on this spectrum, because there's an IPI search hit
                        # that's just as good.
                        SkippingFlag = 1
                        FilterCount += 1
            if PValue < 0.12 and not SkippingFlag:
                OutputLineCount += 1
                Bits[14] = OldPeptideInfo
                OutFile.write(string.join(Bits, "\t"))
        InFile.close()
        OutFile.close()
        print "-->Parsed %s spectra from %s"%(SpectrumCount, FilePath)
        print "Filtered %s spectra out"%FilterCount
        print "Wrote %s/%s lines"%(OutputLineCount, LineCount)
    def ParseIPIResults(self, FilePath):
        """
        Parse IPI search output, and remember the results in self.IPIResults
        """
        self.IPIResults = {}
        try:
            File = open(FilePath, "rb")
        except:
            print "* IPI results not seen at '%s'"%FilePath
            raise
        OldSpectrum = None
        for FileLine in File.xreadlines():
            Bits = FileLine.split("\t")
            try:
                MQScore = float(Bits[5])
                PValue = float(Bits[10])
                ScanNumber = int(Bits[1])
            except:
                continue
            Spectrum = (Bits[0], Bits[1])
            if Spectrum == OldSpectrum:
                continue
            OldSpectrum = Spectrum
            self.IPIResults[ScanNumber] = (PValue, MQScore, Bits)
        File.close()

if __name__ == "__main__":
    import psyco
    psyco.full()
    Filter = ExonGraphFilter()
    Filter.HandleBatch("e:\\ms\\briggs\\resultsXFixed","e:\\ms\\briggs\\resultsIPIFixed", "e:\\ms\\briggs\\resultsxx")
    Filter.HandleBatch("e:\\ms\\PeptideAtlas\\resultsXFixed","e:\\ms\\PeptideAtlas\\resultsIPIFixed", "e:\\ms\\PeptideAtlas\\resultsxx")
