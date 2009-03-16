"""
Given ordinary search results and search results from a splice-db, look for peptides that were found in the
splice-db but NOT in the standard search.  Also, look for overlaps between these new peptides and the old
search results.
"""
import os
import sys
import traceback


File = open("Database\\DCLens.fasta", "r")
OriginalStuff = File.read().replace("I", "L").replace("Q","K").replace("\n","").replace(" ","")

class SpliceAnalyzer:
    def __init__(self):
        pass
    def Analyze(self, OldResultsFile, SpliceResultsFile):
        PeptideDict =  {}
        OriginalHits = {}
        LineCount = 0
        OldSpectrum = None
        for FileLine in OldResultsFile.xreadlines():
            Bits = FileLine.split("\t")
            Spectrum = Bits[0]
            MQScore = float(Bits[5])
            if Spectrum!=OldSpectrum:
                OldSpectrum = Spectrum
                OriginalHits[Spectrum] = (MQScore, Bits[3], Bits[9])
            if MQScore < -.55:
                continue
            Pep = Bits[1].replace("I","L").replace("Q","K")
            PeptideDict[Pep] = PeptideDict.get(Pep, 0) + 1
            LineCount += 1
        print "Original results file: Read %d peptides from %d lines."%(len(PeptideDict.keys()), LineCount)
        sys.stdout.flush()
        # Now examine the new stuff:
        SplicePepCount = 0
        FoundPepCount = 0
        PeptidesToGo = PeptideDict.copy()
        LineNumber = 0
        for FileLine in SpliceResultsFile.xreadlines():
            LineNumber += 1
            if LineNumber%100 == 0:
                sys.stdout.flush()
            SplicePepCount += 1
            Bits = FileLine.split("\t")
            Peptide = Bits[1].replace("I", "L").replace("Q","K")
            if PeptideDict.get(Peptide, None):
                FoundPepCount += 1
                try:
                    del PeptidesToGo[Bits[1]]
                except:
                    pass
                continue
            MQScore=  float(Bits[5])
            if MQScore < -.45:
                continue
            # Ok, it's a new peptide.  Check whether it overlaps:
            Len = len(Peptide)
            Half = Len/2
            LeftAA = Peptide[:Half]
            RightAA = Peptide[Half:]
            PossibleOverlaps = {}
            OverlapCount = 0
            Containment = 0
            if OriginalStuff.find(Peptide) == -1:
                for OldPeptide in PeptideDict.keys():
                    OverLen = self.IsOverlap(Peptide, LeftAA, RightAA, OldPeptide)
                    if OverLen < 0:
                        # Containment!
                        Containment = 1
                        break
                    if OverLen:
                        PossibleOverlaps[OldPeptide] = OverLen
                        OverlapCount = 1
                        if OverlapCount > 5:
                            break
                if not Containment:
                    print Peptide, OverlapCount, PossibleOverlaps
                    Str = "%s\t%s\t%s\t%s\t%s\t"%(Bits[0], Bits[9], Bits[3],Bits[4], Bits[5])
                    Str += "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t"%(Bits[6], Bits[8], str(OriginalHits.get(Bits[0], None)), Bits[12], Bits[13], Bits[14], Bits[15], PossibleOverlaps)
                    print Str
        print "Peptide known in %d of %d cases."%(SplicePepCount, FoundPepCount)
        print "Old peptides not seen: %d"%(len(PeptidesToGo.keys()))
    def IsOverlap(self, Peptide, LeftAA, RightAA, OldPeptide):
        if OldPeptide.find(Peptide)!=-1:
            return -1        
        LeftPos = OldPeptide.find(LeftAA)
        if LeftPos!=-1:
            SideStuff = OldPeptide[LeftPos+len(LeftAA):]
            #print "L:", LeftAA, RightAA, OldPeptide
            if SideStuff.find(RightAA) == 0 or RightAA.find(SideStuff)==0:
                #print "YES."
                OverLen = len(LeftAA) + min(len(RightAA), len(SideStuff))
                return OverLen
        RightPos = OldPeptide.find(RightAA)
        if RightPos!=-1:
            SideStuff = OldPeptide[:RightPos]
            #print "R:", LeftAA, RightAA, OldPeptide
            if LeftAA.find(SideStuff)!=-1 or SideStuff.find(LeftAA)!=-1:
                #print "YES."
                OverLen = len(RightAA) + min(len(LeftAA), len(SideStuff))
                return OverLen
        return 0
if __name__ == "__main__":
    # Arguments:
    # PlainResultsFileName and SplicedResultsFileName
    OldResultsFile = open(sys.argv[1], "r")
    SpliceResultsFile = open(sys.argv[2], "r")
    Analyzer = SpliceAnalyzer()
    Analyzer.Analyze(OldResultsFile, SpliceResultsFile)
    