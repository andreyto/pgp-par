"""
The splice-tolerant search of lens spectra recovered most of the annotations that were found
in a standard database search.  However, SOME high-scoring, unmodified peptide annotations
may have been lost!  That may be due to problems with tagging, or problems with false matches
'crowding out' the true matches and blocking them from being rescored.  But, it also may be
due to the peptide simply not being present in the splice-tolerant database.

For each canonical protein sequence, we'd like to verify that the
sequence - or at least, most of it - can be matched within our splice-tolerant database.
The procedure is as follows:
- Open and read the gene-database
- For each exon: Call a recursive matching-function to match, starting with that
 exon and the start of the sequence.  As soon as you cover most of the protein with such
 a call, stop.
"""
import os
import psyco
import sys
import struct

class SplicedExon:
    def __init__(self):
        self.Sequence = ""
        self.BackLinkExon = []
        self.BackLinkAA = []
        self.BackLinkPower = []
        self.BackLinkIndex = []
        self.ForwardLinkExon = []
        self.ForwardLinkAA = []
        self.ForwardLinkPower = []
        self.ForwardLinkIndex = []

class SplicedGene:
    """
    One record in the splice-tolerant database.  A gene is a graph where the nodes
    are exons, and each exon has forward and backward links to other exons.
    """
    def __init__(self):
        self.Exons = []
    def DebugPrint(self):
        print ">> Debug print gene of %s exons: '%s'"%(len(self.Exons), self.Name)
        for Exon in self.Exons:
            print "Exon %s: length %s (%s-%s)"%(Exon.Index, len(Exon.Sequence), Exon.Start, Exon.End)
            print Exon.Sequence
            for Index in range(len(Exon.BackLinkExon)):
                print " <<< back %s by '%s' to exon %s (...%s)"%(Exon.BackLinkPower[Index], Exon.BackLinkAA[Index], Exon.BackLinkExon[Index].Index,
                                                              Exon.BackLinkExon[Index].Sequence[-10:])
            for Index in range(len(Exon.ForwardLinkExon)):
                print " >>> forw %s by '%s' to exon %s (%s...)"%(Exon.ForwardLinkPower[Index], Exon.ForwardLinkAA[Index], Exon.ForwardLinkExon[Index].Index,
                                                              Exon.ForwardLinkExon[Index].Sequence[:10])
    def LoadFromFile(self, FileName):
        File = open(FileName, "rb")
        self.Name = File.read(256)
        File.read(256) # dummy
        self.Chromosome = struct.unpack("<i", File.read(4))[0]
        self.ExonCount = struct.unpack("<i", File.read(4))[0]
        for ExonIndex in range(self.ExonCount):
            print "Read exon %s..."%ExonIndex
            Exon = SplicedExon()
            Exon.Index = len(self.Exons)
            self.Exons.append(Exon)
            Exon.Start = struct.unpack("<i", File.read(4))[0]
            Exon.End = struct.unpack("<i", File.read(4))[0]
            Exon.Length = struct.unpack("<i", File.read(4))[0]
            Exon.Occurrences = struct.unpack("<i", File.read(4))[0]
            if Exon.Length:
                Exon.Sequence = File.read(Exon.Length)
            Exon.Prefix = File.read(2)
            Exon.Suffix = File.read(2)
            Exon.BackLinkCount = struct.unpack("<i", File.read(4))[0]
            Exon.ForwardLinkCount = struct.unpack("<i", File.read(4))[0]
            for LinkIndex in range(Exon.BackLinkCount):
                Exon.BackLinkIndex.append(struct.unpack("<i", File.read(4))[0])
                Exon.BackLinkPower.append(struct.unpack("<i", File.read(4))[0])
                Exon.BackLinkAA.append(File.read(1))
        self.FixExonLinks()
    def FixExonLinks(self):
        print "Fixup exon links..."
        for Exon in self.Exons:
            for Index in range(len(Exon.BackLinkIndex)):
                ExIndex = Exon.BackLinkIndex[Index]
                LinkedExon = self.Exons[ExIndex]
                Exon.BackLinkExon.append(LinkedExon)
                LinkedExon.ForwardLinkPower.append(Exon.BackLinkPower[Index])
                LinkedExon.ForwardLinkAA.append(Exon.BackLinkAA[Index])
                LinkedExon.ForwardLinkExon.append(Exon)
                LinkedExon.ForwardLinkIndex.append(Exon.Index)
        print "Exon links verified."
    def PrintBushiestBranches(self):
        """
        Find nodes in our graph of high in-degree or high out-degree.  These are positions
        where alternative splicing is likely.
        """
        
def GetOriginalProtienSequence(ProteinName):
    "Read a (standard) protein sequence from the lens FASTA-file."
    FASTAFile = open("Database\\Lens.fasta", "r")
    GettingSequence = 0
    Sequence = ""
    for FileLine in FASTAFile.xreadlines():
        if FileLine[0] == ">":
            if GettingSequence:
                break
            if FileLine.find(ProteinName)==1:
                GettingSequence = 1
                continue
        elif GettingSequence:
            Sequence += FileLine.strip()
    return Sequence

class ProteinMatcher:
    def MatchProtein(self, ProteinName, ProteinSequence):
        """
        Open the splicedb for this protein, and attempt to match the full true
        sequence of the protein by walking through the graph.
        """
        GeneFileName = os.path.join("LensSpliceDB", "%s.dat"%ProteinName)
        Gene = SplicedGene()
        Gene.LoadFromFile(GeneFileName)
        Gene.DebugPrint()
        if len(Gene.Exons) == 0:
            print "No genes -> no match :("
            return
        # Ok, now try to find this protein sequence within our exons:
        ProteinLen = len(ProteinSequence)
        if ProteinLen < 10:
            print "!?? protein sequence looks too short!", ProteinSequence
            return
        print "Protein sequence:"
        print ProteinSequence
        BestMatchLength = -1
        for ExonIndex in range(len(Gene.Exons)):
            print "X:", ExonIndex, Gene.Exons[ExonIndex].Sequence
            (MatchedLength, Trail) = self.MatchSequence(Gene.Exons[ExonIndex], ProteinSequence, 0)
            if (MatchedLength > BestMatchLength):
                BestMatchLength = MatchedLength
                BestTrail = Trail
                BestExonIndex = ExonIndex
        print "**Matched %s of %s bases, starting in exon %s"%(BestMatchLength, ProteinLen, BestExonIndex)
        print "Exons used:", len(BestTrail)
        BestTrail.reverse()
        for Exon in BestTrail:
            print Exon.Index, Exon.Length, Exon.Sequence
        return BestMatchLength
    def MatchSequence(self, Exon, Sequence, IncomingEdgeUsed):
        MaxMatch = 0
        MaxMatchTrail = [Exon,]
        if not Sequence:
            return (0, [])
        if Exon.Sequence:
            print "MS: X %s len %s starts %s (sequence len %s starts %s)"%(Exon.Index, len(Exon.Sequence), Exon.Sequence[:5],
                len(Sequence), Sequence[:5])
        else:
            print "MS: X %s len %s (sequence len %s starts %s)"%(Exon.Index, len(Exon.Sequence),
                len(Sequence), Sequence[:5])
            
        if IncomingEdgeUsed:
            MatchBreakPoint = None
            # match one amino acid at a time:
            for Index in range(min(len(Exon.Sequence), len(Sequence))):
                if Sequence[Index] != Exon.Sequence[Index]:
                    MatchBreakPoint = Index
                    break
            print "MatchBreakPoint:",MatchBreakPoint
            if MatchBreakPoint == None:
                # We matched the full length of the possible alignment:
                if len(Exon.Sequence) >= len(Sequence):
                    # We're done!  We matched everything.
                    print "Match complete!"
                    return (len(Sequence), [Exon,])
                # We finished this exon, but we'll need to use more exons in order
                # to cover everything:
                print "Exon %s covered from start to finish (len %s)."%(Exon.Index, len(Exon.Sequence))
                NextAAIndex = len(Exon.Sequence)
                MaxMatch = NextAAIndex
                MaxMatchTrail = [Exon,]                
                for EdgeIndex in range(Exon.ForwardLinkCount):
                    if Exon.ForwardLinkAA[EdgeIndex] == '\0':
                        (MatchLen, MatchTrail) = self.MatchSequence(Exon.ForwardLinkExon[EdgeIndex], Sequence[NextAAIndex:], 1)
                        MatchLen += NextAAIndex
                        if MatchLen > MaxMatch:
                            MaxMatch = MatchLen
                            MaxMatchTrail = MatchTrail
                            MaxMatchTrail.append(Exon)
                        if MaxMatch == len(Sequence):
                            return (MaxMatch, MaxMatchTrail)
                    elif Exon.ForwardLinkAA[EdgeIndex] == Sequence[NextAAIndex]:
                        (MatchLen, MatchTrail) = self.MatchSequence(Exon.ForwardLinkExon[EdgeIndex], Sequence[NextAAIndex+1:], 1)
                        MatchLen += NextAAIndex
                        if MatchLen > MaxMatch:
                            MaxMatch = MatchLen + 1
                            MaxMatchTrail = MatchTrail
                            MaxMatchTrail.append(Exon)
                        ##MaxMatch = max(MaxMatch, NextAAIndex + self.MatchSequence(Exon, Sequence[NextAAIndex+1:], 1))
                        if MaxMatch == len(Sequence):
                            return (MaxMatch, MaxMatchTrail)
                return (MaxMatch, MaxMatchTrail)
            else:
                # We matched part of the exon, but not all:
                return (MatchBreakPoint, [Exon,])
        else:
            # We haven't used the incoming edge yet, so we can start the match anywhere we want.
            # First, try using incoming edges:
            for EdgeIndex in range(Exon.BackLinkCount):
                if Exon.BackLinkAA[EdgeIndex] == Sequence[0]:
                    (MatchLen, MatchTrail) = self.MatchSequence(Exon, Sequence[1:], 1)
                    if (MatchLen > MaxMatch):
                        MaxMatch = MatchLen
                        MaxMatchTrail = MatchTrail
                    if MaxMatch == len(Sequence):
                        return (MaxMatch, MaxMatchTrail)
            for StartIndex in range(len(Exon.Sequence)):
                MaxLen = min(len(Sequence), len(Exon.Sequence) - StartIndex)
                MatchBreakPoint = None
                for Index in range(MaxLen):
                    if Exon.Sequence[Index + StartIndex] != Sequence[Index]:
                        MatchBreakPoint = Index
                        break
                if MatchBreakPoint > 0:
                    print "SI %s breaks at %s"%(StartIndex, MatchBreakPoint)
                if MatchBreakPoint == None:
                    if MaxLen == len(Sequence):
                        # We're done!  We matched everything.
                        print "MATCH COMPLETE."
                        return (len(Sequence), [Exon,])
                    NextAAIndex = MaxLen
                    print "Exon %s covered from %s to end (len %s) %s."%(Exon.Index, StartIndex, len(Exon.Sequence), MaxLen)
                    MaxMatchLen = MaxLen
                    MaxMatchTrail = [Exon,]
                    for EdgeIndex in range(Exon.ForwardLinkCount):
                        print "Edge %s: Forward by '%s' to exon %s"%(EdgeIndex, Exon.ForwardLinkAA[EdgeIndex], Exon.ForwardLinkExon[EdgeIndex].Index)
                        if Exon.ForwardLinkAA[EdgeIndex] == '\0':
                            (MatchLen, MatchTrail) = self.MatchSequence(Exon.ForwardLinkExon[EdgeIndex], Sequence[NextAAIndex:], 1)
                            MatchLen += MaxLen
                            if (MatchLen > MaxMatch):
                                MaxMatch = MatchLen
                                MaxMatchTrail = MatchTrail
                                MaxMatchTrail.append(Exon)
                            if MaxMatch == len(Sequence):
                                return (MaxMatch, MaxMatchTrail)
                        elif Exon.ForwardLinkAA[EdgeIndex] == Sequence[NextAAIndex]:
                            (MatchLen, MatchTrail) = self.MatchSequence(Exon.ForwardLinkExon[EdgeIndex], Sequence[NextAAIndex+1:], 1)
                            MatchLen += MaxLen+1
                            if (MatchLen > MaxMatch):
                                MaxMatch = MatchLen
                                MaxMatchTrail = MatchTrail
                                MaxMatchTrail.append(Exon)
                            if MaxMatch == len(Sequence):
                                return (MaxMatch, MaxMatchTrail)
                else:
                    if MatchBreakPoint > MaxMatch:
                        MaxMatch = MatchBreakPoint
                        MaxMatchTrail = [Exon]
            print "End up with:", MaxMatch, MaxMatchTrail
            return (MaxMatch, MaxMatchTrail)
##        MaxMatch = 0
##        if not IncomingEdgeUsed:
##            for EdgeIndex in range(len(Exon.BackLinkAA)):
##                if Exon.BackLinkAA[EdgeIndex] == Sequence[0]:
##                    MaxMatch = max(self.MatchSequence(Exon, Sequence[1:], 1))
##                    break
        
                                   
                
if __name__ == "__main__":
    #os.chdir(r"C:\ftproot\swt\Inspect")
    #ProteinName = sys.argv[1].upper()
    psyco.full()
    if len(sys.argv)>1:
        ProteinName = sys.argv[1]
    else:
        ProteinName = "CRYAB_HUMAN"
    
    ProteinSequence = GetOriginalProtienSequence(ProteinName)
    #File = open("Q8N519_HUMAN.txt", "r")
    ##File = open("Isoform2.txt", "r")
    ##ProteinSequence = File.read().replace(" ","").replace("\r","").replace("\n","")
    #ProteinSequence = ProteinSequence[15:] 
    Matcher = ProteinMatcher()
    Length = Matcher.MatchProtein(ProteinName, ProteinSequence)
##    if Length < len(ProteinSequence) / 2:
##        print "Try dropping a prefix..."
##        Length = Matcher.MatchProtein(ProteinName, ProteinSequence[2:])