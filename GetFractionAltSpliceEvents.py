"""GetFractionAltSpliceEvents.py
This script takes in a set of annotation for some spectra.  Then we find
the chromosomal coordinates for each peptide.  Marking spliced peptides
and then figuring out which splicing events dileniate between alternate isoforms


"""
UsageInfo = """GetFractionAltSpliceEvents.py
Takes GFFL results and see which peptides help distinguish alternative
splicing

Required Options
 -r [Directory] Directory where the GFFL results files are kept.
 -d [Trie file] TAIR database
 -g [GFF file] the GFF file for TAIR

"""

import sys
import os
import getopt
import traceback
import ResultsParser
import SelectProteins
import GFF
from Utils import *
Initialize()


class CompileClass(ResultsParser.ResultsParser):
    def __init__(self):
        self.ReferenceResults = None
        self.GFFFile = None
        self.OutputPath = None # the training set
        self.ReferenceDatabasePath = [] #possibly multiple
        ResultsParser.ResultsParser.__init__(self)
        #hopefully all programs will annotate scans with the same peptide, but it's no guarentee
        self.AllReferenceAnnotations = {} # (Aminos) => (Score, Line)
        self.AllGenomicLocations = {} # (Aminos, Chr, ChrStart, ChrEnd) = 1 # to keep from printing multiple
                                    # exons for the same location (because multiple transcripts annotate it)
        self.GlobalCount = 0
        self.GlobalSpliceCount = 0
        self.SpectrumCount = 0
        self.PepToFind = "SVSQNSQAPVAVQENGLVK"
        self.LocusToFind = "AT1G17720"
        self.LocusList = {} # Name->LocusClass Object
        self.LociWithTwoDiffIntrons = {} # a dict to keep track of transcripts which differ at
        self.LociWithOneDiffIntrons = {}
        self.SpliceSitesForAltDonorAcceptor = {} #key = (chr, start, stop) value = observed or not
        self.ObservedIntrons = {} # just to keep stuff. key = (chr, start, stop)
        self.ConsensusIntrons = {} #key = (chr, start, stop) value = observed or not
        self.AllIntrons = {}

        
    def Main(self):
        self.ProteinPicker = SelectProteins.ProteinSelector()
        self.ProteinPicker.LoadMultipleDB(self.ReferenceDatabasePath)
        self.ProcessResultsFiles(self.ReferenceResults, self.ParseReferenceFileCallback)
        self.GFFParser = GFF.GFFClass()
        self.Transcripts = self.GFFParser.ParseGFF(self.GFFFile)
        self.OrganizeTranscripts()
        self.GetAltSplicingEvents()
        print "I have %d sites for alt d/a "%len(self.SpliceSitesForAltDonorAcceptor)
        print "I have %d consensus introns"%len(self.ConsensusIntrons)
        print "I have %d ALL introns"%len(self.AllIntrons)
        print "I input %d introns from parsing"%len(self.ObservedIntrons)

    def GetAltSplicingEvents(self):
        """This function takes the loci from the GFF parser, and looks to
        see if we have alternate protein isoforms.  Then it attempts to find the
        locations that differ, the first one.  Then we see if we observe these locations
        
        """
        ObservedEvents = 0
        ObservedConsensusEvents = 0
        ObservedAllEvents = 0
        for Locus in self.LocusList.values():
            self.MakeAllIntronList(Locus)
            if len (Locus.Transcripts) > 1:
                self.CompareIntrons(Locus)
        ## now compare the observed introns to their alt splice ones
        for Intron in self.ObservedIntrons:
            #print "Checking ObservedIntron", Intron
            if self.SpliceSitesForAltDonorAcceptor.has_key(Intron):
                ObservedEvents += 1
            if self.ConsensusIntrons.has_key(Intron):
                ObservedConsensusEvents += 1
            if self.AllIntrons.has_key(Intron):
                ObservedAllEvents += 1
                
        print "Observing %d introns"%ObservedEvents
        print "observed %d consensus introns"%ObservedConsensusEvents
        print "Observed %d of all introns"%ObservedAllEvents


    def MakeAllIntronList(self, Locus):
        AllIntrons = {} # just a counter
        for Transcript in Locus.Transcripts:
            ## 1. Now look at the introns and count occurances
            for Intron in Locus.IntronsByTranscript[Transcript.Name].keys():
                if not AllIntrons.has_key(Intron):
                    AllIntrons[Intron] = 0
                AllIntrons[Intron] += 1
        ## Now cycle through introns and find those that are not consensus
        for (Intron, Count) in AllIntrons.items():
            Key = (Locus.Chromosome, Intron[0], Intron[1])
            self.AllIntrons[Key] = 0

    def CompareIntrons(self, Locus):
        """Looks at the IntronsByTranscript variable and sees which introns are not shared by all transcripts"""
        #if not Locus.Name == self.LocusToFind:
        #    return
        NumTranscripts = len(Locus.Transcripts)
        AllIntrons = {} # just a counter
        UnSharedIntrons = {} # those not shared by all transcripts
        for Transcript in Locus.Transcripts:
            ## 1. Now look at the introns and count occurances
            for Intron in Locus.IntronsByTranscript[Transcript.Name].keys():
                if not AllIntrons.has_key(Intron):
                    AllIntrons[Intron] = 0
                AllIntrons[Intron] += 1
        ## Now cycle through introns and find those that are not consensus
        for (Intron, Count) in AllIntrons.items():
            #Key = (Locus.Chromosome, Intron[0], Intron[1])
            #self.AllIntrons[Key] = 0
            if not Count == NumTranscripts:
                #print Intron, Count
                UnSharedIntrons[Intron] = Count
            else:
                #Start = Intron[0], stop = Intron[1]
                Key = (Locus.Chromosome, Intron[0], Intron[1])
                self.ConsensusIntrons[Key] = 0
                #print Intron
        ## special hack for me
        if len(UnSharedIntrons) < 2:
            ##single splice diff.  looks like
            ## p1       QWETPOIQTQWET*
            ## p2       QWETPOIQTQW^^^POSGIG
            ## OR            
            ## p1       QWETPOIQTQWET*
            ## p2       QWETP^^QTQWET*
            self.LociWithOneDiffIntrons[Locus.Name] = UnSharedIntrons
        if len(UnSharedIntrons) == 2:
            ## could be alt donor/acceptor thingy.
            ## p1       QWETP^^^GFSDI
            ## p2       QWET^^^^GFSDI
            AltDonorAcceptorType = self.IsAltDonorAcceptor(UnSharedIntrons)
            if AltDonorAcceptorType:
                self.LociWithTwoDiffIntrons[Locus.Name] = UnSharedIntrons
                #print "Alt Donor and acceptor %s"%Locus.Name
                #print UnSharedIntrons
                for (Start, Stop) in UnSharedIntrons.keys():
                    Key = (Locus.Chromosome, Start, Stop)
                    self.SpliceSitesForAltDonorAcceptor[Key] = 0 # observed or not
                    #print "Adding Splice Site ", Key
                    #if Start > Stop:
                    #    print "Not strictly numerical ordering"
        if len(UnSharedIntrons) > 2:
            ## let's try to recover some more of the alt donor acceptor things.  Those are
            ## good things to have
            self.RecoverAltDonorAcceptor(UnSharedIntrons, Locus)

    def RecoverAltDonorAcceptor(self, UnSharedIntrons, Locus):
        """This function is called when we have more than two unshared
        introns.  But we would like to check pairs for being alt donor/acceptors
        """

        ## have to figure out how to distinguish better the exon skipping (where two from one transcript
        ## overlap one from another, versus just multiple donor acceptor.
        return 
        ##1. convert to list, so that I can use indexing to check all pairs
        AllIntrons = UnSharedIntrons.keys()
        NumIntrons = len(AllIntrons)
        PairsFound = {} # Intron1 -> intron2, Intron2-> intron1
        # this has will allow us to see if we've used the same intron multiple times
        # which is a sign of exon skipping, not alt donor acceptor
        
        for Index in range(NumIntrons):
            for Jndex in range(Index + 1, NumIntrons):
                Dictionary = {}
                Dictionary[AllIntrons[Index]] = 0
                Dictionary[AllIntrons[Jndex]] = 0
                Found = self.IsAltDonorAcceptor(Dictionary)
                if Found:
                    # check here to see if it really is recovery
                    SNAFU =0
                    if not PairsFound.has_key(AllIntrons[Index]):
                        PairsFound[AllIntrons[Index]] = AllIntrons[Jndex]
                    else:
                        #possible exon skipping.  Print it out
                        SNAFU = 1
                        print "possible exon skipping at  %s"%Locus.Name
                        print AllIntrons[Index]
                        print AllIntrons[Jndex]
                        print PairsFound[AllIntrons[Index]]
                    if not PairsFound.has_key(AllIntrons[Jndex]):
                        PairsFound[AllIntrons[Jndex]] = AllIntrons[Index]
                    else:
                        #possible exon skipping.  Print it out
                        SNAFU =1
                        print "possible exon skipping at  %s"%Locus.Name
                        print AllIntrons[Jndex]
                        print AllIntrons[Index]
                        print PairsFound[AllIntrons[Jndex]]
                    if not SNAFU:
                        print "recovered another pair on %s"%Locus.Name
                    
        

    def IsAltDonorAcceptor(self, UnSharedIntrons):
        """Alternate donor acceptor have the following look
            ## p1       QWETP^^^GFSDI
            ## p2       QWET^^^^GFSDI
        so we want to check these unshared introns to see if that's the case.
        These should share some intronic nucleotides.  So we make a list of the
        region, and check to see if they share.  It's that easy.
        Here I'm guaranteed that there are only two items.
        """
        Range1 = None
        Range2 = None
        for (Start, Stop) in UnSharedIntrons:
            if not Range1:
                Range1 = range(Start, Stop)
            else:
                Range2 = range(Start, Stop)
        #now I go and see whether these overlap
        for Base in Range1:
            if Base in Range2:
                return 1
        return 0

    def OrganizeTranscripts(self):
        "this is to sort and group transcripts into locus and then collect locus information"
        ##1. group transcripts into locus
        self.GroupIntoLocus()
        self.RemoveRedundantTranscripts()
        self.SetLocusIntronsByTranscript()


    def GroupIntoLocus(self):
        "populate th4e self.LocusList"
        for TranscriptID in self.Transcripts.keys():
            Transcript = self.Transcripts[TranscriptID]
            Locus = Transcript.ParentLocus
            if not self.LocusList.has_key(Locus):
                self.LocusList[Locus]= GFF.LocusClass()
                self.LocusList[Locus].Name = Locus
                self.LocusList[Locus].Chromosome = Transcript.Chromosome
            self.LocusList[Locus].Transcripts.append(Transcript)

    def RemoveRedundantTranscripts(self):
        """In the database there are many transcripts that differ only by their
        un-translated regions.  From the proteomics point of view, they are
        redundant, and indistinguishable.  So let's find them and remove them
        """
        for (LocusID, Locus) in self.LocusList.items():
            if len(Locus.Transcripts) < 2:
                continue
            OldLen =  len (Locus.Transcripts)
            Strings = {}
            RedundantTranscripts = []
            Count = 0
            for Transcript in Locus.Transcripts:
                String = ""
                for Exon in Transcript.Exons:
                    String += "%s..%s,"%(Exon.Start, Exon.Stop)
                if Strings.has_key(String):
                    RedundantTranscripts.append(Count)
                else:
                    Strings[String] = Transcript.Name
                Count +=1
            RedundantTranscripts.sort()
            RedundantTranscripts.reverse()  # must reverse so that we pop in descending order, so the indicies remain valid
            for TranscriptIndex in RedundantTranscripts:
                Locus.Transcripts.pop(TranscriptIndex)
            #print Locus.Name, len(Locus.Transcripts), OldLen


    def SetLocusIntronsByTranscript(self):
        """Set the self.IntronsByTranscript variable for all loci, 
        """
        for (LocusID, Locus) in self.LocusList.items():
            for Transcript in Locus.Transcripts:
                Locus.IntronsByTranscript[Transcript.Name] = {}
                ## THis is a MAJOR cludge, but I can't deal with all this reverse strand crap any more
                ## for some reason, I sometimes want exons ordered according to coding order, and
                ## sometimes by numerical order on the DNA.  Here I want option2, which is not the default
                ##option.  So I am going to quickly reverse them, and pretend that this big kludge does
                ## not exist
                if Transcript.Strand == "-":
                    Transcript.Exons.reverse()
                #Transcript.PrintMe()
                for ExonIndex in range(len(Transcript.Exons) - 1):
                    IntronStart = Transcript.Exons[ExonIndex].Stop
                    IntronStop = Transcript.Exons[ExonIndex + 1].Start
                    if not Locus.Introns.has_key((IntronStart, IntronStop)):
                        Locus.Introns[(IntronStart, IntronStop)] = 0 #not yet observed
                    ## add all introns to the "bytranscript" list
                    Locus.IntronsByTranscript[Transcript.Name][(IntronStart, IntronStop)] = 0



    def ParseReferenceFileCallback(self, FilePath):
        """Called by the ResultsParser.ProcessResultsFiles
        Here I parse out lines get the spliced peptides, keeping just the intron boundaries
        ONLY WORKS with MSMS as the input track
        """
        PeptidesWithASplice = 0
        Handle = open(FilePath, "rb")
        for Line in Handle.xreadlines():
            #print "Parsing line %s"%Line
            if Line[:9] == "reference":
                Line = Line.strip()
                Chromosome =Line[10:] # don't get the return carriage
            if not Line.strip():
                continue
            if not Line[:4] == "MSMS":
                continue
            Bits = list(Line.split("\t"))
            try:
                Coordinates = Bits[2]
            except:
                traceback.print_exc()
                continue # SNAFU
            #print Spectrum
            if Coordinates.find(",") == -1:
                # not a splicing peptide
                continue
            ##now here we deal with the aweful gffl.  the splice boundaries will be the
            ## numbers near the comma, but may be in the wrong order.  we want strictly numerical
            ##ordering for our purposes
            ##Observed Sequences are written in the (aweful) GFFL differently for items belonging to
            ##different strands.  For Example, on the reverse strand we write
            ##525..521,424..388
            ##See how the numbers are totally backwards.  both within a unit (525..521) and between?
            PeptidesWithASplice += 1
            Exons = Coordinates.split(",")
            # here it's possible to have multiple introns, so let's get them all
            for Index in range(len(Exons)-1):
                (Start1, Stop1) = Exons[Index].split("..")
                (Start2, Stop2) = Exons[Index + 1].split("..")
                #stupid formatting.  now convert to Ints
                Boundary = [int(Stop1), int(Start2)]
                Boundary.sort() # always numberical lowest, highest
                if Boundary[0] > Boundary[1]:
                    print "DANG IT.  Need numerical ordering"
                else:
                    self.ObservedIntrons[(Chromosome, Boundary[0], Boundary[1])] = 1

        Handle.close()
        print "I observe %s peptides with a splice"%PeptidesWithASplice


    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:g:d:t:s:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-r":
                # -r results file(s)
                if not os.path.exists(Value):
                    print "** Error: couldn't find results file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.ReferenceResults = Value
            if Option == "-g":
                # -r results file(s)
                if not os.path.exists(Value):
                    print "** Error: couldn't find results file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.GFFFile = Value
            if Option == "-d":
                if not os.path.exists(Value):
                    print "** Error: couldn't find database file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.ReferenceDatabasePath.append( Value)
            if Option == "-t":
                self.TrackName = Value
            if Option == "-s":
                self.SharedUniqueFlag = int(Value)
        if not OptionsSeen.has_key("-r")  or not OptionsSeen.has_key("-g") or not OptionsSeen.has_key("-d"):
            print UsageInfo
            sys.exit(1)


if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(psyco not found - running in non-optimized mode)"
    Gather = CompileClass()
    Gather.ParseCommandLine(sys.argv[1:])
    Gather.Main()