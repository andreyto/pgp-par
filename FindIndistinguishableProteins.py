"""FindIndistinguishable.py
This script takes in a set of annotation for some spectra.  Then we find sets of
proteins which have the exact same set of peptides
"""
UsageInfo = """FindIndistinguishable.py
This script takes in a set of annotation for some spectra.  Then we find sets
of proteins which have the exact same set of peptides.  We compare this list
of MS/MS identical proteins to a list of genomically identically proteins.
Required Options
 -r [Directory] Directory where the reference results files are kept.
 -d [Trie file] TAIR database
 -w [FileName]  Output file


"""

import sys
import os
import getopt
import traceback
import ResultsParser
import SelectProteins
from Utils import *
Initialize()


class CompileClass(ResultsParser.ResultsParser):
    def __init__(self):
        self.ReferenceResults = None
        self.GFFFile = None
        self.OutputPath = "RenameYourOutput.txt"
        self.ReferenceDatabasePath = [] #possibly multiple
        ResultsParser.ResultsParser.__init__(self)
        #hopefully all programs will annotate scans with the same peptide, but it's no guarentee
        self.AllReferenceAnnotations = {} # aminos -> count
        self.SpectrumCount = 0
        self.MSProteinPeptideList = {} # ProteinID -> [Pep1, Pep2, ..]
        self.MSMSIdenticalProteins = {} #Group->[ID2, ID2, ..]
        self.GenomicIdentityListFile = None
        self.GenomicIdentityList = {} #p1 -> [p1, p2, p3, ..]
        self.ListsAreTheSame ={}
        self.ListsAreTheSame["Same"] = 0
        self.ListsAreTheSame["GenomicBigger"] = 0
        self.ListsAreTheSame["MSMSBigger"] = 0
        self.SameLocusCounter =0
        
    def Main(self):
        self.ProteinPicker = SelectProteins.ProteinSelector()
        self.ProteinPicker.LoadMultipleDB(self.ReferenceDatabasePath)
        self.GenerateGenomeIdentityList()
        self.ProcessResultsFiles(self.ReferenceResults, self.ParseReferenceFileCallback)
        self.PutPeptidesInProteins()
        self.RemoveOneHitWonders()
        self.FindMSMSIdenticalProteins()
        self.WriteResults()

    def GenerateGenomeIdentityList(self):
        """simple double loop to check for identities"""
        DatabaseLen = len(self.ProteinPicker.ProteinSequences)
        for Index in range(DatabaseLen):
            Seq1 = self.ProteinPicker.ProteinSequences[Index]
            for Jndex in range(Index + 1, DatabaseLen):
                Seq2 = self.ProteinPicker.ProteinSequences[Jndex]
                if Seq1 == Seq2:
                    #identical.  See if we have a list started
                    PName1 = self.ProteinPicker.ProteinNames[Index].split()[0]
                    PName2 = self.ProteinPicker.ProteinNames[Jndex].split()[0]
                    if not self.GenomicIdentityList.has_key(PName1):
                        self.GenomicIdentityList[PName1] = []
                        self.GenomicIdentityList[PName1].append(PName1)
                    self.GenomicIdentityList[PName1].append(PName2)


    def ParseGenomeIdentityList(self):
        """Given the input database, generate a list of identical loci.  this
        really pisses me off, because I can't get a good list of identities
        """
        if not self.GenomicIdentityListFile:
            return
        Handle = open(self.GenomicIdentityListFile, "rb")
        LinesSeen = {} # hack because my stupid program put out multiple lines
        for Line in Handle.xreadlines():
            if LinesSeen.has_key(Line):
                continue
            LinesSeen[Line] = 1
            Bits = Line.strip().split(" ")
            P1 = Bits[2]
            P2 = Bits[4]
            if not self.GenomicIdentityList.has_key(P1):
                self.GenomicIdentityList[P1] = []
                self.GenomicIdentityList[P1].append(P1)
            self.GenomicIdentityList[P1].append(P2)
        Handle.close()


    def RemoveOneHitWonders(self):
        """Remove all the proteins that have only one peptide in them.
        They don't count towards our end numbers anyway
        """
        ToRemove = []
        for (Key, List) in self.MSProteinPeptideList.items():
            if len(List) < 2:
                ToRemove.append(Key)
        for Key in ToRemove:
            del self.MSProteinPeptideList[Key]

    def WriteResults(self):
        """Here we are writing the the identical protein list with the
        help of an identity list. basically, if the proteins are not
        identical according to my list, then we are going to ignore them
        """
        TotalProteinCounter = 0
        GroupCounter = 0
        Handle = open(self.OutputFile, "wb")
        for GroupID in self.MSMSIdenticalProteins.keys():
            ProteinIDs = self.MSMSIdenticalProteins[GroupID]
            ProteinNames = self.GetProteinNames(ProteinIDs )
            IdentityKey = self.AreTheseAlsoGenomeIdentical(ProteinNames)
            ## now check to see if the accessions exist in our list
            if IdentityKey:
                IdentityValue = self.GenomicIdentityList[IdentityKey]
                Handle.write("Genomic Indistinguishable Proteins: %s\n"%IdentityValue)
                Handle.write( "MSMS Indistinguishable Proteins: %s\n"%ProteinNames)
                self.KeepTrackOfSameLocus(ProteinNames)
                #print  "Proteins: %s\n"%ProteinNames
                Handle.write( "Peptides: %s\n\n"%self.MSProteinPeptideList[ProteinIDs[0]])
                #print "Peptides: %s\n"%self.MSProteinPeptideList[ProteinIDs[0]]
                GroupCounter += 1
                TotalProteinCounter += len(IdentityValue)
        print "%d clusters comprising %d proteins"%(GroupCounter, TotalProteinCounter)
        print "%d of these are from a single genome locus"%self.SameLocusCounter
        Handle.close()

    def KeepTrackOfSameLocus(self, ListOfNames):
        """See if the list of names are all from the same locus or not
        """
        FirstFullAccession = ListOfNames[0]
        Dot = FirstFullAccession.find(".")
        Standard = FirstFullAccession[:Dot]
        for Name in ListOfNames:
            if not Name[:Dot] == Standard:
                return 0 #someone does not match
        self.SameLocusCounter += 1

    def AreTheseAlsoGenomeIdentical(self, MSMSIdenticalProteins):
        """A very ugly method to see if a list of proteins, which are MS/MS identical, are
        also genome identical.  Only return if the list are exactly identical
        """
        ## now see if any of the proteins are in the list
        FoundInGenomeIdentity = None
        for Protein in MSMSIdenticalProteins:
            if self.GenomicIdentityList.has_key(Protein):
                FoundInGenomeIdentity = Protein
                break

        if not FoundInGenomeIdentity:
            ## None of these proteins make it into our dictionary of identital proteins
            return 0
        GenomeIdenticalProteins = self.GenomicIdentityList[FoundInGenomeIdentity]
        # this function was mostly for debugging
        #self.SameSizeIdentities(MSMSIdenticalProteins, GenomeIdenticalProteins)
        if len(MSMSIdenticalProteins) == len(GenomeIdenticalProteins):
            return FoundInGenomeIdentity
        return 0

    def SameSizeIdentities(self, MSMSIdenticalProteins, GenomeIdenticalProteins):
        """Another method aimed to discover whether the two lists are the same
        Here I just see if the sizes are the same, meaning whether the size of the
        MS/MS identity list is the same number of proteins as the genomicIdentityList
        """
        MSMSLen = len(MSMSIdenticalProteins)
        GenomicLen = len(GenomeIdenticalProteins)
        if MSMSLen == GenomicLen:
            # good sign.
            self.ListsAreTheSame["Same"] += 1
        elif MSMSLen > GenomicLen:
            # okay.  there are not enough peptides to exclude some of the
            #proteins that are not truly part of this cluster
            self.ListsAreTheSame["MSMSBigger"] += 1
            #print "MSMSBigger"
            #print MSMSIdenticalProteins
            #print GenomeIdenticalProteins
        else:
            # a strange situation depending on your version of genomic identity
            self.ListsAreTheSame["GenomicBigger"] += 1
            #print "GenomeBigger"
            #print MSMSIdenticalProteins
            #print GenomeIdenticalProteins

    def GetProteinNames(self, IDList):
        """returns only the accession of the protein (AT1G123456.1)
        """
        ToReturn = []
        for ID in IDList:
            Protein = self.ProteinPicker.ProteinNames[ID]
            Accession = Protein.split()[0]
            ToReturn.append(Accession)
        return ToReturn

    def FindMSMSIdenticalProteins(self):
        """ go through the MSProteinPeptideList and find which proteins
        have an identical list of peptides.
        """
        Counter = 0
        SkipID = {} # for IDs that have already been found to be indistinguishable
        ProteinIDs = self.MSProteinPeptideList.keys()
        for Index in range(len(ProteinIDs)):
            ProteinIDi = ProteinIDs[Index]
            if SkipID.has_key(ProteinIDi):
                continue # we've already looked at this thing before and know who it is
            Group = [ProteinIDi,]
            ILen = len(self.MSProteinPeptideList[ProteinIDi])
            #use length for a quick hack at comparing
            for Jndex in range(Index + 1, len(ProteinIDs)):
                ProteinIDj = ProteinIDs[Jndex]
                if SkipID.has_key(ProteinIDj):
                    continue # we've already looked at this thing before and know who it is
                JLen = len(self.MSProteinPeptideList[ProteinIDj])
                if not ILen == JLen:
                    continue
                # at least the same length of peptide sequences attributed to these
                # two proteins, so lets to a good compare
                if self.MSProteinPeptideList[ProteinIDj] == self.MSProteinPeptideList[ProteinIDi]:
                    #these are indistinguishable
                    Group.append(ProteinIDj)
            # now we've searched through our inner loop and possibly have a list in the
            #Group variable.  Add the group, if it exists, and update the skip var
            if len(Group) > 1:
                self.MSMSIdenticalProteins[Counter]=Group
                Counter += 1
                for ID in Group:
                    SkipID[ID] = 1


    def PutPeptidesInProteins(self):
        """Go through the list in self.AllReferenceAnnotations and put
        each peptide wherever it matches.
        """
        for Aminos in self.AllReferenceAnnotations.keys():
            #print Aminos
            Locations = self.ProteinPicker.FindPeptideLocations(Aminos)
            for (ID, Residue) in Locations:
                if not self.MSProteinPeptideList.has_key(ID):
                    self.MSProteinPeptideList[ID] = []
                self.MSProteinPeptideList[ID].append(Aminos)
                #print self.ProteinPicker.ProteinNames[ID]
                #print self.MSProteinPeptideList[ID]
        self.AllReferenceAnnotations = {}  # set to null just for the memory savings
        
    def ParseReferenceFileCallback(self, FilePath):
        """Called by the ResultsParser.ProcessResultsFiles
        Here I parse out lines and get annotations, Putting them in a hash for safe keeping
        """
        Handle = open(FilePath, "rb")
        for Line in Handle.xreadlines():
            #print "Parsing line %s"%Line
            if Line[0] == "#":
                continue # comment line
            if not Line.strip():
                continue
            Bits = list(Line.split("\t"))
            try:
                Annotation = Bits[self.Columns.Annotation]
                Peptide = GetPeptideFromModdedName(Annotation)
                Aminos = Peptide.Aminos
            except:
                traceback.print_exc()
                continue # SNAFU
            
            self.SpectrumCount += 1
            if not self.AllReferenceAnnotations.has_key(Aminos):
                self.AllReferenceAnnotations[Aminos] = 0 
            self.AllReferenceAnnotations[Aminos] += 1
        Handle.close()


    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:d:w:i:")
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
            if Option == "-d":
                if not os.path.exists(Value):
                    print "** Error: couldn't find database file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.ReferenceDatabasePath.append( Value)
            if Option == "-w":
                self.OutputFile = Value
            if Option == "-i":
                self.GenomicIdentityListFile = Value
        if not OptionsSeen.has_key("-r")  or not OptionsSeen.has_key("-w") or not OptionsSeen.has_key("-d"):
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