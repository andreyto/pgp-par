UsageInfo = """CompareCorrectionsInOrthologs.py
Given a set of corrections from a proteogenomics analysis, we
look in neighboring genomes (through ortholog clusters) to
see how they should be altered 

Required Parameters
 -c [FileName] File containing the list of cluster of orthologs
 -b [FileName] Bad genes removed from the original genome
 -n [FileName] Novel genes added in the original genome
 -f [FileName] First observed amino acid by gene
 -m [FileName] Mapping file that converts between accessions
 -w [FileName] Output from this program
"""

import sys
import os
import getopt
import traceback
from Bio import AlignIO

class FinderClass:
    def __init__(self):
        self.BadGenesPath = None
        self.NewGenesPath = None
        self.MappingPath = None
        self.ClusterListPath = None
        self.FirstObservedPath = None
        self.OutputPath = "RenameYourOutput.txt"
        self.GenomeDict = {}
        self.GenomeDict["ntyp10"] = "NC_010159 Yersinia pestis Angola"
        self.GenomeDict["ntyp11"] = "NC_008150 Yersinia pestis Antiqua"
        self.GenomeDict["ntyp12"] = "NC_005810 Yersinia pestis biovar Microtus str. 91001"
        self.GenomeDict["ntyp13"] = "NC_003143 Yersinia pestis CO92"
        self.GenomeDict["ntyp14"] = "NC_004088 Yersinia pestis KIM"
        self.GenomeDict["ntyp15"] = "NC_008149 Yersinia pestis Nepal516"
        self.GenomeDict["ntyp16"] = "NC_009381 Yersinia pestis Pestoides F"
        self.GenomeDict["ntyp17"] = "NC_009708 Yersinia pseudotuberculosis IP 31758"
        self.GenomeDict["ntyp18"] = "NC_006155 Yersinia pseudotuberculosis IP 32953"
        self.GenomeDict["ntyp19"] = "NC_010634 Yersinia pseudotuberculosis PB1/+"
        self.GenomeDict["ntyp20"] = "NC_010465 Yersinia pseudotuberculosis YPIII"
        self.GenomeDict["ntye02"] = "NC_008800 Yersinia enterocolitica subsp. enterocolitica 8081"
        
        
    def Main(self):
        """Parameters: None
        Return: NOne
        Description: the main, start, head dude
        """
        #first we need to read in all the changes in the original genome
        if self.BadGenesPath:
            BadGenes = self.ParseTBLChanges(self.BadGenesPath)
            BadErgatisGenes = self.ParseMapping(self.MappingPath, BadGenes)
            #now roll through the list of clusters and find stuff
            self.WrapClusterParsing(self.ClusterListPath, BadErgatisGenes)
        if self.NewGenesPath:
            NewGenes = self.ParseTBLChanges(self.NewGenesPath)
            NewErgatisGenes = self.ParseMapping(self.MappingPath, NewGenes)
            #now roll through the list of clusters and find stuff
            #self.WrapClusterParsing(self.ClusterListPath, NewErgatisGenes)
            self.ParseClusterForGeneOfInterest(self.ClusterListPath, NewErgatisGenes, self.AnalyzeNovelGeneCluster)
        if self.FirstObservedPath:
            FirstObservation = self.ParseFirstObserved(self.FirstObservedPath)
            ErgatisMapping = self.ParseMapping(self.MappingPath, FirstObservation.keys())
            #now go through and parse all the clusters to find out whether
            #the orthologs are correctly predicted in extant genomes
            self.WrapFindStart(self.ClusterListPath, FirstObservation, ErgatisMapping)

    def ParseClusterForGeneOfInterest(self, ListOfClusterPaths, GenesOfInterest, CallForward):
        """Parameters: a path to a file which contains a list of paths to clusters (clustalw currently)
                       A list of genes that I care about
                       A function to call when I match a cluster with a gene
        Return: None
        Description: the wrapper and sorter.
        """
        Handle = open (ListOfClusterPaths, "rb")
        print GenesOfInterest
        for Line in Handle.xreadlines():
            ClusterPath = Line.strip() # path to the clustalw alignment file
            Alignment = AlignIO.read(open(ClusterPath), "clustal") # now parsed out
            SeqObjectsArray = Alignment.get_all_seqs() # all the sequence objects
            #now check if any of these sequence objects are in our mapping dictionary, which we care about
            for SeqObj in SeqObjectsArray:
                SeqID = SeqObj.id
                if SeqID in GenesOfInterest:
                    # hey I found someone I care about.  Let's do what they suggest.
                    ArgsList = [SeqID, ClusterPath]
                    apply(CallForward, ArgsList)
                    break #out of the SeqObj loop

    def AnalyzeNovelGeneCluster(self, ID,  ClusterPath):
        """Parameters: the ID of the gene of Interest
                        the parsed BioPython AlignmentObject
        Return: None
        Description: Look at the cluster of novel genes.  Because I care about them.
        """
        print "I found a novel gene %s in %s"%(ID, ClusterPath)
        

            
    def ParseFirstObserved(self, FilePath):
        """Parameters: a path to the file with the index of the first observed amino acid
        Return: a dictionary with RefSeq -> zero-based index
        Description : file parser, created elsewhere
        """
        Dict = {}
        Handle = open(FilePath, "rb")
        for Line in Handle.xreadlines():
            Bits = Line.strip().split("\t")
            FastaLine = Bits[0]
            Observation = int (Bits[1])
            #now parse out the fasta crap
            FastaBits = FastaLine.split("|")
            RefSeq = FastaBits[3]
            Dict[RefSeq] = Observation
        return Dict
        
    def WrapClusterParsing(self, ListFile, BadGeneDict):
        """Parameters: file with list of clusters, genes to look for
        Return: none
        Description: wrap the file searching
        """
        Handle = open (ListFile, "rb")
        print BadGeneDict
        for Line in Handle.xreadlines():
            #this is a new file to go out and parse
            ClusterPath = Line.strip()
            self.ParseCluster(ClusterPath, BadGeneDict)
            #break #debug breakout

    def WrapFindStart (self, ListFile, ObservationDict, MappingDict):
        """Parameters: file with list of Clusters, Dictionary with first observed amino acid
                dictionary with mappings of gene names
        Return: none
        Description: Given observed peptides in our target genome, we want to find out
        if we have problems in our extant genomes.  This wraps on to a single file
        """
        Handle = open (ListFile, "rb")
        for Line in Handle.xreadlines():
            ClusterPath = Line.strip()
            Alignment = AlignIO.read(open(ClusterPath), "clustal")
            SeqObjectsArray = Alignment.get_all_seqs()
            #now check if any of these sequence objects are in our mapping dictionary, which we care about
            InObservedList =None
            ObservedList = MappingDict.keys() #the ergatis ID. STupid ergatis mapping through multiple id schemes
            for SeqObj in SeqObjectsArray:
                SeqID = SeqObj.id
                #first I check for the
                if SeqID in ObservedList:
                    InObservedList = SeqID
                    break
            #now if we care about this one, then let's send it out to be looked at
            if not InObservedList:
                continue #next round here I come
            #first let's check to see if we have the right number of proteins in the alignment
            #now we send it out to checkers, I know the ID that I care about
            #so let's get the real ID for it, what we can find stuff with.
            RefSeqFound = MappingDict[InObservedList]
            ObservedAA = ObservationDict[RefSeqFound]
            #MissingAGene = self.RequireAllGenomes(Alignment, RefSeqFound, ClusterPath)
            self.FindStartAgreement(Alignment, InObservedList, ObservedAA, ClusterPath, MappingDict)
            
    def RequireAllGenomes(self, AlignmentObject, RefSeqFound, ClusterPath):
        """
        """
        #first check to see that there are at least 12
        SeqArray = AlignmentObject.get_all_seqs()
        if len(SeqArray) >= 12:
            return 0
        ListOfGenomes = self.GenomeDict.keys()
        for Seq in SeqArray:
            #get the genome name
            Gene = Seq.id
            GeneBits = Gene.split(".")
            Database = GeneBits[0]
            if not Database in ListOfGenomes:
                    #SNAFU parsing, or potentially paralogs clustering together
                    pass
            else:
                ListOfGenomes.remove(Database)
        #now here we should have just the missing ones left in the list
        print "Alignment for %s (%s) had only %s sequences"%(RefSeqFound, ClusterPath, len(SeqArray))
        for Genome in ListOfGenomes:
            print "\tMissing genome %s, %s"%(Genome, self.GenomeDict[Genome])
        return -1

    def FindStartAgreement(self, AlignmentObject, TargetGene, ObservedAA, ClusterPath, MappingDict):
        """Parameters: A BioPython Alignment object, Accession (ergatis) of the observed gene
                The index (zero based) of the first observed aa
                dictionary with mappings of gene names
        Return: none
        Description: Given observed peptides in our target genome, we want to find out
        if we have problems in our extant genomes.  This method looks into the alignment
        and finds those that have a blank or - with nothing upstream.  We also are looking
        for whether each genome has a representative.
        """
        #first we have to find the column within the alignment of the first observed in target
        ObservedColumn = self.FindObservationColumn(AlignmentObject, TargetGene, ObservedAA)
        RefSeqFound = MappingDict[TargetGene]
        BadNews = 0
        if ObservedColumn < 0:
            #SNAFU. Abort
            return -1
        #now we look into whether the observed column is found in all of the alignments
        FirstPrint = 0
        for SeqObj in AlignmentObject.get_all_seqs():
            AlnSeq = SeqObj.seq
            BadNews = 0 #new sequence, new day.
            #first check for a letter at this position
            if ObservedColumn > len(AlnSeq):
                #this is a problem.  We have an observed column that is smaller than our length
                #we should print out something here
                if not FirstPrint:
                    print "mispredictions %s, %s, %s"%(RefSeqFound, ObservedAA, ClusterPath)
                    FirstPrint = 1
                #now print off the specifics of this problem
                print "\t%s does not have enough columns to test the observation"%(SeqObj.id)
                continue #skip over any other computations for this sequence in the alignment
            #we have a letter here, so let's check it
            AlignColumn = AlnSeq[ObservedColumn]
            if not AlignColumn == "-":
                pass #so long as it's not a dash we're okay, because that means it's a letter
            else:
                #it was a dash, let's check the upstream slice and see if it's all dashes, which would be bad = sequence missing the 5' region
                AlnSlice = AlnSeq[:ObservedColumn]
                BadNews = 1 #there's a distinct possibility that we may have some bad news
                for Index in range(len(AlnSlice)):
                    if not AlnSlice[Index] == "-":
                        #whew!, we found something besides the dash.
                        BadNews =0 
                        break #out of the inner for over the slice
                #now we may have bad news
                if BadNews:
                    #let's print off what we know
                    if not FirstPrint:
                        print "mispredictions %s, %s, %s"%(RefSeqFound, ObservedAA, ClusterPath)
                        FirstPrint = 1
                    #now print off the specifics of this problem
                    print "\t%s was too short"%(SeqObj.id)
        
        


    def FindObservationColumn(self, AlignmentObject, TargetGene, ObservedAA):
        """Parameters: A BioPython Alignment object, Accession (ergatis) of the observed gene
                The index (zero based) of the first observed aa
        Return: index, -1 for SNAFU
        Description: find in the alignment of the target gene, the ith amino acid, which is the
        first one observed.
        """
        AlnSequence = None
        for SeqObj in AlignmentObject.get_all_seqs():
            if SeqObj.id == TargetGene:
                #this is the one we want
                AlnSequence = SeqObj.seq
                break
        #now we have a sequence, we'll check to make sure
        if not AlnSequence:
            return -1 # SNAFU
        #just going to do a dummy for loop to count out the columns of nothing up front
        for Index in range (len(AlnSequence)):
            if AlnSequence[Index] == "-":
                pass
            else:
                break
        #now Index should be the column of the starting met (or whatever)
        ObservedColumn = Index + ObservedAA
        return ObservedColumn


            
    def ParseCluster(self, FilePath, GenesToFind):
        """Parameters: Filepath of a clustal alignment, list of genes to find
        Return: a hash. key = genetofind, value = accessions in the cluster.
        Description: we are looking for a set of genes, and if we find one of
        them, then we want to know what else was in the same cluster.
        """
        #print "opening %s"%FilePath
        Handle = open(FilePath, "rb")
        AllAccessionsInCluster = []
        FoundAccession = None 
        for Line in Handle.xreadlines():
            #print Line
            #read all the things with a Name in front, and quit after //
            if Line[:2] == "//":
                break
            if Line[:5] == " Name":
                Bits = Line.strip().split(" ")
                ErgatisAcc = Bits[1]
                AllAccessionsInCluster.append(ErgatisAcc)
                if GenesToFind.has_key(ErgatisAcc):
                    FoundAccession = ErgatisAcc
                    print "I found an accession %s, %s\n%s"%(ErgatisAcc, GenesToFind[ErgatisAcc], FilePath)
        #done with our foor
        
        if not FoundAccession:
            return None
        #now we build our return
        Dictionary = {}
        Dictionary[FoundAccession] = AllAccessionsInCluster
        print "I found %s genes in my cluster"%len(AllAccessionsInCluster)
        print
        #note that this includes the found accession in the value list
        return Dictionary
            
        
    def ParseMapping(self, FilePath, LimitList):
        """Parameters: filepath, list of limited refseq accessions to care about
        Return: dict of ergatis accessions->refseq
        Description: go throug the mapping file and get back ergatis accessions
        """
        ToReturn = {}
        Handle = open (FilePath, "rb")
        for Line in Handle.xreadlines():
            Bits = Line.strip().split("\t")
            Ergatis = Bits[0]
            RefSeq = Bits[2]
            if RefSeq in LimitList:
                ToReturn[Ergatis] = RefSeq
        return ToReturn
        


    def ParseTBLChanges(self, FilePath):
        """Parameters: path to a file
        Return: list of refseq accessions
        Description: During the course of our studies, we may come across bogus
        gene calls, which we have removed in the original genomes.  This method
        parses out that list, from a table file.  We just look for
        protein_id     ref|1234556
        """
        BadGenes = []
        Handle = open (FilePath, "rb")
        for Line in Handle.xreadlines():
            if Line.find("protein_id") > -1:
                Bits = Line.strip().split("\t")
                #print Bits
                #now just strip off the bars
                RefSeqBits = Bits[-1].split("|")
                #print RefSeqBits
                RefSeqAcc = RefSeqBits[1]
                BadGenes.append(RefSeqAcc)
        return BadGenes
            
    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "c:b:w:m:n:f:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-c":
                # -r results file(s)
                if not os.path.exists(Value):
                    print "** Error: couldn't find results file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.ClusterListPath = Value
            elif Option == "-b":
                if not os.path.exists(Value):
                    print "** Error: couldn't find results file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.BadGenesPath = Value
            elif Option == "-n":
                if not os.path.exists(Value):
                    print "** Error: couldn't find results file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.NewGenesPath = Value
            elif Option == "-f":
                if not os.path.exists(Value):
                    print "** Error: couldn't find results file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.FirstObservedPath = Value
            elif Option == "-m":
                if not os.path.exists(Value):
                    print "** Error: couldn't find results file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.MappingPath = Value
            elif Option == "-w":
                self.OutputPath = Value
            else:
                print "Option %s not recognized.  Fatal Error."%Option
                print UsageInfo
                sys.exit(1)
        if not OptionsSeen.has_key("-c") or not OptionsSeen.has_key("-m"):
            print UsageInfo
            sys.exit(1)


if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(psyco not found - running in non-optimized mode)"
    Gumshoe = FinderClass()
    Gumshoe.ParseCommandLine(sys.argv[1:])
    Gumshoe.Main()
