UsageInfo="""WrapClustalW.py
This script is a wrapper and parser for the ClustalW program.
It takes an orthology list, and a set of databases makes a
temp.fasta and has clustal align it.

Required Options:
 -o [FileName] the orthology file 
 -c [TrieFile] Database of S. cerevisiae protein sequences
 -p [TrieFile] Database of S. paradoxus protein sequences
 -b [TrieFile] Database of S. bayanus protein sequences
 -m [TrieFile] Database of S. mikatae protein seqeunces
"""

import os
import sys
import getopt
import SelectProteins

class WrapperClass:
    def __init__ (self):
        self.OrthologyFilePath = None
        self.CerevisiaeDBPath = []
        self.ParadoxusDBPath = []
        self.BayanusDBPath = []
        self.MikataeDBPath = []
        self.ValidSpecies = ["S. paradoxus", "S. bayanus", "S. mikatae"]
        self.ClustalPath = "C:\source\ClustalW\clustalw2.exe"
        self.ClustalFilesStub = "C:\source\ClustalW\OutputDump"
        self.Orthology = {} # key = S.cervisiae locus, value = [(species, orf), (species, orf), ..]
        
    def Main(self):
        self.Cerevisiae = SelectProteins.ProteinSelector()
        self.Cerevisiae.LoadMultipleDB(self.CerevisiaeDBPath)
        self.Paradoxus = SelectProteins.ProteinSelector()
        self.Paradoxus.LoadMultipleDB(self.ParadoxusDBPath)
        self.Bayanus = SelectProteins.ProteinSelector()
        self.Bayanus.LoadMultipleDB(self.BayanusDBPath)
        self.Mikatae = SelectProteins.ProteinSelector()
        self.Mikatae.LoadMultipleDB(self.MikataeDBPath)
        self.ParseOrthology()
        self.ClustalStuff()


    def ParseOrthology(self):
        """go through and put genes in the orthology block
        """
        ORFColumn = 1
        SpeciesColumn = 3
        SourceColumn = 4
        TypeColumn = 5
        CerevisiaeColumn = 0
        Handle = open(self.OrthologyFilePath, "rb")
        for Line in Handle.xreadlines():
            if Line[0] == "#":
                continue
            if not Line.strip():
                continue
            Bits = Line.strip().split("\t")
            CerevisiaeLocus = Bits[CerevisiaeColumn]
            Species = Bits[SpeciesColumn]
            ORF = Bits[ORFColumn]
            Type = Bits[TypeColumn]
            DataSource = Bits[SourceColumn]
            ##now we weed out stuff we don't want
            ## 1. stuff from washu
            ## 2. stuff from species we don't care about
            ## 3. stuff thats not an ortholog
            if not DataSource == "MIT":
                continue
            if not Species in self.ValidSpecies:
                #print Species
                continue
            if not Type in ["ortholog", "paralog"]:
                #some paralogs and unresolved really mess up the alignment
                continue
            if not self.Orthology.has_key(CerevisiaeLocus):
                self.Orthology[CerevisiaeLocus] = []
            Tuple = (Species, ORF)
            self.Orthology[CerevisiaeLocus].append(Tuple)
            
        

    def ClustalStuff(self):        
        """Given the orthology groups that have already been parsed
        find sequences and clustal them
        """
        for CerevisiaeLocus in self.Orthology.keys():
            Sequences = {} # (species, ORF) =>Sequence
            CKey = ("S. cerevisiae", CerevisiaeLocus)
            Sequences[CKey] = self.GetCerevisiaeSequence(CKey)
            ListOfOrthologs = self.Orthology[CerevisiaeLocus]
            if len(ListOfOrthologs) == 0:
                ## should not get here, but it's a security catch
                continue
            for Key in ListOfOrthologs: #(species, ORF)
                Sequences[Key] = self.GetProteinSequence(Key)
            FilePath = self.MakeTempFastaFile(Sequences)
            Command = "%s %s"%(self.ClustalPath, FilePath)
            print "Submitting: %s"%Command
            os.system(Command)
            

    def GetCerevisiaeSequence(self, Key):
        """Cerevisiae database have a different format
        >YAL001C TFC3 SGDID:S000000001, Chr I from 151168-151099,151008-147596, reverse complement, Verified ORF, "Largest of six subunits of the RNA polymerase III transcription initiation factor complex (TFIIIC); part of the TauB domain of TFIIIC that binds DNA at the BoxB promoter sites of tRNA and similar genes; cooperates with Tfc6p in DNA binding"

        So just pull out the orf part of the key        
        """
        (Species, ORFName) = Key
        for (ID, Name) in self.Cerevisiae.ProteinNames.items():
            Bits = Name.split()
            if Bits[0] == ORFName:
                #print ORFName, Name
                return self.Cerevisiae.ProteinSequences[ID]

    def GetProteinSequence(self, Key):
        """Given a key (species, ORF), go through and find the sequence
        all of the sensu stricto MIT stuff has a fast line like these
        >ORFP:5209 YDR378C, Contig c454 4043-4303 reverse complement
        or
        >ORFP:15158 Multiple, Contig c955 2963-4480 reverse complement

        They all have ORFP and then the orf number.  so we use that to parse out.    
        """
        (Species, ORF) = Key
        ORFName = "ORFP:%s"%ORF
        if Species == "S. paradoxus":
            for (ID, Name) in self.Paradoxus.ProteinNames.items():
                Bits = Name.split()
                if Bits[0] == ORFName:
                    #print ORFName, Name
                    return self.Paradoxus.ProteinSequences[ID]
        if Species == "S. bayanus":
            for (ID, Name) in self.Bayanus.ProteinNames.items():
                Bits = Name.split()
                if Bits[0] == ORFName:
                    #print ORFName, Name
                    return self.Bayanus.ProteinSequences[ID]
        if Species == "S. mikatae":
            for (ID, Name) in self.Mikatae.ProteinNames.items():
                Bits = Name.split()
                if Bits[0] == ORFName:
                    #print ORFName, Name
                    return self.Mikatae.ProteinSequences[ID]
                

    def MakeTempFastaFile(self, Sequences):
        """Makes a fasta file from the sequences given
        The first sequence is from Cerevisiae, and it is going to be our key for
        making the filename
        """
        cORF = None
        for (Species, ORF) in Sequences.keys():
            if Species == "S. cerevisiae":
                cORF = ORF
                break
        if not cORF:
            print "No cerevisiae orf found for Sequences"
            print Sequences
            return
        cORFFile = "%s.txt"%cORF            
        FilePath = os.path.join(self.ClustalFilesStub, cORFFile)
        Handle = open(FilePath, "wb")
        for (Species, ORF) in Sequences.keys():
            ## species name all appear like "S. cerevisiae" or "S. paradoxus"
            ## clustal takes the first space delimited thing as the name
            ## and thus everything has the same name according to it
            ## we correct this by removing the space and cat the ORF
            SpeciesNoSpace = Species.replace(" ", "")
            FastaHeader = ">%s.%s"%(SpeciesNoSpace, ORF)
            ThisSequence = Sequences[(Species, ORF)]
            if not ThisSequence:
                ## sometimes there are weird things throw in for an orf name
                ## like Smik_Contig1865.3, which is not in our sequence database
                ## so we have to skip these
                continue
            Handle.write("%s\n%s\n"%(FastaHeader, ThisSequence ))

        Handle.close()
        return FilePath

    def ParseAlignmentResults(self):
        """basically counts the '*' in a file, returns
        MAJOR Assumption: you have no '*' in the protein name or sequence.  Doing such would be DUMB, and make the method fail
        """
        AlignmentResult = self.TempFileName.replace(".fasta", ".aln")
        Handle = open(AlignmentResult, "rb")
        Alignment = Handle.read()
        Handle.close()
        Count = Alignment.count("*")
        return Count


    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "o:c:p:b:m:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-o":
                self.OrthologyFilePath = Value
            if Option == "-c":
                # -r results file(s)
                if not os.path.exists(Value):
                    print "** Error: couldn't find database file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.CerevisiaeDBPath.append(Value)
            if Option == "-p":
                # -r results file(s)
                if not os.path.exists(Value):
                    print "** Error: couldn't find database file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.ParadoxusDBPath.append(Value)
            if Option == "-b":
                # -r results file(s)
                if not os.path.exists(Value):
                    print "** Error: couldn't find database file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.BayanusDBPath.append(Value)
            if Option == "-m":
                # -r results file(s)
                if not os.path.exists(Value):
                    print "** Error: couldn't find database file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.MikataeDBPath.append(Value)
                
        if not OptionsSeen.has_key("-o")  :
            print UsageInfo
            sys.exit(1)


if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(psyco not found - running in non-optimized mode)"
    DoStuff = WrapperClass()
    DoStuff.ParseCommandLine(sys.argv[1:])
    DoStuff.Main()        