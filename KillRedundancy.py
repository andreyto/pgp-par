UsageInfo = """ KillRedundancy.py
Takes an input file(s) and returns the non-redundant peptide
list.  For this it is assumed that you have already run the
PhosphateLocalization, and are concerned about the best
localized version of the peptide, so we keep the spectrum
example with the highest PLS.


Required Options
 -r [FileName] Filename or directory of annotations
 -w [FileName] stub output filename (e.g. S.Cerevisiae)
 -d [TrieFile] Database for the Search
"""


import os
import getopt
import sys
import ResultsParser
import SelectProteins
from Utils import *
Initialize()

class AssassinClass(ResultsParser.ResultsParser):
    def __init__(self):
        self.InputFile = None
        self.OutputFile = None
        self.Header = None
        self.UniquePeptides = {} # Peptide => (MQScore, Bits)
        self.UniqueSites = {} #site in DB => (PLS, Bits)
        self.DBPath = []  # array for possible multiple files
        ResultsParser.ResultsParser.__init__(self)
        
    def Main(self):
        self.ProteinPicker = SelectProteins.ProteinSelector()
        self.ProteinPicker.LoadMultipleDB(self.DBPath)
        self.ProcessResultsFiles(self.InputFile, self.ParseAnnotations)
        self.WritePeptideOutput()
        self.DetermineSiteRedundancy()
        self.WriteNonRedundantSites()

    def DetermineSiteRedundancy(self):
        """Here, after we have written out a simple list of
        non-redundant peptides (which should all have been assigned
        to the most confident site) we will now have to tackle the
        more difficult problem of non-redundant sites.  I don't think
        that there is a quick computational idea to make this easy, so
        first I will just do some simple stuff, and then perhaps revamp
        this for something more complicated.  We find sites location of
        the phosphate in a protein.  (FLAG non-unique peptide strings)
        and then keep a key of that.  Much like for peptides, we keep
        the best scoring site, here the criteria is the PLS and not MQScore
        """
        NonUniquePeptideCount = 0
        NonUniqueFilePath = "%s.NonUniquePeptides.txt"%self.OutputFile
        NonUniqueHandle = open( NonUniqueFilePath, "wb")
        NonUniqueHandle.write("%s\n"%self.Header)
        for Annotation in self.UniquePeptides.keys():
            (MQScore, Bits) = self.UniquePeptides[Annotation]
            #first get just the aminos
            Peptide = GetPeptideFromModdedName(Annotation)
            Aminos = Peptide.Aminos
            Locations = self.ProteinPicker.FindPeptideLocations(Aminos)
            if len(Locations) == 0:
                #bad.
                print "peptide %s not found in database"%Annotation
            if len(Locations) > 1:
                NonUniquePeptideCount += 1
                #multiple locations just means we don't have confidence
                # in the actual phosphorylation site.  We may do these
                #later manually, but for now we are leaving them out
                print "Annotation :%s: has %s locations"%(Annotation, len(Locations))
                (MQScore, Bits) = self.UniquePeptides[Annotation]
                ToPrint = "\t".join(Bits)
                NonUniqueHandle.write("%s\n"%ToPrint)
                continue
            #now we process the rest of it.  Peptides which have a unique protein match
            #get the phosphorylated residue, and the PLS
            (ProteinID, StartResidue) = Locations[0] # StartResidue is ZERO based, most people think 1 based
            StartResidue += 1                        # just remember that for printing out and stuff
            for Offset in Peptide.Modifications.keys():
                PTMClassList = Peptide.Modifications[Offset]
                for PTMClass in PTMClassList:
                    #print PTMClass.Name
                    if PTMClass.Name == "Phosphorylation":
                        #found a phosphorylation, and an offset in the peptide
                        PhosphorylatedResidue = StartResidue + Offset
                        PLS = Bits[-1] # not a stardard column.  Better be the last.!
                        Key = (ProteinID, PhosphorylatedResidue)
                        if self.UniqueSites.has_key(Key):
                            (OldPLS, OldBits) = self.UniqueSites[Key]
                            if PLS > OldPLS:
                                self.UniqueSites[Key] = (PLS, Bits)
                        else:
                            self.UniqueSites[Key] = (PLS, Bits)
                            
                        #print "Protein %s"%self.ProteinPicker.ProteinNames[ProteinID]
                        #print "Annotation %s"%Annotation
                        #print "StartOfPeptide %s, Phos %s"%(StartResidue, PhosphorylatedResidue)
            
            
        print "%d peptides appear in multiple locations"%NonUniquePeptideCount
        NonUniqueHandle.close()

    def WriteNonRedundantSites(self):
       """Given that we've already determine which sites there are, we will now print them out.
       """
       NonRedundantSitesPath = "%s.NonRedundantSites.txt"%self.OutputFile
       NonRedundantSitesHandle = open(NonRedundantSitesPath, "wb")
       NonRedundantSitesHandle.write("%s\t%s\n"%(self.Header, "PhosphorylatedResidue"))
       for Key in self.UniqueSites.keys():#Key = (ProteinID, PhosphorylatedResidue)
           (PLS, Bits) = self.UniqueSites[Key]
           (ProteinID, PhosphorylatedResidue) = Key
           #Sequence = self.ProteinPicker.ProteinSequences[ProteinID]
           #print "I'm looking for index %s of %s" %(PhosphorylatedResidue, len(Sequence))
           PhosphorylatedLetter = self.ProteinPicker.ProteinSequences[ProteinID][PhosphorylatedResidue-1] # remember to subtract one for the zero based
           NewData = "%s,%s"%(PhosphorylatedResidue, PhosphorylatedLetter)# cast to string, so I can string.join
           Bits.append(NewData) 
           ToPrint = "\t".join(Bits)
           NonRedundantSitesHandle.write("%s\n"%ToPrint)
        

        
    def ParseAnnotations(self, FileName):
        Handle = open(FileName, "rb")
        ResidueCol = self.Columns.FileOffset + 2 # ascore, then residue
        ExpectedLen = ResidueCol + 1 # zero based
        for Line in Handle.xreadlines():
            if not Line.strip():
                continue
            Line = Line.strip()
            if Line[0] == "#":
                self.Header = Line
                continue
            Bits = Line.split("\t")
            MQScore = float(Bits[self.Columns.MQScore])
            Annotation = Bits[self.Columns.Annotation]
            if not Annotation.find("phos") > 0:
                continue
            ##1. get things for the peptide list
            if not self.UniquePeptides.has_key(Annotation):
                self.UniquePeptides[Annotation] = (MQScore, Bits)
            else:
                (TopMQScore, TopBits) = self.UniquePeptides[Annotation]
                if MQScore > TopMQScore:
                    self.UniquePeptides[Annotation] = (MQScore, Bits)

        Handle.close()
                
    def WritePeptideOutput(self):
        OutName = "%s.NonRedundantPeptides.txt"%self.OutputFile
        Handle = open(OutName, "wb")
        Handle.write("%s\n"%self.Header)
        for Annotation in self.UniquePeptides.keys():
            (MQScore, Bits) = self.UniquePeptides[Annotation]
            String = "\t".join(Bits)
            Handle.write("%s\n"%String)
        Handle.close()

    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:w:d:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-r":
                # -r results file(s)
                if not os.path.exists(Value):
                    print "** Error: couldn't find results file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.InputFile = Value
            elif Option == "-w":
                self.OutputFile = Value
            elif Option == "-d":
                self.DBPath.append(Value)
        if not OptionsSeen.has_key("-r") or not OptionsSeen.has_key("-w"):
            print UsageInfo
            sys.exit(1)
            

if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(psyco not found - running in non-optimized mode)"
    Guido = AssassinClass()
    Guido.ParseCommandLine(sys.argv[1:])
    Guido.Main()                