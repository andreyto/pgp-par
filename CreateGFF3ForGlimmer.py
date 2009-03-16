"""CreateGFF3ForGlimmer.py
This makes a GFF file out of some PROKARYOTIC results and makes a
GFF3 for Glimmer.

"""
UsageInfo = """CreateGFF3ForGlimmer.py

Required Options
 -r [Directory/File] Inspect Results
 -w [GFF file] the GFF file for Glimmer
 -d [TrieFile] the 6 frame translation of the bacterial genome

Additional Options
 -m MultiLocus Flag.  If set, peptides with multiple locations on the 
     database will be included in the output. Default is to output only 
     uniquely located peptides.
"""

import sys
import os
import getopt
import traceback
import ResultsParser
import SelectProteins
import ProkaryoticGFF3
from Utils import *
Initialize()


class ConvertFormatClass(ResultsParser.ResultsParser):
    def __init__(self):
        self.InspectResults = None
        self.OutputPath = "RenameYourResults.txt"
        self.DBPath = [] # so that you can input more than one, but they all better have the right format!!
        self.ObservedPeptides = {} # simple list of peptide sequences (stripped of anything else), value = best pvalue
        self.AllowMultiLocusPeptides = 0
        ResultsParser.ResultsParser.__init__(self)

        
    def Main(self):
        self.GFF3Writer = ProkaryoticGFF3.GFFClass()
        self.ProteinPicker = SelectProteins.ProteinSelector()
        self.ProteinPicker.LoadMultipleDB(self.DBPath)
        self.ProcessResultsFiles(self.InspectResults, self.ParseInspectCallback) # method inhereted from ResultsParser
        self.WriteToGFF3(self.ObservedPeptides)
        
        
    def WriteToGFF3(self, ObservedPeptides):
        """With the input parameter of a dictionary of observed peptides, I now find the genomic location
        and then send this out to the self.GFF3Writer.  Note that I have a parameter input, not defaulting tothe self.Observedpeptides
        in case anyone else ever wants to use this function
        """
        OutHandle = open(self.OutputPath, "wb")
        LineCounter = 0
        PeptideCounter = 0
        NotInDatabaseCounter =0 
        for (Peptide, FDR) in ObservedPeptides.items():
            # as noted during parsing, Inspect uses a FDR and not a pvalue.  here we switch
            #print "The Inspect PValue is %s"%FDR
            PValue = 1-FDR
            Locations = self.ProteinPicker.FindPeptideLocations(Peptide) #already stripped down to simple Amino string
            if not self.AllowMultiLocusPeptides:
                #if they don't want to allow them, we check, otherwise we don't worry about it.
                if len(Locations) > 1:
                    continue
            if len(Locations) == 0:
                #print "WARNING: no location found for peptide:%s, perhaps in CommonContaminants?"%Peptide
                NotInDatabaseCounter += 1
                continue
            PeptideCounter += 1
            for (ProteinID, StartResidue) in Locations:
                ProteinName = self.ProteinPicker.ProteinNames[ProteinID]
                Line = self.GFF3Writer.WriteProteomicEvidenceLine(Peptide, ProteinName, StartResidue, PValue)
                LineCounter += 1
                OutHandle.write("%s\n"%Line)
        OutHandle.close()
        print "wrote %s locations for %s peptides"%(LineCounter, PeptideCounter)
        print "There were %d peptides not in the supplied database"%NotInDatabaseCounter


    def ParseInspectCallback(self, FilePath):
        """Called by the ResultsParser.ProcessResultsFiles
        Here I parse out lines and get a big list of peptides
        """
        Handle = open(FilePath, "rb")
        for Line in Handle.xreadlines():
            if not Line.strip():
                continue # dont' know why this would happen, but Stephen had it in all of his
            Line = Line.strip()
            if Line[0] == "#": # the header character.  see how its a comment.
                continue
            Bits = Line.split("\t")
            Annotation = Bits[self.Columns.Annotation]
            PValue = float(Bits[self.Columns.PValue]) #note that Inspect uses an FDR (small number good) instead of a pvalue (big number good)
            Peptide = GetPeptideFromModdedName (Annotation) # in Utils.py
            Aminos = Peptide.Aminos # raw amino acid sequence, stripped of PTMs or flanking AAs
            if not self.ObservedPeptides.has_key(Aminos):
                self.ObservedPeptides[Aminos] = 1.0 # crappy default pvalue
            if PValue < self.ObservedPeptides[Aminos]:
                self.ObservedPeptides[Aminos] = PValue
            

    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:w:d:m")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-r":
                # -r results file(s)
                if not os.path.exists(Value):
                    print "** Error: couldn't find results file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.InspectResults = Value
            elif Option == "-d":
                # -r results file(s)
                if not os.path.exists(Value):
                    print "** Error: couldn't find database file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.DBPath.append( Value) # so that you can input more than one

            elif Option == "-w":
                self.OutputPath = Value
            elif Option == "-m":
                self.AllowMultiLocusPeptides = 1
            else:
                print "Option %s not understood"%Option
                print UsageInfo
                sys.exit(1)
        if not OptionsSeen.has_key("-r")  or not OptionsSeen.has_key("-d"):
            print UsageInfo
            sys.exit(1)



if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(psyco not found - running in non-optimized mode)"
    Boredom = ConvertFormatClass()
    Boredom.ParseCommandLine(sys.argv[1:])
    Boredom.Main()