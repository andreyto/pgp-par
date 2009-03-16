"""ExportPeptidesInGFF3.py


"""
UsageInfo = """CreateGFF3ForGlimmer.py
Given some peptide results, we create a GFF3 file for
visualization in GBrowse.  SingleExonPeptides only 


Required Options
 -r [Directory/File] Inspect Results
 -w [GFF file] the GFF3 output
 -d [TrieFile] the 6 frame translation of the bacterial genome

Additional Options
 -m MultiLocus Flag.  If set, peptides with multiple locations on the 
     database will be included in the output. Default is to output only 
     uniquely located peptides.
 -p [float] The pvalue cutoff for inclusion.  Default is 0.05
"""

import sys
import os
import getopt
import traceback
import ResultsParser
import PeptideMapper
from Utils import *
Initialize()


class ConvertFormatClass(ResultsParser.ResultsParser):
    def __init__(self):
        self.InspectResults = None
        self.OutputPath = "RenameYourResults.txt"
        self.DBPath = [] # so that you can input more than one, but they all better have the right format!!
        self.ObservedPeptides = {} # simple list of peptide sequences (stripped of anything else), aminos = best pvalue
        self.AllowMultiLocusPeptides = 0
        self.ChromosomeName = None  # used to overwrite the chromosome name from the database, if needed
        self.PValueLimit = 0.05
        ResultsParser.ResultsParser.__init__(self)

        
    def Main(self):
        self.ProcessResultsFiles(self.InspectResults, self.ParseInspectCallback) # method inhereted from ResultsParser
        self.WriteToGFF3(self.ObservedPeptides)
        
        
    def WriteToGFF3(self, ObservedPeptides):
        """With the input parameter of a dictionary of observed peptides, I now find the genomic location
        and then send this out to the self.GFF3Writer.  Note that I have a parameter input, not defaulting tothe self.Observedpeptides
        in case anyone else ever wants to use this function
        """
        OutHandle = open(self.OutputPath, "wb")
        LocationCounter = 0
        PeptideCounter = 0
        NotInDatabaseCounter =0 
        MapMan = PeptideMapper.PeptideMappingClass()
        MapMan.LoadDatabases(self.DBPath)
        print "ExportPeptidesInGFF3::WriteToGFF3, %s to process"%len(ObservedPeptides)
        for (Peptide, FDR) in ObservedPeptides.items():
            Locations = MapMan.MapMe(Peptide) #already stripped down to simple Amino string
            if not self.AllowMultiLocusPeptides:
                #if they don't want to allow them, we check, otherwise we don't worry about it.
                if len(Locations) > 1:
                    continue
            if len(Locations) == 0:
                #print "WARNING: no location found for peptide:%s, perhaps in CommonContaminants?"%Peptide
                NotInDatabaseCounter += 1
                continue
            PeptideCounter += 1
            for Location in Locations:
                if self.ChromosomeName:
                    #overwrite the chromosomeName with what the user input
                    Location.Chromosome = self.ChromosomeName
                Line = Location.GetGFF3Line(FDR, LocationCounter)
                LocationCounter += 1
                OutHandle.write("%s\n"%Line)
            if (PeptideCounter % 500) == 0:
                print "WriteToGFF3, %s peptides of %s"%(PeptideCounter, len(ObservedPeptides))
        OutHandle.close()
        print "wrote %s locations for %s peptides"%(LocationCounter, PeptideCounter)
        print "There were %d peptides not in the supplied database"%NotInDatabaseCounter


    def ParseInspectCallback(self, FilePath):
        """Called by the ResultsParser.ProcessResultsFiles
        Here I parse out lines and get a big list of peptides and pvalues
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
            ## HACK BELOW. looking to see if people have done the LFDR assignment, if so, we prefer it
            if len(Bits) < self.Columns.LFDR: #meaning that I don't have the column in question
                    continue
            else:
                LFDR = float(Bits[self.Columns.LFDR])
                PValue = LFDR
            if PValue > self.PValueLimit:
                continue
            if not self.ObservedPeptides.has_key(Aminos):
                self.ObservedPeptides[Aminos] = 1.0 # crappy default pvalue
            if PValue < self.ObservedPeptides[Aminos]:
                self.ObservedPeptides[Aminos] = PValue
            

    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:w:d:mc:p:")
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
            elif Option == "-c":
                self.ChromosomeName = Value
            elif Option == "-p":
                self.PValueLimit = float(Value)
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