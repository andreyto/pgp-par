"""PeptideMapper.py
This auxiliary set of classes/functions is to be used in mapping peptides to 
their genomic location.  hopefully, this will be the central place that we 
do all of this activity, regardless of the genome.

Usage: You can call the following two methods.  nothing else.
LoadDatabases(DBPaths)
MapMe(AminoAcidSequencs)

NOTE: because this is genomic mapping, the database must have the standard
notion of a genomic context. you can get this by using the SixFrameFasta.py
script to create your database.

NOTE: this is a utility, and not executable from the command line

In the past, scripts like GetGenomeCoordinatesForTAIR did something like this
but that was pretty specialized, and we hope to just make it simple here
"""

import sys
import os
import traceback
import SelectProteins
import GenomicLocations


class PeptideMappingClass:
    def __init__(self):
        self.DatabasePaths = [] #possibly multiple
        self.CurrentAminos = ""
        
    def LoadDatabases(self, DBPaths):
        self.DatabasePaths = DBPaths
        self.ProteinPicker = SelectProteins.ProteinSelector()
        self.ProteinPicker.LoadMultipleDB(self.DatabasePaths)

    def MapMe(self, Aminos, PValue, WarnNoMatch = 0):
        """Given an amino acid sequence, return (potentially multiple)
        GenomicLocationForPeptide object(s).
        This method parses information out of the protein header 
        if it was created using the SixFrameFasta.py (or has the exact 
        same format). 
        
        WARNING: Don't include a period '.' in your name. we use that to split
        
        Protein256370.Chr:Chr1.Frame1.StartNuc831372.Strand-
        XXX.Protein157841.Chr:Chr1.Frame1.StartNuc4033167.Strand+
        InfoBits[0] = protein unique id
        InfoBits[1] = chromosome
        InfoBits[2] = Frame
        InfoBits[3] = StartNucleotide of the ORF
        InfoBits[4] = Strand
        """
        ReturnLocations = []
        self.CurrentAminos = Aminos #only used for printing in case or error
        DBLocations = self.ProteinPicker.FindPeptideLocations(Aminos)
        if (len(DBLocations) < 1) and WarnNoMatch:
            print "Peptide %s was not found in the database"%Aminos
            
        for (ProteinID, PeptideStartAA) in DBLocations:
            ## PeptideStartAA is the start position within the amino acid sequence
            #1. parse out the features of the protein name (fasta header)
            Location = GenomicLocations.GenomicLocationForPeptide()
            ProteinName = self.ProteinPicker.ProteinNames[ProteinID]
            ### CRAP ! the fame proteins are XXX.Protein1, so the . screws up the splitting
            ## Total hack below.  Find a better splitter than . and redo the sixframefasta.py
            if ProteinName.find("XXX.") == 0:
                ProteinName = ProteinName.replace("XXX.", "XXX")
            InfoBits = ProteinName.split(".")
            #chromosome, strand
            Location.ProteinName = InfoBits[0]
            Location.Chromosome = InfoBits[1].replace("Chr:", "")
            Location.Strand = InfoBits[4].replace("Strand", "")
            Location.Aminos = Aminos
            Location.Score = PValue
            ORFStartNucleotide = int(InfoBits[3].replace("StartNuc", ""))
            ##can get complicated below, because of spliced stuff, which I am not
            ##going to do now.  Not in the least. But the method structure is
            ##here, so it can be expanded upon without really killing the code
            self.SetStartStop(Location, Aminos, PeptideStartAA, ORFStartNucleotide)
            self.SetFrame(Location, InfoBits[2])
            if (len(DBLocations) == 1):
                Location.Unique = 1 #default set to zero in the constructor, so we set to 1 if unique
            #append to list
            ReturnLocations.append(Location)
            
        return ReturnLocations
            
    def SetFrame(self, Location, FrameString):
        """Currently this only works for the single exon version.  I have not figured out
        how to work this for multi, so pretty small method
        """
        Frame = int(FrameString.replace("Frame", ""))
        Location.StartFrame = Frame
        Location.StopFrame = Frame
            
            
    def SetStartStop(self, Location, Aminos, PeptideStartAA, ORFStartNucleotide):
        """This is we start having some seemingly duplicate code.  It was important
        for me to separate the positive and negative strand stuff.  It just got really confusing
        WARNING: Just got single exon mapping working!!
        """
        PeptideStartNucleotide = -1
        PeptideStopNucleotide = -1
        if (Location.Strand == "-"):
            (PeptideStartNucleotide, PeptideStopNucleotide) = self.GetCoordsRevStrand(Aminos, PeptideStartAA, ORFStartNucleotide)
        else:
            #now the position.  We first get the protein start nuc, and then offset to the peptide start
            PeptideStartOffset = PeptideStartAA * 3 # because StartLocation is the start in amino acid space
            PeptideStartNucleotide = PeptideStartOffset + ORFStartNucleotide
            PeptideStopNucleotide = PeptideStartNucleotide + (len(Aminos) * 3) - 1 # we do a minus one because the bases are inclusive
            
        if (PeptideStartNucleotide < 0) or (PeptideStopNucleotide < 0):
            ##SNAFU!!!!!!
            print "ERROR: PeptideMapper:GetStartStop"
            print "##### Can't find peptide start or stop for %s"%self.CurrentAminos
            return
        Location.StartNucleotide = PeptideStartNucleotide
        Location.StopNucleotide = PeptideStopNucleotide

    def GetCoordsRevStrand(self, Aminos, PeptideStartAA, ORFStartNucleotide):
        """Separate method because doing math on the reverse strand is difficult.
        For the visually minded:  codons on the reverse listed below
        1   4   7   10  13  16  19
        AAA TTT CCC GGG AAA TTT CCC
        TTT AAA GGG CCC TTT AAA GGG AGT
         F   K   G   P   F   K   G   *
        Peptide GKFPG starts at position 21 (last C) and includes all bases up to 7, so we should do 
        start = 7, stop = 21
        Peptide KFPG starts at position 18.  So start protein (21) minus AA offset (3)
        NOTE: this gets the start<stop version, not the 5' and 3' version.
        """
        PeptideStartOffset = PeptideStartAA *3
        PeptideStartNucleotide = ORFStartNucleotide - PeptideStartOffset
        PeptideStopNucleotide = PeptideStartNucleotide - (len(Aminos)*3) + 1
        #now because this is on the reverse, and we requires start< stop, we send
        #them back in the reverse order
        return (PeptideStopNucleotide, PeptideStartNucleotide)

    def ClusterPeptidesOnGenome(self, PeptideList, MaxInterpeptideLen):
        """Given some set of peptides which are GenomicLocation objects (could be in a ORF, or random 
        group, or whatever), we find those which are clustered together (according to the criteria).  
        We return the the peptide list, dictionary of [clusterpos]->num peps in cluster.
        NOTE MaxInterpeptideLen is given in nucleotides
        """