"""PeptideMapper.py
This auxiliary set of classes/functions is to be used in mapping peptides to 
their genomic location.  hopefully, this will be the central place that we 
do all of this activity, regardless of the genome.

Usage: You can call the following two methods.  nothing else.
LoadDatabases(DBPaths)
MapMe(AminoAcidSequence)

NOTE: because this is genomic mapping, the database must have the standard
notion of a genomic context. you can get this by using the SixFrameFasta.py
script to create your database.

NOTE: this is a utility, and not executable from the command line

NOTE: this is strictly for bacterial peptides that map unspliced.  There
is a separate program eukPeptideMapper that deals with the headache of splicing
"""

import sys
import os
import traceback
import SelectProteins
import GenomicLocations
import PGPeptide


class PeptideMappingClass:
    def __init__(self):
        self.DatabasePaths = [] #possibly multiple
        self.CurrentAminos = ""
        
    def LoadDatabases(self, DBPaths):
        """
        Parameters: a list of paths to databases.  these should be 6 frame translations
        Return: none
        Description: load up databases in preparation for searching
        """
        self.DatabasePaths = DBPaths
        self.ProteinPicker = SelectProteins.ProteinSelector()
        self.ProteinPicker.LoadMultipleDB(self.DatabasePaths)

    def MapMe(self, Aminos, PValue, WarnNoMatch = 0):
        """
        Parameters: an amino acid string, the best score (pvalue)
        Return: a list of LocatedPeptideObjects (possible list of len 1)
        Description: This is the method that takes amino acids and maps
        them into dna sequence space.
        NOTE: This method parses information out of the protein header 
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
        ReturnList = [] # the list of PGPeptide.LocatedPeptide objects
        self.CurrentAminos = Aminos #only used for printing in case or error
        #1. First find the location(s) of the aminos in the ORF database
        ProteinLocations = self.ProteinPicker.FindPeptideLocations(Aminos)
        if (len(ProteinLocations) < 1) and WarnNoMatch:
            #sometimes we don't care that there's no match.  Like for trypsin.  it's 
            #not part of our 6frame bacteria, but it was in the inspect search
            print "Peptide %s was not found in the database"%Aminos
        #2. now go through the process of converting protein space to nucleotide space    
        for (ProteinID, PeptideStartAA) in ProteinLocations:
            ## PeptideStartAA is the start position within the amino acid sequence
            #1. parse out the features of the protein name (fasta header)
            ProteinName = self.ProteinPicker.ProteinNames[ProteinID]
            ### CRAP ! the fame proteins are XXX.Protein1, so the . screws up the splitting
            ## Total hack below.  Find a better splitter than . and redo the sixframefasta.py
            if ProteinName.find("XXX.") == 0:
                ProteinName = ProteinName.replace("XXX.", "XXX")
            InfoBits = ProteinName.split(".")
            
            ProteinName = InfoBits[0]
            Chromosome = InfoBits[1].replace("Chr:", "")
            Strand = InfoBits[4].replace("Strand", "")
            Frame = int(InfoBits[2].replace("Frame", ""))
            
            ORFStartNucleotide = int(InfoBits[3].replace("StartNuc", "")) # start of the open reading frame, not my peptide
            (Start, Stop) = self.GetStartStop(Aminos, PeptideStartAA, ORFStartNucleotide, Strand)
            SimpleLocation = PGPeptide.GenomicLocation(Start, Stop, Strand)
            SimpleLocation.chromosome = Chromosome
            SimpleLocation.frame = Frame
            
            #now that we have a Location, let's get our Located Peptide Object up and running
            Peptide = PGPeptide.LocatedPeptide(Aminos, SimpleLocation)
            Peptide.bestScore = PValue
            #now we check the letter before us to see if it's tryptic
            ORFSequence = self.ProteinPicker.ProteinSequences[ProteinID]
            PrefixOfPeptide = ORFSequence[PeptideStartAA -1]
            if PrefixOfPeptide in ["R", "K"]:
                Peptide.TrypticNTerm = 1
            if Aminos[-1] in ["R", "K"]:
                Peptide.TrypticCTerm = 1
            if len(ProteinLocations) == 1:
                Peptide.isUnique = 1

            #append to list
            ReturnList.append(Peptide)
            
        return ReturnList


    def GetStartStop(self, Aminos, PeptideStartAA, ORFStartNucleotide, Strand):
        """
        Parameters: amino acid string, index of aminos within protein, index of ORF within dna
        Return: start and stop of aminos on the dna
        Description: Using the supplied offsets, we get nucleotide coordinates for the amino
        acid sequence.  It is important to reiterate, that start < stop.  ALWAYS.  The words start
        and stop do not in any way refer to 5' or 3'.
        """
        PeptideStartNucleotide = -1 #set as impossible values for sanity check below
        PeptideStopNucleotide = -1
        if (Strand == "-"):
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
        return(PeptideStartNucleotide, PeptideStopNucleotide) 


            
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

