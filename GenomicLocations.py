"""GenomicLocations.py
This auxiliary set of classes/functions is to be used in dealing
with genomic locations of any type

NOTE: because this is genomic mapping, the database must have the standard
notion of a genomic context. you can get this by using the SixFrameFasta.py
script to create your database.

NOTE: this is a utility, and not executable from the command line

"""

import sys
import os
import traceback
import GFFIO

class GenomicLocationClass:
    """The most basic genomic location
    """
    def __init__(self):
        self.Chromosome = "UNDEFINED"
        self.StartNucleotide = -1 
        self.StopNucleotide = -1
        self.Strand =  "UNDEFINED" #string '-' or '+'
        
            
class GenomicLocationForORF(GenomicLocationClass):
    """This is for a genomic open reading frame, not a generically named protein
    the startnucloetide will not corresponde to an 'a' of 'atg'.  
    NOTE: it is assumed that a single ORF corresponds to a single protein. Don't abuse that
    Note: ORFs only have one frame (StartFrame).  Don't try to access StopFrame.  You would be stupid.
    NOTE: self.StopNucleotide is upto and INCLUDING the stop codon
    """
    def __init__(self, FastaLine, AALen):
        self.PeptideLocationList = [] # list of GenomicLocationForPeptide objects
        self.ProteinPrediction = None # GenomicLocationForPeptide object
        self.ProteinPredictionName = None #simply the name given in the genome annotation
        GenomicLocationClass.__init__(self)
        self.ORFName = None
        self.AminoAcidLength = AALen #handy to keep around.
        self.ParseFastaHeader(FastaLine)
        self.MisPrecition = None # types include 'MissingProtein' or 'PeptidesBeforeStart' or 'PeptidesAfterStop'

    def AmIMispredicted(self, MinPeptides = 2):
        """This is the main entry point for proteogenomic analysis of mispredicted proteinss.
        This includes: No prediction but supporting peptides, and peptides which fall outside 
        of predicted boundaries. 
        NOTE: return values.  zero means we have no evidence supporting a new prediction
        """
        MinORFLen = 20 #minimum number of amino acids in an orf before we pay attention
        if len(self.PeptideLocationList) < MinPeptides:
            return 0
        #so we have some peptides
        if not self.ProteinPrediction:
            #we have no protein, but we do have peptides.  Let's see if this is a viable ORF
            if self.AminoAcidLength >= MinORFLen:
                print "######UNPREDICTED ORF###########"
                self.PrintMe(1,0)
                print "\n\n"
                self.MisPrediction = "MissingProtein"
                return 1
        else:
            Miscalled = self.ArePeptidesOutsideMyProtein(MinPeptides)
            if Miscalled:
                print "############MisPredicted ORF################"
                self.PrintMe(0,1)
                print "\n\n"
            return Miscalled
        return 0
    
    def ArePeptidesOutsideMyProtein(self, MinPeptides):
        """Looking to see if the peptides map to within my predicted protein boundaries
        Because this is a simple ORF, we assume no splicing.  Because start and stop are numerical order
        (and NOT 5' or 3'), we can actually do a very simple boundary test for peps before start or after stop
        
                   PROTEINSEQUENCEISTHISLINEANDLETTERS
             BEFOREPROTEIN                     LETTERSAFTER
        """
        #print "GenomicLocationForORF::ArePeptidesOutsideMyProtein"
        #self.PrintMe(0,1)
        #print "\nGoing through all of the peptides"
        for PepLocation in self.PeptideLocationList:
            if PepLocation.StartNucleotide < self.ProteinPrediction.StartNucleotide:
                self.MisPrecition = "PeptidesBeforeStart"
                print "I determined that this peptide is before the start"
                PepLocation.PrintMe(1)
                print "\n\n"
                return 1
            if PepLocation.StopNucleotide > self.ProteinPrediction.StopNucleotide:
                self.MisPrecition = "PeptidesAfterStop"
                print "I determined that this peptide is after the stop"
                PepLocation.PrintMe(1)
                print "\n\n"
                return 1
        return 0
            
    def PrintMe(self, PrintPeptides = 0, PrintProtein = 0):
        """pretty printer"""
        print "GenomicLocationForORF object %s"%self.ORFName
        print "Chr %s, Start %s"%(self.Chromosome, self.StartNucleotide)
        if PrintPeptides:
            print "Printing Peptides"
            for Peptide in self.PeptideLocationList:
                Peptide.PrintMe()
        if PrintProtein:
            print "Protein within ORF %s"%self.ProteinPredictionName
            self.ProteinPrediction.PrintMe()
        
    def ParseFastaHeader(self, FastaLine):
        """This takes the fasta line (from sixframefasta.py) and parses out the start and stop
        of the orf.  This may look a lot like some code from PeptideMappingClass.  that's because
        I copied it from there.  it was close, but could not be munged, or abstracted. sadly
        now we have two copies.

        Protein256370.Chr:Chr1.Frame1.StartNuc831372.Strand-
        XXX.Protein157841.Chr:Chr1.Frame1.StartNuc4033167.Strand+
        InfoBits[0] = protein unique id
        InfoBits[1] = chromosome
        InfoBits[2] = Frame
        InfoBits[3] = StartNucleotide of the ORF
        InfoBits[4] = Strand
        
        """
        if FastaLine.find("XXX.") == 0:
            FastaLine = FastaLine.replace("XXX.", "XXX")
        InfoBits = FastaLine.split(".")
        #chromosome, strand
        self.ORFName = InfoBits[0]
        self.Chromosome = InfoBits[1].replace("Chr:", "")
        self.StartFrame = int(InfoBits[2].replace("Frame", ""))   #no stop frame, see class note
        self.Strand = InfoBits[4].replace("Strand", "")
        self.StartNucleotide = int(InfoBits[3].replace("StartNuc", ""))
        #now some math to figure out the stop nuc
        if (self.Strand == "-"):
            CodingStopNucleotide = self.StartNucleotide - (self.AminoAcidLength * 3) + 1
            self.StopNucleotide = CodingStopNucleotide - 3 #three bases of the stop codon
        else:
            #now the position.  We first get the protein start nuc, and then offset to the peptide start
            CodingStopNucleotide = self.StartNucleotide + (self.AminoAcidLength * 3) - 1 # we do a minus one because the bases are inclusive
            self.StopNucleotide = CodingStopNucleotide + 3 # three bases of the stop codon are INCLUDED
        
        
        


class GenomicLocationForPeptide(GenomicLocationClass):
    """This class is to store information about the genomic location of a peptide.
    NOTE: This is an inclusive nucleotide system.  Thus, if a di-amino peptide starts on
    nucleotide 4 of a sequence, it occupies bases 4,5,6,7,8,9.  Start =4, stop=9
    NOTE: Start < Stop.  Numerically.  If you want start to be 5' you should
    ask for the 5' with a method
    """
    def __init__(self):
        self.StartFrame = 0 #only positive numbers 1,2,3 even for the reverse strand
        self.StopFrame = 0
        self.ProteinName = "UNDEFINED"
        self.Aminos = "UNDEFINED"
        self.Score = -1 # this is meant to be the pvalue of a peptide, or other meaningful score.
        GenomicLocationClass.__init__(self)
        
    def GetFivePrimeNucelotide(self):
        """Just checks for the strand, and returns the 5' nucelotide
        """
        if self.Strand == "-":
            return self.StopNucleotide
        return self.StartNucleotide
    def GetThreePrimeNucleotide(self):
        """Just checks for the strand, and returns the 5' nucelotide
        """
        if self.Strand == "-":
            return self.StartNucleotide
        return self.StopNucleotide
        
    def PrintMe(self, Verbose = 0):
        print "GenomicLocationForPeptide object, %s"%self.Aminos
        print "found in protein %s"%self.ProteinName
        print "Chr:%s, Start %s, Stop %s"%(self.Chromosome, self.StartNucleotide, self.StopNucleotide)
        if Verbose:
            print "Strand %s, StartFrame %s, Stop Frame, %s"%(self.Strand, self.StartFrame, self.StopFrame)
            
    def GetGFF3Line(self, PValue, GlobalPeptideCount):
        """output the GFF3 line for this peptide"""
        """ DEPRECATE THIS SOON!!!!!!!!!!!!!"""
        Values = [0] * 9 #clean out before use
        Values[GFFColumns.Chr] = "%s"%self.Chromosome
        Values[GFFColumns.Source] = "Proteomics"
        Values[GFFColumns.Type] = "CDS"
        Values[GFFColumns.Phase] = "0" # all unspliced peptides are phase 0
        Values[GFFColumns.Score] = "%s"%PValue
        PeptideAttributeID = "Peptide%s"%GlobalPeptideCount
        Values[GFFColumns.Attributes] = "ID=%s;Name=%s"%(PeptideAttributeID, self.Aminos)
        Values[GFFColumns.Start] = "%s"%self.StartNucleotide
        Values[GFFColumns.Stop] = "%s"%self.StopNucleotide
        Values[GFFColumns.Strand] = "%s"%self.Strand
        GFF3FormatedLine = "\t".join(Values)
        return GFF3FormatedLine

    def GetGFF3LineNew(self, GlobalPeptideCount):
        """Utilizes the GFFIO functions to get the standard GFF line.  
        First we create the proper dictionary container and then send it out.
        """
        Dictionary = {}
        Dictionary["seqid"] = self.Chromosome
        Dictionary["source"] = "Proteomics"
        Dictionary["type"] = "CDS"  # hack we need to change this when we figure it out.  This should be polypeptide
        Dictionary["start"] = self.StartNucleotide
        Dictionary["end"] = self.StopNucleotide
        Dictionary["score"] = self.Score
        Dictionary["strand"] = self.Strand
        Dictionary["phase"] = 0 # all unspliced peptides are phase 0
        Dictionary["attributes"] = {}
        Dictionary["attributes"]["ID"] = "Peptide%s"%GlobalPeptideCount
        Dictionary["attributes"]["Name"] = self.Aminos
        Dictionary["attributes"]["Parent"] = self.ProteinName
        return GFFIO.MakeGFFLineFromDictionary(Dictionary)
    
    def FillFromGFF(self, Dictionary):
        """Parameters: dictionary of key-value pairs made when parsing gff files
        Return: none
        Description: take the dictionary and fill in values of myself.  You'll notice
        that there are some parts of the dictionary that I dont' care about.  that's life
        """    
        self.Chromosome = Dictionary["seqid"]
        self.StartNucleotide = Dictionary["start"]
        self.StopNucleotide = Dictionary["end"]
        self.Score = Dictionary["score"]
        self.Strand = Dictionary["strand"]
        self.Aminos = Dictionary["attributes"]["Name"] 
        self.ProteinName = Dictionary["attributes"]["Parent"]

    def AddOneAminoAcidFivePrime(self):
        """This is a very very ugly kludge.  to map proteins on to the ORFs, we remove the first amino acid
        because if it's an alternate start site, the ORF has V, and the protein has M.  So now we need to add
        the ability to add back that one amino acid to our position.  Remember, start and stop are numerical order
        and NOT associated with 5' or 3'.  
        """
        if self.Strand == "-":
            self.StopNucleotide += 3
        else:
            self.StartNucleotide -= 3
   
   
class GFFColumns:  # see big comment below, these are zero indexed.
    Chr = 0
    Source = 1 #"Proteomics" or "Inspect" or "6frame"
    Type = 2
    Start = 3
    Stop = 4
    Score = 5
    Strand = 6
    Phase = 7
    Attributes = 8
        
     