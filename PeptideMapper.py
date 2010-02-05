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


import SelectProteins
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
        self.ORFDB = SelectProteins.ProteinSelector()
        self.ORFDB.LoadMultipleDB(self.DatabasePaths)

    def MapPeptide(self, Aminos, PValue, WarnNoMatch = 0):
        """
        Parameters: an amino acid string, the best score (pvalue)
        Return: a list of LocatedPeptide objects (possible list of len 1)
        Description: This is the method that takes amino acids and maps
        them into dna sequence space.  I have tried very hard to label
        things as ORF___ if they are from open reading frames, and protein___
        if they are on actual predicted proteins.
        """
        ReturnList = [] # the list of PGPeptide.LocatedPeptide objects
        self.CurrentAminos = Aminos #only used for printing in case or error
        #1. First find the location(s) of the aminos in the ORF database
        LocationsInORFDB = self.ORFDB.FindPeptideLocations(Aminos)
        if (len(LocationsInORFDB) < 1) and WarnNoMatch:
            #sometimes we don't care that there's no match.  Like for trypsin.  it's 
            #not part of our 6frame bacteria, but it was in the inspect search
            print "Peptide %s was not found in the database"%Aminos
        #2. now go through the process of converting protein space to nucleotide space    
        for (ORFID, PeptideStartAA) in LocationsInORFDB:
            ## PeptideStartAA is the start position within the amino acid sequence
            #1. parse out the features of the protein name (fasta header)
            ORFFastaLine = self.ORFDB.ProteinNames[ORFID]
            ParsedORFInfo = PGPeptide.ORFFastaHeader(ORFFastaLine)
            SimpleLocation = self.MapNucleotideLocation(ParsedORFInfo, PeptideStartAA, len(Aminos))
            #now that we have a Location, let's get our Located Peptide Object up and running
            Peptide = PGPeptide.LocatedPeptide(Aminos, SimpleLocation)
            Peptide.bestScore = PValue
            Peptide.ORFName = ParsedORFInfo.ORFName
            #now we check the letter before us to see if it's tryptic
            ORFSequence = self.ORFDB.ProteinSequences[ORFID]
            PrefixOfPeptide = ORFSequence[PeptideStartAA -1]
            if PrefixOfPeptide in ["R", "K"]:
                Peptide.TrypticNTerm = 1
            if Aminos[-1] in ["R", "K"]:
                Peptide.TrypticCTerm = 1
            if len(LocationsInORFDB) == 1:
                Peptide.isUnique = 1

            #append to list
            ReturnList.append(Peptide)
            
        return ReturnList


    def MapProtein(self, Sequence, ProteinName, WarnNonSingleMatch = 1):
        """
        Parameters: the protein sequence, name of the protein (PROTEIN, not ORF)
        Return: LocatedProtein object
        Description: this method takes a protein sequence, and puts it onto the
        DNA within the bounds of a single Open reading frame.
        
        WARNING: not meant for any mutli-orf proteins, like ribosomal frame shifts
        or self splicing introns.
        
        NOTE: because of the alternate start sites, we have a bit of a problem, in that most proteins
        don't map to the six frame translation, specifically the initial M.  See the following example
            RDVLNRVMYYIILARFINYRLISLSCRSKRMRIFQGVVCGMALFLA  (six frame translation)
                  MMYYIILARFINYRLISLSCRSKRMRIFQGVVCGMALFLA  (protein sequence)
        to remedy this, we are going to axe off the initial M, and then kludge back
        the start site of the mapping

        """
        SearchableSequence = Sequence[1:] # see the problem in the NOTE
        self.CurrentAminos = Sequence #for error printing
        LocationsInORFDB = self.ORFDB.FindPeptideLocations(SearchableSequence)
        if len(LocationsInORFDB) == 0:
            if WarnNonSingleMatch:
                print "WARNING: protein %s does not map to any ORF"%ProteinName
            return None #can't find it, why waste time
        #now here we wonder what to do with proteins that map to multiple ORFs.  I think today
        # that I will just say with an iron fist, that these suck, and should be treated badly
        #we will take the first location. HACK
        (ORFID, ProteinStartAA) = LocationsInORFDB[0]
        ORFFastaLine = self.ORFDB.ProteinNames[ORFID]
        ParsedORFInfo = PGPeptide.ORFFastaHeader(ORFFastaLine)
        SimpleLocation = self.MapNucleotideLocation(ParsedORFInfo, ProteinStartAA, len(SearchableSequence))
        #now we make the protein
        Protein = PGPeptide.LocatedProtein(SimpleLocation)
        Protein.name = ProteinName
        Protein.ORFName = ParsedORFInfo.ORFName
        Protein.AddOneAminoAcidFivePrime() #adding back what we took off because of the NOTE
        Protein.AddStopCodon()
        return Protein
        
        


    def MapNucleotideLocation(self, ParsedORFInfo, StartInProteinSpace, LenAminos):
        """
        Parameters: the fasta line from a 6frame ORF, the index of the match start in protein space, the len of match
        Return: a PGPeptide.GenomicLocation object
        Description: we bust out the mapping of a sequence on to the dna in a single method, because
        both proteins and peptides want to use this method, and they need other stuff different.
        so shared code goes in one place, no?
        """
        (Start, Stop) = self.GetStartStop(LenAminos, StartInProteinSpace, ParsedORFInfo.Start, ParsedORFInfo.Strand)
        SimpleLocation = PGPeptide.GenomicLocation(Start, Stop, ParsedORFInfo.Strand)
        SimpleLocation.chromosome = ParsedORFInfo.Chromosome
        SimpleLocation.frame = ParsedORFInfo.Frame
        return SimpleLocation


    def GetStartStop(self, LenAminos, PeptideStartAA, ORFStartNucleotide, Strand):
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
            (PeptideStartNucleotide, PeptideStopNucleotide) = self.GetCoordsRevStrand(LenAminos, PeptideStartAA, ORFStartNucleotide)
        else:
            #now the position.  We first get the protein start nuc, and then offset to the peptide start
            PeptideStartOffset = PeptideStartAA * 3 # because StartLocation is the start in amino acid space
            PeptideStartNucleotide = PeptideStartOffset + ORFStartNucleotide
            PeptideStopNucleotide = PeptideStartNucleotide + (LenAminos * 3) - 1 # we do a minus one because the bases are inclusive
            
        if (PeptideStartNucleotide < 0) or (PeptideStopNucleotide < 0):
            ##SNAFU!!!!!!
            print "ERROR: PeptideMapper:GetStartStop"
            print "##### Can't find peptide start or stop for %s"%self.CurrentAminos
            return
        return(PeptideStartNucleotide, PeptideStopNucleotide) 

    def GetCoordsRevStrand(self, LenAminos, PeptideStartAA, ORFStartNucleotide):
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
        PeptideStopNucleotide = PeptideStartNucleotide - (LenAminos*3) + 1
        #now because this is on the reverse, and we requires start< stop, we send
        #them back in the reverse order
        return (PeptideStopNucleotide, PeptideStartNucleotide)

