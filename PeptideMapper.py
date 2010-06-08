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

import PGPeptide
import bioseq


class PeptideMappingClass:
    def __init__(self):
        self.DatabasePaths = [] #possibly multiple
        self.CurrentAminos = ""
        self.UniquePeptideCount =0

    def LoadDatabases(self, DBPaths):
        """
        Parameters: a list of paths to databases.  these should be 6 frame translations
        Return: none
        Description: load up databases in preparation for searching
        """
        self.DatabasePaths = DBPaths
        #self.orfIndex = bioseq.TrieIndexSeqs( DBPaths )
        self.orfIndex = bioseq.QGramIndex( DBPaths )
        self.orfIndex.index()

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
        LocationsInORFDB = self.orfIndex.accessionsWhereSeqFound(Aminos)
        if (len(LocationsInORFDB) < 1) and WarnNoMatch:
            #sometimes we don't care that there's no match.  Like for trypsin.  it's 
            #not part of our 6frame bacteria, but it was in the inspect search
            print "Peptide %s was not found in the database"%Aminos
        #2. now go through the process of converting protein space to nucleotide space    
        for (ORFFastaLine, ORFID, PeptideStartAA) in LocationsInORFDB:
            ## PeptideStartAA is the start position within the amino acid sequence
            #1. parse out the features of the protein name (fasta header)
            ParsedORFInfo = PGPeptide.ORFFastaHeader(ORFFastaLine)
            SimpleLocation = self.MapNucleotideLocation(ParsedORFInfo, PeptideStartAA, len(Aminos))
            #now that we have a Location, let's get our Located Peptide Object up and running
            Peptide = PGPeptide.LocatedPeptide(Aminos, SimpleLocation)
            Peptide.bestScore = PValue
            Peptide.ORFName = ParsedORFInfo.ORFName
            Peptide.name = "Peptide%s"%self.UniquePeptideCount
            self.UniquePeptideCount += 1
            #now we check the letter before us to see if it's tryptic
            ORFSequence = self.orfIndex.seqs[ORFID]
            PrefixOfPeptide = ORFSequence[PeptideStartAA -1]
            PeptideEndAA = PeptideStartAA + len(Aminos)
            if PeptideEndAA == len(ORFSequence):
                #this is at the end, so it is really c-term tryptic, even though there
                #may not be an R, K
                Peptide.SetTryptic(PrefixOfPeptide, 1) #1 is for overwrite c term
            else:
                Peptide.SetTryptic(PrefixOfPeptide)
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
                print "WARNING: protein %s does not map to any ORF. No LocatedProtein created"%ProteinName
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
        SimpleLocation = PGPeptide.GenomicLocation.FromHeader(
            ParsedORFInfo,
            LenAminos,
            offsetInAA=StartInProteinSpace
        )
        return SimpleLocation
