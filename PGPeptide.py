from Bio import SeqIO
import GFFIO
import bioseq


###############################################################################

class GenomicLocation(object):
    def __init__(self,start,stop,strand,chromosome=None):
        """Start and stop of a region on a sequence.
        Stop must be > then start.
        """
        if start > stop or stop < 0 or start < 0:
            raise ValueError("Start %d > Stop %d" % (start,stop))

        if not strand in ['+','-']:
            raise ValueError("strand must be + or -")

        self.__start = start
        self.__stop = stop
        self.__frame = None
        self.__strand = strand
        self.chromosome = chromosome

    @classmethod
    def FromHeader(self, ParsedHeader, AALength, offsetInAA=0, addStop=False):
        """
        Parameters: ORFFastaHeader object, length of the amino acids in the ORF
        addStop: is a boolean. On true it includes 1 extra codon at the 3' end.
        offsetInAA: offset location by this many peptides from the AA seq start

        Return: GenomicLocation instance
        Description: Make a GenomicLocation from the information in the header
        """
        FivePrime = ParsedHeader.Start
        ThreePrime = -1 #set below
        Strand = ParsedHeader.Strand
        Chrom = ParsedHeader.Chromosome
        Frame = ParsedHeader.Frame
        offset3Prime = 0
        if addStop:
            offset3Prime = 3

        # now some math to figure out the stop nuc
        # ie on minus strand subtract from start on + strand add
        if Strand == "-":
            '''For the visually minded:  codons on the reverse listed below
        1   4   7   10  13  16  19
        AAA TTT CCC GGG AAA TTT CCC
        TTT AAA GGG CCC TTT AAA GGG AGT
         F   K   G   P   F   K   G   *
        Peptide GKFPG starts at position 21 (last C) and includes all bases up
        to 7, so we should do start = 7, stop = 21
        Peptide KFPG starts at position 18. So start protein 21 minus AA offset 3
        '''
            FivePrime -= offsetInAA * 3
            CodingThreePrime = FivePrime - (AALength * 3) + 1
            ThreePrime = CodingThreePrime - offset3Prime
            if ThreePrime < 1:
                print "Warning setting negative start to 1 %d" % ThreePrime
                ThreePrime = 1
            SimpleLocation = GenomicLocation(ThreePrime, FivePrime, Strand, Chrom)
        else:
            # now the position.  We first get the protein start nuc, and then
            # offset to the peptide start
            # we do a minus one because the bases are inclusive
            FivePrime += offsetInAA * 3
            CodingThreePrime = FivePrime + (AALength * 3) - 1
            ThreePrime = CodingThreePrime + offset3Prime
            SimpleLocation = GenomicLocation(FivePrime, ThreePrime, Strand, Chrom)

        SimpleLocation.frame = Frame
        return SimpleLocation

    def __str__(self):
        """Parameters: None
        Return: the string version of this object
        Description: This method is used to stringify a GenomicLocation object
        that means that you can do
        print "%s"%GenomicLocationObject
        and get a meaningful print out
        """
        #return "%d,%d %s%d" % (self.start, self.stop, self.strand, self.frame)
        return "%s:%d,%d %s" % (self.chromosome, self.start, self.stop, self.strand)

    def __cmp__(self, other):
        """Parameters: Another GenomicLocation object
        Return: -1, 0, 1 following the cmp standards
        Description: this does a strict NUMERICAL sort of the objects.  If
        starts are equal, then return the smaller stop. By defining this
        function, we can use a simple call to sort() on a list of
        GenomicLocation objects and it will happen. This has nothing to
        do with 5' or 3' sorting. NOTHING.
        """
        Value = cmp(self.start, other.start)
        if Value == 0: #this is checked out just in case
            return cmp(self.stop, other.stop)
        return Value

    def SortFivePrime(self, other):
        """Parameters: antother GenomeLocation object
        Return: -1, 0, 1 following the cmp standards
        Description: This does a biological sort of the objects, returning
        such that a 5' peptide should appear first. It is designed to be

        NOTE: We ASSUME that the two locations are on the same strand,
        because otherwise 5' would be meaning less and you would be
        stupid for calling this function
        """
        #here's the switch on strandedness,
        if self.strand == "+":
            #this can just be a simple thing
            return self.__cmp__(other)
        #now here' the hard part
        #return negative if self is 5' of other
        if self.stop > other.stop:
            return -1
        if self.stop == other.stop:
            #they are equal, so let's check the 3' end. The convention will mirror
            #what a numerical sort does on + strand stuff.
            if self.start > other.start:
                return -1
            if self.start == other.start:
                return 0
        return 1

    def AddOneAminoAcidFivePrime(self):
        'This is accessed to add a single amino acid to the front of a protein.'
        #because this location was created, and initially passed the start < stop
        #we can be sure that this adjustment will not change that
        #start < stop ==> start -3 < stop.  Start only gets smaller
        if self.strand == "-":
            self.__stop +=  3
        else:
            self.__start -= 3

    def AddOneAminoAcidThreePrime(self):
        'This is accessed to add the stop codon for a protein'
        if self.strand == "-":
            self.__start -= 3
        else:
            self.__stop += 3

    def GetThreePrime(self):
        if self.strand == "+":
            return self.__stop
        return self.__start

    @property
    def start(self):
        'Lesser coordinate on the sequence.'
        return self.__start

    @property
    def stop(self):
        'Greater coordinate on the sequence.'
        return self.__stop

    @property
    def strand(self):
        'Strand of location + or -'
        return self.__strand

    @property
    def frame(self):
        'Frame of the translation.'
        return self.__frame

    @frame.setter
    def frame(self,value):
        if not value in [1,2,3]:
            raise ValueError("Frame must be 1, 2 or 3.")
        self.__frame = value

    def overlap(self,otherLocation):
        """Returns None if the locations don't overlap.
        Otherwise, it returns a tuple of the overlapping region.
        """
        (myStart,myStop) = (self.start, self.stop)
        (otherStart,otherStop) = (otherLocation.start, otherLocation.stop)
        if myStop < otherStart or myStart > otherStop:
            return None
        elif myStart <= otherStart:
            if otherStop < myStop:
                # other is contained sub sequenc
                return (otherStart,otherStop)
            else:
                return (otherStart, myStop)
        else: # mSta otherStart
            if myStop < otherStop:
            # self is contained sub sequence
                return (myStart,myStop)
            else:
                return (myStart, otherStop)


###############################################################################

class LocatedPeptide(object):
    def __init__(self, Aminos, location=None):
        self.location = location
        self.name = None
        self.aminos = Aminos
        self.bestScore = None
        self.spectrumCount = 0
        self.isUnique = None
        self.TrypticNTerm = None
        self.TrypticCTerm = None
        self.ORFName = None #this is actually the ORF Name.  as we put them into ORFs

    @property
    def chromosome(self):
        return self.location.chromosome

    def SetTryptic(self, Prefix):
        """Parmeters: the prefix of the peptide (letter immediately before)
        Return: None
        Description: This function sets the two tryptic variables, n and c term
        """
        TrypticLetters = ["R", "K"]
        if Prefix in TrypticLetters:
            self.TrypticNTerm = 1
        if self.aminos[-1] in TrypticLetters:
            self.TrypticCTerm = 1


    def isTryptic(self):
        """Parameters: None
        Return: true/falst on whether a peptide is tryptic
        Description: current implementation is that any tryptic
        endpoint will do.
        """
        if self.TrypticCTerm:
            return True
        if self.TrypticNTerm:
            return True
        return False

    def GetFivePrimeNucleotide(self):
        """Just checks for the strand, and returns the 5' nucelotide
        """
        if self.location.strand == "-":
            return self.location.stop
        return self.location.start

    def GetThreePrimeNucleotide(self):
        """Just checks for the strand, and returns the 5' nucelotide
        """
        if self.location.strand == "-":
            return self.location.start
        return self.location.stop

    def __str__(self):
        """Parameters: None
        Return: the string version of this object
        Description: This method is used to stringify a GenomicLocation object
        that means that you can do
        print "%s"%GenomicLocationObject
        and get a meaningful print out
        """
        return "%s in %s, unique=%s, %s"%(self.aminos, self.ORFName, self.isUnique, self.location)

    def __cmp__(self, other):
        """Parameters: Another GenomicLocation object
        Return: -1, 0, 1 following the cmp standards
        Description: this does a strict NUMERICAL sort of the objects.  If
        starts are equal, then return the smaller stop. By defining this
        function, we can use a simple call to sort() on a list of
        GenomicLocation objects and it will happen. This has nothing to
        do with 5' or 3' sorting. NOTHING.
        """
        return cmp(self.location, other.location)

    def GetStart(self):
        return self.location.start

    def GetStop(self):
        return self.location.stop

    def Strand(self):
        return self.location.strand


###############################################################################

class LocatedProtein(object):
    def __init__(self,location):
        self.location = location
        self.name = None #like asparagine synthase
        self.ORFName = None #like Protein1233

    def AddOneAminoAcidFivePrime(self):
        """This is a very very ugly kludge.  to map proteins on to the ORFs, we remove the first amino acid
        because if it's an alternate start site, the ORF has V, and the protein has M.  So now we need to add
        the ability to add back that one amino acid to our position.  Remember, start and stop are numerical order
        and NOT associated with 5' or 3'.
        """
        self.location.AddOneAminoAcidFivePrime()

    def AddStopCodon(self):
        """Add the three nucs to the location object representing my stop codon, which is part of me
        as far as NCBI is concerned
        """
        self.location.AddOneAminoAcidThreePrime()

    def GetFivePrimeNucleotide(self):
        """Return the 5' nucleotide"""
        if self.location.strand == "-":
            return self.location.stop
        return self.location.start

    def GetStart(self):
        return self.location.start

    def GetStop(self):
        return self.location.stop
    def GetORFName(self):
        return self.ORFName

    def __str__(self):
        return "%s in %s, %s"%(self.name, self.ORFName, self.location)


###############################################################################

class OpenReadingFrame(object):
    def __init__(self, FastaHeader=None, AASequence=None, name=None):
        self.location = None
        self.__peptides = [] #LocatedPeptide objects
        self.annotatedProtein = None # a LocatedProtein Object, can't have more than one. If there's a fight, longest one wins
        self.name = None #the unique identifier of the open reading frame, e.g. Protein12345
        self.aaseq = AASequence
        self.naseq = None
        if FastaHeader:
            ParsedHeader = ORFFastaHeader(FastaHeader)
            self.name = ParsedHeader.ORFName
            self.location = GenomicLocation.FromHeader(ParsedHeader,
                                                       len(AASequence),
                                                       addStop=True)
        if name:
            self.name = name

        self.GCWholeORF = None
        self.GCPredictedProtein = None
        self.GCObservedRegion = None
        self.CDS = None # biopython SeqFeature for the CDS record from a genbank file

    def __str__(self):
        NumUniquePeptides = 0
        for Peptide in self.peptideIter():
            if Peptide.isUnique:
                NumUniquePeptides += 1
        return "%s as %s, %s (%s) peptides, %s"%(self.name, self.annotatedProtein, self.numPeptides(), NumUniquePeptides, self.location)

    @property
    def chromosome(self):
        return self.location.chromosome

    @chromosome.setter
    def chromosome(self,chrom):
        self.location.chromosome = chrom

    def GetTranslation(self):
        """Parameters: none
        Return: Amino acid translation of the ORF
        Description: gets the translation of the entire open reading frame
        """
        return self.aaseq

    def GetStrand(self):
        return self.location.strand

    def GetObservedDNACoords(self):
        """Parameters: NOne
        Return: tuple (start, stop) where start<stop
        Description: get the DNA coordinates of the observed sequence,
        from the first peptide to the stop codon
        """
        FirstPeptide = self.GetFivePrimePeptide() #note here that we are not using the uniqueness filter
        FivePrime = FirstPeptide.GetFivePrimeNucleotide()
        ThreePrime = self.location.GetThreePrime()
        if self.location.strand == "+":
            return (FivePrime, ThreePrime)
        return (ThreePrime, FivePrime)

    def GetObservedSequence(self):
        """Parameters: none
        Return: amino acid sequence
        Description: We want to get the sequence that we have evidence as
        being translated.  This is from the first peptide through the stop
        """
        #1. get the first peptide, in order
        FirstPeptide = self.GetFivePrimePeptide() #note here that we are not using the Uniqueness filter
        #2. map it to the aaseq
        FirstResidue = self.aaseq.find(FirstPeptide.aminos)
        #3. return the slice
        return self.aaseq[FirstResidue:]

    def GetFivePrimePeptide(self, Unique = 0):
        """Parameters: None required.  Optional Unique parameter
        Return: the five prime most peptide
        Description: cycle through the peptides and return the one closest
        to the start. If you specify Unique = 1, then you are specifying
        that you want the 5' most UNIQUELY MAPPING peptide.  Non-uniquely
        mapping peptides will be passed over
        """
        #1. Sort the list of peptides
        #self.__peptides.sort(LocatedPeptide.SortScore)
        self.__peptides.sort(cmp=lambda x,y: x.location.SortFivePrime(y.location))
        #2. walk through and find the first one that is unique
        for Peptide in self.peptideIter():
            if not Unique: #you don't care about uniqueness
                return Peptide #just return right off the bat
            if Peptide.isUnique:
                return Peptide

    def GetNucleotideStartOfTranslation(self):
        """Returns the start of translation"""
        return self.annotatedProtein.GetFivePrimeNucleotide()

    def numPeptides(self):
        'Returns the number of peptides in the ORF.'
        return len(self.__peptides)

    def addLocatedProtein(self, Protein):
        self.annotatedProtein = Protein

    def GetLocatedProtein(self):
        if self.annotatedProtein:
            return self.annotatedProtein
        return None

    def addLocatedPeptide(self, Peptide):
        'Adds a single LocatedPeptide objects to the ORF.'
        if self.chromosome != Peptide.location.chromosome:
            raise ValueError("Adding peptide %s to ORF %s" % (Peptide, self))
        self.__peptides.append( Peptide )


    def addLocatedPeptides(self,peptideList):
        'Adds a list of LocatedPeptide objects to the ORF.'
        self.__peptides.extend( peptideList )

    def addRawPeptides(self,peptideList):
        'Adds a list of peptide sequences without location.'
        for pep in peptideList:
            pepObj = LocatedPeptide(pep)
            self.__peptides.append( pepObj )

    def filterPeptides( self, filterFunc ):
        'Takes a filter function returning True or False for a LocatedPeptide object'
        self.__peptides = filter( filterFunc, self.__peptides)

    def peptideIter( self ):
        'An iterator for the peptides. Usage: for pep in orf.peptideIter():'
        for pep in self.__peptides:
            yield pep


###############################################################################

class ORFFastaHeader(object):
    """This tiny class just holds and deals with the information encoded in an ORF fasta header
    as created by the SixFrameTranslation script.

    WARNING: Don't include a period '.' in your name. we use that to split

    Protein256370.Chr:Chr1.Frame1.StartNuc831372.Strand-
    XXX.Protein157841.Chr:Chr1.Frame1.StartNuc4033167.Strand+
    InfoBits[0] = protein unique id
    InfoBits[1] = chromosome
    InfoBits[2] = Frame
    InfoBits[3] = StartNucleotide of the ORF
    InfoBits[4] = Strand

    """
    def __init__(self, Header):
        self.Chromosome = None
        self.ORFName = None
        self.Start = None
        self.Frame = None
        self.Strand = None
        self.String = Header
        self.Parse(Header)

    def __str__(self):
        return '%s.Chr:%s.Frame%d.StartNuc%d.Strand%s' % (
            self.ORFName, self.Chromosome, self.Frame,self.Start,self.Strand )

    def Parse(self, String):
        """Parameters: the fasta header line
        Return: none
        Description: parse it out into components
        """
        ### CRAP ! the fame proteins are XXX.Protein1, so the . screws up the splitting
        ## Total hack below.  Find a better splitter than . and redo the sixframefasta.py
        if String.find("XXX.") == 0:
            String = String.replace("XXX.", "XXX")
        InfoBits = String.split(".")

        self.ORFName = InfoBits[0]
        self.Chromosome = InfoBits[1].replace("Chr:", "")
        self.Strand = InfoBits[4].replace("Strand", "")
        self.Frame = int(InfoBits[2].replace("Frame", ""))
        self.Start = int(InfoBits[3].replace("StartNuc", "")) # start of the open reading frame, not my peptide


###############################################################################

class GFFPeptide(GFFIO.File):
    'This class reads/writes GFF stuff that pertains specifically to our PGPeptide Implementation.'

    ### Inherits the constructor of the GFFIO.File ###

    def generateORFs(self, sequenceFile, genome):
        '''Parameters: A sequence file supported by SequenceIO, a Genome() object
        to populate with OpenReadingFrame objects and their LocatedPeptides,
        and optionally a class for parsing the sequence accession.
        Description: Reads the peptides from the GFF, and the ORFs from the sequence file
        '''
        # Read in the peptides from the GFF file, creating ORFs as needed
        for gffRec in self:
            protein = gffRec.attributes['Parent']
            chrom = genome.chromosomes[ gffRec.seqid ]
            orf = chrom.getOrf( protein )
            if not orf:
                # ORF doesn't exist in chromosome, so add as an pepOnlyOrf
                orf = OpenReadingFrame(name=protein)
                # Real location is read in from 6frame file in addSeqToPepOnlyOrfs
                orf.location = GenomicLocation(0,0,'+',gffRec.seqid)
                chrom.addOrf( orf, 'PepOnly' )

            location = GenomicLocation(gffRec.start, gffRec.end, gffRec.strand, gffRec.seqid)
            # Peptide is encoded as the name, since it's generally short
            peptide = LocatedPeptide( gffRec.attributes['Name'], location)
            peptide.name = gffRec.attributes['ID']
            peptide.bestScore = gffRec.score
            orf.addLocatedPeptide( peptide )

        genome.addSeqToPepOnlyOrfs( sequenceFile )

    def writeORFPeptides(self, orf):
        gffRec = GFFIO.Record() #create empty record.  add values below then write each one
        gffRec.source = 'Proteomics'
        gffRec.type = 'polypeptide'
        gffRec.seqid = orf.chromosome
        gffRec.attributes['Parent'] = orf.name

        for peptide in orf.peptideIter():
            gffRec.start = peptide.GetStart()
            gffRec.end   = peptide.GetStop()
            gffRec.score = peptide.bestScore
            gffRec.strand= peptide.Strand()
            gffRec.attributes['Name'] = peptide.aminos
            gffRec.attributes['ID'] = peptide.name
            self.write( gffRec )


###############################################################################

class Chromosome(object):
    """A collection of ORFs and annotated proteins across a single large NA sequence.
    Currently, ORFs are stored in three separate collections, simpleOrfs,
    complexOrfs and pepOnlyOrfs.
    Simple ORFs correspond directly to a contiguous annotated protein.
    Complex ORFs are noncontiguous,
    pepOnly ORFs are supported only by peptides.
    """
    def __init__(self,accession=None,seq=None):
        self.accession = accession
        self.sequence = seq  # NA sequence of chromosome, currently a biopython Seq object
        self.simpleOrfs = {} # ORFs corresponding to a contiguous annotated protein
        self.complexOrfs= {} # complex ORFs  have no direct end to end protein map
        self.pepOnlyOrfs= {} # ORFs created only via peptides
        self.endToCDS   = {} # maps the 3' end of an protein to its SeqFeature object

    def addOrf(self,orf,orfType):
        orfHash = self.simpleOrfs
        if orfType == 'PepOnly':
            orfHash = self.pepOnlyOrfs
        elif orfType == 'Complex':
            orfHash = self.complexOrfs

        if orfHash.has_key( orf.name ):
            raise KeyError("Duplicate orf %s in chrom %s"%(orf.name,self.accession))
        else:
            orfHash[ orf.name ] = orf

    def getOrf(self,protName):
        "Given an ORFName returns the orf regardless of simple or complex membership"
        if self.simpleOrfs.has_key( protName ):
            return self.simpleOrfs[ protName ]
        elif self.complexOrfs.has_key( protName ):
            return self.complexOrfs[ protName ]
        elif self.pepOnlyOrfs.has_key( protName ):
            return self.pepOnlyOrfs[ protName ]
        else:
            return None

    def numORFsWithPeptides(self):
        i = len(self.pepOnlyOrfs)
        for orf in self.simpleOrfs.values() + self.complexOrfs.values():
            if orf.numPeptides() > 0:
                i += 1
        return i

###############################################################################

class Genome(object):
    """A collection of Chromosome objects indexed by accession."""
    def __init__(self,taxon=None):
        self.taxon = taxon
        self.chromosomes = {}

    def addOrf(self, orf, orfType):
        self.chromosomes[ orf.chromosome ].addOrf( orf, orfType )

    def numChromosomes(self):
        return len(self.chromosomes)

    def numOrfs(self,orfType='All'):
        count = 0
        for (name,chrom) in self.chromosomes.items():
            if orfType in ['All','Simple']:
                count += len(chrom.simpleOrfs)
            if orfType in ['All','PepOnly']:
                count += len(chrom.pepOnlyOrfs)
            if orfType in ['All','Complex']:
                count += len(chrom.complexOrfs)
        return count

    def makeChromosome(self,accession,seq=None):
        """Given an accession, creates a chromosome for it in the Genome and
        returns the Chromosome object"""
        if self.chromosomes.has_key( accession ):
            raise RuntimeError("Trying to create Chromosome %s again."%accession)

        chrom = Chromosome(accession,seq)
        self.chromosomes[accession] = chrom
        return chrom

    def getOrf(self,protName,chromName):
        """Given an ORFName and chromName returns the orf regardless of
        simple or complex membership.
        """
        return self.chromosomes[chromName].getOrf( protName )


    def filterORFs(self, filterList):
        "Takes a PGORFFilter.FilterList object and applies it to all simple ORFs"
        for (name,chrom) in self.chromosomes.items():
            filteredOrfDict = filterList.ApplyAllFilters( chrom.simpleOrfs )
            chrom.simpleOrfs = filteredOrfDict

            filteredOrfDict = filterList.ApplyAllFilters( chrom.pepOnlyOrfs )
            chrom.pepOnlyOrfs = filteredOrfDict

            filteredOrfDict = filterList.ApplyAllFilters( chrom.complexOrfs )
            chrom.complexOrfs = filteredOrfDict

    def addSeqToPepOnlyOrfs( self, sequenceFile, definitionParser=ORFFastaHeader):
        seqReader = bioseq.SequenceIO( sequenceFile )
        # Read in only the needed ORFs from the sequence file
        for seq in seqReader:
            seqDef = definitionParser( seq.acc )
            chromName = seqDef.Chromosome

            if self.chromosomes.has_key( chromName ):
                chrom = self.chromosomes[ chromName ]
                if chrom.pepOnlyOrfs.has_key( seqDef.ORFName ):
                    orf = chrom.pepOnlyOrfs[ seqDef.ORFName ]
                    orf.aaseq = seq.seq
                    orf.location = GenomicLocation.FromHeader( seqDef,
                                                              len(seq.seq),
                                                              addStop=True)
                else:
                    pass # ORF is not a pepOnly, so skip for now
            else:
                print "WARNING! unknown Chromosome %s" % chromName

###############################################################################

class GenbankGenomeReader(bioseq.FlatFileIO):
    """Returns a single Genome object populated with Chromosomes from each
    sequence in a genbank file, and locates ORFs from a six frame sequence file
    onto their annotated proteins. Uses biopython to parse the Genbank file.
    """
    def __init__(self, gbFile, sixFrameFile):
        bioseq.FlatFileIO.__init__(self,gbFile)
        self.orfReader = bioseq.SequenceIO(sixFrameFile)

    def makeGenomeWithProteinORFs(self):
        genome = Genome()
        # First read in all the CDS features from the genbank file
        for gb_rec in SeqIO.parse(self.io, 'genbank'):
            # Each gb record becomes it's own chromsome
            chrom = genome.makeChromosome( gb_rec.name, gb_rec.seq )

            # Store the cds SeqFeatures in the chromsome indexed by 3' end
            for feat in gb_rec.features:
                if feat.type == 'CDS':
                    if feat.strand == 1:
                        chrom.endToCDS[ feat.location.end.position ] = feat
                    else:
                        # biopython 1.53 seems to use 0, or space based coords
                        # so start is 1 less then what's in the genbank file
                        chrom.endToCDS[ feat.location.start.position + 1 ] = feat

        # Now read in the ORFs from the 6Frame file and match their
        # ends with the ends of the annotated proteins
        unusedOrfs = 0
        for orfSeq in self.orfReader:

            tmpOrf = OpenReadingFrame( orfSeq.acc, orfSeq.seq )

            if tmpOrf.name.startswith('XXX'): # Skip the decoy ORFs for now
                continue

            chrom = genome.chromosomes[ tmpOrf.chromosome ]
            orfThreePrime = tmpOrf.location.GetThreePrime()

            if chrom.endToCDS.has_key( orfThreePrime ):
                cds = chrom.endToCDS.pop( orfThreePrime )
                tmpOrf.CDS = cds # Not sure if we really need this but keep around for a little while
                orfStart = tmpOrf.location.start
                orfStop  = tmpOrf.location.stop
                # +1 is back to 1 based
                prot5Prime = cds.location.start.position + 1
                if cds.strand == -1:
                    prot5Prime = cds.location.end.position

                # separate simple ORF from complex ORFs for now
                # not sure if we'll need to further separate complex ORFs
                if prot5Prime >= orfStart and prot5Prime <= orfStop:
                    # Create a LocatedProtein object for this protein
                    locProt = LocatedProtein( GenomicLocation(
                        cds.location.start.position + 1,
                        cds.location.end.position,
                        cds.strand == 1 and '+' or '-',
                        tmpOrf.chromosome
                    ))
                    locProt.name = cds.qualifiers['product'][0]
                    locProt.ORFName = tmpOrf.name
                    # and add it to the ORF
                    tmpOrf.addLocatedProtein( locProt )
                    chrom.addOrf( tmpOrf, 'Simple' )

                elif len(cds.sub_features) > 0:
                    chrom.addOrf( tmpOrf, 'Complex' )
                else:
                    chrom.addOrf( tmpOrf, 'Complex' )
            else:
                # ORF without a 3' mapping to a protein
                unusedOrfs += 1

        print "Read %d chromosomes, with %d Simple and %d Complex ORFs" % (
            genome.numChromosomes(),
            genome.numOrfs('Simple'),
            genome.numOrfs('Complex') )
        print "%d unused ORFs from 6frame fasta." % unusedOrfs
        for acc,chrom in genome.chromosomes.items():
            for cds in chrom.endToCDS.values():
                print "Warning, unmapped protein %s on chrom %s" % (
                        cds.qualifiers['locus_tag'][0], acc )

        return genome
