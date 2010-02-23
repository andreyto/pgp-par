import GFFIO
import bioseq

class GenomicLocation(object):
    def __init__(self,start,stop,strand):
        """Start and stop of a region on a sequence.
        Stop must be > then start.
        """
        if start > stop:
            print "you tried %s <? %s"%(start, stop)
            raise ValueError("Start > Stop")

        if not strand in ['+','-']:
            raise ValueError("strand must be + or -")

        self.__start = start
        self.__stop = stop
        self.__frame = None
        self.chromosome = None
        self.strand = strand
        
    def __str__(self):
        return "%d,%d %s%d" % (self.start,self.stop,self.strand,self.frame)

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

    @property
    def start(self):
        'Lesser coordinate on the sequence.'
        return self.__start
    

    @property
    def stop(self):
        'Greater coordinate on the sequence.'
        return self.__stop

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
        return "%s in %s, unique=%s, %s"%(self.aminos, self.ORFName, self.isUnique, self.location)

    def GetStart(self):
        return self.location.start

    def GetStop(self):
        return self.location.stop

    def Strand(self):
        return self.location.strand

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
            self.SetLocation(ParsedHeader, len(AASequence))
        if name:
            self.name = name

        self.GCWholeORF = None
        self.GCPredictedProtein = None
        self.GCObservedRegion = None
        
    def __str__(self):
        NumUniquePeptides = 0
        for Peptide in self.peptideIter():
            if Peptide.isUnique:
                NumUniquePeptides += 1
        return "%s as %s, %s (%s) peptides, %s"%(self.name, self.annotatedProtein, self.numPeptides(), NumUniquePeptides, self.location)

    def Chromosome(self):
        return self.location.chromosome

    def SetLocation(self, ParsedHeader, AALength):
        """
        Parameters: ORFFastaHeader object, length of the amino acids in the ORF
        Return: none
        Description: Fill in the information needed for the location object, create 
        and assign object
        """
        FivePrime = ParsedHeader.Start
        ThreePrime = -1 #set below
        Strand = ParsedHeader.Strand
        Chromosome = ParsedHeader.Chromosome
        Frame = ParsedHeader.Frame
        #now some math to figure out the stop nuc
        if (Strand == "-"):
            CodingThreePrime = FivePrime - (AALength * 3) + 1
            ThreePrime = CodingThreePrime - 3 #three bases of the stop codon
        else:
            #now the position.  We first get the protein start nuc, and then offset to the peptide start
            CodingThreePrime = FivePrime + (AALength * 3) - 1 # we do a minus one because the bases are inclusive
            ThreePrime = CodingThreePrime + 3 # three bases of the stop codon are INCLUDED
        if Strand == "+":
            SimpleLocation = GenomicLocation(FivePrime, ThreePrime, Strand)
        else:
            SimpleLocation = GenomicLocation(ThreePrime, FivePrime, Strand)
        SimpleLocation.chromosome = Chromosome
        SimpleLocation.frame = Frame
        self.location = SimpleLocation

    def GetTranslation(self):
        """Parameters: none
        Return: Amino acid translation of the ORF
        Description: gets the translation of the entire open reading frame
        """
        return self.aaseq

    def GetStrand(self):
        return self.location.strand
    
    def GetObservedSequence(self):
        """Parameters: none
        Return: amino acid sequence 
        Description: We want to get the sequence that we have evidence as
        being translated.  This is from the first peptide through the stop
        """
        #1. get the first peptide, in order
        FirstPeptide = self.GetFivePrimePeptide()
        #2. map it to the aaseq
        FirstResidue = self.aaseq.find(FirstPeptide.aminos)
        #3. return the slice
        return self.aaseq[FirstResidue:]

    def GetFivePrimePeptide(self):
        """Parameters: None
        Return: the five prime most peptide
        Description: cycle through the peptides and return the one closest
        to the start
        """
        PeptideWinner = None
        FivePrimeWinner = None
        for Peptide in self.peptideIter():
            if not PeptideWinner:
                PeptideWinner = Peptide
                FivePrimeWinner = Peptide.GetFivePrimeNucleotide()
                continue
            if self.location.strand == "+":
                #if + strand, get the smallest
                if Peptide.GetFivePrimeNucleotide() < FivePrimeWinner:
                    PeptideWinner = Peptide
                    FivePrimeWinner = Peptide.GetFivePrimeNucleotide()
            else: # if - strand, get the biggest 
                if Peptide.GetFivePrimeNucleotide() > FivePrimeWinner:
                    PeptideWinner = Peptide
                    FivePrimeWinner = Peptide.GetFivePrimeNucleotide()
        return PeptideWinner
                
                
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
       

class GFF(GFFIO.File):
    'This class creates OpenReadingFrame objects from GFF files.'

    def generateORFs(self, sequenceFile, definitionParser=ORFFastaHeader):
        '''Parameters: A sequence file supported by SequenceIO and
        optionally a Class for parsing the sequence accession.
        Return: Dictionary of OpenReadingFrame objects, and their LocatedPeptides
        Description: Reads the peptides from the GFF, and the ORFs from the sequence file
        '''
        seqReader = bioseq.SequenceIO( sequenceFile )
        observedORFs = {}
        
        for gffRec in self:
            protein = gffRec.attributes['Parent']
            if not observedORFs.has_key( protein ):
                observedORFs[ protein ] = OpenReadingFrame(name=protein)

            location = GenomicLocation(gffRec.start, gffRec.end, gffRec.strand)
            # Peptide is encoded as the name, since it's generally short 
            peptide = LocatedPeptide( gffRec.attributes['Name'], location) 
            peptide.name = gffRec.attributes['ID']
            peptide.bestScore = gffRec.score

            orf = observedORFs[ protein ]
            orf.addLocatedPeptide( peptide )

        for seq in seqReader:
            seqDef = definitionParser( seq.acc )
            orfName = seqDef.ORFName
            if observedORFs.has_key( orfName ):
                orf = observedORFs[ orfName ]
                orf.aaseq = seq.seq
                orf.SetLocation( seqDef, len(seq.seq))

        return observedORFs

    def writeORFPeptides(self, orf):
        gffRec = GFFIO.Record()
        gffRec.source = 'Proteomics'
        gffRec.type = 'polypeptide'
        gffRec.seqid = orf.Chromosome()
        gffRec.attributes['Parent'] = orf.name

        for peptide in orf.peptideIter():
            gffRec.start = peptide.GetStart()
            gffRec.end   = peptide.GetStop()
            gffRec.score = peptide.bestScore
            gffRec.strand= peptide.Strand()
            gffRec.attributes['Name'] = peptide.aminos
            gffRec.attributes['ID'] = peptide.name
            self.write( gffRec )
