
class GenomicLocation(object):
    def __init__(self,start,stop,strand):
        """Start and stop of a region on a sequence.
        Stop must be > then start.
        """
        if start > stop:
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
        self.proteinName = None #this is actually the ORF Name.  as we put them into ORFs

    def isTryptic(self):
        # TBD
        return False

    def GetFivePrimeNucelotide(self):
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
        return "%s in %s, unique=%s, %s"%(self.aminos, self.proteinName, self.isUnique, self.location)


class LocatedProtein(object):
    def __init__(self,location):
        self.location = location
        self.name = None
        self.ORFName = None
        
    def AddOneAminoAcidFivePrime(self):
        """This is a very very ugly kludge.  to map proteins on to the ORFs, we remove the first amino acid
        because if it's an alternate start site, the ORF has V, and the protein has M.  So now we need to add
        the ability to add back that one amino acid to our position.  Remember, start and stop are numerical order
        and NOT associated with 5' or 3'.  
        """
        self.location.AddOneAminoAcidFivePrime()
        return
        
        
class OpenReadingFrame(object):
    def __init__(self):
        self.location = None
        self.__peptides = []
        self.annotatedProtein = None
        self.name = None
        self.aaseq = None
        self.naseq = None

    def numPeptides(self):
        'Returns the number of peptides in the ORF.'
        return len(self.__peptides)

    def addLocatedPeptides(self,peptideList):
        'Adds a list of LocatedPeptide objects to the ORF.'
        self.__peptides.extend( peptideList )

    def addRawPeptides(self,peptideList):
        'Adds a list of peptide sequences without location.'
        for pep in peptideList:
            pepObj = LocatedPeptide()
            pepObj.aminos = pep
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
       