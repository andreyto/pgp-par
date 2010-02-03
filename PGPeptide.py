
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
                # other is contained sub sequenc              return (otherStart,otherStop)
            else:
                return (otherStart, myStop)
        else: # mSta otherStart
            if myStop < otherStop:
            # self is contained sub sequence
                return (myStart,myStop)
            else:
                return (myStart, otherStop)

class LocatedPeptide(object):
    def __init__(self,Aminos, location=None):
        self.location = location
        self.name = None
        self.aminos = Aminos
        self.bestScore = None
        self.spectrumCount = 0
        self.isUnique = None
        self.TrypticNTerm = None
        self.TrypticCTerm = None
        self.proteinName = None

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
