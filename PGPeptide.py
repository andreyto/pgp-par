
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
                # other is contained sub sequence
                return (otherStart,otherStop)
            else:
                return (otherStart, myStop)
        else: # myStart > otherStart
            if myStop < otherStop:
                # self is contained sub sequence
                return (myStart,myStop)
            else:
                return (myStart, otherStop)

class LocatedPeptide(object):
    def __init__(self,location):
        self.location = location
        self.name = None
        self.amino = None
        self.bestScore = None
        self.spectrumCount = 0
        self.unique = None

    def isTryptic(self):
        # TBD
        return False


