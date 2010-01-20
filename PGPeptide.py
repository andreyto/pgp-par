
class GenomicLocation:
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
        self.__strand = strand
        
    def getStart(self):
        return self.__start

    def getStop(self):
        return self.__stop

    def overlap(self,otherLocation):
        """Returns None if the locations don't overlap.
        Otherwise, it returns a tuple of the overlapping region.
        """
        (myStart,myStop) = (self.getStart(), self.getStop())
        (otherStart,otherStop) = (otherLocation.getStart(), otherLocation.getStop())
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

class LocatedPeptide:
    def __init__(self,location):
        self.location = location

    def isTryptic(self):
        # TBD
        return False


