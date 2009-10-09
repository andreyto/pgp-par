'''
@author: Eli Venter

A generic sequence class for passing sequences around.
Also support classes for doing sequence IO.
'''

class Sequence():
    '''
    A generic sequence class
    '''

    def __init__(self,id=-1,seq='',desc=''):
        '''
        id is a sequence id
        seq is the sequence string itself
        desc is a description of the sequence ie defline
        '''
        self.id = id  
        self.seq = seq
        self.desc = desc
        
class FastaReader():
    '''
    A class that reads fasta files and produces Sequence objects.
    Should be iterable as well.
    '''

    def __init__(self,file):
        if type(file) is file:
            self.io = file

        elif type(file) is str:
            self.name = str
            self.io = open( file, 'r')
        else:
            raise TypeError

    def __iter__(self):
        seq = None
        seqList = None
        for line in self.io.xreadlines():
            line = line.rstrip()
            if line[0] == '>':
                # Next defline, now done with previous seq, so yield it
                if seqList:
                    seq.seq = ''.join(seqList)
                    yield seq

                space = line.index(' ')
                seq = Sequence( id=line[1:space] )
                seq.desc = line[space+1:]
                seqList = []
            else:
                seqList.append(line)

        if seqList:
           seq.seq = ''.join(seqList)
           yield seq

        raise StopIteration
