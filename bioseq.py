'''
@author: Eli Venter

A generic sequence class for passing sequences around.
Also support classes for doing sequence IO.
'''

class Sequence():
    '''
    A generic sequence class
    '''

    def __init__(self,acc=-1,seq='',desc=''):
        '''
        acc is a sequence id
        seq is the sequence string itself
        desc is a description of the sequence ie defline
        '''
        self.acc = acc  
        self.seq = seq
        self.desc = desc
        
class FlatFileReader():
    '''
    A parent class for all the classes that read flat files to
    inherit from. Children should always be iterable.
    '''        
    def __init__(self,infile):
        if type(infile) is file:
            self.io = infile

        elif type(infile) is str:
            self.name = str
            self.io = open( infile, 'r')
        else:
            raise TypeError

class FlatFileWriter():
    '''
    A parent class for all the classes that write flat files to
    inherit from.
    '''        
    def __init__(self,infile):
        if type(infile) is file:
            self.io = infile

        elif type(infile) is str:
            self.name = str
            self.io = open( infile, 'w')
        else:
            raise TypeError
        
    def close(self):
        self.io.close()

class FastaReader(FlatFileReader):
    '''
    A class that iterates through fasta files and produces Sequence objects.
    '''
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
                seq = Sequence( line[1:space] )
                seq.desc = line[space+1:]
                seqList = []
            else:
                seqList.append(line)
        # Don't forget to yield the final seq in the file
        if seqList:
            seq.seq = ''.join(seqList)
            yield seq

        raise StopIteration
    
class FastaOut(FlatFileWriter):
    '''
    A class that takes Sequence objects and writes them out in fasta format 
    to the given file. Includes a linesize accessor.
    '''
    def __init__(self,out):
        FlatFileWriter.__init__(self, out)
        self.linesize = 60
        
    def write(self,seq):
        self.io.write(">%s %s\n" % (seq.acc,seq.desc))
        i = 0
        l = len(seq.seq)
        while i < l:
            self.io.write("%s\n" % seq.seq[i:i+self.linesize])
            i += self.linesize
            