'''
@author: Eli Venter

A generic sequence class for passing sequences around.
Also support classes for doing sequence IO.
'''

import os, re

class Sequence(object):
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
            
class FlatFileReader(object):
    '''
    A parent class for all the classes that read flat files to
    inherit from. Children should always be iterable.
    '''        
    def __init__(self,infile):
        if type(infile) is file:
            self.io = infile

        elif type(infile) is str:
            self.name = infile
            self.io = open( infile, 'r')
        else:
            raise TypeError

    def __del__(self):
        self.io.close()

class FlatFileWriter(object):
    '''
    A parent class for all the classes that write flat files to
    inherit from.
    '''        
    def __init__(self,outfile):
        if type(outfile) is file:
            self.io = outfile

        elif type(outfile) is str:
            self.name = outfile
            self.io = open( outfile, 'w')
        else:
            raise TypeError
        
    def __del__(self):
        self.io.close()

    def close(self):
        self.io.close()

    def write(self,seq):
        pass

    def writeMulti(self,seqList):
        for i in seqList:
            self.write(i)

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

                defline = re.split('\s+',line[1:],1)
                seq = Sequence( defline[0] )
                if len(defline) == 1:
                    seq.desc = ''
                else:
                    seq.desc = defline[1]

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
        self.linesize = 80
        
    def write(self,seq):
        if seq.desc:
            self.io.write(">%s %s\n" % (seq.acc,seq.desc))
        else:
            self.io.write(">%s\n" % seq.acc )

        for i in xrange(0, len(seq.seq), self.linesize):
            self.io.write("%s\n" % seq.seq[i:i+self.linesize])

class SequenceIO(object):
    FormatTable = {'fasta': (FastaReader,FastaOut)}
    
    def __init__(self,fileName,mode='r'):
        '''
        Wrapper for all types of Sequence IO.
        fileName is the name of the file to read or write
        mode defaults to r for reading and use w for writing
        '''
        self.fileName = fileName
        self.mode = mode == 'w' and 1 or 0
        
        ext = os.path.splitext(fileName)[1].lower()
        if ext in ['.fa','.fasta','.fsa']:
            self.handle = SequenceIO.FormatTable['fasta'][self.mode](fileName)
    
    def __iter__(self):
        return self.handle.__iter__()
        
    def write(self,seqData):
        if type(seqData) is list:
            self.handle.multiWrite(seqData)
        else:
            self.handle.write(seqData)                

    def close(self):
        self.handle.close()

    def set(self, name, value):
        'Method to set properties on the implementation classes.'
        setattr(self.handle, name, value)