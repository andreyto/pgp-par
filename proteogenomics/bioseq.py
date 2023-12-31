###############################################################################
#                                                                             # 
#       Copyright (c) 2009 J. Craig Venter Institute.                         #     
#       All rights reserved.                                                  #
#                                                                             #
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.    #
#                                                                             #
###############################################################################



'''
@author: Eli Venter

A generic sequence class for passing sequences around.
Also support classes for doing sequence IO.
'''

import os, re, fileinput, struct, bisect
import gzip, bz2

from StringIO import StringIO

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

class FlatFileIO(object):
    '''
    A parent class for all the classes that do IO to flat files.
    Children should always be iterable.
    '''

    # Support transparent compression based on filename extension
    openTable = { '.bz2': bz2.BZ2File, '.gz': gzip.open }

    def __init__(self, fileForIO, mode='r'):
        self.io = None
        if type(fileForIO) is file:
            self.io = fileForIO

        elif type(fileForIO) is str: #assumes this is the path for a file to open .  mode default is read.  
            self.name = fileForIO
            self.io = self.open( fileForIO, mode)

        elif isinstance(fileForIO, list): # combines a list of files together into 1
            if len(fileForIO) == 1:
                # Treat just like a single file
                self.name = fileForIO[0]
                self.io = self.open( fileForIO[0], mode)
            else:
                self.io = fileinput.input(fileForIO)

        elif isinstance(fileForIO, StringIO):
            self.io = fileForIO
        else:
            raise TypeError("Need string or io handle, got %s" % type(fileForIO))

    def __del__(self):
        if self.io:
            self.io.close()

    def open(self, fileName, mode):
        (base, extension) = os.path.splitext(fileName)
        openCall = None
        if FlatFileIO.openTable.has_key( extension ):
            openCall = FlatFileIO.openTable[ extension ]
        else:
            openCall = open
        # Open the file using the correct compression method as needed 
        return openCall( fileName, mode )

    def close(self):
        self.io.close()

    def write(self,data):
        self.io.write(data)

class TrieReader(FlatFileIO):
    def __iter__(self):
        seqs = self.io.read().split('*')
        self.close()

        IndexPath = os.path.splitext(self.name)[0] + ".index"
        IndexFile = open(IndexPath, "rb")
        BlockSize = struct.calcsize("<qi80s")
        i = 0
        while (1):
            Block = IndexFile.read(BlockSize)
            if not Block:
                break
            Info = struct.unpack("<qi80s", Block)
            Name = Info[2]
            NullPos = Name.find("\0")
            if NullPos != -1:
                Name = Name[:NullPos]
            seq = Sequence( Name, seqs[i], Info[1] )
            yield seq
            i += 1

        IndexFile.close()

class TrieIndexSeqs(object):
    def __init__(self,trieFile):
        self.reader = TrieReader(trieFile)
        self.positions = []
        self.ids = []
        self.seqs = []
        self.trie = ''

    def index(self):
        offset = 0
        for seq in self.reader:
            self.ids.append( seq.acc )
            self.seqs.append( seq.seq )
            # Store begin of sequence
            self.positions.append( offset );
            offset += len(seq.seq) + 1 # extra 1 for the * in the trie

        self.trie = '*'.join( self.seqs )

    def indexAtOffset(self, beginInTrie):
        # bisect right will always give us the sequence after 
        # the one we want, so subtract 1
        return bisect.bisect_right( self.positions, beginInTrie) - 1

    def accessionsWhereSeqFound( self, seqToFind):
        '''Given a sequence returns a list of tuples of locations in the index.
        The tuple contains the sequence accession, sequence number, and 
        start of the search sequence relative to it's parents begin.
        '''
        accIndexPairs = []
        offset = 0
        while 1:
            offset = self.trie.find( seqToFind, offset )
            if offset == -1:
                break
            index = self.indexAtOffset( offset )
            offsetInSeq = offset - self.positions[ index ]
            accIndexPairs.append( (self.ids[index], index, offsetInSeq) )
            offset += 1

        return accIndexPairs

import IndexSearch

class QGramIndex(TrieIndexSeqs):
    def index(self):
        TrieIndexSeqs.index(self)
        self.search = IndexSearch.IndexSearch( self.trie )

    def accessionsWhereSeqFound( self, seqToFind):
        accIndexPairs = []
        for offset in  self.search.find( seqToFind ):
            index = self.indexAtOffset( offset )
            offsetInSeq = offset - self.positions[ index ]
            accIndexPairs.append( (self.ids[index], index, offsetInSeq) )

        return accIndexPairs


class FastaReader(FlatFileIO):
    '''
    A class that iterates through fasta files and produces Sequence objects.
    '''
    def __iter__(self):
        seq = None
        seqList = None
        for line in self.io:
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

class FastaOut(FlatFileIO):
    '''
    A class that takes Sequence objects and writes them out in fasta format 
    to the given file. Includes a linesize accessor.
    '''
    def __init__(self,out):
        FlatFileIO.__init__(self, out, mode='w')
        self.linesize = 80

    def write(self,seq):
        if seq.desc:
            self.io.write(">%s %s\n" % (seq.acc,seq.desc))
        else:
            self.io.write(">%s\n" % seq.acc )

        for i in xrange(0, len(seq.seq), self.linesize):
            self.io.write("%s\n" % seq.seq[i:i+self.linesize])

class SequenceIO(object):
    FormatTable = {'fasta': (FastaReader,FastaOut),
                   '.trie': (TrieReader,None) }

    def __init__(self,fileName,mode='r'):
        '''
        Wrapper for all types of Sequence IO.
        fileName is the name of the file to read or write
        mode defaults to r for reading and use w for writing
        '''
        self.fileName = fileName
        if isinstance(fileName, list):
            self.fileName = fileName[0]
        self.mode = mode == 'w' and 1 or 0

        ext = os.path.splitext(self.fileName)[1].lower()
        if ext in ['.fa','.fasta','.fsa', '.fna', '.faa']:
            self.handle = SequenceIO.FormatTable['fasta'][self.mode](fileName)
        else:
            self.handle = SequenceIO.FormatTable[ext][self.mode](fileName)

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
