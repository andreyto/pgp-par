#!/usr/bin/env python

"""
Shuffle all the records in a database.  Useful for generating a database of all bogus records,
to get an idea of the rate at which matches of a given quality are generated when there's
nothing real to match.  Or a database of half valid, half invalid records.
"""
import os
import sys
import string
import getopt
import struct
import random
random.seed(1)

UsageInfo = """
ShuffleDB - Produce a 'red herring' database of erroneous peptide records.
Options:
 -r [Trie file name]: Path to input database
 -w [FileName]: Path to output database
 -s: If set, proteins will be REVERSED.  Default behavior is to SHUFFLE.
 -b: If set, ONLY the scrambled proteins are written out.  Default
     behavior is to write both forward and scrambled proteins.
 -p: In shuffled mode, avoid repeating peptides of length 8 or more
     in the shuffled database.  (Treat I and L as identical. 
     Don't treat Q and K as identical)  Requires a little longer to run;
     some "bad words" (repeated 8mers) will still be seen for repetitive
     records.
 -t: Number of shuffled copies to write out (defaults to 1)
 
Example:
 ShuffleDB.py -r database\Shewanella.trie -w database\ShewanellaHalf.trie -p
"""

##def EncodeQuad(Str):
##    ValA = ord(Str[0]) - 65
##    ValB = ord(Str[1]) - 65
##    ValC = ord(Str[2]) - 65
##    ValD = ord(Str[3]) - 65
##    return ValA + 32*ValB + 32*32*ValC + 32*32*32*ValD
##def EncodeOct(Str):
##    return struct.pack("<ii", EncodeQuad(Str[:4]), EncodeQuad(Str[4:]))
##def DecodeOct(Str):
##    (ValA, ValB) = struct.unpack("<ii", Str)
##    return DecodeQuad(ValA) + DecodeQuad(ValB)
##def DecodeQuad(Value):
##    ValA = Value % 32
##    Value -= ValA
##    Value /= 32
##    ValB = Value % 32
##    Value -= ValB
##    Value /= 32
##    ValC = Value % 32
##    Value -= ValC
##    Value /= 32
##    ValD = Value
##    print ValA, ValB, ValC, ValD
##    return chr(ValA + 65) + chr(ValB + 65) + chr(ValC + 65) + chr(ValD + 65)

class Shuffler:
    MAX_SHUFFLE_ATTEMPTS = 100
    PARANOIA_PEPTIDE_LENGTH = 8
    def __init__(self):
        self.ShuffleFlag = 1
        self.WriteBothFlag = 1
        self.TrieFileName = None
        self.OutputFileName = None
        self.BogusTimes = 1
        self.BogusProteins = 0
        self.ParanoiaFlag = 0
        self.TotalBadWordCount = 0
    def LoadProteinNames(self, IndexPath):
        File = open(IndexPath, "rb")
        self.ProteinNames = []
        BlockSize = struct.calcsize("<qi80s")
        while (1):
            Block = File.read(BlockSize)
            if not Block:
                break
            Tuple = struct.unpack("<qi80s", Block)
            Name = Tuple[-1]
            NullPos = Name.find('\0')
            if NullPos != -1:
                Name = Name[:NullPos]
            self.ProteinNames.append(Name)
        File.close()
    def Main(self):
        IndexPath = os.path.splitext(self.TrieFileName)[0] + ".index"
        self.LoadProteinNames(IndexPath)
        TrieFile = open(self.TrieFileName, "rb")
        if self.ParanoiaFlag:
            self.ForbiddenPeptides = {}
            DB = TrieFile.read()
            for Pos in range(len(DB) - self.PARANOIA_PEPTIDE_LENGTH):
                if Pos % 10000 == 0:
                    print "%s (%s)..."%(Pos, len(self.ForbiddenPeptides.keys()))
                Peptide = DB[Pos:Pos + self.PARANOIA_PEPTIDE_LENGTH].replace("I", "L")
                if Peptide.find("X")!=-1:
                    # Peptides containing X need not be forbidden, because they will
                    # never be matched!
                    continue
                self.ForbiddenPeptides[Peptide] = 1
            TrieFile.seek(0)
            print "(Note: %s peptide words are forbidden)"%len(self.ForbiddenPeptides.keys())
        NewIndexPath = os.path.splitext(self.OutputFileName)[0] + ".index"
        self.OutputTrieFile = open(self.OutputFileName, "wb")
        self.OutputIndexFile = open(NewIndexPath, "wb")
        Sequence = ""
        ProteinIndex = 0
        while (1):
            Data = TrieFile.read(1024)
            if not Data:
                break
            Sequence += Data
            Pos = Sequence.find("*")
            while (Pos != -1):
                self.WriteProtein(Sequence[:Pos], ProteinIndex)
                ProteinIndex += 1
                Sequence = Sequence[Pos+1:]
                Pos = Sequence.find("*")
        if (Sequence):
            self.WriteProtein(Sequence, ProteinIndex)
            ProteinIndex += 1
            #List = list(Sequence)
            #List.reverse()
            #Protein = string.join(List,"")
            #ReversedTrieFile.write(Protein)
            #ReversedTrieFile.write("*")
        self.OutputTrieFile.close()
        self.OutputIndexFile.close()
        print "Wrote out %s proteins."%ProteinIndex
        print "Wrote out %d bogus proteins."%self.BogusProteins
        print "Total bad words:", self.TotalBadWordCount
    def ShuffleProtein(self, Sequence):
        """
        Produce the invalid (shuffled) version of a protein.
        """
        Residues = list(Sequence)
        if not self.ShuffleFlag:
            Residues.reverse()
            return string.join(Residues, "")
        if not self.ParanoiaFlag:
            random.shuffle(Residues)
            return string.join(Residues, "")
        # And now, the interesting case: We shall shuffle the protein, and we shall apply some
        # heuristics along the way to minimize the number of shared k-mers.
        BestBadWordCount = 9999
        BestPeptideString = None
        for AttemptIndex in range(10):
            random.shuffle(Residues)
            BadWordCount = 0
            for Pos in range(len(Residues) - self.PARANOIA_PEPTIDE_LENGTH):
                WordResidues = Residues[Pos:Pos + self.PARANOIA_PEPTIDE_LENGTH]
                Word = string.join(WordResidues, "").replace("I", "L")
                if self.ForbiddenPeptides.has_key(Word):
                    # Damn, this shuffling shares a word!  Maybe we can re-shuffle this
                    # word and solve the problem:
                    FixedFlag = 0
                    for WordShuffleIndex in range(10):
                        random.shuffle(WordResidues)
                        FixedWord = string.join(WordResidues, "").replace("I", "L")
                        if self.ForbiddenPeptides.has_key(FixedWord):
                            # The shuffled word is no good!
                            continue
                        # We shuffled a word, and in so doing, we changed the preceding
                        # words.  Let's check whether they're now forbidden:
                        BrokeOldWord = 0
                        for StepsBack in range(1, self.PARANOIA_PEPTIDE_LENGTH):
                            PrevPos = Pos - StepsBack
                            if PrevPos < 0:
                                break
                            PrevWord = Residues[PrevPos:Pos]
                            PrevWord.extend(WordResidues[:-StepsBack])
                            PrevWord = string.join(PrevWord, "").replace("I", "L")
                            if self.ForbiddenPeptides.has_key(PrevWord):
                                BrokeOldWord = 1
                                #print "Preceding word %s is '%s': FORBIDDEN!"%(StepsBack, PrevWord)
                                break
                            #print "Preceding word %s is '%s'"%(StepsBack, PrevWord)
                        if not BrokeOldWord:
                            FixedFlag = 1
                            break
                    if FixedFlag:
                        # This word (and the previous words that overlap it) is now ok.
                        Residues[Pos:Pos + self.PARANOIA_PEPTIDE_LENGTH] = WordResidues
                    else:
                        # We couldn't fix the word by shuffling it.  Increment the bad word count:
                        BadWordCount += 1
            if BadWordCount == 0:
                #print "Protein '%s...' shuffled with no bad words"%(Sequence[:20])
                return string.join(Residues, "")
            if BadWordCount < BestBadWordCount:
                BestBadWordCount = BadWordCount
                BestPeptideString = string.join(Residues, "")
        print "Protein '%s...' shuffled with %s bad words"%(Sequence[:20], BestBadWordCount)
        self.TotalBadWordCount += BestBadWordCount
        return BestPeptideString
    def WriteProtein(self, Sequence, ProteinIndex):
        """
        Given a protein sequence, and protein index number (for looking up the name),
        write a scrambled or reversed record to the output database.  (And write the
        original, if the -b flag was specified)
        """
        for ShuffleIndex in range(self.BogusTimes):
            ShuffledProtein = self.ShuffleProtein(Sequence)
            OutputFilePos = self.OutputTrieFile.tell()
            self.OutputTrieFile.write(ShuffledProtein)
            self.OutputTrieFile.write("*")
            if self.BogusTimes > 1:
                ShuffledName = "XXX.%d.%s"%(ShuffleIndex, self.ProteinNames[ProteinIndex])
            else:
                ShuffledName = "XXX.%s"%self.ProteinNames[ProteinIndex]
            Block = struct.pack("<qi80s", 0, OutputFilePos, ShuffledName)
            self.OutputIndexFile.write(Block)
            self.BogusProteins +=1
        # If we're writing both the red herrings and the originals,
        # then write the original protein as well now:
        if self.WriteBothFlag:
            OutputFilePos = self.OutputTrieFile.tell()
            self.OutputTrieFile.write(Sequence)
            self.OutputTrieFile.write("*")
            Block = struct.pack("<qi80s", 0, OutputFilePos, self.ProteinNames[ProteinIndex])
            self.OutputIndexFile.write(Block)
    def ParseCommandLine(self, Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:w:sbt:p")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-r":
                self.TrieFileName = Value
            elif Option == "-w":
                self.OutputFileName = Value
            elif Option == "-s":
                self.ShuffleFlag = 0
            elif Option == "-b":
                self.WriteBothFlag = 0
            elif Option == "-t":
                self.BogusTimes = int(Value)
            elif Option == "-p":
                self.ParanoiaFlag = 1
            else:
                print "** Warning: Option %s not understood"%Option
        
def main(args):
    App = Shuffler()
    App.ParseCommandLine(args)
    if not App.TrieFileName or not App.OutputFileName:
        raise Exception(UsageInfo)
    App.Main()

if __name__ =="__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "* Warning: psyco not found"

    main(sys.argv[1:])
