#!/usr/bin/env python

UsageInfo = """ SixFrameFasta.py
This program takes a fasta of genome sequence, and turns out the six
frame protein sequence in fasta format. > Protein.Chr.startnucleutide.frame
It appears that this produces GFF compliant 1 based indexing.
 -r [FileName] Genomic Sequence
 -w [FileName] Output filename
 -c [Name] Chromosome name, e.g. Chr3
"""


import os
import getopt
import sys

import bioseq

class ProteinTranslationClass:
    def __init__(self):
        """Simply the translation table, so it can be accessed
        """
        self.Table = {}
        self.Table["ATT"] = "I" #Isoleucine      I    ATT, ATC, ATA
        self.Table["ATC"] = "I"
        self.Table["ATA"] = "I"
        self.Table["CTT"] = "L" #Leucine         L    CTT, CTC, CTA, CTG, TTA, TTG
        self.Table["CTC"] = "L"
        self.Table["CTA"] = "L"
        self.Table["CTG"] = "L"
        self.Table["TTA"] = "L"
        self.Table["TTG"] = "L"
        self.Table["GTT"] = "V" #Valine          V    GTT, GTC, GTA, GTG
        self.Table["GTC"] = "V"
        self.Table["GTA"] = "V"
        self.Table["GTG"] = "V"
        self.Table["TTT"] = "F" #Phenylalanine   F    TTT, TTC
        self.Table["TTC"] = "F"
        self.Table["ATG"] = "M" #Methionine      M    ATG
        self.Table["TGT"] = "C" #Cysteine        C    TGT, TGC
        self.Table["TGC"] = "C"
        self.Table["GCT"] = "A" #Alanine         A    GCT, GCC, GCA, GCG
        self.Table["GCC"] = "A"
        self.Table["GCA"] = "A"
        self.Table["GCG"] = "A"
        self.Table["GGT"] = "G" #Glycine         G    GGT, GGC, GGA, GGG
        self.Table["GGC"] = "G"
        self.Table["GGA"] = "G"
        self.Table["GGG"] = "G"
        self.Table["CCT"] = "P" #Proline         P    CCT, CCC, CCA, CCG
        self.Table["CCC"] = "P"
        self.Table["CCA"] = "P"
        self.Table["CCG"] = "P"
        self.Table["ACT"] = "T" #Threonine       T    ACT, ACC, ACA, ACG
        self.Table["ACC"] = "T"
        self.Table["ACA"] = "T"
        self.Table["ACG"] = "T"
        self.Table["TCT"] = "S" #Serine          S    TCT, TCC, TCA, TCG, AGT, AGC
        self.Table["TCC"] = "S"
        self.Table["TCA"] = "S"
        self.Table["TCG"] = "S"
        self.Table["AGT"] = "S"
        self.Table["AGC"] = "S"
        self.Table["TAT"] = "Y" #Tyrosine        Y    TAT, TAC
        self.Table["TAC"] = "Y"
        self.Table["TGG"] = "W" #Tryptophan      W    TGG
        self.Table["CAA"] = "Q" #Glutamine       Q    CAA, CAG
        self.Table["CAG"] = "Q"
        self.Table["AAT"] = "N" #Asparagine      N    AAT, AAC
        self.Table["AAC"] = "N"
        self.Table["CAT"] = "H" #Histidine       H    CAT, CAC
        self.Table["CAC"] = "H"
        self.Table["GAA"] = "E" #Glutamic acid   E    GAA, GAG
        self.Table["GAG"] = "E"
        self.Table["GAT"] = "D" #Aspartic acid   D    GAT, GAC
        self.Table["GAC"] = "D"
        self.Table["AAA"] = "K" #Lysine          K    AAA, AAG
        self.Table["AAG"] = "K"
        self.Table["CGT"] = "R" #Arginine        R    CGT, CGC, CGA, CGG, AGA, AGG
        self.Table["CGC"] = "R"
        self.Table["CGA"] = "R"
        self.Table["CGG"] = "R"
        self.Table["AGA"] = "R"
        self.Table["AGG"] = "R"
        self.Table["TAA"] = "*" #Stop codons     *    TAA, TAG, TGA
        self.Table["TAG"] = "*" 
        self.Table["TGA"] = "*" 


class AbacusClass(ProteinTranslationClass):
    def __init__(self):
        self.InputFile = None
        self.OutputFile = "RenameYourOutput.txt"
        self.ChromosomeName = "Chr1" #you should rename this
        self.FrameStrings = ['','','']
        self.ProteinCount = 0
        self.FrameStarts = [1,2,3] #start nucleotide. 1 based
        self.MinProteinLength = 7 # if it's a small ORF, don't bother printing it out
        self.TempRTFileName = "Temp.RT.fasta"
        self.TotalDNAIndex = 0 #start at base 0, increment whenever you consume a letter (including the first time, which then starts us at base1)
        ProteinTranslationClass.__init__(self)

    def Main(self):
        self.OutHandle = open(self.OutputFile, "wb")
        self.TranslateChromosome(self.InputFile)
        self.ReverseTranscribeChromosome(self.InputFile)
        self.TranslateChromosomeOnReverse(self.TempRTFileName)
        self.OutHandle.close()
        
    def ReverseTranscribeChromosome2(self, FileName):
        """Here's my attempt to do this smartly if there is a cache overload
        break into multiple pieces, reverse each piece, then reverse the order
        
        0.1.2.3.4.5..6.7.8.9 - > (0.1.2), (3.4.5), (6,7,8), 9 -> (2.1.0), (5.4.3), (8,7,6), (9) -> (9), (8.7.6), (5.4.3). (2.1.0)
        
        """
        Handle = open(FileName, "rb")
        TempHandle = open(self.TempRTFileName, "wb")
        #RTString = ""
        MaxLength = 10000# 10k, should fit in cache, right?
        #string.reverse()
        #string.replace (A, t), etc.
        #string.toUpper()
        Buffer = ""
        KeyCounter = 0
        LetterCounter = 0
        SonicatedDNA = {} #a dictionary meant to hold the DNA in pieces, to make reversing quick
        for Line in Handle.xreadlines():
            Line = Line.strip()
            #print ":%s:"%Line
            if Line[0] == ">": # fasta header                
                continue
            #past all the protection, time for real work
            LetterCounter += len(Line)
            if len(Buffer) < MaxLength:                #add to string
                Buffer += Line
            else:            #it's big enough.  Store current, start over
                SonicatedDNA[KeyCounter] = Buffer
                Buffer = Line # repopulate with the string that would have made us go over
                KeyCounter += 1
        #now add the last Buffer
        SonicatedDNA[KeyCounter] = Buffer
        Handle.close()
        #now RT all pieces
        for (Key,DNA) in SonicatedDNA.items():
            Reversed = DNA[::-1]# some crazy hacking for reverse that I found online. apparently there is no built in reverse for python
            Reversed = Reversed.replace("A", "t")
            Reversed = Reversed.replace("T", "a")
            Reversed = Reversed.replace("C", "g")
            Reversed = Reversed.replace("G", "c")
            RT = Reversed.upper()
            SonicatedDNA[Key] = RT
        Keys = SonicatedDNA.keys()
        Keys.sort()
        Keys.reverse() 
        FullRTChromosome = ""
        for Key in Keys:
            FullRTChromosome += SonicatedDNA[Key]
        ## to make this a bit prettyier, and easier to parse on the other side, I'm going to put in \n every 70 chars as I write to file
        TempHandle.write(">ReversedChromosome\n")
        for Index in range(0, len(FullRTChromosome), 70):
            ToPrint = FullRTChromosome[Index:Index+70]
            TempHandle.write("%s\n"%ToPrint)
        
        TempHandle.close()

    def ReverseTranscribeChromosome(self, FileName):
        """I'm assuming a single-feature fasta formatted dna file
        #string.reverse()
        #string.replace (A, t), etc.
        #string.upper()
        """
        chromReader = bioseq.FastaReader(FileName)
        chromRevOut = bioseq.FastaOut(self.TempRTFileName)

        chromosome = bioseq.Sequence("ReversedChromosome")

        # For a multi fasta it will just concatenate all the seqs together
        for fastaEntry in chromReader:
            chromosome.seq += fastaEntry.seq

        #now Reverse and then transcribe
        ChromosomeSequenceInArray = list(chromosome.seq)
        ChromosomeSequenceInArray.reverse()
        Reversed = "".join(ChromosomeSequenceInArray)
        Reversed = Reversed.replace("A", "t")
        Reversed = Reversed.replace("T", "a")
        Reversed = Reversed.replace("C", "g")
        Reversed = Reversed.replace("G", "c")
        chromosome.seq = Reversed.upper()

        chromRevOut.write( chromosome )

        
    def TranslateChromosomeOnReverse(self, FileName ):
        """The reason for doing as a separate method is to get the numbering straight,
        and the numbering of the reverse is such a hassle.  So if the DNA looked like (on both strands
        i  123.456.789.10  13  16  19  22
        5' TTA.ATG.AAA.TTT.CCC.GGG.TAA.CAT
           AAT.TAC.TTT.AAA.GGG.CCC.ATT.GTA  5'
           
        The hard part is getting the numbering correct on the start base.  we are inclusive, meaning
        that on the reverse strand, first A of ATG is at index 24, and thus our starting nucleotide is 24
        - handy is the self.TotalDNAIndex, which was built on the forward strand, which has total number of
        nucleotides used in that. but it never counted the last two nucleotides (because they were part of the
        last readable codon.  so the total number of bases is self.TotalDNAIndex + 2, which is the base we start
        with when we translate back in frame 1
        """    
        Handle = open(FileName, "rb")
        MaxLength = 200# number of dna letters in the buffer before we send out
        #set here before we loop through any dna.  Don't reset while looping.  then you would kill all genes overlapping a buffer
        self.FrameStrings = ['','','']
        #put these here?
        self.FrameStarts = [ self.TotalDNAIndex + 2, self.TotalDNAIndex + 1, self.TotalDNAIndex ]
        
        self.TotalDNAIndex += 3
        
        self.translateStrand( Handle, '-', MaxLength )

        # now after the last buffer's been processed, if there are some framebuffers that have yet to be written out, let's do that too
        for i in 0,1,2:
            if len(self.FrameStrings[i]) >= self.MinProteinLength:
                FastaLine = self.fastaLine(i+1, '-')
                self.OutHandle.write("%s\n%s\n"%(FastaLine, self.FrameStrings[i]))
                self.ProteinCount += 1

        Handle.close()

    def TranslateDNAOnReverse(self, Buffer, StartInFrame):
        """Again, here the reason for doing a separate method is to get the counting straight
        """
        CurrentFrame = StartInFrame
        WriteToFile = 0
        #self.TotalDNAIndex += 3 # start one more than the real, so that the decrement at the beginning of the loop works right
        for Index in range(len(Buffer)-2):
            self.TotalDNAIndex -= 1 # I've consumed a letter, account for it
            #print "TotalDNAIndex %d, and letter %s"%(self.TotalDNAIndex, Buffer[Index])
            Codon = Buffer[Index:Index+3]
            AminoAcid = self.Table.get(Codon, "*") # if the codon is not there then return stop (most likely the NNNNN sequences)
            # check for end of protein, if so we print to file
            if AminoAcid == "*": # the terminator
                WriteToFile = 1
                
            frame = CurrentFrame % 3
            if frame == 0: # Frame 3 is last, an even multiple of 3
                frame = 3

            #print Index, self.TotalDNAIndes, Codon, AminoAcid, (CurrentFrame % 3) 
            if WriteToFile:
                self.writeFrame( frame, '-' )
                WriteToFile = 0 #I'm done writing, reset
            else:                     
                #simply add on
                self.FrameStrings[frame-1] += AminoAcid

            CurrentFrame += 1
        #now we want to return the last translated frame.  Since we increment at the end of the loop
        # we need to subtract one to get back to the frame that last worked
        return (CurrentFrame - 1) % 3
        
    def callTranslate( self, buffer, startInFrame, strand):
        lastTranslatedFrame = None
        if strand == '+':
            lastTranslatedFrame = self.TranslateDNA(buffer, startInFrame)
        elif strand == '-':
            lastTranslatedFrame = self.TranslateDNAOnReverse(buffer, startInFrame)
        else:
            raise ValueError
        return lastTranslatedFrame

    def translateStrand( self, handle, strand, maxLength ):
        Buffer = ''
        ## Using 1 based stuff!, because people talk about frames 1,2,3 and not the zero based
        startInFrame = 1
        for Line in handle.xreadlines():
            if not Line.strip():
                continue
            Line = Line.strip()
            #print ":%s:"%Line
            if Line[0] == ">": 
                # fasta header
                continue
            if len(Buffer) < maxLength:           #add to string
                Buffer += Line
            else:          #it's big enough.  Process it, strip it back down, then add the line
                LastTranslatedFrame = self.callTranslate( Buffer, startInFrame, strand )

                Buffer = Buffer[-2:] # get only the last two, because we translated the all up till that
                #now add the line onto the buffer 
                Buffer += Line
                """ yes I know that I subtracted one in the return and now I'm about to add one
                but this is the way that we can keep the terminology straight and know what
                we are talking about.
                """
                startInFrame = LastTranslatedFrame + 1  

        #now at the end of the file, we send in the last buffer
        self.callTranslate( Buffer, startInFrame, strand )

    def TranslateChromosome(self, FileName ):    
        Handle = open(FileName, "rb")
        MaxLength = 10000# number of dna letters in the buffer before we send out
                            
        self.translateStrand( Handle, '+', MaxLength )

        #after this last buffer's been sent, output any latent ORFs
        for i in 0,1,2:
            if len(self.FrameStrings[i]) >= self.MinProteinLength:
                self.OutHandle.write("%s\n%s\n"%(self.fastaLine(i+1,'+'), self.FrameStrings[i]))
                self.ProteinCount += 1

        #print "I consumed %d letters of DNA"%self.TotalDNAIndex
        Handle.close()

    def TranslateDNA(self, Buffer, StartInFrame):
        """Given some buffer, and the frame that I should start recording in, I translate the DNA
        The StartInFrame is used to tell me which frame I am already in, NOT used to offset things
        """
        CurrentFrame = StartInFrame
        WriteToFile = 0
        for Index in range(len(Buffer)-2):
            self.TotalDNAIndex += 1 # I've consumed a letter, account for it
            Codon = Buffer[Index:Index+3]
            AminoAcid = self.Table.get(Codon, "*") #if we don't find a codon, then return *, likely caused by NNNNN
            # check for end of protein, if so we print to file
            if AminoAcid == "*": # the terminator
                WriteToFile = 1
                
            #print Index, self.TotalDNAIndes, Codon, AminoAcid, (CurrentFrame % 3) 

            frame = CurrentFrame % 3
            if frame == 0: # Frame 3 is last, an even multiple of 3
                frame = 3

            #print Index, self.TotalDNAIndex, Codon, AminoAcid, (CurrentFrame % 3)
            if WriteToFile:
                self.writeFrame(frame,'+')
                WriteToFile = 0 #I'm done writing, reset
            else:                     
                #simply add on
                self.FrameStrings[frame-1] += AminoAcid

            CurrentFrame += 1
        #now we want to return the last translated frame.  Since we increment at the end of the loop
        # we need to subtract one to get back to the frame that last worked
        return (CurrentFrame - 1) % 3
            
    

    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:w:c:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-r":
                # -r results file(s)
                if not os.path.exists(Value):
                    print "** Error: couldn't find results file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.InputFile = Value
            elif Option == "-w":
                self.OutputFile = Value
            elif Option == "-c":
                self.ChromosomeName = Value
        if not OptionsSeen.has_key("-r") or not OptionsSeen.has_key("-w"):
            print UsageInfo
            sys.exit(1)

    def fastaLine(self,frame,strand):
        return ">Protein%s.Chr:%s.Frame%s.StartNuc%s.Strand%s" % ( self.ProteinCount,
                self.ChromosomeName, frame, self.FrameStarts[frame-1], strand )

    def writeFrame(self,frame,strand):
        idx = frame - 1
        trans = self.FrameStrings[idx]
        if len(trans) >= self.MinProteinLength:
            fastaLine = self.fastaLine( frame, strand )
            self.OutHandle.write("%s\n%s\n" % (fastaLine, trans))
            self.ProteinCount += 1
        #now reset, even if we don't print it out, we still reset
        
        self.FrameStrings[idx] = ""
        if strand == '+':
            self.FrameStarts[idx] = self.TotalDNAIndex + 3 # the next ORF from this frame will start at +3
        elif strand == '-':
            self.FrameStarts[idx] = self.TotalDNAIndex - 3 # the next ORF from this frame will start at -3
        else:
            raise ValueError("Invalid strand %s" % strand)


if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(psyco not found - running in non-optimized mode)"
    BeanCounter = AbacusClass()
    BeanCounter.ParseCommandLine(sys.argv[1:])
    BeanCounter.Main()                
