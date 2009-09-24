UsageInfo = """ FabricateProteome.py
Makes a database of fabricated proteins, by emitting proteins of a random 
sequences (according to some average amino acid frequency) of length 100 to 
1000 amino acids It also makes a database of proteins with no-hydrophobic 
amino acids.

Required Options
 -d [TrieFile] Database(s) as seed 
"""


import os
import getopt
import sys
import random
import struct
import SelectProteins

class LiarClass:
    def __init__(self):
        self.OutputFileName = "RenameYourOutput.txt"
        self.DBPath = [] # for multiple databases
        self.AACountTable = {} 
        self.AAFrequencyTable = {}
        self.AllAminoAcids = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
        self.NonAminoAcids = ["B", "J", "O", "U", "X", "Z"] # just for error checking
        self.CumulativeProbabilities = [] #holds the amino acid frequencies, so we can get letters easily
        for AminoAcid in self.AllAminoAcids:
            self.AACountTable[AminoAcid] = 0
            self.AAFrequencyTable[AminoAcid] = 0.0 
        
    def Main(self):
        self.ProteinPicker = SelectProteins.ProteinSelector()
        self.ProteinPicker.LoadMultipleDB(self.DBPath)
        #now I've got a big long dict of multiple databases (potentially).
        # Let's fabricate one.  First a composite, then one without hydrophobic residues
        # I think that 20,000 proteins should do, with somewhere between 100 and 1000 amino acids
        self.CalculateAminoAcidFrequencies()
        self.OutputFileName = "AverageProteinFrequencies.trie"
        self.FabricateDatabase(self.OutputFileName)
        print "Now for the no-hydrophobic"
        for AminoAcid in self.AllAminoAcids:
            self.AACountTable[AminoAcid] = 0
            self.AAFrequencyTable[AminoAcid] = 0.0 
        HydrophobicResidues = ["G", "A", "P", "V", "L", "I", "M", "F", "Y", "W"]
        self.CalculateAminoAcidFrequencies(HydrophobicResidues)
        self.OutputFileName = "NoHydrophobic.trie"
        self.FabricateDatabase(self.OutputFileName)
        print "now for the low-hydrophobic"
        for AminoAcid in self.AllAminoAcids:
            self.AACountTable[AminoAcid] = 0
            self.AAFrequencyTable[AminoAcid] = 0.0 
        self.CalculateAminoAcidFrequencies(HydrophobicResidues, 0.5)
        self.OutputFileName = "LowHydrophobic.trie"
        self.FabricateDatabase(self.OutputFileName)
        
        
        
    def FabricateDatabase(self, FileName):
        """Going to make 20,000 proteins under this name
        """
        print "Making output %s"%FileName
        NewIndexPath = os.path.splitext(FileName)[0] + ".index"
        self.OutputTrieFile = open(FileName, "wb")
        self.OutputIndexFile = open(NewIndexPath, "wb")
        NumProteins = 20000
        for ProteinIndex in range(NumProteins):
            NewSequence = self.FabricateSingleProtein()
            self.WriteProtein(NewSequence, ProteinIndex)
        self.OutputTrieFile.close()
        self.OutputIndexFile.close()
        

    def WriteProtein(self, Sequence, ProteinIndex):
        """
        Given a protein sequence, and protein index number (for looking up the name),
        write a scrambled or reversed record to the output database.  (And write the
        original, if the -b flag was specified)
        """
        ShuffledProtein = Sequence
        OutputFilePos = self.OutputTrieFile.tell()
        self.OutputTrieFile.write(ShuffledProtein)
        self.OutputTrieFile.write("*")
        ShuffledName = "XXX.%s"%ProteinIndex
        Block = struct.pack("<qi80s", 0, OutputFilePos, ShuffledName)
        self.OutputIndexFile.write(Block)

        
    def FabricateSingleProtein(self):
        """Using the amino acid frequencies (actually the cumulativeProbabilities table) we
        create totally fake amino acid sequences
        """
        random.seed()
        RandomNum = random.random()
        ProteinLength = int(RandomNum * 1000 )# to getsequences between 100 and 1000 amino acids long
        ProteinSequence = ""
        for Index in range (ProteinLength):
            RandomNum = random.random()
            # now go through to cum probabilitiy dist and get the letter
            
            for AAIndex in range(len(self.AllAminoAcids)):
                if RandomNum < self.CumulativeProbabilities[AAIndex]:
                    #emit letter
                    ProteinSequence += self.AllAminoAcids[AAIndex]
                    break
        return ProteinSequence
        
        
        
    def CalculateAminoAcidFrequencies(self, BlackoutLetters = [], BlackoutLevel = 0.0):
        """from the database sequence we count frequencies (blacking out letters as desired)
        You have the option of probabilistically excluding blackout letters (like say half 
        of the time) making them more of greyout than blackout
        """
        TotalLetters = 0
        for Sequence in self.ProteinPicker.ProteinSequences.values():
            for Letter in Sequence:
                if Letter in BlackoutLetters: 
                    #here we test to see if we want to include this letter in the count
                    #this is designed to exclude it XX% of the time.  If you want these
                    #letters only 20% of the time, then have blackoutlevel be 0.2
                    RandomNum = random.random()
                    if RandomNum > BlackoutLevel:
                        continue
                if Letter in self.NonAminoAcids:
                    continue
                self.AACountTable[Letter] += 1
                TotalLetters += 1
        #Now get Frequencies
        RunningTotal = 0.0 # probabilities
        Iterater =0
        for Letter in self.AllAminoAcids:
            Frequency = self.AACountTable[Letter] / float(TotalLetters)
            self.AAFrequencyTable[Letter] = Frequency
            print "AA: %s \t %.3f"%(Letter, Frequency)
            #print "%.3f"%Frequency
            RunningTotal += Frequency
            self.CumulativeProbabilities.append( RunningTotal)
            Iterater += 1

    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "d:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-d":
                if not os.path.exists(Value):
                    print "** Error: couldn't find results file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.DBPath.append(Value)
        if not OptionsSeen.has_key("-d"):
            print UsageInfo
            sys.exit(1)
            

if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(psyco not found - running in non-optimized mode)"
    Abagnale = LiarClass()
    Abagnale.ParseCommandLine(sys.argv[1:])
    Abagnale.Main()                