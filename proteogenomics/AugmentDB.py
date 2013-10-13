"""
Augment an existing trie-database with more sequences from a .fasta file.
Useful for 'needle-in-haystack' search where we add any missing sequences
to a large database like sprot.  Also useful if we need to add contaminants,
like porcine trypsin, to a standard species database.
"""
import os
import struct
import sys
import getopt


class Augmentor:
    def __init__(self):
        self.SourceDatabasePath = None
        self.AppendRecordsPath = None
        self.FinalOutputPath = None
    def LoadSourceDB(self, OriginalDBName):
        (Dir, FileName) = os.path.split(OriginalDBName)
        (Stub, Extension) = os.path.splitext(FileName)
        if not os.path.exists(OriginalDBName):
            Path = os.path.join("Database", OriginalDBName)
            if os.path.exists(Path):
                OriginalDBName = Path
            else:
                Path = os.path.join("Database", "%s.trie"%Stub)
                if os.path.exists(Path):
                    OriginalDBName = Path
                else:
                    print "* Error: Original database '%s' not found"%OriginalDBName
                    return
        print "Reading DB..."
        OriginalFile = open(OriginalDBName, "r")
        self.Text = OriginalFile.read()
        OriginalFile.close()
        print "Reading index..."
        OriginalIndexFile = open(os.path.join(Dir, "%s.index"%Stub), "rb")
        self.IndexText = OriginalIndexFile.read()
        OriginalIndexFile.close()
    def AugmentFromDir(self, Dir):
        for FileName in os.listdir(Dir):
            RootName = os.path.splitext(FileName)[0]
            Bits = RootName.split("_")
            Aminos = Bits[-1].replace("#","").replace("+","")
            print Aminos
            Pos = self.Text.find(Aminos)
            if Pos == -1:
                NewPos = len(Text)
                self.Text += "%s*"%Aminos
                print "Appended new peptide %s"%Aminos
                PepName = FileName[:40]
                self.IndexText += struct.pack("<qi80s", 0, NewPos, PepName)
            else:
                print "Already got %s"%Aminos
    def Augment(self, Sequence, Name):
        Pos = self.Text.find(Sequence)
        if Pos == -1:
            NewPos = len(self.Text)
            self.Text += "%s*"%Sequence
            print "Appended new peptide %s"%Name #Sequence
            self.IndexText += struct.pack("<qi80s", 0, NewPos, Name)
        else:
            #print "Already got %s"%Name
            pass
    def WriteAugmentedDB(self, Name):
        (Dir, FileName) = os.path.split(Name)
        if not Dir:
            Dir = "Database"
        (Stub, Extension) = os.path.splitext(FileName)
        if not Extension:
            FileName = "%s.trie"%FileName
        IndexFileName = "%s.index"%Stub
        NewFile = open(os.path.join(Dir, FileName), "w")
        NewFile.write(self.Text)
        NewFile.close()
        NewIndexFile = open(os.path.join(Dir, IndexFileName), "wb")
        NewIndexFile.write(self.IndexText)
        NewIndexFile.close()
    def AugmentWithFasta(self, FASTAName):
        File = open(FASTAName, "r")
        Sequence = ""
        Name = ""
        for FileLine in File.xreadlines():
            FileLine = FileLine.strip()
            if not FileLine:
                continue
            if FileLine[0] == ">":
                # finish any pending sequence
                if Sequence:
                    self.Augment(Sequence, Name)
                Sequence = ""
                Name = FileLine[1:]
            else:
                Sequence += FileLine
        if Sequence:
            self.Augment(Sequence, Name)
    def ParseCommandLine(self, Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:a:w:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-r":
                self.SourceDatabasePath = Value
            elif Option == "-a":
                self.AppendRecordsPath = Value
            elif Option == "-w":
                self.FinalOutputPath = Value
    def Main(self):
        self.LoadSourceDB(self.SourceDatabasePath)
        self.AugmentWithFasta(self.AppendRecordsPath)
        self.WriteAugmentedDB(self.FinalOutputPath)
        
if __name__ == "__main__":
    Bob = Augmentor()
    Bob.ParseCommandLine(sys.argv[1:])
    Bob.Main()
