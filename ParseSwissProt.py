"""
ParseSwissProt.py
This script runs through a .dat file from swiss sprot and looks for PTM annotations
Swiss Prot is actually a pain to process, because the mass deltas for each PTM are
stored in a separate file. Additionally, there are PTMs stored under a variety of
keywords (parse hassle).

Swiss Prot style Annotations
1. feature (FT), modification (MOD_RES), crosslink (CROSSLNK), or lipid (LIPID)
2. shows residue number (start, end)
3. Does NOT have mol weights listed in this file!!!
#eg.
FT   MOD_RES      53     53       Phosphotyrosine.
FT   MOD_RES     130    130       N6,N6-dimethyllysine (Potential).
FT   MOD_RES      83     83       Sulfotyrosine (Potential).
FT   CROSSLNK    119    119       Glycyl lysine isopeptide (Lys-Gly)  ### ubiquitin
"""
import getopt
import os
import sys



class FinderClass:
    def __init__ (self):
        "constructor"
        self.SwissProtFileName = None
        self.SwissProtPTMFileName = None
        self.OutputFileName = None
        self.Count =0
        self.DictyFlag =0
        self.ModMasses = {}

    def ParsePTMList(self):
        """ parses out the useful information in the swiss prot PTM file.
        PTMName should match the string at the end of the the "FT" feature line.
        """
        Handle = open(self.SwissProtPTMFileName, "rb")
        for Line in Handle.xreadlines():
            if Line[0] == '#':
                continue # skip commented lines
            Bits = Line.strip().split("\t")
            try:
                PTMName = Bits[0] #not really the nice name, but at least it's unique
                DeltaMassFloat = float(Bits[5])
                DeltaMass = int(round(DeltaMassFloat))
            except:
                continue
            self.ModMasses[PTMName] = DeltaMass
        
    def ParseFile (self):
        """
        This function takes the whole file and starts to parse it.
        This one only parses out a record, and then sends that on
        for further work
        """
        print "Opening %s"%self.SwissProtFileName
        FileHandle = open(self.SwissProtFileName, "rb")
        OutHandle = open(self.OutputFileName, "wb")
        MEG = 1024*1024
        Text = ""
    
        ## Process a block of data on each pass
        Block = FileHandle.read(MEG)
        if not Block:
            print "no text in input file"
            return
        Text = Block
        Pos = Text.find("ID   ",0) # get the first start

        while 1:
            NextPos = Text.find("\nID   ", Pos+1) 
            if NextPos == -1:
                TextRemainder = Text[Pos:]
                Block = FileHandle.read(MEG)
                if not Block: # reached the end of the file
                    Record = Text[Pos:]
                    self.ParseRecordForMod(Record, OutHandle)
                    break
                Text = TextRemainder + Block
                Pos =0
                continue
            #I've found the start of the next one, which means that I
            #can get all of this one
            Record = Text[Pos:NextPos]
            self.ParseRecordForMod(Record, OutHandle)
            Pos = NextPos
        FileHandle.close()
        OutHandle.close()
        print "I got %d records"%self.Count

    def ParseRecordForMod (self,Text, FileHandle):
        """
        Looks at a single record.  Determine whether it contains a PTM I care about
        FileHandle is an output file
        """
        HasPTM =0
        if Text.find("MOD_RES") >= 0:
            HasPTM =1
        if Text.find("interchain with G-Cter in") >= 0: #ubiquitin, sumo containing crosslinks
            HasPTM =1
        if Text.find("FT   LIPID") >= 0:
            HasPTM =1
        if not HasPTM:
            return
        if self.DictyFlag:
            if Text.find("_DICDI") < 0:
                return
        #now process it
        #1. Get Basic Info
        self.Count +=1
        FirstLineStart = Text.find("ID   ")
        FirstLineEnd = Text.find("\n",FirstLineStart+1)
        FirstLine = Text[FirstLineStart:FirstLineEnd]
        AcessionLineStart = Text.find("\nAC   ")
        AcessionLineEnd = Text.find("\n",AcessionLineStart +1)
        AcessionLine = Text[AcessionLineStart+1: AcessionLineEnd]
        NameLineStart = Text.find("\nDE   ")
        NameLineEnd = Text.find("\n", NameLineStart+1)
        NameLine = Text[NameLineStart+1 :NameLineEnd] #plus one because of the \n at the start of the find
        SequenceStart = Text.find("\nSQ   SEQUENCE")
        if SequenceStart == -1:
            print "no sequence with record."
            return
        #SequenceEnd
        Sequence = Text[SequenceStart + 1:] #assume that sequence ends at the end of the record
        #2. Get PTM Lines
        ModResPos = Text.find("FT   MOD_RES")#one type of PTM keyword
        ModResString = ""
        while not (ModResPos == -1):
            ModResLineEnd = Text.find("\n", ModResPos +1)
            ModResLine = Text[ModResPos:ModResLineEnd]
            #print ("Line: %s"%ModResLine)
            ModResString += "%s\n"%ModResLine
            ModResPos = Text.find("FT   MOD_RES", ModResPos +1)
        #unbiquitin and Sumoylation are annotated as crosslinks, a category which
        #is largely useless to me.  So I do some hackish parsing
        UbiquitinString = ""
        UbiquitinPos = Text.find ("interchain with G-Cter in ") #unique string for ubiquitin and sumo
        while not (UbiquitinPos == -1):
            UbiquitinLineEnd = Text.find("\n", UbiquitinPos +1)
            UbiquitinLineStart = Text.rfind("FT   CROSSLNK", 0, UbiquitinPos) #find last before the ubiquitin
            UbiquitinLine = Text[UbiquitinLineStart:UbiquitinLineEnd]
            #if UbiquitinLine.find("SUMO") == -1 and UbiquitinLine.find("ubiquitin") == -1:
            #    print "This is a rogue line \n%s\n\n"%UbiquitinLine
            UbiquitinPos = Text.find ("interchain with G-Cter in ", UbiquitinPos +1) 
            UbiquitinString += "%s\n"%UbiquitinLine
            #print "This is a ubiquitin line\n%s\n\n"%UbiquitinLine
        #processing Lipid additions
        LipidPos = Text.find ("FT   LIPID") #final kind of PTM i care about
        LipidString = ""
        while not (LipidPos == -1):
            LipidLineEnd = Text.find("\n", LipidPos +1)
            LipidLine = Text[LipidPos:LipidLineEnd]
            LipidString += "%s\n"%LipidLine
            LipidPos = Text.find("FT   LIPID", LipidPos +1)
        #printing out all the data gathered. 
        FileHandle.write("%s\n"%FirstLine)
        FileHandle.write("%s\n"%AcessionLine)
        FileHandle.write("%s\n"%NameLine)
        FileHandle.write("%s"%ModResString)
        FileHandle.write("%s"%UbiquitinString)
        FileHandle.write("%s"%LipidString)
        FileHandle.write("%s\n"%Sequence)

        ###############################################################
        ### instead of printing to a file, we would just store things
        ### in whatever kind of struct you want
        ################################################################

    def ParseCommandLine(self, Arguments):
        (Options, Args) = getopt.getopt(sys.argv[1:], "f:g:o:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-f":
                # -f SwissProtFileName file(s)
                if not os.path.exists(Value):
                    print "** Error: couldn't find swissprot file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.SwissProtFileName = Value
            elif Option == "-g":
                if not os.path.exists(Value):
                    print "** Error: couldn't find swissprot file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.SwissProtPTMFileName = Value
            elif Option == "-o":
                #-o output
                self.OutputFileName = Value
            elif Option == "-d":
                self.DictyFlag =1
        # Error out, if we didn't see required options:
        if not OptionsSeen.has_key("-f") or not OptionsSeen.has_key("-o") or not OptionsSeen.has_key("-g"):
            print "*******************************************************************"
            print "** Error: Please specify swiss prot database file (-f) and ptm file (-g) and output file (-o)\n"
            sys.exit(1)


if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(psyco not found - running in non-optimized mode)"
    Finder = FinderClass()
    Finder.ParseCommandLine(sys.argv[1:])
    Finder.ParsePTMList()
    Finder.ParseFile()        
    