UsageInfo = """ AddFixedMod.py
Takes an input file(s) and adds a fixed modification
explicitely to the annotation column.  This is most
useful for things like itraq, or silac labeling


Required Options
 -r [FileName] Filename or directory of annotations
 -a [AminoAcid] the amino acid with a modification
 -m [Mass] Mass of the modification
"""


import os
import getopt
import sys
import ResultsParser
from Utils import *
Initialize()

class AssassinClass(ResultsParser.ResultsParser):
    def __init__(self):
        self.InputFile = None
        self.OutputFile = None
        self.Mass = 0
        self.AminoToAlter = None
        ResultsParser.ResultsParser.__init__(self)
        
    def Main(self):
        self.ProcessResultsFiles(self.InputFile, self.ParseAnnotations)

        
    def ParseAnnotations(self, FileName):
        Handle = open(FileName, "rb")
        OutFilePath = "%s.fixed"%FileName
        OutHandle = open (OutFilePath, "wb")
        for Line in Handle.xreadlines():
            if not Line.strip():
                continue
            Line = Line.strip()
            if Line[0] == "#":
                OutHandle.write("%s\n"%Line)
                continue
            Bits = Line.split("\t")
            Annotation = Bits[self.Columns.Annotation]
            NewAnnotation = self.MakeNewAnnotation(Annotation)
            Bits[self.Columns.Annotation] = NewAnnotation
            ToPrint = "\t".join(Bits)
            OutHandle.write("%s\n"%ToPrint)
            #break
        Handle.close()
        OutHandle.close()


    def MakeNewAnnotation(self, OldAnnotation):
        """Given an annotation and some modification to make
        we change it
        """
        Peptide = GetPeptideFromModdedName(OldAnnotation)
        Length = len(Peptide.Aminos)
        for Index in range(Length):
            #print "checking %s, %s"%(Index, Peptide.Aminos[Index])
            if Peptide.Aminos[Index] == self.AminoToAlter:
                ## make the mod object
                if self.Mass > 0:
                    ModName = "+%s"%self.Mass
                else:
                    ModName = "-%s"%self.Mass
                Mod = PTModClass(ModName)
                Mod.Mass = self.Mass
                
                ##put it in the peptide
                if not Peptide.Modifications.has_key(Index):
                    Peptide.Modifications[Index] = []
                Peptide.Modifications[Index].append(Mod)
        ## now a new loop to check things out
        #print "Checking things out after"
        #for AminoIndex in Peptide.Modifications.keys():
        #    ModArray = Peptide.Modifications[AminoIndex]
        #    print "I'm a mod at %s"%AminoIndex
        #    for Mod in ModArray:
        #        print "My Name is :%s:"%Mod.Name
        #        print "My Mass is :%s:"%Mod.Mass
        
        return Peptide.GetFullModdedName()

                

    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:a:m:")
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
            elif Option == "-m":
                self.Mass = int (Value)
            elif Option == "-a":
                self.AminoToAlter = Value
        if not OptionsSeen.has_key("-r") :
            print UsageInfo
            sys.exit(1)
            

if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(psyco not found - running in non-optimized mode)"
    Guido = AssassinClass()
    Guido.ParseCommandLine(sys.argv[1:])
    Guido.Main()                