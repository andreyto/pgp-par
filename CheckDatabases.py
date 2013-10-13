UsageInfo = """CheckDatabases.py
This is a quicky script to check the validity of an abbreviated database

Required Options
 -a [Trie] Database made from Summary.py
 -d [Trie file] Database

"""

import sys
import os
import getopt
import traceback
import SelectProteins


class CompileClass():
    def __init__(self):
        self.ReferenceDatabasePath = [] #possibly multiple
        self.AbbreviatedDatabasePath = []

        
    def Main(self):
        self.AbbreviatedProteins = SelectProteins.ProteinSelector()
        self.AbbreviatedProteins.LoadMultipleDB(self.AbbreviatedDatabasePath)
        # now get full database
        self.OriginalProteins = SelectProteins.ProteinSelector()
        self.OriginalProteins.LoadMultipleDB(self.ReferenceDatabasePath)
        self.VerifyConsistency()

    def VerifyConsistency(self):
        """Now we have two databases.  Let's check one against the other, in a full double loop
        """
        AllAbbrev =0
        SNAFUAbbrev = 0
        for (AbbrevID, AbbrevName) in self.AbbreviatedProteins.ProteinNames.items():
            AllAbbrev += 1
            AbbrevSequence = self.AbbreviatedProteins.ProteinSequences[AbbrevID]
            AbbrevAccession = self.GetAccessionFromFasta(AbbrevName)
            if AbbrevAccession[:3] == "XXX":
                continue
            OriginalMatches = self.OriginalProteins.FindPeptideLocations(AbbrevSequence)
            #print "AbbrevSequence %s has %d matches"%(AbbrevAccession, len(OriginalMatches))
            Verified = 0
            for (MatchID, ResidueNumber) in OriginalMatches:
                # Match is an ID in the OriginalProteins object
                OriginalName = self.OriginalProteins.ProteinNames[MatchID]
                OriginalAccession = self.GetAccessionFromFasta(OriginalName)
                if OriginalAccession == AbbrevAccession:
                    Verified = 1
                    break

            if not Verified:
                SNAFUAbbrev += 1
                print "Abbrev database entry %s did not find a match"%AbbrevAccession
        print "of %d total sequences, %s are screwed up"%(AllAbbrev, SNAFUAbbrev)

    def GetAccessionFromFasta(self, Header):
        """just get the thing up to the first space"""
        Bits = Header.strip().split(" ")
        return Bits[0]


    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "a:d:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-a":
                # -r results file(s)
                if not os.path.exists(Value):
                    print "** Error: couldn't find database file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.AbbreviatedDatabasePath.append(Value)
            if Option == "-d":
                if not os.path.exists(Value):
                    print "** Error: couldn't find database file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.ReferenceDatabasePath.append( Value)
        if not OptionsSeen.has_key("-a") or not OptionsSeen.has_key("-d"):
            print UsageInfo
            sys.exit(1)


if __name__ == "__main__":
    try:
        import psyco
        psyco.full()
    except:
        print "(psyco not found - running in non-optimized mode)"
    Gather = CompileClass()
    Gather.ParseCommandLine(sys.argv[1:])
    Gather.Main()
