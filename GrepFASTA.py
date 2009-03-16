"""
Look for a peptide in the IPI database.
"""
import sys
def Find(DesiredPeptide):
    HitCount = 0
    #File = open("Database\\IPIv315.fasta", "rb")
    File = open("Database\\IPI.fasta", "rb")
    CurrentSequence = ""
    CurrentProtein = ""
    for FileLine in File.xreadlines():
        if FileLine[0] == ">":
            Pos = CurrentSequence.find(DesiredPeptide)
            if Pos != -1:
                print "%s found at pos %s in %s"%(DesiredPeptide, Pos, CurrentProtein)
                Prefix = CurrentSequence[max(0, Pos - 10):Pos]
                Suffix = CurrentSequence[Pos + len(DesiredPeptide):Pos + len(DesiredPeptide) + 10]
                print "%s.%s.%s"%(Prefix, DesiredPeptide, Suffix)
                HitCount += 1
            CurrentProtein = FileLine[1:].strip()
            CurrentSequence = ""
        else:
            CurrentSequence += FileLine.strip()
    Pos = CurrentSequence.find(DesiredPeptide)
    if Pos != -1:
        print "%s found at pos %s in %s"%(DesiredPeptide, Pos, CurrentProtein)
        Prefix = CurrentSequence[max(0, Pos - 10):Pos]
        Suffix = CurrentSequence[Pos + len(DesiredPeptide):Pos + len(DesiredPeptide) + 10]
        print "%s.%s.%s"%(Prefix, DesiredPeptide, Suffix)
        HitCount += 1
    File.close()
    return HitCount

if __name__ == "__main__":
    DesiredPeptide = sys.argv[1]
    Find(DesiredPeptide)