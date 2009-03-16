import sys
import string
File = open(sys.argv[1], "r")
for FileLine in File.xreadlines():
    NewBits = [""]*80
    Bits = list(FileLine.strip().split())
    NewBits[0] = Bits[0]
    for Index in range(1, len(Bits)):
        ColonPos = Bits[Index].find(":")
        try:
            FeatureNumber = int(Bits[Index][:ColonPos])
            NewBits[FeatureNumber] = Bits[Index][ColonPos+1:]
        except:
            continue
        
    print string.join(NewBits, "\t")
    