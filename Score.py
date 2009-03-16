"""
Score.py: Simple wrapper for inspect scoring
"""
import sys
import string
import PyInspect

def FormatTuple(Tuple):
    Str = "("
    for Entry in Tuple:
        Str += "%.4g, "%Entry
    Str = Str[:-2]
    Str += ")"
    return Str

ColonBits = sys.argv[1].split(":")
try:
    FileOffset = int(ColonBits[-1])
    FileName = string.join(ColonBits[:-1], ":")
except:
    FileName = sys.argv[1]
    FileOffset = 0
    
Spectrum = PyInspect.Spectrum(FileName, FileOffset)
#Result = Spectrum.ScorePeptideDetailed(sys.argv[2])
Result = Spectrum.ScorePeptideDetailed(sys.argv[2])
Str = "MQ %.4f %s"%(Result[0], FormatTuple(Result[1:]))
print Str