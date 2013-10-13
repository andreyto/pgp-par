import traceback
import sys

try:
    import PyInspect
    Spectrum = PyInspect.Spectrum("ShewPTM\\Consensus\\K.DFSQIDNAP+16EER.E.mgf")
    Spectrum.ScorePeptide("K.DFSQIDNAP+16EER.E")
except:
    traceback.print_exc()

print "** Test complete.  Press >>>ENTER<<<"
sys.stdin.readline()