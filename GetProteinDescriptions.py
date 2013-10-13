"""
GetProteinDecriptions:
Because the IPI database records have uninformative names, we go to the IPI website and retrieve their
descriptions.
"""
import os
import sys
import urllib

File = open(sys.argv[1], "r")
for FileLine in File.xreadlines():
    Bits = FileLine.strip().split("\t")
    if len(Bits)<4:
        continue
    try:
        ProteinID = int(Bits[0])
    except:
        continue
    ProteinName = Bits[1]
    IPIID = ProteinName.split("|")[0][4:].split(".")[0] # barftacular.
    URL = "http://srs.ebi.ac.uk/srsbin/cgi-bin/wgetz?-e+[IPI-acc:%s]+-vn+2"%IPIID
    TargetFile = "IPI\\%s.txt"%IPIID
    urllib.urlretrieve(URL, TargetFile)
    IPIFile = open(TargetFile, "r")
    for IPILine in IPIFile.xreadlines():
        if IPILine[:3] == "DE ":
            print "%s\t%s\t%s"%(ProteinName, IPIID, IPILine[3:].strip())
    IPIFile.close()
    