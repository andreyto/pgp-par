UsageInfo = """PepXMLToInspect.py
This script takes the output of Mascot in XML and converts it into
the output format of inspect.  For the MQScore column, it puts in the ion score from Mascot.
Where Inspect has values that Mascot doesn't have, it puts "Dummy".
Only the top protein of each scan is included in the output.

FileName\tScan\tAnnotation\tProtein\tCharge\tMQScore .... \tP-value ... byteOffset

Required Parameters
-r [FileName] - Mascot output file in XML
-w [FileName] - Output of this program
-m [FileName] - File mgf (for scan number and byte offset)

"""

import sys
import os
import getopt
import re
import GetByteOffset
from xml.sax import make_parser
from xml.sax.handler import feature_namespaces
from xml.sax import ContentHandler

class ParseMascotXMLClass(ContentHandler):
    def __init__(self):
        self.InputFilePath = None
        self.OutputFilePath = None
        self.MGFFile = None
        self.InputFile = None
        self.OutputFile = None
        self.ScanOffset = {}
        self.ScanDict= {}
        # for parsing
        self.SpectrumFile = None
        self.ScanNum = None
        self.Annotation = None
        self.CompleteAnnotation = None      # with prefix and suffix
        self.Description = None     # spectrum description
        self.Prefix = None
        self.Suffix = None
        self.Protein = None
        self.Charge = None
        self.Score = None
        self.PValue = None
        self.TopProteinp = False
        self.ByteOffset = 'Dummy'
        self.ScanText = ""
        
    def Start(self):
        self.InputFile = open(self.InputFilePath, "r")
        self.OutputFile = open (self.OutputFilePath, "w")
        
    def startElement(self, name, attrs):
        if name == 'parameter' and attrs.get('name', None) == "FILE":
            self.SpectrumFile = attrs.get('value', None)

        if name == 'spectrum_query':
            # reset variables
            self.ScanNum = None
            self.Annotation = None
            self.CompleteAnnotation = None      # with prefix and suffix
            self.Description = None     # spectrum description
            self.Prefix = None
            self.Suffix = None
            self.Protein = None
            self.Charge = None
            self.Score = None
            self.PValue = None
            self.TopProteinp = False
            self.ByteOffset = 'Dummy'
            self.ScanText = ""
            
            self.Description = attrs.get('spectrum', None).strip()            
            self.Charge = attrs.get('assumed_charge', None)
            # get scan number from dictionary using spectrum description
            self.ScanNum = self.ScanDict[self.Description]
            
        if name == 'search_hit':
            if attrs.get('hit_rank', None) == "1":    # top ranking protein 
                self.TopProteinp = True
                self.Annotation = attrs.get('peptide', None)
                self.Prefix = attrs.get('peptide_prev_aa', None)
                self.Suffix = attrs.get('peptide_next_aa', None)
                self.Protein = attrs.get('protein', None)

                if self.Prefix == '-':
                    self.Prefix = '*'
                if self.Suffix == '-':
                    self.Suffix = '*'
                self.CompleteAnnotation = self.Prefix + "." + self.Annotation + "." + self.Suffix

        if name == 'search_score' and self.TopProteinp:
            if attrs.get('name', None) == 'ionscore':
                self.Score = attrs.get('value', None)
            elif attrs.get('name', None) == 'expect':
                self.PValue = attrs.get('value', None)
         
    def endElement(self, name):
        if name == 'spectrum_query' and self.Protein != None:    # this query has a hit
            # write hit to output file
            self.WriteLine()

        if name == 'search_hit':
            self.TopProteinp = False
              
    def ReadMGF(self):
        """Parses an individual MGF file and saves the scan num and description
        into an dictionary called self.ScanDict
        """

        FileName = self.MGFFile
        self.ScanDict = {} #reset for the new file
        print "Opening MGF file %s"%FileName
        File = open(FileName, "r")

        # read one line at a time
        for line in File:
            p = re.compile('^TITLE=([^\n]*)')
            m = p.match(line)
            if m != None:   # start with TITLE=
                title = m.group(1)
                
            p = re.compile('^SCANS=(\d*)')
            m = p.match(line)
            if m != None:   # start with SCANS=
                scan = m.group(1)
                # store title,scan pair to dictionary
                self.ScanDict[title.strip()] = scan
                #print "adding %s, %s to dict" %(title,scan)
        
        print "Searching %d scans from %s"%(len(self.ScanDict), FileName)
        File.close()

    def WriteLine(self):
        if self.ScanOffset != {}:
            self.ByteOffset = self.ScanOffset[int(self.ScanNum)]
        # write one line of Inspect format to output file
        self.OutputFile.write("%s\t%s\t%s\t%s\t%s\t%s\tDummy\tDummy\tDummy\tDummy\tDummy\tDummy\tDummy\t%s\tDummy\tDummy\tDummy\tDummy\tDummy\t%s\n"
            %(self.SpectrumFile,self.ScanNum,self.CompleteAnnotation,self.Protein,self.Charge,self.Score,self.PValue,self.ByteOffset))

    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "r:w:m:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-r":
                # -r results file(s)
                if not os.path.exists(Value):
                    print "** Error: couldn't find results file '%s'\n\n"%Value
                    print UsageInfo
                    sys.exit(1)
                self.InputFilePath = Value
            if Option == "-w":
                self.OutputFilePath = Value
            if Option == "-m":
                self.MGFFile = Value
        if not OptionsSeen.has_key("-r") or not OptionsSeen.has_key("-w") or not OptionsSeen.has_key("-m"):
            print UsageInfo
            sys.exit(1)

    def Finish(self):   
        self.InputFile.close()
        self.OutputFile.close()

if __name__ == '__main__':
    try:
        import psyco
        psyco.full()
    except:
        print "(psyco not found - running in non-optimized mode)"
    Fix = ParseMascotXMLClass()
    Fix.ParseCommandLine(sys.argv[1:])
    if Fix.MGFFile != None:
        Fix.ReadMGF() # get scan #
        # get byte offsets from MGF
        Abacus = GetByteOffset.Abacus()
        Fix.ScanOffset = Abacus.GetByteOffset(Fix.MGFFile)
    
    Fix.Start()
    
    # Create a parser
    parser = make_parser()

    # Tell the parser we are not interested in XML namespaces
    parser.setFeature(feature_namespaces, 0)

    # Tell the parser to use our handler
    parser.setContentHandler(Fix)

    # Parse the input
    parser.parse(Fix.InputFile)

    # finish up
    Fix.Finish()
