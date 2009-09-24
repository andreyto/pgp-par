UsageInfo = """InspectToPepXML.py
This script takes the output of Inspect and converts it into
the format of PepXML (based on Mascot XML output format).  

Required Parameters
-r [Directory] - Directory containing the outputs file of Inspect 
-w [Directory] - Directory of converted file in PepXML

Optional
-m [Directory] - Directory of spectra files
    (for extracting spectrum titles from mgf files)
"""

import sys
import os
import getopt
import re
import GetByteOffset
import ResultsParser
from Utils import *

Initialize()

class InspectToPepXMLClass(ResultsParser.ResultsParser):
    def __init__(self):
        self.InputFileDir = None
        self.OutputFileDir = None
        self.SpectraDir = None
        self.ScanOffset = {}
        self.ScanDict= {}
        self.SpectrumFileType = ""
        self.SpectrumFileBase = ""
        ResultsParser.ResultsParser.__init__(self)
        
    def Main(self):
        self.ProcessResultsFiles(self.InputFileDir, self.ConvertInspectToPepXML)

    def ConvertInspectToPepXML(self, FilePath):
        """ convert a single raw Inpsect output file to PepXML """

        InspectHandle = open(FilePath, "rb")
        FileName = os.path.split(FilePath)[1]
        FileName = FileName.replace(".txt", ".xml")
        NewPath = os.path.join(self.OutputFileDir, FileName)
        PepXMLHandle = open(NewPath, "wb")

        NumQueries = 0
        NumHits = 0
        LastScanNumber = -1
        for Line in InspectHandle.xreadlines():
            if Line[0] == "#":  # comments
                continue
            Bits = list(Line.split("\t"))

            # fields of the entry
            SpectrumFilePath = Bits[self.Columns.SpectrumFile]
            SpectrumFileName = os.path.split(SpectrumFilePath)[1]
            SpectrumFilePath = os.path.join(self.SpectraDir, SpectrumFileName)
            self.SpectrumFileBase = SpectrumFilePath.replace(os.path.splitext(SpectrumFilePath)[1], "")
            self.SpectrumFileType = os.path.splitext(SpectrumFilePath)[1]
            FileOffset = int(Bits[self.Columns.FileOffset])

            # get info from spectrum file
            SpectrumTitle = ""
            PepMass = ""
            if os.path.exists(SpectrumFilePath):
                (Stub, Ext) = os.path.splitext(SpectrumFilePath)
                if Ext.lower() == ".mzxml":
                    self.SpectrumFileType = ".mzXML"
                elif Ext.lower() == ".mgf":
                    self.SpectrumFileType = ".mgf"
                    (SpectrumTitle, PepMass) = self.GetSpectrumInfoFromMGF(SpectrumFilePath, FileOffset)

            if LastScanNumber == -1: # first line
                self.WritePepXMLOpening(PepXMLHandle, NewPath)
                
            ScanNumber = Bits[self.Columns.ScanNumber]
            Annotation = Bits[self.Columns.Annotation]
            Prefix = Annotation[0]
            Peptide = Annotation[2:-2]
            Suffix = Annotation[-1]
            Protein = Bits[self.Columns.ProteinName]
            Charge = Bits[self.Columns.Charge]
            MQScore = Bits[self.Columns.MQScore]
            FScore = Bits[self.Columns.FScore]
            DeltaScore = Bits[self.Columns.DeltaScore]
            PValue = Bits[self.Columns.PValue]
            ProteinID = Bits[self.Columns.ProteinID]

            if LastScanNumber != ScanNumber:    # result for a new spectrum
                if LastScanNumber != -1:            # write closing tags for results for the last spectrum
                    PepXMLHandle.write('</search_result>\n')
                    PepXMLHandle.write('</spectrum_query>\n')
                NumHits = 0
                NumQueries = NumQueries + 1
                Query = '<spectrum_query spectrum="%s" start_scan="%s" end_scan="%s" precursor_neutral_mass="%s" assumed_charge="%s" index="%s">\n'%(SpectrumTitle,ScanNumber,ScanNumber,"",Charge,NumQueries)
                PepXMLHandle.write(Query)
                PepXMLHandle.write('<search_result search_id="1">\n')

            LastScanNumber = ScanNumber
            NumHits = NumHits + 1
            CalcMass = ""   #GetMass(Peptide)
            MassDiff = ""   #CalcMass - float(PepMass)
            Hit = '<search_hit hit_rank="%s" peptide="%s" peptide_prev_aa="%s" peptide_next_aa="%s" protein="%s" num_tot_proteins="%s" num_matched_ions="%s" calc_neutral_pep_mass="%s" massdiff="%s" num_missed_cleavages="%s" is_rejected="%s" protein_descr="%s" protein_mw="%s">\n'%(NumHits,Peptide,Prefix,Suffix,Protein,'','',CalcMass,MassDiff,'','','','')
            PepXMLHandle.write(Hit)
            PepXMLHandle.write('<search_score name="mqscore" value="%s"/>\n'%MQScore)
            PepXMLHandle.write('<search_score name="expect" value="%s"/>\n'%PValue)
            PepXMLHandle.write('<search_score name="fscore" value="%s"/>\n'%FScore)
            PepXMLHandle.write('<search_score name="deltascore" value="%s"/>\n'%DeltaScore)
            PepXMLHandle.write('</search_hit>\n')

        PepXMLHandle.write('</search_result>\n')
        PepXMLHandle.write('</spectrum_query>\n')
        self.WritePepXMLClosing(PepXMLHandle)
        InspectHandle.close()
        PepXMLHandle.close()

    def WritePepXMLOpening(self, PepXMLHandle, PepXMLFilePath):
        PepXMLHandle.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        PepXMLHandle.write('<msms_pipeline_analysis xmlns="http://regis-web.systemsbiology.net/pepXML" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://regis-web.systemsbiology.net/pepXML http://mascot1/mascot/xmlns/schema/pepXML_v18/pepXML_v18.xsd" date="" summary_xml="%s">\n'%PepXMLFilePath)
        PepXMLHandle.write('<msms_run_summary base_name="%s" raw_data_type="raw" raw_data="%s">\n'%(self.SpectrumFileBase, self.SpectrumFileType))
        PepXMLHandle.write('<sample_enzyme name="trypsin" description="" fidelity="" independent="">\n')
        PepXMLHandle.write('<specificity sense="" min_spacing="" cut="KR" no_cut="P" sense="C"/>\n')
        PepXMLHandle.write('</sample_enzyme>\n')
        PepXMLHandle.write('<search_summary base_name="" search_engine="Inspect" precursor_mass_type="monoisotopic" fragment_mass_type="monoisotopic" out_data_type="out" out_data=".txt" search_id="1">\n')
        PepXMLHandle.write('<search_database local_path="" database_name="" database_release_identifier="" size_in_db_entries="" size_of_residues="" type="AA"/>\n')
        PepXMLHandle.write('<enzymatic_search_constraint enzyme="Trypsin" max_num_internal_cleavages="0" min_number_termini="2"/>\n')
        PepXMLHandle.write('<parameter name="CHARGE" value="2+ and 3+"/>\n')
        PepXMLHandle.write('<parameter name="CLE" value="Trypsin"/>\n')
        PepXMLHandle.write('<parameter name="DB" value=""/>\n')
        PepXMLHandle.write('<parameter name="FILE" value=""/>\n')
        PepXMLHandle.write('<parameter name="FORMAT" value=""/>\n')
        PepXMLHandle.write('<parameter name="FORMVER" value=""/>\n')
        PepXMLHandle.write('<parameter name="INSTRUMENT" value="ESI-QUAD-TOF"/>\n')
        PepXMLHandle.write('<parameter name="ITOL" value=""/>\n')
        PepXMLHandle.write('<parameter name="ITOLU" value="Da"/>\n')
        PepXMLHandle.write('<parameter name="MASS" value="Monoisotopic"/>\n')
        PepXMLHandle.write('<parameter name="REPORT" value=""/>\n')
        PepXMLHandle.write('<parameter name="REPTYPE" value="Peptide"/>\n')
        PepXMLHandle.write('<parameter name="RULES" value=""/>\n')
        PepXMLHandle.write('<parameter name="SEARCH" value=""/>\n')
        PepXMLHandle.write('<parameter name="TAXONOMY" value=""/>\n')
        PepXMLHandle.write('<parameter name="TOL" value=""/>\n')
        PepXMLHandle.write('<parameter name="TOLU" value="Da"/>\n')
        PepXMLHandle.write('</search_summary>\n')

    def WritePepXMLClosing(self, PepXMLHandle):
        PepXMLHandle.write('</msms_run_summary>\n')
        PepXMLHandle.write('</msms_pipeline_analysis>\n')

    def GetSpectrumInfoFromMGF(self, FilePath, FileOffset):
        """ returns the spectrum title and peptide mass corresponding to
            the spectrum at the given file offset in the given mgf file
        """
        File = open(FilePath, "r")
        File.seek(FileOffset)
        Title = None
        MatchTitle = re.compile('^TITLE=([^\n]*)')
        MatchMass = re.compile('^PEPMASS=([^\n]*)')
        # read one line at a time
        for Line in File:
            Match = MatchTitle.match(Line)
            if Match != None:   # start with TITLE=
                Title = Match.group(1)
                continue         
           
            Match = MatchMass.match(Line)
            if Match != None:   # start with PEPMASS
                Mass = Match.group(1)
                File.close()
                return (Title,Mass)

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
                self.InputFileDir = Value
            if Option == "-w":
                self.OutputFileDir = Value
            if Option == "-m":
                self.SpectraDir = Value
        if not OptionsSeen.has_key("-r") or not OptionsSeen.has_key("-w"):
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
    Fix = InspectToPepXMLClass()
    Fix.ParseCommandLine(sys.argv[1:])  
    Fix.Main()
