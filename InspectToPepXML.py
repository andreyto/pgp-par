#!/usr/bin/python

UsageInfo = \
"""InspectToPepXML.py: Converts output of InsPecT search engine
to PepXML format.  Written by Samuel Payne, Venter Institute,
and Terry Farrah, Institute for Systems Biology  October 2008

Required Parameters
-i [Filename] - InsPecT results file from search against mzXML (input)
-o [Filename] - Converted file in PepXML (output)

Optional Parameters
-p [Filename] - InsPecT input (parameter) file
                  default: inspect.params
-m [Dirname] - Dir containing .mzXML spectrum file
                  default: current working directory
-d N   - write at most N hits per assumed charge

Assumes InsPecT results file is TSV containing header line and
 one record per peptide prediction sorted by scan #, then by rank.
User must manually edit PepXML file and insert correct information
  near top of file for precursor and fragment mass types --
  either average or monoisotopic.
If database file mentioned in parameter file is not in fasta
  format (.fasta or .fa), you must create a fasta format file of
  the same base name in the same dir. Use TrieToFASTA.py.
This script, InspectToPepXML.py, must reside in the same directory
  as the rest of the InsPecT code.
"""

import sys
import os
import glob
import getopt
import re
import time
import GetByteOffset
import ResultsParser
import Utils
import Global
from xml.sax import saxutils          #xml.sax is for reading mzXML
from xml.sax import ContentHandler
from xml.sax import make_parser
from xml.sax.handler import feature_namespaces

global initial_dir
global spectrum_query_count

# ========================================================
# Read tables with standard data such as amino acid masses
# ========================================================

# chdir() is Hack to make pgm invokable from any dir
#  (Utils makes use of auxiliary files in same dir as code)
initial_dir = os.getcwd()
os.chdir(sys.path[0])
Utils.Initialize()  
os.chdir(initial_dir)

# ===========================================================
# Define classes to hold spectra (scans) and peptides (hits)
# ===========================================================

class InspectSpectrumClass:
    """Stores the relevant InsPecT output file data for a spectrum"""
    def __init__(self):
        self.ScanNumber = -1
        self.PrecursorMz = -1.0
        self.RetentionTime = -1.0
        self.HitList = []   # store a hit list for each charge state
        for i in range (1,6): self.HitList = self.HitList + [[]]

    def WriteSpectrumQueries(self, PepXMLHandle, SpectrumFileName, enzyme,
            MaxHitsPerCharge):
        """ Write <spectrum_query> tags for this spectrum.

            There is one tag for each assumed charge that has any hits.
        """
        global spectrum_query_count

        for charge in range(1,5):     # for each charge state
            if len(self.HitList[charge]) > 0: # if any hits
                spectrum_query_count = spectrum_query_count + 1
                SpectrumTitle="%s.%05d.%05d.%s" % \
                   (os.path.splitext(SpectrumFileName)[0],
                     self.ScanNumber,self.ScanNumber,
                     charge)
                _proton_mass = 1.007276
                PrecursorNeutralMass =  \
                    (self.PrecursorMz * charge) - \
                    (charge * _proton_mass)
                Query = '<spectrum_query\n' + \
                    '    spectrum="%s"\n' % SpectrumTitle + \
                    '    start_scan="%s"\n' % self.ScanNumber + \
                    '    end_scan="%s"\n' % self.ScanNumber + \
                    '    precursor_neutral_mass="%.4f"\n' % \
                                  PrecursorNeutralMass + \
                    '    assumed_charge="%s"\n' % charge + \
                    '    index="%s"\n' % spectrum_query_count + \
                    '    retention_time_sec="%.2f"\n' % \
                                  self.RetentionTime + \
                    '>\n'
                PepXMLHandle.write(Query)
                PepXMLHandle.write('<search_result search_id="1">\n')
                for i in range(min(MaxHitsPerCharge,
                         len(self.HitList[charge]))):
                    self.HitList[charge][i].PrecursorNeutralMass = \
                          PrecursorNeutralMass
                    self.HitList[charge][i].WriteSearchHit(PepXMLHandle,
                          i+1, enzyme)
                PepXMLHandle.write('</search_result>\n')
                PepXMLHandle.write('</spectrum_query>\n')

class InspectOutputRecordClass:
    """Stores the relevant data from a single line of InsPecT output.

    Each line represents a search hit--a predicted peptide for a spectrum.
    """
    def __init__(self):
        self.Spectrum = None
        self.FileOffset = -1
        self.Protein = ""
        self.Charge = -1
        self.MQScore = ""
        self.FScore = ""
        self.DeltaScore = ""
        self.PValue = ""
        self.ProteinID = ""
        self.Prefix = ""
        self.Peptide = ""
        self.Suffix = ""
        self.OptModList = []
        self.PrecursorNeutralMass = -1.0

    def WriteSearchHit(self, PepXMLHandle, rank, enzyme):
        """ Write <search_hit> tag for this this line of InsPecT output
        """           
        global initial_dir
        os.chdir(sys.path[0]) # hack to make pgm invokable from any dir
        # GetMass adds on fixed modifications, but not optional ones
        CalcMass = Utils.GetMass(self.Peptide) + 18.01528 #add h2o mass
        for mod in self.OptModList:
            CalcMass = CalcMass + float(mod[2])
        os.chdir(initial_dir)
        MassDiff = self.PrecursorNeutralMass - CalcMass
        # If the enzyme is trypsin, count all KR except
        # final one, and except when followed by P (proline).
        if enzyme.lower() == "trypsin":
            MissedCleavages = self.Peptide[:-1].count("K") + \
               self.Peptide[:-1].count("R") - \
               (self.Peptide[:-1].count("KP") + self.Peptide[:-1].count("RP"))
        elif enzyme.lower() == "none":
            MissedCleavages = 0
        else: MissedCleavages = -1
        # Break up Protein into accession # and description
        first_space = self.Protein.find(' ')
        if first_space >= 0:
          Protein = self.Protein[:first_space]
          ProteinDescr = self.Protein[first_space+1:]
          ProteinDescr = ProteinDescr.replace(">","&gt;")
          ProteinDescr = ProteinDescr.replace("<","&lt;")
          ProteinDescr = ProteinDescr.replace("&","&amp;")
        else:
          Protein = self.Protein
          ProteinDescr = Protein
        Hit = '<search_hit\n' + \
              '    hit_rank="%s"\n' % (rank) + \
              '    peptide="%s"\n' % (self.Peptide) + \
              '    peptide_prev_aa="%s"\n' % (self.Prefix) + \
              '    peptide_next_aa="%s"\n' % (self.Suffix) + \
              '    protein="%s"\n' % (Protein) + \
              '    protein_descr="%s"\n' % (ProteinDescr) + \
              '    num_tot_proteins="0"\n' + \
              '    num_matched_ions="0"\n'  + \
              '    tot_num_ions="0"\n' + \
              '    calc_neutral_pep_mass="%s"\n' % (CalcMass) + \
              '    massdiff="%s"\n' % (MassDiff) + \
              '    num_tol_term="%s"\n' % "2" + \
              '    num_missed_cleavages="%d"\n'%(MissedCleavages) + \
              '    is_rejected="0"\n' + \
              '>\n'
        PepXMLHandle.write(Hit)
        # Create a dictionary of masses of all amino acids that
        # are modified, indexed by peptide position.
        # First, add to the dictionary all aa's that have optional mods.
        # Then, add fixed mods. Use monoisotopic mass for basic AA.
        ModMassDict = {}
        for mod in self.OptModList:
           aa = mod[0]
           pos = mod[1]
           mod_mass = mod[2]
           if pos in ModMassDict:
             ModMassDict[pos] += float(mod_mass)
           else:
             ModMassDict[pos] = float(mod_mass) + Global.AminoMass[aa]
        for i in range(len(self.Peptide)):
           aa = self.Peptide[i]
           pos = i + 1
           if aa in Global.FixedMods:
               mod_mass = Global.FixedMods[aa]
               if pos in ModMassDict:
                   ModMassDict[pos] += float(mod_mass)
               else:
                   ModMassDict[pos] = float(mod_mass) + \
                                            Global.AminoMass[aa]
        # Now, create a pepXML string with an element for each modified AA.
        ModString = ''
        for i in range(len(self.Peptide)):
            pos = i + 1
            if pos in ModMassDict:
                ModString = ModString + '<mod_aminoacid_mass ' + \
                  'position="%d" ' % pos + \
                  'mass="%.4f" />' % ModMassDict[pos]
        if len(ModString) > 0:
            ModInfo = '<modification_info>%s</modification_info>\n' % \
                   ModString
            PepXMLHandle.write(ModInfo)
        PepXMLHandle.write(
          '<search_score name="mqscore" value="%s"/>\n'%self.MQScore)
        PepXMLHandle.write(
          '<search_score name="expect" value="%s"/>\n'%self.PValue)
        PepXMLHandle.write(
          '<search_score name="fscore" value="%s"/>\n'%self.FScore)
        PepXMLHandle.write(
          '<search_score name="deltascore" value="%s"/>\n'%self.DeltaScore)
        PepXMLHandle.write('</search_hit>\n')


# ======================================================================
# Virtually all the code is contained within class InspectToPepXMLClass
# ======================================================================

class InspectToPepXMLClass(ResultsParser.ResultsParser):
    def __init__(self):
        """Initialize fields of InspectToPepXMLClass instance to null values
        """
        self.InputFilePath = None
        self.OutputFilePath = None
        self.SpectraDir = os.getcwd()
        self.MaxHitsPerCharge = 10000   #effectively maxint
        self.ParamFilePath = os.path.join(os.getcwd(), "inspect.params")
        self.ScanOffset = {}
        self.ScanDict= {}
        self.SpectrumFileType = ""
        self.SpectrumFileBase = ""
        ResultsParser.ResultsParser.__init__(self)
        
    #---------------------------------------------------------------------

    def Main(self):
        """Convert raw InsPecT output file to PepXML

        Initially designed to handle entire directories of files.
        """
        try:
            import psyco
            psyco.full()
        except:
            print "(psyco not found - running in non-optimized mode)"
        # Line directly below needed only if we want to handle directories
        #self.self.ProcessResultsFiles(self.InputFilePath,
        #    self.ConvertInspectToPepXML)
        self.ConvertInspectToPepXML(self.InputFilePath)

    #---------------------------------------------------------------------

    def ConvertInspectToPepXML(self, FilePath):
        """ Convert a single raw InsPecT output file to PepXML
        """

        global spectrum_query_count

        # ------------------------------------------------------------
        # Open input/output files and gather info from auxiliary files
        # ------------------------------------------------------------

        # Get input filename; open output file handle
        FileName = os.path.split(FilePath)[1]
        #FileName = FileName.replace(".txt", ".xml")
        #NewPath = os.path.join(self.OutputFilePath, FileName)
        #PepXMLHandle = open(NewPath, "wb")
        PepXMLHandle = open(self.OutputFilePath, "wb")

        # Glean info from inspect input file
        if not os.path.exists(self.ParamFilePath):
            print >> sys.stderr, "Inspect input file %s does not exist" % \
                  self.ParamFilePath
            sys.exit()
        ParamFile = open(self.ParamFilePath, "r")
        nmods_allowed_per_spectrum = 0
        nmods_in_params = 0
        self.mod_weight = []
        self.mod_aa = []
        self.mod_type = []
        self.mod_name = []
        # reset Global.FixedMods to empty; Global.py initializes it to
        # {"C":57.0518}, but this is a hack we don't want
        Global.FixedMods = {}
        self.instrument = "UNKNOWN"
        self.protease = "trypsin"
        self.search_db = "UNKNOWN"
        for Line in ParamFile.readlines():
           Line = Line.strip()  #remove leading and trailing whitespace
           if Line.lower().startswith("mods,"):
               nmods_allowed_per_spectrum = int(Line[len("mods,"):])
           elif Line.lower().startswith("mod,"):
               tokens = Line.split(",")
               this_mod_weight =  float(tokens[1])
               this_mod_aa_string = tokens[2].strip()
               if this_mod_aa_string == "*":
                 this_mod_aa_string = "ACDEFGHIKLMNPQRSTVWY"
               this_mod_type = "opt"
               if len(tokens) > 3:
                   this_mod_type = tokens[3].strip()
               if len(tokens) > 4:
                   this_mod_name = tokens[4].strip()
               else: this_mod_name = None
               for this_mod_aa in this_mod_aa_string:
                   self.mod_weight = self.mod_weight + [this_mod_weight]
                   self.mod_aa = self.mod_aa + [this_mod_aa]
                   self.mod_type = self.mod_type + [this_mod_type]
                   if this_mod_type == "fix":
                       Global.FixedMods[this_mod_aa] = this_mod_weight
                   self.mod_name = self.mod_name + [this_mod_name]
                   nmods_in_params = nmods_in_params + 1
           elif Line.lower().startswith("instrument,"):
               self.instrument = Line[len("instrument,"):].strip()
           elif Line.lower().startswith("protease,"):
               self.protease = Line[len("protease,"):].strip()
           elif Line.lower().startswith("db,"):
               search_db = Line[len("db,"):].strip()
               search_db_ext = os.path.splitext(search_db)[1]
               # Find the Fasta format of the .trie file used
               if search_db_ext not in [".fa", ".fasta"]:
                   search_db_root = os.path.splitext(search_db)[0]
                   search_db_file_list = set(
                           glob.glob("%s.*" % search_db_root))
                   #print search_db_file_list
                   ext_list = [os.path.splitext(f)[1]
                           for f in search_db_file_list]
                   #print ext_list
                   try: ext_list.remove(".index")
                   except: pass
                   try: ext_list.remove(".trie")
                   except: pass
                   #print ext_list
                   if len(ext_list) == 1:
                       search_db = search_db_root + ext_list[0]
                   elif ".fasta" in ext_list:
                       search_db = search_db_root + ".fasta"
                   elif ".fa" in ext_list:
                       search_db = search_db_root + ".fa"
                   else:
                       print >> sys.stderr, \
           "Can't find a RefreshParser compatible database " + \
           "file corresponding to %s " % search_db + \
           "(such as a .fasta or .fa file with same root)"
                       sys.exit(1)
               self.search_db = search_db
        self.nmods = nmods_in_params

        # Read just first line of inspect output to get spectrum filename
        InspectHandle = open(FilePath, "r")
        for Line in InspectHandle.xreadlines():
            if Line[0] == "#":  # comments
                continue
            Bits = list(Line.split("\t"))
            break
        InspectHandle.close()
        InspectHandle = open(FilePath, "rb")

        # Glean info from spectrum file
        SpectrumFilePath = Bits[self.Columns.SpectrumFile]
        SpectrumFileName = os.path.split(SpectrumFilePath)[1]
        SpectrumFilePath = os.path.join(self.SpectraDir, SpectrumFileName)
        self.SpectrumFileBase = \
          SpectrumFilePath.replace(os.path.splitext(SpectrumFilePath)[1], "")
        self.SpectrumFileType = os.path.splitext(SpectrumFilePath)[1]

        if not os.path.exists(SpectrumFilePath):
            print >> sys.stderr, "Spectrum file %s does not exist" % \
                  SpectrumFilePath
            sys.exit()
        (Stub, Ext) = os.path.splitext(SpectrumFilePath)
        if Ext.lower() == ".mzxml":
            self.SpectrumFileType = ".mzXML"
            (retentionTimeDict, precursorMzDict) =  \
                self.GetSpectrumInfoFromMzXML(SpectrumFilePath)
        elif Ext.lower() == ".mgf":
            self.SpectrumFileType = ".mgf"
        else:
            print >> sys.stderr, \
               "Spectrum file %s lacks .mzXML or .mgf extension" % \
                  SpectrumFilePath
            sys.exit()

        # ------------------------------------------------------------
        # - Write opening info to PepXML
        # - Process InsPecT output file line by line and write to PepXML
        # - Write closing info to PepXML
        # ------------------------------------------------------------

        self.WritePepXMLOpening(PepXMLHandle, self.OutputFilePath)

        LastScanNumber = -1
        spectrum_query_count = 0

        # Each line represents a predicted peptide for a spectrum (scan).
        # A scan can have multiple predicted peptides (hits).
        #  All hits for a scan are grouped together in the file.
        for Line in InspectHandle.xreadlines():
            if Line[0] == "#": continue  # skip comments
            # create a record for this line and read the fields into Bits
            this_rec = InspectOutputRecordClass()
            Bits = list(Line.split("\t"))

            this_rec.FileOffset = int(Bits[self.Columns.FileOffset])
            ScanNumber = int(Bits[self.Columns.ScanNumber])

            if (LastScanNumber != ScanNumber): 
                if (LastScanNumber != -1):
                    # write results for last spectrum
                    this_scan.WriteSpectrumQueries(PepXMLHandle,
                        SpectrumFileName, self.protease,
                        self.MaxHitsPerCharge)
                # initialize new spectrum
                this_scan = InspectSpectrumClass()
                this_scan.ScanNumber = ScanNumber
                # get info about spectrum from spectrum file
                # (not fully implemented for .mgf files)
                if self.SpectrumFileType == ".mgf":
                    print "Sorry, can't handle Inspect searches on MGF files."
                    print "Please perform your search on an mzXML file."
                    sys.exit(1)
                    #(MgfSpectrumTitle, MgfPepMass) = \
                    #   self.GetSpectrumInfoFromMGF(SpectrumFilePath,
                    #   this_rec.FileOffset)
                    #PrecursorNeutralMass = float(MgfPepMass)
                elif self.SpectrumFileType == ".mzXML":
                  if not retentionTimeDict.has_key(this_scan.ScanNumber):
                    print "WARNING: RT for scan %s not found in spectrum file; using RT of -1"
                    this_scan.RetentionTime = "-1.0"
                  else:
                    this_scan.RetentionTime = \
                          retentionTimeDict[this_scan.ScanNumber]
                  if not precursorMzDict.has_key(this_scan.ScanNumber):
                    print "WARNING: m/z for scan %s not found in spectrum file; using m/z of -1"
                    this_scan.PrecursorMz = "-1.0"
                  else:
                    this_scan.PrecursorMz = \
                           precursorMzDict[this_scan.ScanNumber]

            this_rec.Spectrum = this_scan
            LastScanNumber = ScanNumber

            # ---------------------------
            # Process data about this hit
            # ---------------------------
            Annotation = Bits[self.Columns.Annotation]
            Peptide = Annotation[2:-2]

            # process peptide string --TMF
            # I think there is already code to do this in Utils.py
            # Sam, you may want to replace my code with a call to that.
            def ExtractAAModifications(search, peptide):
              '''Given peptide like TVAM+16GGKYphosLV, extract the numbers
                 and other modification symbols.

                 Return (a) the peptide without the mods, and
                 (b) a list of (aa, aa-pos, number) tuples -- 
                 aa/aa-pos describe the aa posessing the mod.
              '''
              i = 0
              mod_list = []
              stripped_peptide = ""
              while i < len(peptide):
                if peptide[i].isupper():
                  stripped_peptide = stripped_peptide + peptide[i]
                  i = i + 1
                  continue
                j = i + 1
                while j < len(peptide) and not peptide[j].isupper():
                  j = j + 1
                aa = peptide[i-1]
                added_mod = peptide[i:j]
                added_mod_pos = len(stripped_peptide) #counting starts at 1
                # modifications with names in the param file
                # will be represented by their names embedded in the
                # peptide. Look up their weights.
                for k in range(nmods_in_params):
                  if search.mod_name[k]:
                     truncated_name = search.mod_name[k][:4]
                     this_weight = int(search.mod_weight[k])
                     if this_weight > 0:
                       weight_string = "+" + str(this_weight)
                     added_mod = added_mod.replace(truncated_name,
                         weight_string)
                # added_mod could be a concatenation of several mods,
                #  as in AEQDNLGKSVM-5+16IPTK;
                # store each one as a separate mod.
                this_mod = ""
                for i in range(len(added_mod)):
                  c = added_mod[i]
                  if (c == "+" or c == "-"):
                    # store the previous mod
                    if len(this_mod) > 0:
                      mod_list = mod_list + [(aa, added_mod_pos, this_mod)]
                    # start a new mod
                    this_mod = c
                  else:
                    this_mod = this_mod + c
                # store the last mod
                mod_list = mod_list + [(aa, added_mod_pos, this_mod)]
                i = j
              return (stripped_peptide, mod_list)

            (this_rec.Peptide, this_rec.OptModList) = \
                 ExtractAAModifications(self, Peptide)

            # done processing peptide string
   
            this_rec.Prefix = Annotation[0]
            this_rec.Suffix = Annotation[-1]
            this_rec.Protein = Bits[self.Columns.ProteinName]
            this_rec.Charge = int(Bits[self.Columns.Charge])
            this_rec.MQScore = Bits[self.Columns.MQScore]
            this_rec.FScore = Bits[self.Columns.FScore]
            this_rec.DeltaScore = Bits[self.Columns.DeltaScore]
            this_rec.PValue = Bits[self.Columns.PValue]
            this_rec.ProteinID = Bits[self.Columns.ProteinID]

            this_scan.HitList[this_rec.Charge] = \
                this_scan.HitList[this_rec.Charge] + [this_rec]
             
            # done processing a single line of InsPecT output file

        self.WritePepXMLClosing(PepXMLHandle)
        InspectHandle.close()
        PepXMLHandle.close()

    #---------------------------------------------------------------------

    def WritePepXMLOpening(self, PepXMLHandle, PepXMLFilePath):
        """Write stuff that belongs at the top of the pepXML file"""
        PepXMLHandle.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        PepXMLHandle.write('<?xml-stylesheet type="text/xsl" href="pepXML_std.xsl"?>\n')
        datestr = time.strftime('%Y-%m-%dT%H:%M:%S')
        PepXMLHandle.write(
           '<msms_pipeline_analysis ' + 
             'date="%s" ' % datestr  +
             'summary_xml="%s" ' %PepXMLFilePath +
             'xmlns="http://regis-web.systemsbiology.net/pepXML" ' +
             'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" ' +
             'xsi:schemaLocation="http://regis-web.systemsbiology.net/pepXML /tools/bin/TPP/tpp/schema/pepXML_v112.xsd" ' +
           '>\n' )
           #'xsi:schemaLocation="http://regis-web.systemsbiology.net/pepXML http://mascot1/mascot/xmlns/schema/pepXML_v18/pepXML_v18.xsd" ' +
        PepXMLHandle.write(
           '<msms_run_summary ' +
             'base_name="%s" ' % self.SpectrumFileBase +
             'search_engine="InsPecT" ' +
             'msManufacturer="UNKNOWN" ' +
             'msModel="%s" ' % self.instrument +
             'msIonization="UNKNOWN" ' +
             'msMassAnalyzer="UNKNOWN" ' +
             'msDetector="UNKNOWN" ' +
             'raw_data_type="raw" ' +
             'raw_data="%s" ' % self.SpectrumFileType +
           '>\n')
        PepXMLHandle.write('<sample_enzyme name="%s">\n' % self.protease)
        PepXMLHandle.write('<specificity cut="KR" no_cut="P" sense="C"/>\n')
        PepXMLHandle.write('</sample_enzyme>\n')
        PepXMLHandle.write(
           '<search_summary ' +
             'base_name="%s" ' % self.SpectrumFileBase +
             'search_engine="InsPecT" ' +
             'precursor_mass_type="monoisotopic" ' +
             'fragment_mass_type="monoisotopic" ' +
             'search_id="1" ' +
             'out_data_type="out" ' +
             'out_data=".txt" ' +
           '>\n')
        PepXMLHandle.write(
           '<search_database ' +
             'local_path="%s" ' % self.search_db +
             'type="AA" ' +
           '/>\n')
        #'database_name="" ' +
        #'database_release_identifier="" ' +
        #'size_in_db_entries="" ' +
        #'size_of_residues="" ' +
        PepXMLHandle.write(
           '<enzymatic_search_constraint ' +
              'enzyme="%s" ' % self.protease +
              'max_num_internal_cleavages="2" ' +
              'min_number_termini="2" ' +
           '/>\n')
        for i in range(self.nmods):
            mod_aa = self.mod_aa[i]
            mod_weight = self.mod_weight[i]
            mass = mod_weight + Global.AminoMass[mod_aa]
            if self.mod_type[i] == "opt": mod_variable="Y"
            elif self.mod_type[i]=="fix": mod_variable="N"
            else: mod_variable="UNKNOWN"   # are there other types?
            PepXMLHandle.write(
               '<aminoacid_modification ' +
                'aminoacid="%s" ' % mod_aa +
                'massdiff="%.4f" ' % mod_weight +
                'mass="%.4f" ' % mass +
                'variable="%s" ' % mod_variable +
               '/>\n')
        PepXMLHandle.write('<parameter name="CHARGE" value="2+ and 3+"/>\n')
        PepXMLHandle.write('<parameter name="CLE" value="Trypsin"/>\n')
        PepXMLHandle.write('<parameter name="DB" value=""/>\n')
        PepXMLHandle.write('<parameter name="FILE" value=""/>\n')
        PepXMLHandle.write('<parameter name="FORMAT" value=""/>\n')
        PepXMLHandle.write('<parameter name="FORMVER" value=""/>\n')
        PepXMLHandle.write('<parameter name="INSTRUMENT" value="%s"/>\n' % \
          self.instrument)
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

    #---------------------------------------------------------------------

    def WritePepXMLClosing(self, PepXMLHandle):
        """Write stuff that belongs at the end of the pepXML file"""
        PepXMLHandle.write('</msms_run_summary>\n')
        PepXMLHandle.write('</msms_pipeline_analysis>\n')

    #---------------------------------------------------------------------

    def GetSpectrumInfoFromMGF(self, FilePath, FileOffset):
        """ returns the spectrum title and peptide mass corresponding to
            the spectrum at the given file offset in the given mgf file
        """
        File = open(FilePath, "r")
        File.seek(FileOffset)
        Title = None
        MatchTitle = re.compile('^TITLE=([^\n]*)')
        MatchMass = re.compile('^PEPMASS=([^\n])\s*')
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

    #--------------------------------------------------------------------

    def GetSpectrumInfoFromMzXML(self, FilePath):
        """ compiles dictionaries of the precursorMz and retentionTime
            for each spectrum in an mzXML file
        """

	def normalize_whitespace(text):
	    "Remove redundant whitespace from a string"
	    return ' '.join(text.split())

	class MzXMLHandler(ContentHandler):

	    def __init__(self):
                self.this_scan = None
                self.this_precursorMz = None
                self.precursorMz = dict()
                self.retentionTime = dict()
                self.inPrecursorMzContent = 0

	    def startElement(self, name, attrs):
		# If it's not a comic element, ignore it
		if name == 'scan':
		    # Look for the title and number attributes
		    num = int(normalize_whitespace(attrs.get('num', None)))
		    retentionTime = normalize_whitespace(
				     attrs.get('retentionTime', None))
		    self.this_scan = int(num)
		    self.retentionTime[self.this_scan] = \
                        float(retentionTime[2:-1])
		elif name == 'precursorMz':
		    self.inPrecursorMzContent = 1
		    self.thisprecursorMz = ""

	    def characters(self, ch):
		if self.inPrecursorMzContent:
		    self.thisprecursorMz = self.thisprecursorMz + ch

	    def endElement(self, name):
		if name == 'precursorMz':
		    self.inPrecursorMzContent = 0
                    i = self.this_scan
                    self.precursorMz[i] = float(self.thisprecursorMz)
		elif name == 'scan':
                    pass

	# Create an XML parser and tell it
	# we are not interested in XML namespaces
	MzXMLparser = make_parser()
	MzXMLparser.setFeature(feature_namespaces, 0)

	# Create a handler and tell the parser to use it
	mh = MzXMLHandler()
	MzXMLparser.setContentHandler(mh)

	# Parse the file
        File = open(FilePath, "r")
	MzXMLparser.parse(File)

        return (mh.retentionTime, mh.precursorMz)

    #---------------------------------------------------------------------

    def ParseCommandLine(self,Arguments):
        (Options, Args) = getopt.getopt(Arguments, "i:o:m:p:d:")
        OptionsSeen = {}
        for (Option, Value) in Options:
            OptionsSeen[Option] = 1
            if Option == "-i":
                if not os.path.exists(Value):
                  print "** Error: couldn't find results file '%s'\n\n"%Value
                  print UsageInfo
                  sys.exit(1)
                self.InputFilePath = Value
            if Option == "-o":
                self.OutputFilePath = Value
            if Option == "-m":
                self.SpectraDir = Value
            if Option == "-p":
                self.ParamFilePath = Value
            if Option == "-d":
                self.MaxHitsPerCharge = int(Value)
        if not OptionsSeen.has_key("-i") or not OptionsSeen.has_key("-o"):
            print UsageInfo
            sys.exit(1)

    def Finish(self):   
        self.InputFile.close()
        self.OutputFile.close()

#-------------------------------------------------------------------------

if __name__ == '__main__':
    Fix = InspectToPepXMLClass()
    Fix.ParseCommandLine(sys.argv[1:])  
    Fix.Main()
