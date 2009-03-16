"""
MZXML and mzData peak parsing
"""

import os
import sys
import struct
import xml.sax
import xml.sax.handler
import base64
import MSSpectrum

if hasattr(base64, "b64decode"):
    B64Decode = base64.b64decode
    B64Encode = base64.b64encode
else:
    B64Decode = base64.decodestring
    B64Encode = base64.encodestring

def GetSpectrumPeaksMZXML(Spectrum, File):
    Spectrum.Peaks = []
    SAXParser = xml.sax.make_parser()
    Handler = MZXMLPeakParser(Spectrum)
    Handler.Parser = SAXParser
    SAXParser.setContentHandler(Handler)
    try:
        SAXParser.parse(File)
    except xml.sax.SAXException, XMLException:
        Message = XMLException.getMessage()
        # If there are no peaks, then all exceptions are raised:
        if not len(Spectrum.Peaks):
            raise
        # If we did succeed in getting peaks, then the error likely arose
        # after the end of the peaks tag.
        if Message == "junk after document element":
            pass
        elif Message == "not well-formed (invalid token)":
            pass
        else:
            raise

def GetSpectrumPeaksMZData(Spectrum, File):
    Spectrum.Peaks = []
    SAXParser = xml.sax.make_parser()
    Handler = MZDataPeakParser(Spectrum)
    Handler.Parser = SAXParser
    SAXParser.setContentHandler(Handler)
    try:
        SAXParser.parse(File)
    except xml.sax.SAXException, XMLException:
        Message = XMLException.getMessage()
        if Message == "junk after document element":
            pass
        elif Message == "not well-formed (invalid token)":
            pass
        else:
            raise

class MZXMLParseStates:
    SpectrumComplete = -1
    Skipping = 0
    Peaks = 1
    PrecursorMZ = 2

class XMLDictionaryHandler(xml.sax.handler.ContentHandler):
    """
    A simple wrapper for the skeletal ContentHandler class.  Fixes broken API names, and
    supports the use of "triage dictionaries" self.StartHandlers and self.EndHandlers
    to find the handlers for tags.
    """
    def __init__(self):
        # Repair names:
        self.startElement = self.StartElement
        self.endElement = self.EndElement
        self.characters = self.HandleCharacters
        #
        self.VerboseFlag = 0
    def StartElement(self, Name, Attributes):
        #print "Start <%s>@%s"%(Name, self.Parser._parser.CurrentByteIndex)
        Handler = self.StartHandlers.get(Name, None)
        if self.VerboseFlag:
            print "<%s> %s"%(Name, Handler)
        if Handler:
            Handler(Attributes)
    def EndElement(self, Name):
        #print "  End <%s>"%Name
        Handler = self.EndHandlers.get(Name, None)
        if self.VerboseFlag:
            print "</%s> %s"%(Name, Handler)
        if Handler:
            Handler()
    def HandleCharacters(self, String):
        pass
    
class MZXMLPeakParser(XMLDictionaryHandler):
    def __init__(self, Spectrum):
        self.State = MZXMLParseStates.Skipping
        self.StartHandlers = {"peaks":self.StartPeaks,
                              "precursorMz":self.StartPrecursorMZ
                              }
        self.EndHandlers = {"peaks":self.EndPeaks,
                            "scan":self.EndScan,
                            "precursorMz":self.EndPrecursorMZ
                            }
        self.Spectrum = Spectrum
        self.PeakBuffer = ""
        XMLDictionaryHandler.__init__(self)
    def HandleCharacters(self, String):
        if self.State == MZXMLParseStates.PrecursorMZ:
            self.PrecursorMZBuffer += String
            return
        if self.State == MZXMLParseStates.Peaks:
            self.PeakBuffer += String
            return
    def StartPrecursorMZ(self, Attributes):
        if self.State == MZXMLParseStates.SpectrumComplete:
            return
        self.State = MZXMLParseStates.PrecursorMZ
        self.PrecursorMZBuffer = ""
    def EndPrecursorMZ(self):
        if self.State == MZXMLParseStates.SpectrumComplete:
            return
        #print "Precursor MZ -> %s"%self.PrecursorMZBuffer
        self.Spectrum.PrecursorMZ = float(self.PrecursorMZBuffer)
        self.State = MZXMLParseStates.Skipping
    def StartPeaks(self, Attributes):
        if self.State == MZXMLParseStates.SpectrumComplete:
            return
        self.State = MZXMLParseStates.Peaks
        self.PeakBuffer = ""
        ByteOrder = Attributes.get("byteOrder", "network")
        if ByteOrder == "little" or ByteOrder == "little-endian":
            self.ByteOrder = "little"
        else:
            self.ByteOrder = "big"
    def EndScan(self):
        self.State = MZXMLParseStates.SpectrumComplete
    def EndPeaks(self):
        if self.State == MZXMLParseStates.SpectrumComplete:
            return
        DecodedPeaks = B64Decode(self.PeakBuffer)
        StringPos = 0
        self.Peaks = []
        while StringPos < len(DecodedPeaks):
            if self.ByteOrder == sys.byteorder:
                Mass = struct.unpack("f", DecodedPeaks[StringPos:StringPos+4])[0]
                Intensity = struct.unpack("f", DecodedPeaks[StringPos+4:StringPos+8])[0]                
            else:
                Mass = struct.unpack("!f", DecodedPeaks[StringPos:StringPos+4])[0]
                Intensity = struct.unpack("!f", DecodedPeaks[StringPos+4:StringPos+8])[0]
            Peak = MSSpectrum.PeakClass(Mass, Intensity)
            StringPos += 8
            #print Peak.Mass, Peak.Intensity
            self.Spectrum.Peaks.append(Peak)

class MZDataParseStates:
    SpectrumComplete = -1
    Skipping = 0
    MZArray = 1
    MZArrayData = 2
    IntensityArray = 3
    IntensityArrayData = 4

class MZDataPeakParser(XMLDictionaryHandler):
    def __init__(self, Spectrum):
        self.State = MZDataParseStates.Skipping
        self.StartHandlers = {"data":self.StartData,
                              "mzArrayBinary":self.StartMZArrayBinary,
                              "intenArrayBinary":self.StartIntensityArrayBinary,
                              "cvParam":self.StartCVParam,
                              }
        self.EndHandlers = {"data":self.EndData,
                            "spectrum":self.EndSpectrum,
                            }
        self.Spectrum = Spectrum
        self.PeakBuffer = ""
        XMLDictionaryHandler.__init__(self)
    def StartCVParam(self, Attributes):
        Name = Attributes.get("name", None)
        Value = Attributes.get("value", None)
        if Name == "mz":
            self.Spectrum.PrecursorMZ = float(Value)
    def StartMZArrayBinary(self, Attributes):
        if self.State == MZDataParseStates.SpectrumComplete:
            return
        self.State = MZDataParseStates.MZArray
    def StartIntensityArrayBinary(self, Attributes):
        if self.State == MZDataParseStates.SpectrumComplete:
            return
        self.State = MZDataParseStates.IntensityArray
    def HandleCharacters(self, String):
        if self.State in (MZDataParseStates.MZArrayData, MZDataParseStates.IntensityArrayData):
            self.PeakBuffer += String
    def EndData(self):
        if self.State in (MZDataParseStates.MZArrayData, MZDataParseStates.IntensityArrayData):
            # Parse the float array:
            FloatList = []
            DecodedPeaks = B64Decode(self.PeakBuffer)
            StringPos = 0
            while StringPos < len(DecodedPeaks):
                if self.ByteOrder == sys.byteorder:
                    if self.Precision == 64:
                        Value = struct.unpack("d", DecodedPeaks[StringPos:StringPos + 8])[0]
                    else:
                        Value = struct.unpack("f", DecodedPeaks[StringPos:StringPos + 4])[0]
                else:
                    if self.Precision == 64:
                        Value  = struct.unpack("!d", DecodedPeaks[StringPos:StringPos + 8])[0]
                    else:
                        Value = struct.unpack("!f", DecodedPeaks[StringPos:StringPos + 4])[0]
                #Peak = MSSpectrum.PeakClass(Mass, Intensity)
                FloatList.append(Value)
                if self.Precision == 64:
                    StringPos += 8
                else:
                    StringPos += 4
                #print Peak.Mass, Peak.Intensity
                #self.Spectrum.Peaks.append(Peak)
            if self.State == MZDataParseStates.MZArrayData:
                self.MZList = FloatList
            else:
                self.IntensityList = FloatList
            #print "...parsed %s values!"%len(FloatList)
            self.State = MZDataParseStates.Skipping
    def StartData(self, Attributes):
        if self.State == MZDataParseStates.SpectrumComplete:
            return
        self.Precision = Attributes.get("precision", "32")
        ByteOrder = Attributes.get("endian", "network")
        if ByteOrder == "little" or ByteOrder == "little-endian":
            self.ByteOrder = "little"
        else:
            self.ByteOrder = "big"
        if self.State == MZDataParseStates.MZArray:
            self.State = MZDataParseStates.MZArrayData
            self.PeakBuffer = ""
            return
        if self.State == MZDataParseStates.IntensityArray:
            self.State = MZDataParseStates.IntensityArrayData
            self.PeakBuffer = ""
    def EndSpectrum(self):
        if self.State != MZXMLParseStates.SpectrumComplete:
            self.State = MZXMLParseStates.SpectrumComplete
            for PeakIndex in range(len(self.MZList)):
                Mass = self.MZList[PeakIndex]
                Intensity = self.IntensityList[PeakIndex]
                Peak = MSSpectrum.PeakClass(Mass, Intensity)
                self.Spectrum.Peaks.append(Peak)
