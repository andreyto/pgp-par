"""PGPrimaryStructure.py
THis is a set of classes designed to help analyze the primary
structure of proteins, and whether they agree with the proteomic
evidences

NOTE: this is a utility, and not executable from the command line

"""
import os, DNA, GFFIO

class PrimaryStructure:
    """Class PrimaryStructure: This is an analysis object to help
    us determine whether the predicted protein is in harmony with
    the peptides observed in proteomics experiments
    Variables:
        -
    Functions: CheckStructure(PGPeptide.OpenReadingFrame)
    """
    # Class variable to count start codons
    startCodonCount = 1

    def __init__(self, OutputPath, NucleotideSequence):
        """Parameters: none
        Return: none
        Description: trivial constructor
        """
        #open a file handle to print out stuff
        (Path, Ext) = os.path.splitext(OutputPath)
        self.OutputStub = Path
        NovelFastaPath = "%s.%s"%(self.OutputStub, "novel.faa")
        #Use append mode so output from multiple chromosomes is not overwritten
        #in testing make sure to remove or rename your previous output
        self.NovelFastaHandle = open(NovelFastaPath, "a")
        NovelInfoPath = "%s.%s"%(self.OutputStub, "novel.info")
        self.NovelInfoHandle = open(NovelInfoPath, "a")
        UnderPredictionInfoPath = "%s.%s"%(self.OutputStub, "underprediction.info")
        self.UnderpredictionInfoHandle = open(UnderPredictionInfoPath, "a")
        ShortProteinPath = "%s.%s"%(self.OutputStub, "shortprotein.info")
        self.ShortProteinHandle = open(ShortProteinPath, "a")
        startCodonPath = "%s.%s"%(self.OutputStub, "starts.gff")
        self.startCodonGFF = GFFIO.File(startCodonPath, "a")
        self.HypotheticalCount =0
        self.NamedCount = 0
        self.NovelGC = []
        self.NotNovelGC = []
        self.NovelLen = []
        self.NotNovelLen = []
        self.DNA = NucleotideSequence
        # Start codons to look for
        self.StartCodonTable = [ 'ATG', 'GTG', 'TTG' ]

    def CheckStructure(self, ORF):
        """
        Parameters: a PGPeptide.OpenReadingFrame object
        Return: ?
        Description: various QA/QC checks
        """
        #it is possible that someone called this with an ORF that
        #has no peptide support, so we should not do much for them
        if ORF.numPeptides() == 0:
            return "NO EVIDENCE"


        Novel = self.IsItNovel(ORF)
        if Novel:
            self.FindAllUpstreamStartCodons(ORF)
            return "NOVEL"

        self.ShortProteins(ORF)
        ProteinName = ORF.annotatedProtein.name
        Location = ProteinName.find("hypothetical")
        if Location == -1:
            self.NamedCount += 1
        else:
            self.HypotheticalCount += 1
        UnderPredicted = self.IsItUnderPredicted(ORF)
        if UnderPredicted:
            self.FindAllUpstreamStartCodons(ORF)
            return "UNDERPREDICTED"

    def ShortProteins(self, ORF):
        # I want to do a quick sweep of short proteins.
        ProteinLength = ORF.GetPredictedProteinLength()
        if ProteinLength > 100:
            return
        self.ShortProteinHandle.write(">%s\n"%ORF)
        self.ShortProteinHandle.write("%s\n"%ORF.GetProteinSequence())




    def CalculateGC(self, Sequence):
        """Parameters: DNA sequence
        Return: integer percent GC, eg 45 for 45%
        Description: count g and c
        """
        GCCount =0
        TotalCount = len(Sequence)
        for Letter in Sequence:
            if Letter in ["G", "C", "g", "c"]:
                GCCount += 1
        Fraction = GCCount / float (TotalCount)
        Percent = int (Fraction *100)
        return Percent

    def OutputGCFiles(self):
        """
        Parameters: None
        Return: None
        Description: I have spent some time in the Novel checker to
        save the GC content of the observed sequences.  This Function
        will take these lists and put them out to a file
        """

        NovelPath = "%s.%s"%(self.OutputStub, "GC.novel.txt")
        Handle = open(NovelPath, "wb")
        String = "\t".join(map(str, self.NovelGC))
        Handle.write(String)
        Handle.close()

        NotNovelPath = "%s.%s"%(self.OutputStub, "GC.notnovel.txt")
        Handle = open(NotNovelPath, "wb")
        String = "\t".join(map(str, self.NotNovelGC))
        Handle.write(String)
        Handle.close()


    def OutputLenFiles(self):
        NovelPath = "%s.%s"%(self.OutputStub, "len.novel.txt")
        Handle = open(NovelPath, "wb")
        String = "\t".join(map(str, self.NovelLen))
        Handle.write(String)
        Handle.close()

        NotNovelPath = "%s.%s"%(self.OutputStub, "len.notnovel.txt")
        Handle = open(NotNovelPath, "wb")
        String = "\t".join(map(str, self.NotNovelLen))
        Handle.write(String)
        Handle.close()


    def IsItNovel(self, ORF):
        """
        Parameters: a PGPeptide.OpenReadingFrame object
        Return: true/false
        Description: check to see if this ORF has peptides but not a
        predicted protein

        NOTE: in the old version of this, I had a minimum ORF length.  Should
        I put that in here? I'm not inclined yet, ,mostly because I don't
        remember why I put that in there, and don't have any really solid
        research compelling me to do so.
        """
        #One quick thing here before we start the novelty.
        # I want to get the GC content , and store that separately
        #for the Novel ones and the normal ones
        (Start, Stop) = ORF.GetObservedDNACoords()
        Sequence = self.DNA[Start:Stop]
        if Sequence:
            GC = self.CalculateGC(Sequence)
        else:
            print "WARNING: cannot extract observed sequence for %s (%s, %s)"%(ORF, Start, Stop)
            GC = 0
            for pep in ORF.peptideIter():
                print "A Pep %s" % pep

        #now on with the novelty, I swear
        PredictedProtein = ORF.GetLocatedProtein()
        if PredictedProtein:
            self.NotNovelGC.append(GC)
            self.NotNovelLen.append(len(Sequence)/3) #div by three because it's protein length that I care about
            return False #there is something here, so no it's not novel
        #put in the GC quickly
        self.NovelGC.append(GC)
        self.NovelLen.append(len(Sequence)/3)
        #well, if we get here, then we have something interesting
        #let's try and get the observed sequence.  That's firstpeptide->stop
        ObservedSequence = ORF.GetObservedSequence()
        NovelInfoString = "I got this observed sequence from %s\n"%ORF
        Count =0
        for Peptide in ORF.peptideIter():
            Count += 1
            NovelInfoString += "Peptide %s, %s\n"%(Count, Peptide)
        NovelInfoString += "%s\n\n"%ObservedSequence
        self.NovelInfoHandle.write("%s"%NovelInfoString)
        #print "I got this observed sequence from %s\n%s\n\n"%(ORF, ObservedSequence)
        Fasta = "Observed.%s"%ORF.name
        self.NovelFastaHandle.write(">%s\n%s\n"%(Fasta, ObservedSequence))
        return True

    def IsItUnderPredicted(self, ORF):
        """
        Parameters: a PGPeptide.OpenReadingFrame object
        Return: true/false
        Description: Check to see if this ORF has peptides upstream
        of the currently predicted start (at least two amino acids
        upstream). We should also include whether there is a start
        site upstream to use
        """
        #1. we get the nuc of the 5' peptide
        FirstObservedPeptide = ORF.GetFivePrimePeptide(Unique=1) #here we care about uniqueness
        if not FirstObservedPeptide:
            print "Hey dummy, there are no unique peptides for %s"%ORF
            return False
        FirstObservedNucleotide = FirstObservedPeptide.GetFivePrimeNucleotide()
        StartCodon = ORF.GetNucleotideStartOfTranslation()
        Strand = ORF.GetStrand()
        # we'd also like to get a count on the number of peptides upstream.
        NumPeptidesUpstream = ORF.GetUpstreamPeptideCount()

        #2. Strand switch that we use all over the place
        if Strand == "+":
            if FirstObservedNucleotide + 3 < StartCodon: # do the plus three
                #because we want more than a single amino acid upstream.
                UpstreamExtent = StartCodon - FirstObservedNucleotide
                self.UnderpredictionInfoHandle.write("Peptide %s is %s bases upstream of protein in %s\n"%(FirstObservedPeptide, UpstreamExtent, ORF))
                self.UnderpredictionInfoHandle.write("There are %s peptides upstream\n\n"%NumPeptidesUpstream)
                #print "Peptide %s is %s bases upstream of protein in %s\n\n"%(FirstObservedPeptide, UpstreamExtent, ORF)
                return True
        else:
            if FirstObservedNucleotide > StartCodon + 3:
                UpstreamExtent = FirstObservedNucleotide - StartCodon
                self.UnderpredictionInfoHandle.write("Peptide %s is %s bases upstream of protein in %s\n"%(FirstObservedPeptide, UpstreamExtent, ORF))
                self.UnderpredictionInfoHandle.write("There are %s peptides upstream\n\n"%NumPeptidesUpstream)
                #print "Peptide %s is %s bases upstream of protein in %s\n\n"%(FirstObservedPeptide, UpstreamExtent, ORF)
                return True

        #3. Is the actual start codon observed as a L or V, instead of M
        if FirstObservedNucleotide == StartCodon:
            FirstObservedAminoAcid = FirstObservedPeptide.aminos[0]
            if not FirstObservedAminoAcid == "M":
                self.UnderpredictionInfoHandle.write("Peptide %s is at the start codon of protein in %s\n\n"%(FirstObservedPeptide, ORF))
                #print "Peptide %s is at the start codon of protein %s\n\n"%(FirstObservedPeptide, ORF)
                return True

    def FindAllUpstreamStartCodons(self, ORF):
        #1. we get the nuc of the 5' peptide
        FirstObservedPeptide = ORF.GetFivePrimePeptide()
        if not FirstObservedPeptide:
            print "No peptides for %s"%ORF
            return
        # Subtract one to use zero based coordinates for DNA array access
        firstNucleotide = FirstObservedPeptide.GetFivePrimeNucleotide()
        orfStart = ORF.location.GetFivePrime()
        # Keep track of the peptide position in the ORF
        peptidePosition = 1
        peptideIncr = 1
        # GFF record to output the start codons
        gffRec = GFFIO.Record()
        gffRec.source = 'Proteomics'
        gffRec.type = 'polypeptide'
        gffRec.seqid = ORF.chromosome
        gffRec.attributes['Parent'] = ORF.name
        if not ORF.CDS:
            # Novel ORF output the Observed Coords
            if ORF.location.strand == '+':
                (gffRec.start, gffRec.end) = ORF.GetObservedDNACoords()
            else:
                (gffRec.end, gffRec.start) = ORF.GetObservedDNACoords()
            gffRec.score = 0
            gffRec.strand= ORF.location.strand
            gffRec.attributes['Name'] = "Observed2Stop"
            gffRec.attributes['ID'] = PrimaryStructure.startCodonCount
            PrimaryStructure.startCodonCount += 1
            self.startCodonGFF.write( gffRec )

        # For the reverse strand switch postions so we're still going up
        if ORF.location.strand == '-':
            (orfStart, firstNucleotide) = (firstNucleotide, orfStart)
            # We need to go back one AA to check the 1st AA
            orfStart -= 2 # one AA in genomic base based coordinates
            nucDiff = firstNucleotide - orfStart + 1
            if nucDiff % 3 != 0:
                print "Warning frame is off diff is %d" % nucDiff
            peptidePosition = nucDiff / 3
            peptideIncr = -1

        print "Looking for upstream starts from %d-%d for %s 1st codons %s" % (
            orfStart, firstNucleotide, ORF, self.DNA[orfStart-1:orfStart+8])
        first = True
        while (orfStart < firstNucleotide):
            # Need zero base based coordinates for the array access
            codon = self.DNA[orfStart-1:orfStart+2].upper()
            codonLoc = orfStart
            if ORF.location.strand == '-':
                codon = DNA.ReverseComplement( codon )

            if (first and codon == 'ATG') or (
            not first and codon in self.StartCodonTable):
                gffRec.start = codonLoc
                gffRec.end   = codonLoc+2
                gffRec.score = 0
                gffRec.strand= ORF.location.strand
                gffRec.attributes['Name'] = "%s@%d" % (codon,peptidePosition)
                gffRec.attributes['ID'] = PrimaryStructure.startCodonCount
                if ORF.CDS:
                    gffRec.attributes['locus_tag'] = ORF.CDS.qualifiers['locus_tag'][0]

                PrimaryStructure.startCodonCount += 1
                self.startCodonGFF.write( gffRec )
                print "Start Codon %s found at %d pep %d" % (codon, codonLoc, peptidePosition)
            orfStart += 3
            peptidePosition += peptideIncr
            first = False
