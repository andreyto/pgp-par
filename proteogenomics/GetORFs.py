#Eventually change to have relative pathnames
#inFile = open("/Users/nataliecastellana/ArabadopsisData/AT.trie", 'r') 
#does not create a file that doesn't exist
#outFile = open("/Users/nataliecastellana/ArabadopsisData/ORFS.gff", 'w')

import sys

GeneticCode = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
               "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
               "TAT":"Y", "TAC":"Y", "TAA":"X", "TAG":"X",
               "TGT":"C", "TGC":"C", "TGA":"X", "TGG":"W",
               "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
               "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
               "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
               "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
               "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
               "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
               "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
               "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
               "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
               "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
               "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
               "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
               }

RCDict = {"A":"T", "G":"C", "T":"A", "C":"G",
          "a":"T", "g":"C", "t":"A", "c":"G"}

def Translate(DNA):
    "Returns the peptide translation of a sequence."
    Codon = DNA.upper()
    AA = GeneticCode.get(Codon, "*")
    return AA

def TranslateRev(DNA):

	#print DNA
	revDNA = ""
	for i in range(0,3):
		revDNA += RCDict.get(DNA[i], DNA[i])
        #print revDNA
        
	AA = GeneticCode.get(revDNA.upper(),"*")
        #print AA
	#raw_input()
        return AA

if len(sys.argv) < 4:
	print "Usage: python ORFFinder [trie file] [output fasta file] [min ORF length] (sequence name)"
	sys.exit()

inFile = open(sys.argv[1], 'r')
outFile = open(sys.argv[2], 'w')
cutoff = int(sys.argv[3])

print sys.argv[1]
print sys.argv[2]
print sys.argv[3]

if len(sys.argv) == 5:
	seqName = sys.argv[4]
else:
	seqName = "At"

stopCodonsRev = ["ATC", "ACT", "ATT"]
stopCodons = ["TAG", "TGA", "TAA"]
lineNum = 1

orfCountFor = 0
orfCount = 0
totalLengthFor = 0
totalLength = 0
#Mouse DGenes for IgH locus
#RangeStart = 2463846
#RangeEnd = 2559216

#Mouse JGenes for IgH locus
RangeStart = 2580122
RangeEnd = 2581440

for line in inFile.xreadlines():
    line = line[RangeStart:RangeEnd]
    print "****************"
    print "Line Num: %s"%lineNum
    line = line.strip()
    print "Total characters: %s"%str(len(line))
    lineNum = lineNum + 1
    for phase in range(3):
        start = phase
        end = phase
        sequence = ""
        print "phase " + str(phase)
        while(1):		
            if end < len(line) - 3: #If we have a full codon
                codon = line[end:end+3].upper()
                for a in stopCodons: #Check if it's a stop codon
                    stopFlag = 0
                    if a == codon :
                        stopFlag = 1
                        lenOrf = end-start
                        if lenOrf >= cutoff :
                            totalLength = totalLength + lenOrf
                            totalLengthFor += lenOrf
                            if end%3 == start%3: #Sanity Check
                                
                                orfCount += 1
                                orfCountFor +=1
                                outFile.write(">" + seqName + "(" + str(start +RangeStart) + "-" + str(end + RangeStart) + ")+\n")
                                outFile.write(sequence + "\n")
                                sequence = ""
                                
                            else :
                                print "Fuck"
				
                        start = end + 3
                        end = start
                        sequence = ""
                                
                        break
				
                else :
                    end = end + 3
                    if stopFlag == 0:
                        sequence += Translate(codon)
                                                
				
            else:
                print "End of line\n"
                break
    line = line[::-1].strip()

    orfCountRev = 0
    totalLengthRev = 0

    print "******REVERSE**********"
    
    lineNum = lineNum + 1
	
    for phase in range(3):
	
        start = phase
        end = phase
        sequence = ""
        NSeq = ""
        stopFlag = 0
        print "phase " + str(phase)
        while(1):	
            
            if end+3 <= len(line):
                codon = line[end:end+3].upper()
                for a in stopCodonsRev:
                    stopFlag = 0           
                    if a == codon :
                        stopFlag = 1
                        lenOrf = end-start
                        if lenOrf >= cutoff :
                            totalLength = totalLength + lenOrf
                            totalLengthRev += lenOrf

                            if end%3 == start%3:
                                s = seqName +"_" + str(orfCount) + "\tORFFinder\texon\t" + str(len(line) + RangeStart-end) + "\t" + str(len(line) + RangeStart-start) + "\t1\t-\t0\t" + "phase:" + str(phase) + "\n"				
                                orfCount += 1
                                orfCountRev += 1
                                outFile.write(">" + seqName + "(" + str(len(line)-end) + "-" + str(len(line)-start) + ")-\n")
                                outFile.write(sequence + "\n")

                                sequence = ""
                                NSeq = ""
                            else :
                                print "Fuck"
							
                        start = end + 3            
                        end = start
                        sequence = ""
                        NSeq = ""
                        break			
                else :
                    end = end + 3
                    if stopFlag == 0:
                        NSeq += codon
                        sequence += TranslateRev(codon)
            else:
                print "End of line\n"
                break


print "ORFcount Forward: " + str(orfCountFor)
print "Average Length Forward: " + str(totalLengthFor/orfCountFor)
print "ORFcount Reverse: " + str(orfCountRev)
print "Average Length Reverse: " + str(totalLengthRev/orfCountRev)
print "ORFcount Total: " + str(orfCount)
print "Average Length Total: " + str(totalLength/orfCount)



