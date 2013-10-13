"""PLSUtils.py

These are a set of related functions which help calculate the Phosphate
Localization Score or PLS.  The differences between this and the Ascore
are described in Albuquerque MCP 2008.
1. The winner and runner up are determined by the MQScore, not peptide score.
2. We do not repeatedly filter the peaks to find the optimal
AScore.  That just takes too long.  We use the default Inspect peak filtering which
leaves ~12 peaks / 100 m/z.  This corresponds to a p of 0.12 for the
binomial.
3. We will not do any thing for peptides which contain
more than 2 sites of phosphorylation.  I don't trust those annotations
anyway.
"""

import math
import copy
from Utils import*
Initialize()

class PLSClass:
    def __init__(self):
        self.ChooseTable = {} # (N, k) -> value of NchooseK
        self.Factorial = {}
        self.Factorial[0] = 1

    def ComputeBinomial(self, N, n, p=0.12):
        """ Make sure that you have populated the Choose table and Factorial table
            I have defaulted the p = 0.12 because it's inspect peak density.
        """
        #print "computing binomial with N %s and n %s"%(N,n)
        Sum = 0
        for k in range(n, N+1): #stupid range operator, it's gonna kill me
            # n<=k<=N
            Sum += self.ChooseTable[(N,k)] * pow(p,k) * pow((1-p),(N-k))
        return Sum

    def FillOutChooseTable(self, N):
        """Simply dose all the possible N choose k for k < N
        """
        self.FillOutFactorialTable(N)
        for k in range(N+1):
            if self.ChooseTable.has_key((N,k)):
                continue
            NFac = self.Factorial[N]
            kFac = self.Factorial[k]
            NMinusKFac = self.Factorial[N-k]
            Value = NFac / (kFac * NMinusKFac)
            self.ChooseTable[(N, k)] = Value

    def FillOutFactorialTable(self, X):
        """Get all the factorials of numbers <= X """
        for x in range(1, X+1):
            if self.Factorial.has_key(x):
                continue # already have this
            ##here we calculate it for the value x
            Prod = x
            for num in range(1, x):
                Prod *= num
            self.Factorial[x] = Prod
        #for Pair in self.Factorial.items():
        #    print "Factorials ",Pair

    def ComputePLS(self, N, nWinner, nRunnerUp):
        """From the Ascore Paper paper:
        AScore = ScoreWinner - ScoreRunnerUp
        Score = -10 * Log(P)
        P = Sum_{k=n}^N [N choose k] p^k * (1-p)^{N-k}
        k = iterator
        N = number of distinguishing ions
        n = number of distinguishing ions found
        p = some very poorly explained probability number.  It appears to be the probabilty of
            a peak in 100 m/z.  Inspect filters stuff automatically and I don't want to mess
            with it. so this number will be PeakDensity of 100 m/z units /100
        """
        self.FillOutChooseTable(N) # the N, or number of potential peaks
        TopBinomial = self.ComputeBinomial(N, nWinner)
        RunnerUpBinomial = self.ComputeBinomial(N, nRunnerUp)
        TopScore = -10 * math.log(TopBinomial, 10)
        RunnerUpScore = -10 * math.log(RunnerUpBinomial, 10)
        PLS = TopScore - RunnerUpScore # The AScore
        return PLS

    def GetDistinguishingPeaks(self, Peptide1, Peptide2):
        """ Given two peptides, find all the peaks that
        distinguish between the two phos placements, and return
        those in a list.  e.g. SphosEPTIDE vs. SEPTphosIDE
        distinguishing B fragments = b1: Sphos or S
                                    b2: SphosE or SE
                                    b3: SphosEP or SEP
                                    b4: SphosEPT and SEPTphos have the same mass = NOT DISTINGUISHING
    
        The general case is: if the b fragment has different number of phosphorylations
        between the two annotations, then both b and y are distinguishing.
        """
        Differences = [] #list of indicies that are different [1] would be differences on b1, yn-1
        ModIndex1 = Peptide1.Modifications.keys() # a list of indicies
        #print "Peptide %s"%Peptide1.GetModdedName()
        for Index in ModIndex1:
            #print "This is my mod on %s, %s"%(Index, Peptide1.Modifications[Index])
            PTMList = Peptide1.Modifications[Index]
            FoundPhos =0
            for Item in PTMList:
                if Item.Name == "Phosphorylation":
                    FoundPhos = 1
            if not FoundPhos:
                ModIndex1.remove(Index)
        ModIndex2 = Peptide2.Modifications.keys() # a list of indicies
        #print "Peptide %s"%Peptide2.GetModdedName()
        for Index in ModIndex2:
            #print "This is my mod on %s, %s"%(Index, Peptide2.Modifications[Index])
            PTMList = Peptide2.Modifications[Index]
            FoundPhos =0
            for Item in PTMList:
                if Item.Name == "Phosphorylation":
                    FoundPhos = 1
            if not FoundPhos:
                ModIndex2.remove(Index)
        Count1 =0
        Count2 =0
        for Index in range(len(Peptide1.Aminos)):
            if Index in ModIndex1:
                Count1 +=1
            if Index in ModIndex2:
                Count2 +=1
            if not Count1 == Count2:
                Differences.append(Index+1)
        ## now we have a list of the b indicies.  Let's make the return list
        DistinguishingPeaks = []
        for B in Differences:
            YIndex = len(Peptide1.Aminos) - B
            #print "a distinguishing peak B%s, Y%s"%(B, YIndex)
            DistinguishingPeaks.append("B%s"%B)
            DistinguishingPeaks.append("Y%s"%YIndex)
        return DistinguishingPeaks
            

    def GetAlternateAnnotations(self, Peptide):
        """Given an annotation(SAMPAYphosNE), return all
        alternate annotations.  This version should work correctly
        in the presence of non-phosphorylation modifications,
        e.g. SAM+16PAYphosNE
        """
        #Dummy = GetPeptideFromModdedName("SAM+16PAYphosNESphosT")
        #Peptide = Dummy
        NumPhos = Peptide.GetModdedName().count("phos")
        AllAnnotations = []
        if not NumPhos in [1,2]:
            return AllAnnotations # empty list
        ## now see if the number of phos == number of potential residues
        ## this means that an AScore is impossible, only one possible placement
        Count = Peptide.Aminos.count("S")
        Count += Peptide.Aminos.count("T")
        Count += Peptide.Aminos.count("Y")
        if Count == NumPhos:
            return AllAnnotations #key for "N/A"

        (Dephos, PhosPTM) = self.RemovePhosFromPeptide(Peptide)
        for Index in range(len(Dephos.Aminos)):
            if Dephos.Aminos[Index] in ["S", "T", "Y"]:
                #place a phosphorylation
                CreateNewLevel1 =0
                if not Dephos.Modifications.has_key(Index):
                    CreateNewLevel1 = 1
                    Dephos.Modifications[Index] = []
                Dephos.Modifications[Index].append(PhosPTM)
                #do I need to place a second phosphate?
                if NumPhos == 2:
                    for Jndex in range(Index + 1, len(Dephos.Aminos)):
                        if Dephos.Aminos[Jndex] in ["S", "T", "Y"]:
                            CreateNewLevel2 = 0
                            if not Dephos.Modifications.has_key(Jndex):
                                Dephos.Modifications[Jndex] = []
                                CreateNewLevel2 = 1
                            Dephos.Modifications[Jndex].append(PhosPTM)
                            #add string to list, remove phos and move on.
                            #print Dephos.GetModdedName()
                            AllAnnotations.append(Dephos.GetModdedName())
                            Dephos.Modifications[Jndex].pop()
                            if CreateNewLevel2:
                                del Dephos.Modifications[Jndex]
                else:
                    #only one phos.  Add string to list
                    AllAnnotations.append(Dephos.GetModdedName())
                    #print Dephos.GetModdedName()
                #only add string to the list if it's a single phos peptide
                #regardless, get rid of the modification now that we've done all it's combinations
                Dephos.Modifications[Index].pop() # get rid of the most recently added Phos PTM
                if CreateNewLevel1:
                    del Dephos.Modifications[Index]
        ## now one last thing remains.  We have to remove the original annotation from the list
        try:
            AllAnnotations.remove(Peptide.GetModdedName())
        except:
            print "The original peptide was not in the list.  that's BAAAAAAAADDDDDDDDDDDDDD"
            return []
        return AllAnnotations
                    

    def RemovePhosFromPeptide(self, Peptide):
        """Given a Peptide object, remove the phosphorylation PTMod
        objects and return the neutered copy
        """
        Clone = copy.deepcopy(Peptide)
        PTModPhos = None
        RemoveElement = [] #indicies in Clone.Modifications that become empty
        for AminoIndex in Clone.Modifications:
            ModificationList = Clone.Modifications[AminoIndex]
            for PTMod in ModificationList:
                if PTMod.Name == "Phosphorylation":
                    if not PTModPhos:
                        PTModPhos = PTMod
                    ModificationList.remove(PTMod)
            if len(ModificationList) == 0:
                RemoveElement.append(AminoIndex)
        ## now clone is without phosphorylations, but with other modifications
        ## clean up the empty Modification keys
        for Index in RemoveElement:
            del Clone.Modifications[Index]
        return (Clone, PTModPhos)

    def GetSupportingPeaks(self, FullPeakList, DistinguishingPeakList):
        SupportingPeaks = {} # key = b8 value = 1
        for Tuple in FullPeakList:
            Ion = Tuple[2]
            Index = Tuple[3]
            if Ion in ["B", "Y", "B2", "Y2"]:
                Peak = "%s%s"%(Ion[0], Index)
                if Peak in DistinguishingPeakList:
                    #print "Supporting peak found %s"%Peak
                    SupportingPeaks[Peak]= 1
        return len(SupportingPeaks)