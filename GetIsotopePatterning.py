"""
Scaffolding script: Compute the expected ratio of primary isotopic peaks to secondary isotopic peaks, for
weights up to 1750.  Uses 'averagine' (i.e. average elemental composition)
"""
def Binomial(n, k):
    Result = 1.0
    for X in range(k+1, n+1):
        Result *= X
    for X in range(1, n-k+1):
        Result /= float(X)
    return Result

def GetOdds(AtomCount, Odds):
    if not AtomCount:
        return (1.0, 0, 0, 0, 0, 0)
    AllIonOdds = []
    #print "%s atoms, heavy chance %s"%(AtomCount, Odds)
    for IonCount in range(0, AtomCount):
        Ways = Binomial(AtomCount, IonCount)
        OddsProduct = (Odds)**IonCount * (1.0 - Odds)**(AtomCount - IonCount)
        IonOdds = Ways * OddsProduct
        #print "%s: %s ways, odds %s, %s overall"%(IonCount, Ways, OddsProduct, IonOdds)
        AllIonOdds.append(IonOdds)
    AllIonOdds.append(0)
    AllIonOdds.append(0)
    AllIonOdds.append(0)
    AllIonOdds.append(0)
    return AllIonOdds

def GetOdds2(AtomCount, OddsPlus1, OddsPlus2):
    if AtomCount == 0:
        return (1,0,0,0,0,0)
    OddsPlain = 1.0 - OddsPlus1 - OddsPlus2
    AllIonOdds = [0,0,0,0,0,0,0]
    AllIonOdds[0] = OddsPlain**AtomCount
    AllIonOdds[1] = Binomial(AtomCount, 1) * OddsPlain**(AtomCount-1) * OddsPlus1
    # Plus 2: One +2 isotope or two +1 isotopes:
    AllIonOdds[2] = Binomial(AtomCount, 1) * OddsPlain**(AtomCount-1) * OddsPlus2
    AllIonOdds[2] += Binomial(AtomCount, 2) * OddsPlain**(AtomCount-2) * OddsPlus1**2
    # Plus 3: One +2 and one +1, or three +1:
    AllIonOdds[3] = AtomCount*(AtomCount-1) * OddsPlain**(AtomCount-2) * OddsPlus2 * OddsPlus1
    AllIonOdds[3] += Binomial(AtomCount, 3) * OddsPlain**(AtomCount-3) * OddsPlus1**3
    return AllIonOdds

def CombineOdds(OldOdds, NewOdds):
    while len(OldOdds)<len(NewOdds):
        OldOdds.append(0.0)
    for Index in range(len(NewOdds)):
        OldOdds[Index] += NewOdds[Index]
    return OldOdds

def AdjustForMass2(OxygenOdds):
    NewOxygenOdds = []
    for Index in range(len(OxygenOdds)):
        NewOxygenOdds.append(OxygenOdds[Index])
        NewOxygenOdds.append(0.0)
    return NewOxygenOdds

##print "10 choose 4:", Binomial(10, 4)    
##print "10 choose 5:", Binomial(10, 5)
##print "10 choose 6:", Binomial(10, 6)
##print "10 choose 10:", Binomial(10, 10)
##print "10 choose 0:", Binomial(10, 0)
##print "10 choose 1:", Binomial(10, 1)
##
for Mass in range(0, 1750):
    CarbonCount = int(round(Mass * 0.54 / 12.0))
    # CarbonOdds[n] is the odds of adding n to the base peak
    CarbonOdds = GetOdds(CarbonCount, .011070)
    #OverallIonOdds = CarbonOdds[:]
    #while len(OverallIonOdds)<10:
    #    OverallIonOdds.append(0.0)
    #print "CARBON:", CarbonOdds
    #
    HydrogenCount = int(round(Mass * 0.0666)/1.007805)
    HydrogenOdds = GetOdds(HydrogenCount, .000150)
    #CombineOdds(OverallIonOdds, HydrogenOdds)
    #print "HYDROGEN:", HydrogenOdds
    #
    NitrogenCount = int(round(Mass * 0.170903)/14.00307)
    NitrogenOdds = GetOdds(NitrogenCount, .003663)
    #CombineOdds(OverallIonOdds, NitrogenOdds)
    #
    OxygenCount = int(round(Mass * 0.195212265)/15.99491)
    OxygenOdds = GetOdds2(OxygenCount, .000375, .002035)
    OxygenOdds = AdjustForMass2(OxygenOdds)
    #CombineOdds(OverallIonOdds, OxygenOdds)
    #
    SulfurCount = int(round(Mass * 0.026910932)/31.9721)
    SulfurOdds = GetOdds2(SulfurCount, .00760, .04290)
    SulfurOdds = AdjustForMass2(SulfurOdds)
    #CombineOdds(OverallIonOdds, SulfurOdds)
    Odds = 1.0
    if CarbonCount:
        Odds *= CarbonOdds[0]
    if HydrogenCount:
        Odds *= HydrogenOdds[0]
    if NitrogenCount:
        Odds *= NitrogenOdds[0]
    if OxygenCount:
        Odds *= OxygenOdds[0]
    if SulfurCount:
        Odds *= SulfurOdds[0]
    Odds0 = Odds
    Odds1 = CarbonOdds[1]*HydrogenOdds[0]*NitrogenOdds[0]*OxygenOdds[0]*SulfurOdds[0]
    Odds1 += CarbonOdds[0]*HydrogenOdds[1]*NitrogenOdds[0]*OxygenOdds[0]*SulfurOdds[0]
    Odds1 += CarbonOdds[0]*HydrogenOdds[0]*NitrogenOdds[1]*OxygenOdds[0]*SulfurOdds[0]
    Odds1 += CarbonOdds[0]*HydrogenOdds[0]*NitrogenOdds[0]*OxygenOdds[1]*SulfurOdds[0]
    Odds1 += CarbonOdds[0]*HydrogenOdds[0]*NitrogenOdds[0]*OxygenOdds[0]*SulfurOdds[1]
    ####################
    Odds2 = CarbonOdds[2]*HydrogenOdds[0]*NitrogenOdds[0]*OxygenOdds[0]*SulfurOdds[0]
    Odds2 += CarbonOdds[0]*HydrogenOdds[2]*NitrogenOdds[0]*OxygenOdds[0]*SulfurOdds[0]
    Odds2 += CarbonOdds[0]*HydrogenOdds[0]*NitrogenOdds[2]*OxygenOdds[0]*SulfurOdds[0]
    Odds2 += CarbonOdds[0]*HydrogenOdds[0]*NitrogenOdds[0]*OxygenOdds[2]*SulfurOdds[0]
    Odds2 += CarbonOdds[0]*HydrogenOdds[0]*NitrogenOdds[0]*OxygenOdds[0]*SulfurOdds[2]
    ####
    Odds2 += CarbonOdds[1]*HydrogenOdds[1]*NitrogenOdds[0]*OxygenOdds[0]*SulfurOdds[0]
    Odds2 += CarbonOdds[1]*HydrogenOdds[0]*NitrogenOdds[1]*OxygenOdds[0]*SulfurOdds[0]
    Odds2 += CarbonOdds[1]*HydrogenOdds[0]*NitrogenOdds[0]*OxygenOdds[1]*SulfurOdds[0]
    Odds2 += CarbonOdds[1]*HydrogenOdds[0]*NitrogenOdds[0]*OxygenOdds[0]*SulfurOdds[1]
    Odds2 += CarbonOdds[0]*HydrogenOdds[1]*NitrogenOdds[1]*OxygenOdds[0]*SulfurOdds[0]
    Odds2 += CarbonOdds[0]*HydrogenOdds[1]*NitrogenOdds[0]*OxygenOdds[1]*SulfurOdds[0]
    Odds2 += CarbonOdds[0]*HydrogenOdds[1]*NitrogenOdds[0]*OxygenOdds[0]*SulfurOdds[1]
    Odds2 += CarbonOdds[0]*HydrogenOdds[0]*NitrogenOdds[1]*OxygenOdds[1]*SulfurOdds[0]
    Odds2 += CarbonOdds[0]*HydrogenOdds[0]*NitrogenOdds[1]*OxygenOdds[0]*SulfurOdds[1]
    Odds2 += CarbonOdds[0]*HydrogenOdds[0]*NitrogenOdds[0]*OxygenOdds[1]*SulfurOdds[1]
    
    #print "%s\t%s\t%s\t%s\t"%(Mass, Odds0, Odds1, Odds2)
    print "%s\t%s\t%s\t"%(Mass, Odds1/float(Odds0), Odds2/float(Odds0))