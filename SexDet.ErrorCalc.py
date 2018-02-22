#!/usr/bin/env python3
import sys
from math import sqrt
from collections import OrderedDict

def CalcErrors(AutSnps, XSnps, YSnps, NrAut, NrX, NrY):
    SNPs=[AutSnps, XSnps, YSnps]
    Reads=[NrAut, NrX, NrY]
    p={}
    ErrNr={}
    dp={}
    Errdp={}
    rate={}
    rateErr={}
    
    for Bin,Idx in zip(["Aut", "X", "Y"], range(3)):
        Total=sum(Reads)
        p[Bin]=Reads[Idx]/Total
        ErrNr[Bin]=sqrt(Total*p[Bin])
        dp[Bin]=Reads[Idx]/SNPs[Idx]
        Errdp[Bin]=ErrNr[Bin]/SNPs[Idx]
    
    for Bin in ["X","Y"]:
        rate[Bin]=dp[Bin]/dp["Aut"]
        rateErr[Bin]=sqrt((Errdp[Bin]/dp["Aut"])**2 + (Errdp["Aut"]*dp[Bin]/(dp["Aut"]**2))**2)
    
    return (rate, rateErr)

#### MAIN ####

Names=OrderedDict()
Reads={}
AutSnps=0
YSnps=0
XSnps=0
for line in sys.stdin:
    fields=line.strip().split()
    Chrom=fields[0]
    if fields[0][0]=="#":
        Zip=zip(fields[2:],range(len(fields[2:])))
        for Sample,Index in Zip:
            Names.update({Sample:Index})
        NrAut=[0 for x in range(len(Names))]
        NrX=[0 for x in range(len(Names))]
        NrY=[0 for x in range(len(Names))]
        Totals=[0 for x in range(len(Names))]
        continue
        # print (Names)
        # print (Totals)
        # print (pAut)
    depths=[int(x) for x in fields[2:]]
    for x in Names:
        # Totals[Names[x]]+=depths[Names[x]]
        if Chrom != "Y" and Chrom != "X":
            AutSnps+=1
            NrAut[Names[x]]+=depths[Names[x]]
        if Chrom == "Y":
            YSnps+=1
            NrY[Names[x]]+=depths[Names[x]]
        if Chrom == "X":
            XSnps+=1
            NrX[Names[x]]+=depths[Names[x]]

SortNames=OrderedDict(sorted(Names.items(), key=lambda t: t[1]))
print ("#Sample", "#SnpsAut", "#SNPsX", "#SnpsY", "Nr on target", "NrAut", "NrX", "NrY", "x-rate", "y-rate", "Err(x-rate)", "Err(y-rate)", sep="\t", file=sys.stdout)
for Ind in Names:
    rate,rateErr=CalcErrors(AutSnps, XSnps, YSnps, NrAut[Names[Ind]], NrX[Names[Ind]], NrY[Names[Ind]])
    print (Ind, AutSnps, XSnps, YSnps, NrAut[Names[Ind]], NrX[Names[Ind]], NrY[Names[Ind]], rate["X"], rate["Y"], rateErr["X"], rateErr["Y"], sep="\t", file=sys.stdout)
