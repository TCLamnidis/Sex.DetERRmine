#!/usr/bin/env python3
import sys, argparse
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

parser = argparse.ArgumentParser(description="Calculate the relative X- and Y-chromosome coverage of data, as well as the associated error bars for each.")
parser.add_argument("-I", "--Input", metavar="<INPUT FILE>", type=argparse.FileType('r'), help="The input samtools depth file. Omit to read from stdin.", required=False)
parser.add_argument("-f", "--SampleList", type=argparse.FileType('r'), help="A list of samples/bams that were in the depth file. One per line. Should be in the order of the samtools depth output.")
args = parser.parse_args()

if args.Input == None:
    args.Input = sys.stdin

Names=OrderedDict()
if args.SampleList != None:
    Samples = [line.strip() for line in args.SampleList]
    for idx,Sample in enumerate(Samples):
        Names.update({Sample:idx})
    NrAut  = [0 for x in range(len(Names))]
    NrX    = [0 for x in range(len(Names))]
    NrY    = [0 for x in range(len(Names))]
    # Totals = [0 for x in range(len(Names))]
    
Reads={}
AutSnps=0
YSnps=0
XSnps=0
for line in args.Input:
    fields=line.strip().split()
    Chrom=fields[0]
    if fields[0][0]=="#":
        if args.SampleList==None:
            Zip    = zip(fields[2:],range(len(fields[2:])))
            for Sample,Index in Zip:
                Names.update({Sample:Index})
            NrAut  = [0 for x in range(len(Names))]
            NrX    = [0 for x in range(len(Names))]
            NrY    = [0 for x in range(len(Names))]
            # Totals = [0 for x in range(len(Names))]
            continue
        else:
            continue
    depths=[int(x) for x in fields[2:]]
    if Chrom != "Y" and Chrom != "X":
        AutSnps+=1
    if Chrom == "Y":
        YSnps+=1
    if Chrom == "X":
        XSnps+=1
    for x in Names:
        # Totals[Names[x]]+=depths[Names[x]]
        if Chrom != "Y" and Chrom != "X":
            NrAut[Names[x]]+=depths[Names[x]]
        if Chrom == "Y":
            NrY[Names[x]]+=depths[Names[x]]
        if Chrom == "X":
            NrX[Names[x]]+=depths[Names[x]]

SortNames=OrderedDict(sorted(Names.items(), key=lambda t: t[1]))
print ("#Sample", "#SnpsAut", "#SNPsX", "#SnpsY", "NrAut", "NrX", "NrY", "x-rate", "y-rate", "Err(x-rate)", "Err(y-rate)", sep="\t", file=sys.stdout)
for Ind in Names:
    rate,rateErr=CalcErrors(AutSnps, XSnps, YSnps, NrAut[Names[Ind]], NrX[Names[Ind]], NrY[Names[Ind]])
    print (Ind, AutSnps, XSnps, YSnps, NrAut[Names[Ind]], NrX[Names[Ind]], NrY[Names[Ind]], rate["X"], rate["Y"], rateErr["X"], rateErr["Y"], sep="\t", file=sys.stdout)
