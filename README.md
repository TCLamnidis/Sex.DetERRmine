# Sex.DetERRmine
A python script carry out calculate the relative coverage of X and Y chromosomes, and their associated error bars, out of capture data.


# Instructions
The python script takes a modified output from `samtools depth` as input, via stdin. The samtools depth file should be manually modified to include a header that begins with a `#` and is including the sample names (generic or specific) as column headers, like below:
```
#Chr	Pos	Sample1	Sample2	Sample3	Sample4	Sample5
1	752566	1	0	1	0	1
1	776546	0	0	0	0	0
1	832918	0	1	0	0	0
1	842013	0	1	0	3	1
 ...
```
You can then run the script as below:
```
SexDet.ErrorCalc.py <Input.Depth.file.txt 
```

The script will print out the number of SNPs and the number of reads found on each of Autosomes/X/Y, as well as the relative X/Y coverage and their associated errors.

# Mathematical explanation
We assume that sequenced reads are distributed along the genome randomly and independently from each other. The "genome" here is made up only of positions in the input depth file. 

<p align="center"><img src="https://latex.codecogs.com/gif.latex?\LARGE&space;\textsl{N}=\sum_{i}&space;\textsl{N}_{i}" title="\textsl{N}=\sum_{i} \textsl{N}_{i}" /></p>

_N<sub>i</sub>_ is the number of sequenced reads in a a chunk of the genome _i_, the sum of which is the total number of reads on target, _N_.

We can then calculate:

<p align="center"><img src="https://latex.codecogs.com/gif.latex?\LARGE&space;\textsl{p}_{i}=\frac{\textsl{N}_{i}}{\textsl{N}}" title="\textsl{p}_{i}=\frac{\textsl{N}_{i}}{\textsl{N}}" /></p>

Where _p<sub>i</sub>_ is the proportion of all sequenced reads that map to SNPs in _i_, estimated from the input depths. And: 

<p align="center"><img src="https://latex.codecogs.com/gif.latex?\LARGE&space;\textsl{d}_{i}=\frac{\textsl{N}_{i}}{\textsl{S}_{i}}" title="\textsl{d}_{i}=\frac{\textsl{N}_{i}}{\textsl{S}_{i}}" /></p>

Where _d<sub>i</sub>_ is the average depth on SNPs within _i_, and _S<sub>i</sub>_ is the number of SNPs in _i_.

The relative coverage on the X and Y chromosomes can then be calculated as:

<p align="center"><img src="https://latex.codecogs.com/gif.latex?\LARGE&space;\textsl{rate}=\frac{\textsl{d}_{i}}{\textsl{d}_{Aut}}" title="\textsl{rate}=\frac{\textsl{d}_{X/Y}}{\textsl{d}_{Aut}}" /></p>

