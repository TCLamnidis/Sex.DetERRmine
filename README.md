# Sex.DetERRmine
A python script carry out calculate the relative coverage of X and Y chromosomes, and their associated error bars, from the depth of coverage at specified SNPs.

<br></br>
_Mathematical equations added to README using [this tool](https://www.codecogs.com/latex/eqneditor.php)._

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
Alternatively, a Sample/bam list can be provided using the `-f` option. This list should include 1 name per line, and can be the same list used for the `samtools depth` command.

For instructions on the options available you can try running the script with the `-h` flag:
```
$SexDet.ErrorCalc.py -h

usage: SexDet.ErrorCalc.py [-h] [-I <INPUT FILE>] [-f SAMPLELIST]

Calculate the relative X- and Y-chromosome coverage of data, as well as the
associated error bars for each.

optional arguments:
  -h, --help            show this help message and exit
  -I <INPUT FILE>, --Input <INPUT FILE>
                        The input samtools depth file. Omit to read from
                        stdin.
  -f SAMPLELIST, --SampleList SAMPLELIST
                        A list of samples/bams that were in the depth file.
                        One per line. Should be in the order of the samtools
                        depth output.

```

The script will print out the number of SNPs and the number of reads found on each of Autosomes/X/Y, as well as the relative X/Y coverage and their associated errors.

It is possible to pipe the `samtools depth` output directly to this script:
```
samtools depth -a -q30 -Q30 -b <BED File> -f <BAM file list> | SexDet.ErrCalc.py -f <BAM file list>
```

# Mathematical explanation
We assume that sequenced reads are distributed along the genome randomly and independently from each other. The "genome" here is made up only of positions in the input depth file. 

<p align="center"><img src="https://latex.codecogs.com/gif.latex?\LARGE&space;\textsl{N}=\sum_{i}&space;\textsl{N}_{i}" title="\textsl{N}=\sum_{i} \textsl{N}_{i}" /></p>

_N<sub>i</sub>_ is the number of sequenced reads in a a chunk of the genome _i_, the sum of which is the total number of reads on target, _N_.

We can then calculate:

<p align="center"><img src="https://latex.codecogs.com/gif.latex?\LARGE&space;\textsl{p}_{i}=\frac{\textsl{N}_{i}}{\textsl{N}}" title="\textsl{p}_{i}=\frac{\textsl{N}_{i}}{\textsl{N}}" /></p>

<p align="center"><img src="https://latex.codecogs.com/gif.latex?\LARGE&space;\textsl{Err}\left&space;(&space;N_i&space;\right&space;)=\sqrt{N\times&space;p_i\times&space;(1-p_i)}" title="\LARGE \textsl{Err}\left ( N_i \right )=\sqrt{N\times p_i\times (1-p_i)}" /></p>

Where _p<sub>i</sub>_ is the proportion of all sequenced reads that map to SNPs in _i_, estimated from the input depths. The error around _N<sub>i</sub>_ is the error of the binomial distribution. Then: 

<p align="center"><img src="https://latex.codecogs.com/gif.latex?\LARGE&space;\textsl{d}_{i}=\frac{\textsl{N}_{i}}{\textsl{S}_{i}}" title="\textsl{d}_{i}=\frac{\textsl{N}_{i}}{\textsl{S}_{i}}" /></p>

<p align="center"><img src="https://latex.codecogs.com/gif.latex?\LARGE&space;\textsl{Err(d}_{i})=\frac{\textsl{Err(N}_{i})}{S_{i}}" title="\LARGE \textsl{Err(d}_{i})=\frac{\textsl{Err(N}_{i})}{S_{i}}" /></p>

Where _d<sub>i</sub>_ is the average depth on SNPs within _i_, and _S<sub>i</sub>_ is the number of SNPs in _i_.

The relative coverage on the X and Y chromosomes can then be calculated as:

<p align="center"><img src="https://latex.codecogs.com/gif.latex?\LARGE&space;\textsl{rate}=\frac{\textsl{d}_{i}}{\textsl{d}_{Aut}}" title="\textsl{rate}=\frac{\textsl{d}_{X/Y}}{\textsl{d}_{Aut}}" /></p>

We can then use error propagation to calculate the errors around the relative X and Y coverages:

<p align="center"><img src="https://latex.codecogs.com/gif.latex?\LARGE&space;\textsl{Err(x/y&space;rate)}=\sqrt{\left&space;(&space;Err(d_{x/y})\times&space;\frac{1}{d_{aut}}\right&space;)^{2}&plus;&space;\left&space;(Err(d_{aut})\times&space;\frac{d_{x/y}}{{d_{aut}}^{2}}&space;\right&space;)^{2}}" title="\LARGE \textsl{Err(x/y rate)}=\sqrt{\left ( Err(d_{x/y})\times \frac{1}{d_{aut}}\right )^{2}+ \left (Err(d_{aut})\times \frac{d_{x/y}}{{d_{aut}}^{2}} \right )^{2}}" /></p>
