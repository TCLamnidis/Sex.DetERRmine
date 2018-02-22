# Sex.DetERRmine
A python script carry out calculate the relative coverage of X and Y chromosomes, and their associated error bars, out of capture data.


# Instructions:
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
