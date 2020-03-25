# microTyper
genotyping microhaplotypes

Mainly focused on genotyping using reads from amplicon sequencing.  
  
Reads bam files using the bamtools API, which is bundled with the  
microTyper code for ease of install.
  
  
Installation and example

First, we unpack and compile  
```
  tar -xzvf microTyper.tar.gz
  cd microTyper/
  cmake .
  make
```
  
Now, navigate to the example folder  
```
  cd example/  
```
  
Here, there are three samples, a position file, and a reference file.
We then calculate log-likelihoods for each genotype and sample  
```
  ../mtype -f *.bam -p examplePositionFile.txt -r exampleReference.fasta  
```
  
Then, we call genotypes using a minimum posterior probability of .99 and a 
minimum number of reads that perfectly match one of the alleles in the genotype of 10  
```
  ../genoCaller -f microtyper_llh.mhgenos -c .99 -m 10
```
  

