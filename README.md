# PCPI
A Poisson-based cumulative paternity index analysis system for non-invasive prenatal paternity testing. PCPI provides a novel CPI model by combining traditional and Poisson modelling with a previous judicial CPI model, where users can upload sequencing data of alleged father genomic DNA, mother gDNA and mother cfDNA. In addition, NIPPT-specific CPE model was also developed with improved accuracy. Unique molecular identifiers (UMI) is also supported to reduce sequencing errors in early pregnancy samples.

## Install
	tar -xvf PCPI.tar.gz
	cd ./PCPI
	make install

Dependancies installed by `Makefile`:

+ [fgbio](https://github.com/fulcrumgenomics/fgbio)
+ [bwa](https://github.com/lh3/bwa)
+ [samtools](https://github.com/samtools/samtools)
+ [bcftools](https://github.com/samtools/bcftools)  

And various perl modules.

Please see the respective licence for each before use.

## Usage
For paired-end sequencing:  
`perl pcpi.pl --snp SNPLIST --f1 FASTQ1 --f2 FASTQ2 --m1 FASTQ1 --m2 FASTQ2 --c1 FASTQ1 --c2 FASTQ2 --w 6 [options]`  
For single-end sequencing:  
`perl pcpi.pl --snp SNPLIST --f1 FASTQ1 --m1 FASTQ1 --c1 FASTQ1 --w 6 [options]`  

The options are:
| Option | Effect |
| ---- | ---- |
| f1 | .fastq file 1 of father |
| f2 | .fastq file 2 of father |
| m1 | .fastq file 1 of mather |
| m2 | .fastq file 1 of mather |
| c1 | .fastq file 1 of child |
| c2 | .fastq file 1 of child |
| snp | A list of snp for paternity test |
| w | gestational weeks |
| out | Output Directory |
| n | Threads number. Default is 4. |
| r | Same as --read-structure in fgbio ExtractUmisFromBam. '3M2S+T,3M2S+T' is set as default. |
| t | Same as --molecular-index-tags in fgbio ExtractUmisFromBam. 'ZA,ZB' is set as default. |
| umi | enable the unique molecular identifier sequences analysis |  
| sur | enable the paternity analysis for surrogacy duo cases (conflicts with --umi) |

## Configuration file format
A SNP list is required for the generation and alignment of upstream and downstream sequences. If you need to customize the SNP panel, you need to refer to PCPI/db/nipptRESHAPE.tsv and PCPI/ref/nipptRESHAPE.hg19.fa to modify the database.
Below is an example of a tsv file:  
| SNPid | Chr.hg19 | Pos.hg19 | Chr.GRCh38 | Pos.GRCh38 | Alleles with Frequency |
| ---- | ---- | ---- | ---- | ---- | ---- |
| rs75062661 | chr1 | 69511 | - | - | A | 0.348 | G | 0.652 |
| rs77418980 | chr1 | 91536 | - | - | G | 0.6814 | T | 0.3186 |
| rs151118460 | chr1 | 91581 | - | - | G | 0.6717 | A | 0.3283 |
| rs4951859 | chr1 | 729679 | - | - | C | 0.3407 | G | 0.6593 |
| rs143214544 | chr1 | 748878 | - | - | G | 0.3843 | T | 0.6157 |  

Where location information is not necessary.

In addition, the reference sequence should ensure that the SNP locus is located at position 501.  
