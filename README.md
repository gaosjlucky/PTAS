## PCPI
# Description
	Calculate the CPI for non-invasive prenatal paternity testing
# Parameter
	--f1	.fastq file 1 of father
	--f2	.fastq file 2 of father
	--m1	.fastq file 1 of mather
	--m2	.fastq file 1 of mather
	--c1	.fastq file 1 of child
	--c2	.fastq file 1 of child
	--snp	A list of snp for paternity test
	--w	gestational weeks
	--out	Output Directory
	--n	Threads number. Default is 4.
	--r	ReadStructure, the read name will be formatted '<NAME>+<UMIs1><UMIs2>', and this parameter can be defined twice. '3M2S+T,3M2S+T' is set as default.
		Four kinds of operators are recognized:
  		1. 'T' identifies a template read
  		2. 'B' identifies a sample barcode read
  		3. 'M' identifies a unique molecular index read
  		4. 'S' identifies a set of bases that should be skipped or ignored
	--t	molecular index tags, the number of tags must be same as ReadStructure. 'ZA,ZB' is set as default.
	--umi	enable the unique molecular identifier sequences analysis
	--help	Show this information
# Exmple 
	PE:	perl bcpi.pl --snp snp.lst --f1 test_father_1.fq.gz --f2 test_father_2.fq.gz --m1 test_mather_1.fq.gz --m2 test_mather_2.fq.gz --c1 test_child_1.fq.gz --c2 test_child_2.fq.gz --out /test_dir/ --w 6
	SE:	perl bcpi.pl --snp snp.lst --f1 test_father_1.fq.gz --m1 test_mather_1.fq.gz --c1 test_child_1.fq.gz --out /test_dir/ --w 6
