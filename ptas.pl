#!/usr/bin/env perl
## Name:           ptas.pl
## Date Created:   Thu Aug 8 2019
## Date Modified:  Thu Oct 1 2021
## By:             LBW

## Update 
#  poisson model
#  surrogacy

use strict;
use Getopt::Long;
use FindBin qw($RealBin);
if ($FindBin::VERSION < 1.51) {
        warn "[!]Your Perl is too old, thus there can only be ONE `genlst.pl` file in your PATH. [FindBin Version: $FindBin::VERSION < 1.51]\n\n"
}
FindBin::again();

my $usage=<<USAGE;

Description
	Calculate the CPI for human family
Parameter
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
	--sur	enable the paternity analysis for surrogacy duo cases (conflicts with --umi)
	--help	Show this information
Exmple 
	PE:	perl $0 --snp snp.lst --f1 test_father_1.fq.gz --f2 test_father_2.fq.gz --m1 test_mather_1.fq.gz --m2 test_mather_2.fq.gz --c1 test_child_1.fq.gz --c2 test_child_2.fq.gz --out /test_dir/ --w 8
	SE:	perl $0 --snp snp.lst --f1 test_father_1.fq.gz --m1 test_mather_1.fq.gz --c1 test_child_1.fq.gz --out /test_dir/ --w 8

USAGE

my ($f_1,$f_2,$m_1,$m_2,$c_1,$c_2,$snp_list,$outdir,$help,$umi,$ReadStructure,$IndexTag,$ThreadsNum,$week,$sur);
GetOptions(
	"f1=s"=>\$f_1,
	"f2=s"=>\$f_2,
	"m1=s"=>\$m_1,
	"m2=s"=>\$m_2,
	"c1=s"=>\$c_1,
	"c2=s"=>\$c_2,
	"snp=s"=>\$snp_list,
	"out=s"=>\$outdir,
	"r=s"=>\$ReadStructure,
	"t=s"=>\$IndexTag,
	"n=i"=>\$ThreadsNum,
	"w=i"=>\$week,
	"umi"=>\$umi,
	"sur"=>\$sur,
	"help" => \$help,
);
die "$usage" if(!$f_1 || !$m_1 || !$c_1 || !$outdir || !$snp_list || $help);
#die "! Only paired-end read has been supported in UMI sequences analysis.\n" if(!$c_2 && $umi);

if(!defined $f_2){$f_2="NA";}
if(!defined $m_2){$m_2="NA";}
if(!defined $c_2){$c_2="NA";}
if(!defined $ReadStructure){$ReadStructure="3M2S+T,3M2S+T";}
if(!defined $IndexTag){$IndexTag="ZA,ZB";}
if(!defined $ThreadsNum){$ThreadsNum=4;}
if(!defined $week){$week=8;}

$ReadStructure =~ s/,/ /g;
$IndexTag =~ s/,/ /g;
my $count1 = $ReadStructure =~ tr/,/,/;
my $count2 = $IndexTag =~ tr/,/,/;
die "The number of index tags must be same as ReadStructure.\n" unless ($count1 == $count2);
###############################################################
#dir
my $info = "$outdir/0info";
my $fqout = "$outdir/1bam";
my $vcfout = "$outdir/2vcf";
my $tsvout = "$outdir/3result";
system("rm -rf $outdir") if (-e "$outdir");
system("mkdir -p $info") unless (-e "$info");
system("mkdir -p $fqout") unless (-e "$fqout");
system("mkdir -p $vcfout") unless (-e "$vcfout");
system("mkdir -p $tsvout") unless (-e "$tsvout");
system("echo Mother     Father  Child > $info/family.lst");

open ER,">$outdir/run.err" or die($!);

#Prepare
my $tempDB="$outdir/0info/DBsnp.tsv";
my $tempPOS="$outdir/0info/DBsnp.bed";
my $tempREF="$outdir/0info/REFsnp.fa";
my $database = "$RealBin/db/nipptRESHAPE.tsv";
my $reference = "$RealBin/ref/nipptRESHAPE.hg19.fa";
my $week_thres = "$RealBin/db/week_threshold.txt";

my %need_snp;
my $snp_count = 0;
open SNP,"<$snp_list" or die($!);
while (<SNP>){
	chomp;
	if ($_ =~ /^rs\d+$/){
		$need_snp{$_}++;
		$snp_count++;
	}else{
		print ER "! Reference establishment failed, please check the format of your list.\n";
		die("Step0 failed!\n");
	}
}
close SNP;

if ($snp_count > 10000){
	print ER "! Reference establishment failed, the number of SNP shall not exceed 10000.\n";
	die("Step0 failed!\n");
}

open DB,"<$database" or die($!);
open DOUT,">$tempDB" or die($!);
open DPOS,">$tempPOS" or die($!);
while (<DB>){
	chomp;
	if ($_ =~ /^#/){
		print DOUT "$_\n";
	}else{
		my @data = split /\t/,$_;
		if (defined $need_snp{$data[0]}){
			print DOUT "$_\n";
			print DPOS "$data[0]\t501\n";
		}
	}
}
close DB;
close DOUT;
close DPOS;

local $/ = "\n>";
open REF,"<$reference" or die($!);
open ROUT,">$tempREF" or die($!);
while (<REF>){
	chomp;
	$_ =~ s/^>//;
	my ($seqin,$seq) = (split /\n/,$_,2)[0,1];
	$seq =~ s/\n//g;
	my ($id) = $seqin =~ /^(\S+)/;
	if (defined $need_snp{$id}){
		print ROUT ">$_\n";
	}
}
close REF;
close ROUT;
local $/ = "\n";

system("$RealBin/bin/bwa index $tempREF");
system("java -jar $RealBin/bin/picard.jar CreateSequenceDictionary R=$tempREF");
system("$RealBin/bin/samtools faidx $tempREF");

##############################################################

my @dBWA;
push @dBWA,&rBWA("Father",$f_1,$f_2);
push @dBWA,&rBWA("Mother",$m_1,$m_2);

if ($umi){
	push @dBWA,&uBWA("Child",$c_1,$c_2);
}else{
	push @dBWA,&rBWA("Child",$c_1,$c_2);
}
push @dBWA,"echo step1 done.";

my $dBWA = join " && ",@dBWA;
my $return1 = readpipe("$dBWA");

my $return2;
if ($return1 =~ /step1 done/){
	if ($umi){
		$return2 = readpipe("$RealBin/bin/samtools mpileup $fqout/Mother.bam $fqout/Father.bam $fqout/Child.fmbam -l $tempPOS -d 50000 -Q 10 -f $tempREF -v -t 'DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR' -p -o $vcfout/Family.vcf.gz 2> $vcfout/mpileup.log && $RealBin/bin/bcftools call -Oz -V indels -m $vcfout/Family.vcf.gz -o $vcfout/Family.snp.gz && $RealBin/bin/bcftools index $vcfout/Family.vcf.gz && $RealBin/bin/bcftools index $vcfout/Family.snp.gz && echo step2 done.");
	}else{
		$return2 = readpipe("$RealBin/bin/samtools mpileup $fqout/Mother.bam $fqout/Father.bam $fqout/Child.bam -l $tempPOS -d 4000 -Q 30 -f $tempREF -v -t 'DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR' -p -o $vcfout/Family.vcf.gz 2> $vcfout/mpileup.log && $RealBin/bin/bcftools call -Oz -V indels -m $vcfout/Family.vcf.gz -o $vcfout/Family.snp.gz && $RealBin/bin/bcftools index $vcfout/Family.vcf.gz && $RealBin/bin/bcftools index $vcfout/Family.snp.gz && echo step2 done.");
	}
}else{
	print ER "! BWA failed, please check the log in $fqout.\n";
	die("Step1 failed!\n");
}

my @tags = qw{Father Mother Child};
if ($return2 =~ /step2 done/){
	for (@tags){
		if ($_ eq "Child"){
			system("$RealBin/bin/bcftools query -f '%CHROM\\t%POS\\t%REF,%ALT\\t%QUAL[\\t%TGT;%AD]\\n' -s $_ -i'POS=501 && MIN(FMT/DP)>=100' $vcfout/Family.snp.gz >$tsvout/$_.tsv");
		}else{
			system("$RealBin/bin/bcftools query -f '%CHROM\\t%POS\\t%REF,%ALT\\t%QUAL[\\t%TGT;%AD]\\n' -s $_ -i'POS=501 && MIN(FMT/DP)>=50' $vcfout/Family.snp.gz >$tsvout/$_.tsv");
		}
	}
	if ($umi){
		system("$RealBin/bin/calposin-doublecpi-fornon.pl $tempDB $week_thres --vcf $vcfout/Family.snp.gz --gweek $week --errorate 0.0001 --thresh 0.001 > $tsvout/r01.result.cpibayes");
	}elsif($sur){
		system("$RealBin/bin/calposin-doublecpi-fornon-DUO.pl $tempDB $week_thres --vcf $vcfout/Family.snp.gz --gweek $week > $tsvout/r01.result.cpibayes");
	}else{
		system("$RealBin/bin/calposin-doublecpi-fornon.pl $tempDB $week_thres --vcf $vcfout/Family.snp.gz --gweek $week > $tsvout/r01.result.cpibayes");
	}
	system("$RealBin/bin/heterozygosity.pl $info/family.lst $tsvout $tsvout/r02.heterozygosity.txt");
	system("$RealBin/bin/contamination.MC.pl $info/family.lst $tsvout $tsvout/r03.contamination.MC.txt");
	system("$RealBin/bin/prepare.distribution.pl $info/family.lst $tsvout $info/cffDNA.txt $info/depth.txt");
	my $checkR = system("Rscript -e \"library(ggplot2)\"");
	my $checkP = system("perl -e \"use Statistics::TTest;\"");
	if ($checkP){
		print ER "! Can't locate Statistics/TTest.pm in \@INC. Failed to generate statistical graph.\n";
	}else{
		system("$RealBin/bin/cffDNA_ttest.pl $info/family.lst $tsvout/r04.ttest.txt");
	}
	if ($checkR){
		print ER "! There is no R-package called ‘ggplot2’. Failed to generate statistical graph.\n";
	}else{
		system("Rscript $RealBin/bin/ggplot_bar.cffdna.R $info/cffDNA.txt $tsvout/r05.cffDNA.pdf");
		system("Rscript $RealBin/bin/ggplot_bar.depth.R $info/depth.txt $tsvout/r06.locus_depth.pdf");
	}
}else{
	print ER "! mpileup failed, please check the log in $fqout.\n";
	die("Step2 failed!\n");
}
close ER;

##############################################################
sub rBWA {
	my ($tag,$fq1,$fq2) = @_;
	if ($fq2 eq "NA"){
		return "$RealBin/bin/bwa mem -t $ThreadsNum -Y $tempREF -R \"\@RG\\tID:$tag\\tSM:$tag\" $fq1 2>$fqout/$tag.log | samtools view -bS - | samtools sort -m 2G -T $fqout/$tag.tmp -o $fqout/$tag.bam";
	}else{
		return "$RealBin/bin/bwa mem -t $ThreadsNum -Y $tempREF -R \"\@RG\\tID:$tag\\tSM:$tag\" $fq1 $fq2 2>$fqout/$tag.log | samtools view -bS - | samtools sort -m 2G -T $fqout/$tag.tmp -o $fqout/$tag.bam";
	}
}

sub uBWA {
	my ($tag,$fq1,$fq2) = @_;
	my @tmp_order;
	if ($fq2 eq "NA"){
                push @tmp_order,"java -Xmx8G -jar $RealBin/bin/picard.jar FastqToSam F1=$fq1 O=$fqout/$tag.ubam SM=$tag RG=$tag SORT_ORDER=queryname ALLOW_AND_IGNORE_EMPTY_LINES=false CREATE_INDEX=false TMP_DIR=$fqout/../tmp";
                push @tmp_order,"java -Xmx8G -jar $RealBin/bin/fgbio-1.1.0.jar ExtractUmisFromBam --input=$fqout/$tag.ubam  --output=$fqout/$tag.withUMI.ubam --read-structure=$ReadStructure --molecular-index-tags=$IndexTag --single-tag=RX";
                push @tmp_order,"java -Xmx8G -jar $RealBin/bin/picard.jar SamToFastq I=$fqout/$tag.withUMI.ubam F=/dev/stdout INTERLEAVE=true TMP_DIR=$fqout/../tmp | $RealBin/bin/bwa mem -p -t $ThreadsNum $tempREF /dev/stdin | java -Xmx8G -jar $RealBin/bin/picard.jar MergeBamAlignment ALIGNED=/dev/stdin UNMAPPED=$fqout/$tag.withUMI.ubam O=$fqout/$tag.mbam R=$tempREF SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=true MAX_GAPS=-1 ORIENTATIONS=FR VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true TMP_DIR=$fqout/../tmp";
                push @tmp_order,"java -Xmx8G -Djava.io.tmpdir=$fqout/../tmp -jar $RealBin/bin/fgbio-1.1.0.jar GroupReadsByUmi --input=$fqout/$tag.mbam --output=$fqout/$tag.gbam --strategy=Adjacency --edits=1 --min-map-q=30";
                push @tmp_order,"java -Xmx8G -Djava.io.tmpdir=$fqout/../tmp -jar $RealBin/bin/fgbio-1.1.0.jar CallMolecularConsensusReads --input=$fqout/$tag.gbam --output=$fqout/$tag.cubam --error-rate-post-umi=30 --min-reads=1";
                push @tmp_order,"java -Xmx8G -Djava.io.tmpdir=$fqout/../tmp -jar $RealBin/bin/fgbio-1.1.0.jar FilterConsensusReads --input=$fqout/$tag.cubam --output=$fqout/$tag.fubam --ref=$tempREF --reverse-per-base-tags=true --min-reads=5 -E 0.05 -N 40 -e 0.1 -n 0.1";
                push @tmp_order,"java -Xmx8G -jar $RealBin/bin/picard.jar SortSam I=$fqout/$tag.fubam O=$fqout/$tag.sort.fubam SORT_ORDER=queryname TMP_DIR=$fqout/../tmp";
                push @tmp_order,"java -Xmx8G -jar $RealBin/bin/picard.jar SamToFastq I=$fqout/$tag.sort.fubam  F=/dev/stdout INTERLEAVE=true TMP_DIR=$fqout/../tmp | $RealBin/bin/bwa mem -p -t $ThreadsNum $tempREF /dev/stdin | java -Xmx8G -jar $RealBin/bin/picard.jar MergeBamAlignment ALIGNED=/dev/stdin UNMAPPED=$fqout/$tag.sort.fubam O=$fqout/$tag.sort.fmbam R=$tempREF SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=true MAX_GAPS=-1 ORIENTATIONS=FR VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true TMP_DIR=$fqout/../tmp";		
	}else{
		push @tmp_order,"java -Xmx8G -jar $RealBin/bin/picard.jar FastqToSam F1=$fq1 F2=$fq2 O=$fqout/$tag.ubam SM=$tag RG=$tag SORT_ORDER=queryname ALLOW_AND_IGNORE_EMPTY_LINES=false CREATE_INDEX=false TMP_DIR=$fqout/../tmp";
		push @tmp_order,"java -Xmx8G -jar $RealBin/bin/fgbio-1.1.0.jar ExtractUmisFromBam --input=$fqout/$tag.ubam  --output=$fqout/$tag.withUMI.ubam --read-structure=$ReadStructure --molecular-index-tags=$IndexTag --single-tag=RX";
		push @tmp_order,"java -Xmx8G -jar $RealBin/bin/picard.jar SamToFastq I=$fqout/$tag.withUMI.ubam F=/dev/stdout INTERLEAVE=true TMP_DIR=$fqout/../tmp | $RealBin/bin/bwa mem -p -t $ThreadsNum $tempREF /dev/stdin | java -Xmx8G -jar $RealBin/bin/picard.jar MergeBamAlignment ALIGNED=/dev/stdin UNMAPPED=$fqout/$tag.withUMI.ubam O=$fqout/$tag.mbam R=$tempREF SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=true MAX_GAPS=-1 ORIENTATIONS=FR VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true TMP_DIR=$fqout/../tmp";
		push @tmp_order,"java -Xmx8G -Djava.io.tmpdir=$fqout/../tmp -jar $RealBin/bin/fgbio-1.1.0.jar GroupReadsByUmi --input=$fqout/$tag.mbam --output=$fqout/$tag.gbam --strategy=Adjacency --edits=1 --min-map-q=30";
		push @tmp_order,"java -Xmx8G -Djava.io.tmpdir=$fqout/../tmp -jar $RealBin/bin/fgbio-1.1.0.jar CallMolecularConsensusReads --input=$fqout/$tag.gbam --output=$fqout/$tag.cubam --error-rate-post-umi=30 --min-reads=1";
		push @tmp_order,"java -Xmx8G -Djava.io.tmpdir=$fqout/../tmp -jar $RealBin/bin/fgbio-1.1.0.jar FilterConsensusReads --input=$fqout/$tag.cubam --output=$fqout/$tag.fubam --ref=$tempREF --reverse-per-base-tags=true --min-reads=5 -E 0.05 -N 40 -e 0.1 -n 0.1";
		push @tmp_order,"java -Xmx8G -jar $RealBin/bin/picard.jar SortSam I=$fqout/$tag.fubam O=$fqout/$tag.sort.fubam SORT_ORDER=queryname TMP_DIR=$fqout/../tmp";
		push @tmp_order,"java -Xmx8G -jar $RealBin/bin/picard.jar SamToFastq I=$fqout/$tag.sort.fubam  F=/dev/stdout INTERLEAVE=true TMP_DIR=$fqout/../tmp | $RealBin/bin/bwa mem -p -t $ThreadsNum $tempREF /dev/stdin | java -Xmx8G -jar $RealBin/bin/picard.jar MergeBamAlignment ALIGNED=/dev/stdin UNMAPPED=$fqout/$tag.sort.fubam O=$fqout/$tag.fmbam R=$tempREF SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=true MAX_GAPS=-1 ORIENTATIONS=FR VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true TMP_DIR=$fqout/../tmp";
	}
	return join "&&",@tmp_order;
}
