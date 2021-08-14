#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BGI <huxuesong@genomics.org.cn>
Version: 2.0.0 @ 20190727 update by libowen
=cut
use strict;
use warnings;
use POSIX;
use FindBin qw($RealBin);
use List::Util qw/sum/;
use bignum;
if ($FindBin::VERSION < 1.51) {
        warn "[!]Your Perl is too old, thus there can only be ONE `oykb.pl` file in your PATH. [FindBin Version: $FindBin::VERSION < 1.51]\n\n"
}
FindBin::again();
use lib "$RealBin/../";
require mBCPI::BayesCPI;

my $DBsuffix = shift;

our @Bases;

sub getTPE(@) {
        my @AFs = @_;
        my ($p1,$p2)=(0,0);
        for my $i (0 .. $#AFs) {
		$p1 += $AFs[$i]**2;
		$p2 += $AFs[$i] * ((1-$AFs[$i])**2);
        }
        return $p1*$p2;
}
sub getDPE(@) {
        my @AFs = @_;
        my ($p1,$p2)=(0,0);
        for my $i (0 .. $#AFs) {
		$p1 += $AFs[$i]**2;
		$p2 += ($AFs[$i]**2) * ((1-$AFs[$i])**2);
        }
        return $p1*$p2;
}
our (%Markers,%MarkerAF);
open DB,'<',"$DBsuffix" or die $?;
while (<DB>) {
        next if /^#/;
        chomp;
        my ($rs,$chr,$pos,$chr38,$pos38,@d) = split /\t/,$_;
        $MarkerAF{$rs} = {@d};
        my @AFs = values %{$MarkerAF{$rs}};
	my $TPE = getTPE(@AFs);
        my $DPE = getDPE(@AFs);
        $Markers{$rs} = [$chr,$pos,$DPE,$TPE];
}
close DB;

sub get_homo {
        my ($hash,$file,$tag) = @_;
        my $homo = 0;
	my %temp;
        open TE,"<$file" or die($!);
        while (my $line = <TE>){
                chomp($line);
                my @data = split /\t/,$line;
                next if ($data[3] eq '.' or $data[3] < 20);
		next if ($data[2] =~ /[ACGTN][ACGTN]+/);
                my @tM = splice @data,4;
                my %Dep;
		@Bases = split /,/,$data[2];
        	my $sum = &depcheck(\%Dep,\@tM,50);
        	next if ($sum eq "NA");

                my @alleles = sort {$Dep{$b}<=>$Dep{$a}} keys %Dep;
                if (scalar @alleles > 1){
                        if ($Dep{$alleles[1]} <= $Dep{$alleles[0]} * 0.01){
                                $hash->{$data[0]} = $alleles[0];
				$temp{$data[0]} = join "/",$alleles[0],$alleles[0];
                        }elsif ($Dep{$alleles[1]} >= $Dep{$alleles[0]} * 0.1){
				$temp{$data[0]} = join "/",$alleles[0],$alleles[1];
			}
                }else{
                        $hash->{$data[0]} = $alleles[0];
			$temp{$data[0]} = join "/",$alleles[0],$alleles[0];
                }
        }
        close TE;
	if ($tag){
		return %temp;
	}
}

sub depcheck {
	my ($hash,$array,$cutoff) = @_;
	my $resum = 0;
	my $checkNum = 1;
	for (@$array){
		my $depsum = 0;
		my @Depinfo = split /[;,]/,$_;
		for my $i (1..scalar @Depinfo - 1){
			if ($Depinfo[$i] eq '.'){$Depinfo[$i] = 0;}
			$hash->{$Bases[$i - 1]} += $Depinfo[$i];
			$depsum += $Depinfo[$i];
		}
		if ($depsum < $cutoff){
			$checkNum *= 0;
		}else{
			$checkNum *= 1;
		}
		$resum += $depsum;
	}
	if ($checkNum == 1){
		return $resum;
	}else{
		return "NA";
	}
}

sub printExp($) {
        my $lnV = $_[0]/log('10'); # bignum doesn't properly handle log; Math::BigInt->new(10);
        my $lnInt = floor($lnV);
        my $lnExt = $lnV - $lnInt;
        my $prefix = sprintf("%.3f", exp($lnExt*log('10')));
        my $str = join('e',$prefix,$lnInt);
        return $str;
}

################################################################
my $Quality=shift;
my $ErrP = 10**(0 - $Quality / 10);

my $mother=shift;
my $father=shift;
my $child=shift;
my $outprefix=shift;

my @RefBases = qw{A C G T};

my (%homoM,%genoM,%homoF,%genoF);
%genoM = &get_homo(\%homoM,$mother,1);
%genoF = &get_homo(\%homoF,$father,1);


my @cffDNA;;
my %child;
open FC,"<$child" or die($!);
while (my $line = <FC>){
	chomp($line);
	my @data = split /\t/,$line;
	next if ($data[3] eq '.' or $data[3] < 22);
	next if ($data[2] =~ /[ACGTN][ACGTN]+/);
	my @tmpInfo = splice @data,4;
	my %Dep;
	@Bases = split /,/,$data[2];
	my $sum = &depcheck(\%Dep,\@tmpInfo,200);
        next if ($sum eq "NA");

	for (@RefBases){
		unless (defined $Dep{$_}){
			$Dep{$_} = 0;
		}
	}

	my @refs = sort {$Dep{$b} <=> $Dep{$a}} keys %Dep;

        if (defined $homoM{$data[0]} && defined $homoF{$data[0]} && $homoM{$data[0]} ne $homoF{$data[0]}){
                my $value = $Dep{$homoF{$data[0]}} / $sum * 2;
		push @cffDNA,$value;
	}
	$child{$data[0]} = join("\t",join("\t",@data),join(",",$Dep{A},$Dep{C},$Dep{G},$Dep{T}));
}
close FC;

## useless
my ($Peq,$Pne,$Chomo) = (0,0,0,0);
foreach my $locus (keys %homoM){
	if (defined $homoF{$locus} && defined $child{$locus}){
		if ($homoM{$locus} eq $homoF{$locus}){
			$Peq++;
		}else{
			$Pne++;
			my @Cdata = split /\t/,$child{$locus};
			my @depths = split /,/,$Cdata[4];
			my @sortDep = sort {$b <=> $a} @depths;
			if ($sortDep[1] / $sortDep[0] <= 0.001){
				$Chomo++;
			}
		}
	}
}
my $forceHETE = $Chomo / ($Peq + $Pne);

open OUT,">$outprefix.cpibayes" or die($!);
if (scalar @cffDNA > 0){
	my $sum_dna = sum @cffDNA;
	my $ave = sprintf("%.3f", $sum_dna / scalar @cffDNA);
	my $logcpi = 0;
	my $spe = 0;
	foreach my $marker (sort keys %child){
		my $switch = 0;
		if (defined $homoM{$marker} && defined $genoF{$marker}){
			my @data = split /\t/,$child{$marker};
			my @depths = split /,/,$data[4];
			my @sortDep = sort {$b <=> $a} @depths;
			if (defined $homoF{$marker} && defined $homoM{$marker} && $homoF{$marker} ne $homoM{$marker}){
				if ($sortDep[1] / $sortDep[0] <= 0.001){
					$switch = 1;
				}
			}
			my $TotalDep = sum @depths;
			my %reDep;
			my %need_alleles;
			for my $i (0..scalar @RefBases - 1){
				if ($RefBases[$i] eq $homoM{$marker}){
					$reDep{$RefBases[$i]} = int($depths[$i] - $TotalDep * (1-$ave));
				}else{
					$reDep{$RefBases[$i]} = $depths[$i];
				}
				if ($reDep{$RefBases[$i]} <= 0){
					$reDep{$RefBases[$i]} = 0;
				}else{
					$need_alleles{$RefBases[$i]}++;
				}
			}
			my @alleleF = split /\//,$genoF{$marker};
			for (@alleleF){
				$need_alleles{$_}++;
			}
			my %GT;
			foreach my $allele1 (sort keys %need_alleles){
				foreach my $allele2 (sort keys %need_alleles){
					my @geno = sort ($allele1,$allele2);
					my $genotype = join "/",@geno;
					$GT{$genotype}++;
				}
			}
			my $le_sum = 0;
			my %record;
			my $pe;
			foreach my $genotype (keys %GT){
				my $cret = getcpiT($marker,"NA","A,C,G,T","NA","$homoM{$marker}/$homoM{$marker}",$genoF{$marker},$genotype,$ErrP);
				$record{pi}{$genotype} = $cret->[0];
				$pe = $cret->[1];
				if ($switch == 1){
					$le_sum = 1;
					if ($genotype eq "$homoM{$marker}/$homoM{$marker}"){
						$record{le}{$genotype} = 1 - $forceHETE;
					}elsif ($genotype =~ /$homoM{$marker}/ && $genotype =~ /$homoF{$marker}/){
						$record{le}{$genotype} = $forceHETE;
					}else{
						$record{pi}{$genotype} = $cret->[0];
						$record{le}{$genotype} = 0;
					}
					next;
				}
				$record{pi}{$genotype} = $cret->[0];
				$record{le}{$genotype} = getle(\%reDep,$marker,$genotype,$ErrP);
				$le_sum += $record{le}{$genotype};
			}
			my $corpi = 0;
			my @look;
			foreach my $genotype (sort keys %{$record{pi}}){
				$corpi += $record{pi}{$genotype} * ($record{le}{$genotype} / $le_sum);
				my $freq = $record{le}{$genotype} / $le_sum;
				push @look,"$genotype ($freq)";
			}
			#print "$marker\t$corpi\t";
			#print join(",",@look),"\n";
			$logcpi += log($corpi);
			$spe += log(1-$pe);
		}
	}
	my $sCPI = printExp($logcpi);
	my $sPCPE = printExp($spe);
	print OUT "# cffDNA: $ave\n";
	print OUT "# CPE: 1-$sPCPE\n";
	print OUT "# correct CPI: $sCPI\n";
	#print "# cffDNA: $ave\n";
	#print "# correct CPI: $sCPI\n";
}else{
	print OUT "# abandoned\n";
	#print "# abandoned\n";
}
close OUT;
