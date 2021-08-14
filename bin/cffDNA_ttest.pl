#!/usr/bin/env perl
use strict;
use Statistics::TTest;
use List::Util qw/sum/;

my ($list,$output) = @ARGV;

my @path = split /\//,$list;
$path[$#path] = "../3result";
my $store = join "/",@path;

open LI,"<$list" or die($!);
open OUT,">$output" or die($!);
print OUT "#Paternity inclusion: An inclusion is normally reported with, 1) a parentage index of 1e+04 or more, 2) power of exclusion of 99.99% or more,\n";
print OUT "#Paternity exclusion: An exclusion is normally reported with, 1) a parentage index of 1e-04 or less, 2) power of exclusion of 99.99% or more,\n";
print OUT "#Unknow: An unknow is normally reported as a sample with an ambiguous tendency in paternity test.\n\n";
print OUT "Paternity_Test\tChildID\tcffDNA_concentration\tCPE\tCPI\tHomoNum\tHomoContami\tHomoZero\tHeteNum\tHeteContami\tHeteZero\tT-Score\tPvalue\n";
while (my $line = <LI>){
	chomp($line);
	my @mem = split /\s+/,$line;
	if (-e "$store/$mem[2].tsv"){
		my (@homoCon,@heteCon);
		my ($homoLost,$heteLost) = (0,0);
		my $fileF = "$store/$mem[1].tsv";
		my $fileM = "$store/$mem[0].tsv";
		my $fileC = "$store/$mem[2].tsv";
		my $result = "$store/r01.result.cpibayes";
		my ($recffDNA,$cpi,$cpe) = &read_cpi($result);

		my (%homoM,%homoF);
		&get_homo(\%homoM,$fileM);
		&get_homo(\%homoF,$fileF);

		open CH,"<$fileC" or die($!);
		while (<CH>){
			chomp;
			my @data = split /\t/,$_;
			next unless (defined $homoM{$data[0]} && defined $homoF{$data[0]});
			next if ($data[3] eq '.' or $data[3] < 100);
			my @tM = splice @data,4;
			my %Dep;
			my $depcheck = 1;
			my @bases = split /,/,$data[2];
			my $sum = 0;
			for (@tM){
				my $depsum = 0;
				my @Depinfo = split /[;,]/,$_;
				for my $i (1..scalar @Depinfo - 1){
					if ($Depinfo[$i] eq '.'){$Depinfo[$i] = 0;}
					$Dep{$bases[$i - 1]} += $Depinfo[$i];
					$depsum += $Depinfo[$i];
				}
				$sum+=$depsum;
				if ($depsum <= 50){
					$depcheck *= 0;
				}else{
					$depcheck *= 1;
				}
			}
			next if ($depcheck == 0);

			if ($homoM{$data[0]} eq $homoF{$data[0]}){
				my $totalCon = 0;
				for (@bases){
					unless ($_ eq $homoM{$data[0]}){
						$totalCon += $Dep{$_} / $sum;
					}
				}
				if ($totalCon == 0){
					$homoLost++;
				}
				push @homoCon,$totalCon;
			}else{
				my $totalCon = 0;
                                for (@bases){
                                        unless ($_ eq $homoM{$data[0]}){
                                                $totalCon += $Dep{$_} / $sum;
					}
				}
				if ($totalCon == 0){
					$heteLost++;
				}
				push @heteCon,$totalCon;
			}
		}
		close CH;

		my $homo_sum = sum @homoCon;
		my $hete_sum = sum @heteCon;
		my $homo_num = scalar @homoCon;
		my $hete_num = scalar @heteCon;
		my ($homo_ave,$hete_ave);
		if ($homo_num == 0){
			$homo_ave = "NA";
		}else{
			$homo_ave = sprintf("%.3f", $homo_sum / $homo_num);
		}
		if ($hete_num == 0){
			$hete_ave = "NA";
		}else{
			$hete_ave = sprintf("%.3f", $hete_sum / $hete_num);
		}

                my $ttest = new Statistics::TTest;
                $ttest->load_data(\@homoCon,\@heteCon);
		$ttest->set_significance(99);

		if ($ttest->t_statistic > 5 && $ttest->t_statistic < 10){
			print OUT "# ! Warning: the T-score is between 5 and 10, means a less creditable result in paternity test. This result should be considered with caution.\n\n";
		}

		if ($cpi > 1e+04 && $cpe > 1 - 1e-04){
			if ($ttest->t_statistic < 5){
				print OUT "# ! Error: the T-score is less than 5, means an unreliable result in paternity test.\n\n";
			}
			print OUT "Paternity inclusion\t$mem[2]\t$recffDNA\t$cpe\t$cpi\t$homo_num\t$homo_ave\t$homoLost\t$hete_num\t$hete_ave\t$heteLost\t",$ttest->t_statistic,"\t",$ttest->{t_prob},"\n";
		}elsif ($cpi < 1e-04 && $cpe > 1 - 1e-04){
			if ($ttest->t_statistic < 10){
				print OUT "# ! Error: the T-score is more than 10, means an unreliable result in paternity test.\n\n";
			}
			print OUT "Paternity exclusion\t$mem[2]\t$recffDNA\t$cpe\t$cpi\t$homo_num\t$homo_ave\t$homoLost\t$hete_num\t$hete_ave\t$heteLost\t",$ttest->t_statistic,"\t",$ttest->{t_prob},"\n";
		}else{
			print OUT "Unknown\t$mem[2]\t$recffDNA\t$cpe\t$cpi\t$homo_num\t$homo_ave\t$homoLost\t$hete_num\t$hete_ave\t$heteLost\t",$ttest->t_statistic,"\t",$ttest->{t_prob},"\n";
		}
	}
}
close LI;
close OUT;
##################################################
sub get_homo {
        my ($hash,$tempin) = @_;
        open TE,"<$tempin" or die($!);
        while (my $line = <TE>){
                chomp($line);
                my @data = split /\t/,$line;
                next if ($data[3] eq '.' or $data[3] < 100);
                my @tM = splice @data,4;
                my %Dep;
                my $depcheck = 1;
                my @bases = split /,/,$data[2];
                for (@tM){
                        my $depsum = 0;
                        my @Depinfo = split /[;,]/,$_;
                        for my $i (1..scalar @Depinfo - 1){
                                if ($Depinfo[$i] eq '.'){$Depinfo[$i] = 0;}
                                $Dep{$bases[$i - 1]} += $Depinfo[$i];
                                $depsum += $Depinfo[$i];
                        }
                        if ($depsum <= 50){
                                $depcheck *= 0;
                        }else{
                                $depcheck *= 1;
                        }
                }
                next if ($depcheck == 0);

                my @values = sort {$Dep{$b}<=>$Dep{$a}} keys %Dep;
                if (scalar @values > 1){
                        if ($Dep{$values[1]} <= $Dep{$values[0]} * 0.01){
                                $hash->{$data[0]} = $values[0];
                        }
                }else{
                        $hash->{$data[0]} = $values[0];
                }
        }
        close TE;
}

sub read_cpi {
	my $temp = shift;
	my %pos;
	open TE,"<$temp" or die($!);
	my $head = <TE>;
	chomp($head);
	my @title = split /\t/,$head;
	for my $i (0..$#title){
		$pos{$title[$i]} = $i;
	}
	while (<TE>){
		chomp;
		my @data = split /\t/,$_;
		return ($data[$pos{'fetra'}],$data[$pos{'cpi'}],$data[$pos{'cpe'}]);
	}
	close TE;
}
