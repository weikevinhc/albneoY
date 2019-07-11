#!usr/bin/perl
# vcf.shared.poly.pl
# parse vcfs that have multiple individuals for the number of shared polymorphism at biallelic sites
use strict;
use warnings;
use POSIX;

my $input1 = $ARGV[0]; #sample 1 vcf
my @pop1 = split ",", $ARGV[1]; #comma separated list of group 1
my @pop2 = split ",", $ARGV[2]; #comma separated list of group 2
my $out = $ARGV[3]; # basename of output

open VCF, $input1 or die;

my %pos;

my @pop1dex;
my @pop2dex;

my $win = 50000; #window size;


my %shared_hash;
while (<VCF>) {
	s/[\r\n]+$//;
	my @line = split("\t", $_);
	if ($_ =~ /^#/) {
		if ($_ =~ /^#C/) {
			my %sampledex;
			foreach my $i (9..(scalar(@line)-1)) {
				$sampledex{$line[$i]} = $i;
			}
			foreach my $p (@pop1) {
				if ($sampledex{$p}) {
					push @pop1dex, $sampledex{$p};
				} else {
					die "cannot find $p\n";
				}
			}
			foreach my $p (@pop2) {
				if ($sampledex{$p}) {
					push @pop2dex, $sampledex{$p};
				} else {
					die "cannot find $p\n;"
				}
			}
		}
		next;
	}
	my @format = split(":", $line[8]);
	next if (scalar @format < 5);
	if ($format[0] ne "GT" || $format[1] ne "AD" || $format[3] ne "GQ") {
#		warn "odd GT: $_\n";
		next;
	} ## examine that the formatting field of the genotypes are approrpiate
	my %all_alleles;
	
	my %pop1_alleles;
	foreach my $p (@pop1dex) {
		my @field = split ":", $line[$p];
		my @GT = split "\/", $field[0];
		foreach my $gt (@GT) {
			next if ($gt eq ".");
			$pop1_alleles{$gt} ++;
		}
	}
	
	my %pop2_alleles;
	foreach my $p (@pop2dex) {
		my @field = split ":", $line[$p];
		my @GT = split "\/", $field[0];
		foreach my $gt (@GT) {
			next if ($gt eq ".");
			$pop2_alleles{$gt} ++;
		}
	}
	
	foreach my $gt (keys %pop1_alleles) {
		if ($pop1_alleles{$gt} <= 2) {
			delete $pop1_alleles{$gt}
		} else {
			$all_alleles{$gt} ++;
		}
	} ## remove singletons
	foreach my $gt (keys %pop2_alleles) {
		if ($pop2_alleles{$gt} <= 2) {
			delete $pop2_alleles{$gt}
		} else {
			$all_alleles{$gt} ++;
		}
	} ## remove singletons
	
	next if (scalar keys %all_alleles > 2); ## ignore sites with more than two alleles
	next if (scalar keys %pop1_alleles == 0 || scalar keys %pop2_alleles == 0); ## ignore sites with more than two alleles
	my $shared = 0;
	
	foreach my $p (keys %pop1_alleles) {
		$shared ++ if ($pop2_alleles{$p});
	}
	
	if ($shared >= 2) {
#		print "p1:";
#		for my $i (@pop1dex) {
#			my @field = split ":", $line[$i];
#			print "\t$field[0]";
#		}
#		print "\np2:";
#		for my $i (@pop2dex) {
#			my @field = split ":", $line[$i];
#			print "\t$field[0]";
#		}
		$shared_hash{$line[0]}{$line[1]} = 1; ## 1 = shared polymorphism in both
#		print "\n\n";
	} elsif ($shared == 1) {
		if (scalar keys %pop1_alleles == 1 && scalar keys %pop2_alleles == 1) {
		} elsif (scalar keys %pop1_alleles > 1 && scalar keys %pop2_alleles > 1) {
			$shared_hash{$line[0]}{$line[1]} = 4; ## 4 = unclear
		} elsif (scalar keys %pop1_alleles > 1 && scalar keys %pop2_alleles == 1) {
			$shared_hash{$line[0]}{$line[1]} = 2; ## 2 = fixed in neo-Y
		} elsif (scalar keys %pop1_alleles == 1 && scalar keys %pop2_alleles > 1) {
			$shared_hash{$line[0]}{$line[1]} = 3; ## 3 = fixed in neo-X
		}
	} else {
		if (scalar keys %pop1_alleles == 1 && scalar keys %pop2_alleles == 1) {
			$shared_hash{$line[0]}{$line[1]} = 0; ## 0 = fixed in both
		} else {
			$shared_hash{$line[0]}{$line[1]} = 4; ## 4 = unclear
		}
		
	}
}

open OUT1, ">$out.poly.sites";
my %win_hash;
foreach my $chr (keys %shared_hash) {
	foreach my $pos (sort {$a <=> $b} keys %{$shared_hash{$chr}}) {
		print OUT1 "$chr\t$pos\t$shared_hash{$chr}{$pos}\n";
		my $window = int($pos/$win) * $win;
		
		$win_hash{$chr}{$window}{$shared_hash{$chr}{$pos}} ++;
	}
}

open OUT2, ">$out.poly.window";
print OUT2 "chr\tpos\tfixed_diff\tshared_poly\tfixed_in_B\tfixed_in_A\tunknown\n";
foreach my $chr (keys %shared_hash) {
	foreach my $w (sort {$a <=> $b} keys %{$win_hash{$chr}}) {
		print OUT2 "$chr\t$w";
		foreach my $p (0..4) {
			if (exists $win_hash{$chr}{$w}{$p}) {
				print OUT2 "\t$win_hash{$chr}{$w}{$p}";
			} else {
				print OUT2 "\t0";
			}
		}
		print OUT2 "\n";
	}
}
