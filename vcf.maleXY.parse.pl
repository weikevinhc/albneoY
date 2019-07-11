#!usr/bin/perl
# vcf.maleXY.parse.pl
# determine neoY alleles from VCF with female and male genotypes
# the female must be first and male second in the vcf (column 10 and 11 respectively)
use strict;
use warnings;
use POSIX;

my $input1 = $ARGV[0]; # VCF must be female first male second

open VCF1, $input1 or die;

while (<VCF1>) {
	s/[\r\n]+$//;
	if ($_ =~ /^##/) {
		print "$_\n";
		next;
	}
	my @line = split("\t", $_);
	if ($_ =~ /^#CHROM/) {
		$line[9] = $line[10] . "X";
		$line[10] = $line[10] . "Y";
		print join("\t", @line), "\n";
		next;
	}
	my @format = split(":", $line[8]);
	if ($format[0] ne "GT" || $format[1] ne "AD") {
		warn "odd GT: $_\n";
		next;
		
	}
	my @alt = split(",", $line[4]);
	
	my @female = split(":", $line[9]);
	my @GTfemale = split("\/", $female[0]);
	
	my @male = split(":", $line[10]);
	my @GTmale = split("\/", $male[0]);
#	print "$_\n";
	print join("\t", @line[0..8]), "\t";
	if ($GTmale[0] ne $GTmale[1]) { #heterozygous male genotype
		if ($GTfemale[0] eq $GTfemale[1] && ($GTmale[0] eq $GTfemale[0] || $GTmale[1] eq $GTfemale[0])) { #homozygous female genotype
			my @AD = split(",", $male[1]);
			my @AD_X = ((0) x (scalar(@alt)+1));
			$AD_X[$GTfemale[0]] = $AD[$GTfemale[0]];
			my @AD_Y = ((0) x (scalar(@alt)+1));
			if ($GTmale[0] eq $GTfemale[0]) {
				$AD_Y[$GTmale[1]] = $AD[$GTmale[1]];
				print "$GTfemale[0]\/$GTfemale[0]:", join(",", @AD_X), ":$AD_X[$GTfemale[0]]:$male[3]:$male[4]\t";
				print "$GTmale[1]\/$GTmale[1]:", join(",", @AD_Y), ":$AD_Y[$GTmale[1]]:$male[3]:$male[4]\n";
			} elsif ($GTmale[1] eq $GTfemale[0]) {
				$AD_Y[$GTmale[0]] = $AD[$GTmale[0]];
				print "$GTfemale[1]\/$GTfemale[1]:", join(",", @AD_X), ":$AD_X[$GTfemale[1]]:$male[3]:$male[4]\t";
				print "$GTmale[0]\/$GTmale[0]:", join(",", @AD_Y), ":$AD_Y[$GTmale[0]]:$male[3]:$male[4]\n";
#				<STDIN>;
			}
		} else {
			print join(":", @male), "\t", ".\/.", "\n";
#			print "$_\n";
#			<STDIN>;
		}
	} else {
			print join(":", @male), "\t", join(":", @male), "\n";
	}
}

