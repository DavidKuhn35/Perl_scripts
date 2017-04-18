#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: June 27, 2012
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: makes a composite file of differences between the SNP calls for B,H,S,T from chip and sequencing
# <input>		: four diff files for B,H,S,T 		
# <output>		:    one file based on SNPs
# 			
# Dependencies	: N/A
# Usage			: avo_diff_problems.pl
#                    
##################################################################################
use strict;
use warnings;
use diagnostics;
use Bio::Seq; 
use Bio::SeqIO;
use Bio::SeqUtils;
use Bio::Tools::OddCodes;
use Bio::PrimarySeq;
my %diff = ();
my %counter = ();

#####################################################################################################################
#Read in diff files that have been concatenated avo_bhst_diff 2RRR_chip_seq_compare_out
open (IN,"2RRR_chip_seq_compare_out") || die "Can't open avo_bhst_diff infile $!\n";
	while (<IN>)	
{
	chomp;
		if (m/^SHRSPa/)	{
		my @x = split /\t/;
		my $z = shift @x;
		my $y = join "\t", @x;
		$diff{$z} .= "$y\t";
		$counter{$z}++;
		}
}
####################################################################################################################
open (OUT, ">>2AAA_avo_bhst_diff_out") || die "Can't open avo_bhst_diff outfile $! \n";
	my $snptotal = keys %diff;
	print OUT "The total number of SNPs in the differences between chip and seq is $snptotal.\n";	
#	foreach my $key (sort keys %diff)	{
#		print OUT "$key\t$counter{$key}\t$diff{$key}\n";
#	}
	my @problems = 	sort by_number_and_snp keys %counter;
	my $problemtotal = @problems;
	for (my $i=0; $i < $problemtotal; $i++)	{
	print OUT "$problems[$i]\t$counter{$problems[$i]}\t$diff{$problems[$i]}\n";
	}
###################################################################################################################
#  subroutine to sort counter hash by count and by snp name 
sub by_number_and_snp  {
	$counter{$a} <=> $counter{$b} 
	or
 	$a cmp $b
	}
		
