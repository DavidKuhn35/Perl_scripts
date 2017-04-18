#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: May 9,2013
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	:  Reads in ~/avocado/SNPs/fixit_avo_May9_2013_seq_genotypes_out, which has seq genotypes no Xs; reads in ~/avocado_Illumina_data/Avo1_14f_6kSNP_list_for_Illumina.txt
#		which has SNPnames for SNP chip so seq genotypes can be identified with chip SNPname; reads in ~/avocado_Illumina_data/5050chip_geno_seq_SNP_name.txt which has chip 
#		genotypes for just the 5050 good SNPs on the chip and the seq SNP names to pull seq genotypes out of fixit
# <input>		:  see above		
# <output>		:  chip_seq_comparison_table.txt which has the seq genotypes and het no het info for all chip SNPs and the chip genotypes and het no het info for the 5050 SNPs 
# 			
# Dependencies	: N/A
# Usage			: 
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
my $seqsnpname = ();
my %chipgeno = ();
my %seqgeno = ();
my %chipseqgeno = ();


###################################################################################################################
#####################################################################################################################
#	read in the chip genotypes file from ~/avocado_Illumina_data/5050chip_geno_seq_SNP_name.txt to get chip genotypes and link chip SNP name to seq SNP name	
open (CHIP5050, "/home/david/avocado_Illumina_data/5050chip_geno_seq_SNP_name.txt") || die "Can't open chip genotype data file $!\n";
	while (<CHIP5050>)
{
	chomp;
	if (m/^SHRSPa/)	{
		my @x = split /\t/;
		$chipgeno{$x[1]} = join "\t", @x; #hash has seq name and chip geno info
		}
}
close CHIP5050;

#####################################################################################################################################
#	read in the seq genotypes from ~/avocado/SNPs/fixit_avo_May9_2013_seq_genotypes_out
 open (SEQ, "/home/david/avocado/SNPs/fixit_avo_May9_2013_seq_genotypes_out") || die "Can't open seq genotype data file $!\n";
	while (<SEQ>)
{
	chomp;
	if (m/^g/)	{
		my @y = split /\t/;
		if (exists $chipgeno{$y[1]})	{
			$seqgeno{$y[1]} = join "\t", @y; #hash of 5050 seq genotypes based on seq snp name
		}
	}
}
close SEQ;
					
				
######################################################################################################################
#   calculate het no het info for seqgeno hash
	foreach my $key (sort keys %seqgeno)	{
		my @z = split /\t/, $seqgeno{$key};
			if ($z[2] eq $z[3])	{
				push (@z,"No Het Bacon");
			}else{
				push (@z,"Het Bacon");
			}
			if ($z[5] eq $z[6])	{
				push (@z,"No Het Hass");
			}else{
				push (@z,"Het Hass");
			}
			if ($z[8] eq $z[9])	{
				push (@z,"No Het Simmonds");
			}else{
				push (@z,"Het Simmonds");
			}
			if ($z[11] eq $z[12])	{
				push (@z, "No Het Tonnage");
			}else{
				push (@z, "Het Tonnage");
		}
		$seqgeno{$key} = join "\t", @z;
}

#######################################################################################################################
#  Put the two hashes together 
	foreach my $x (sort keys %chipgeno)	{
			$chipseqgeno{$x} = "$chipgeno{$x}\t$seqgeno{$x}";
		}

#########################################################################################################################
#	print out the joined hashes
	open (OUT, ">>compare_seq_chip_May9_2013_out.txt") || die "Can't open out file $!\n";
		foreach my $key (sort keys %chipseqgeno)	{
			print OUT "$chipseqgeno{$key}\n";
		}
##########################################################################################################################
#

