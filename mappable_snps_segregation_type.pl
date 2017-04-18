#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: November 6, 2012
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: 	Take in list of parents genotypes for a mapping population.  Remove unmappable SNPs and designate segregation type.
# <input>		:  List of parents genotypes with missing data removed	
# <output>		:   List of mappable SNPs and segregation type for each SNP locus
# 			
# Dependencies	: N/A
# Usage			: avo_mappable_snps_segregation_type.pl
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
my $snp = 0;
my $parent1 = 0;
my $parent2 = 0;
my %map_snp = ();
my %parent1_hash = ();
my %parent2_hash = ();
my %segregation_type = ();
###################################################################################################################
#	Define parent 1 and parent 2 from ARGV [0] and [1]
	my $parent1_name = "$ARGV[0]";
	my $parent2_name = "$ARGV[1]";

#####################################################################################################################
#	Read in list of parents genotypes.  parent1_parent2_consensus_genotypes.txt
open (IN,"$ARGV[2]") || die "Can't open parents genotypes file $!\n";
	while (<IN>)		{
		chomp;
		if (m/^SHRSPa/)	{	#contig for cacao, SHRSPa for avocado
			($snp, $parent1, $parent2) = split /\t/; # tab for cacao, comma for avocado
			unless ((($parent1 eq "AA") || ($parent1 eq "CC") || ($parent1 eq "GG") || ($parent1 eq "TT") || ($parent1 eq "--")) && (($parent2 eq "AA")|| ($parent2 eq "CC") || ($parent2 eq "GG") || ($parent2 eq "TT") || ($parent2 eq "--")))	{
		$map_snp{$snp} = $snp;
		$parent2_hash{$snp} = $parent2;
		$parent1_hash{$snp} = $parent1;
		}
	}	
}
						
######################################################################################################################
#  create list of mappable snps with segregation types
	foreach my $key (sort keys %map_snp)	{	
 		if ((($parent1_hash{$key} eq "AA") || ($parent1_hash{$key} eq "CC") || ($parent1_hash{$key} eq "GG") || ($parent1_hash{$key} eq "TT")) && (($parent2_hash{$key} ne 				"AA") && ($parent2_hash{$key} ne "CC") && ($parent2_hash{$key} ne "GG") && ($parent2_hash{$key} ne "TT")))	{
			$segregation_type{$key} = "\<nnxnp\>";  #parent1 x parent2 segregation type
#			$segregation_type{$key} = "\<lmxll\>";	#parent2 x parent1 segregation type
		}elsif ((($parent1_hash{$key} ne "AA") && ($parent1_hash{$key} ne "CC") && ($parent1_hash{$key} ne "GG") && ($parent1_hash{$key} ne "TT")) && (($parent2_hash{$key} eq 				"AA") || ($parent2_hash{$key} eq "CC") || ($parent2_hash{$key} eq "GG") || ($parent2_hash{$key} eq "TT")))	{
			$segregation_type{$key} = "\<lmxll\>";	#parent1 x parent2 segregation type
#			$segregation_type{$key} = "\<nnxnp\>";	#parent2 x parent1 segregation type
		}elsif ((($parent1_hash{$key} ne "AA") && ($parent1_hash{$key} ne "CC") && ($parent1_hash{$key} ne "GG") && ($parent1_hash{$key} ne "TT")) && (($parent2_hash{$key} ne 				"AA") && ($parent2_hash{$key} ne "CC") && ($parent2_hash{$key} ne "GG") && ($parent2_hash{$key} ne "TT")))	{
			$segregation_type{$key} = "\<hkxhk\>";
		}
	}

#######################################################################################################################
#  print out file with mappable snps, segregation types and parent1 x parent2 genotypes parent1 first
open (OUT, ">>$parent1_name\_x_$parent2_name\_mappable_snps_segregation_types.txt") || die "Can't print out parent1 x parent2 mappable snps and segregation types $!\n";
#open (OUT, ">>parent2xparent1_mappable_snps_segregation_types.txt") || die "Can't print out parent2 x parent1 mappable snps and segregation types $!\n";
		print OUT "SNP name\t Segregation Type\t $parent1_name genotype\t $parent2_name genotype\n";
	foreach my $key (sort keys %map_snp)	{
		print OUT "$key\t$segregation_type{$key}\t$parent1_hash{$key}\t$parent2_hash{$key}\n";
		}


