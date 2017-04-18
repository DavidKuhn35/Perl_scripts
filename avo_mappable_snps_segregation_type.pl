#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: October 4, 2012
# Version		: ver. 1.1	December 11,2012
# History		: ver. 1.0
# Description	: 	Take in list of parents genotypes for a mapping population with unmappable SNPs removed and designate segregation type for both reciprocal crosses.
# <input>		:  List of parents genotypes with missing data removed	
# <output>		:   List of mappable SNPs and segregation type for each SNP locus
# 			
# Dependencies	: N/A
# Usage			: avo_mappable_snps_segregation_type.pl, parent1 name, parent2 name
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
my $parent1 = 0;	#as defined by first column in parental genotypes file
my $p1name = ();
my $parent2 = 0;	#as defined by second column in parental genotypes file
my $p2name = ();
my %map_snp = ();
my %parent2_hash = ();
my %parent1_hash = ();
my %segregation_type_p1xp2 = ();
my %segregation_type_p2xp1 = ();
#####################################################################################################################
#	Read in list of parents genotypes with missing data removed.  simmonds_tonnage_consensus_genotypes.txt
open (IN,"$ARGV[0]") || die "Can't open parents genotypes file $!\n";
	while (<IN>)		{
		chomp;
		if (m/^SHRSPa/)	{
			($snp, $parent1, $parent2) = split /\t/;
			unless ((($parent1 eq "AA") || ($parent1 eq "CC") || ($parent1 eq "GG") || ($parent1 eq "TT")) && (($parent2 eq "AA")|| ($parent2 eq "CC") || ($parent2 eq 					"GG") || ($parent2 eq "TT")))	{
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
			$segregation_type_p1xp2{$key} = "\<nnxnp\>";  #parent1 x parent2 segregation type
			$segregation_type_p2xp1{$key} = "\<lmxll\>";	#parent2 x parent1 segregation type
		}elsif ((($parent1_hash{$key} ne "AA") && ($parent1_hash{$key} ne "CC") && ($parent1_hash{$key} ne "GG") && ($parent1_hash{$key} ne "TT")) && (($parent2_hash{$key} eq 				"AA") || ($parent2_hash{$key} eq "CC") || ($parent2_hash{$key} eq "GG") || ($parent2_hash{$key} eq "TT")))	{
			$segregation_type_p1xp2{$key} = "\<lmxll\>";	#parent1 x parent2 segregation type
			$segregation_type_p2xp1{$key} = "\<nnxnp\>";	#parent2 x parent1 segregation type
		}elsif ((($parent1_hash{$key} ne "AA") && ($parent1_hash{$key} ne "CC") && ($parent1_hash{$key} ne "GG") && ($parent1_hash{$key} ne "TT")) && (($parent2_hash{$key} ne 				"AA") && ($parent2_hash{$key} ne "CC") && ($parent2_hash{$key} ne "GG") && ($parent2_hash{$key} ne "TT")))	{
			$segregation_type_p1xp2{$key} = "\<hkxhk\>";	#parent1 x parent2 segregation type
			$segregation_type_p2xp1{$key} = "\<hkxhk\>";	#parent2 x parent1 segregation type
		}
	}

#######################################################################################################################
#  print out file with mappable snps, segregation types and parents genotypes 
	$p1name = $ARGV[1];	#reaD parent1 name based on data in first column
	$p2name = $ARGV[2];	#read in parent2 name based on data in second column

open (P1XP2, ">>$p1name\_x_$p2name\_mappable_snps_segregation_types.txt") || die "Can't print out $p1name x $p2name mappable snps and segregation types $!\n";
		print P1XP2 "SNP name\t Segregation Type\t $p1name genotype\t $p2name genotype\n";
	foreach my $key (sort keys %map_snp)	{
		print P1XP2 "$key\t$segregation_type_p1xp2{$key}\t$parent1_hash{$key}\t$parent2_hash{$key}\n";
		}
close P1XP2;

open (P2XP1, ">>$p2name\_x_$p1name\_mappable_snps_segregation_types.txt") || die "Can't print out $p2name x $p1name mappable snps and segregation types $!\n";
		print P2XP1 "SNP name\t Segregation Type\t $p2name genotype\t $p1name genotype\n";
	foreach my $key (sort keys %map_snp)	{
		print P2XP1 "$key\t$segregation_type_p2xp1{$key}\t$parent2_hash{$key}\t$parent1_hash{$key}\n";
		}
close P2XP1;
#########################################################################################################################
#

##########################################################################################################################
#

