#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: November 6, 2012
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: 	Read in SNP chip data.  Read in cacao parents.  Make consensus genotypes for cacao parents.  Determine segregation patterns for each SNP locus
# <input>		:  SNP chip data, PNG parents WG#		
# <output>		:   consensus genotypes for PNG parents, segregation patterns for each locus.
# 			
# Dependencies	: N/A
# Usage			: PNG_parents_consensus_genotype_Nov6.pl
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
my %hoh = ();
my %snp_hash = ();
my $snp = 0;
my $dna = 0;
my $geno1 = 0;
my $geno2 = 0;
my %WG = ();
my %parent1 = ();	# PNG_K82
my $parent1_header = "SNPname\t";
my %parent2 = ();	#	PNG_KA2-101
my $parent2_header = "SNPname\t";
my %parent1_consensus = ();
my %parent2_consensus = ();

#####################################################################################################################
#	main	Read in SNP chip data.
open (SNPDATA,$ARGV[0]) || die "Can't open SNP chip data $!\n";
	while (<SNPDATA>)	{
		chomp;
		if (m/^contig/)	{	#if line has SNP name in first column
			($snp, $dna, $geno1, $geno2) = split /\t/;	#split on tab
			my @genotype = ($geno1, $geno2);	#put genotypes in an array
			my @genosort = sort (@genotype);	#sort genotypes ASCIIbetically A,C,G,T
			$hoh{$snp}{$dna} = "$genosort[0]$genosort[1]";	#hoh with SNP key and dna key and genotype value
			$snp_hash{$snp} = $snp;	
		}
}	
				
######################################################################################################################
# Read in manifest data and make hashes of parents  
open (DNADATA, $ARGV[1]) || die "Can't open manifest data $!\n";
	while (<DNADATA>)	{
		chomp;
		if (m/^WG/)	{	#if line starts with Well name
		my @dnadata = split /\t/;
		$WG{$dnadata[0]} = $dnadata[1];	#make hash with well numbers as keys and common names as values
		}
	}	
#######################################################################################################################
#  Get all parents and make hashes of WG names
	foreach my $key (sort keys %WG)		{
		if ($WG{$key}=~ /PNG_K82/)	{	#PNG_K82 is parent 1
			$parent1{$key} = $key;
		} elsif ($WG{$key}=~ /PNG_KA2-101/)	{	#PNG_KA2-101 is parent 2
			$parent2{$key} = $key;
		}
	}

		
#########################################################################################################################
#	Go through SNPs for each parent and make a file with all genotype data for each parent
	foreach my $key (sort keys %snp_hash)	{	#go through SNPs one at a time
		foreach my $x (sort keys %parent1)	{  # go through all  parent1 K82
			$parent1_consensus{$key} .= "$hoh{$key}{$x}\t";
			}
		foreach my $x (sort keys %parent2)	{ # go through all parent 2 KA2-101
			$parent2_consensus{$key} .= "$hoh{$key}{$x}\t";
			}
		}

##########################################################################################################################
#	Print outfile with all genotype data for each parent
open (PARENT1_OUT, ">>K82_consensus_parent_genotypes.txt") || die "Can't open parent1 consensus parent genotypes outfile $!\n";
		foreach my $x (sort keys %parent1)	{
			$parent1_header .= "$x\t";
		}
		print PARENT1_OUT "$parent1_header\n";
			foreach my $x (sort keys %parent1_consensus)	{
				print PARENT1_OUT "$x\t$parent1_consensus{$x}\n";
			}
close PARENT1_OUT;				
open (PARENT2_OUT, ">>KA2-101_consensus_parent_genotypes.txt") || die "Can't open parent 2 consensus parent genotypes outfile $!\n";
		foreach my $x (sort keys %parent2)	{
			$parent2_header .= "$x\t";
		}
		print PARENT2_OUT "$parent2_header\n";
			foreach my $x (sort keys %parent2_consensus)	{
				print PARENT2_OUT "$x\t$parent2_consensus{$x}\n";
			}
close PARENT2_OUT;		 

