#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: October 12, 2012
# Version		: ver. 1.0
# History		: ver. 1.0
# Description		: calculate segregation pattern (lm x ll, nnxnp, hkxhk) from parental genotype data 
# <input>		:  SNPchip data and parent names		
# <output>		:   file with segregation pattern for each SNP for each type of population
# 			
# Dependencies	: N/A
# Usage			: avo_segregation_pattern.pl
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
my $snp = 0;
my $dna = 0;
my $geno1 = 0;
my $geno2 = 0;
my $snpname = 0;
my $SNPseq = 0;
my %snptype = ();
my $mother = 0;
my $father = 0;
my %matgen = ();
my %patgen = ();
my %snp_pattern = ();
#####################################################################################################################
#	main  read in SNPchip datafile and make a hoh with SNP => DNA => genotype
	
	open (IN, "$ARGV[0]") || die "Can't open SNPchip datafile $!\n";
		while (<IN>)	{
			chomp;
			if (m/^SHRSPa/)	{	#if line has SNP name in first column
			($snp, $dna, $geno1, $geno2) = split /,/;	#split on comma
			my @genotype = ($geno1, $geno2);	#put genotypes in an array
			my @genosort = sort (@genotype);	#sort genotypes ASCIIbetically A,C,G,T
			$hoh{$snp}{$dna} = "$genosort[0]\t$genosort[1]";	#hoh with SNP key and dna key and genotype value			
		}
	}
	
########################################################################################################################
#	read in SNP sequence data to get SNP type and make hash of SNP type
	
	open (SNP, "$ARGV[1]") || die "Can't open file with SNP sequence data $!\n";
		while (<SNP>)	{
			chomp;
			if (m/^SHRSPa/)	{	#if line has SNP name in first column
			($snpname, $SNPseq) = split /\t/;
			if (m!\[(.)/(.)\])	{
				my @SNPgeno = ($1, $2);
				my @SNPgenosort = sort @SNPgeno;
				$snptype{$snpname} = "$SNPgenosort[0]\t$SNPgenosort[1]";
			}
		}
	}

#######################################################################################################################
#  Compare genotypes of parents at each SNP and assign segregation type to SNP (lmxll, nnxnp, hkxhk)
	$mother = $ARGV[2];
	$father = $ARGV[3];
	foreach my $key (sort keys %hoh)	{
		if ($mother = $hoh{$key})	{
			$matgen{$key}= $hoh{$key}{$mother}; 
		} elsif ($father = $hoh{$key})	{  
			$patgen{$key} = $hoh{$key}{$father};
		}
	}
	
#########################################################################################################################
#	For each SNP compare SNP type to maternal and paternal genotypes to assign segregation type.  Identify which nucleotide is l,m,n,p,h,k
	foreach $x (sort keys %hoh)	{
		my @snp_pattern_nuc = split /\t/, $snptype{$x};
			if (($snptype{$x} eq $matgen{$x}) && ($snptype{$x} eq $patgen{$x}))	{
				$snp_pattern{$x}{type}= "\<hkxhk\>";
				$snptype{$x}{h} = $snp_pattern_nuc[0];
				$snptype{$x}{k} = $snp_pattern_nuc[1];
			} elsif (($snptype{$x} eq $matgen{$x}) && ($snptype{$x} ne $patgen{$x}))	{
				$snp_pattern{$x}{type} = "\<lmxll\>";
					if ($patgen{$x}=~ /$snp_pattern_nuc[0]/)	{
						$snptype{$x}{l} = $snp_pattern_nuc[0];
						$snptype{$x}{m} = $snp_pattern_nuc[1];
					} else {
						$snptype{$x}{l} = $snp_pattern_nuc[1];
						$snptype{$x}{m} = $snp_pattern_nuc[0];
					}
			}elsif (($snptype{$x} ne $matgen{$x}) && ($snptype{$x} eq $patgen{$x}))	{
				$snp_pattern{$x}{type} = "\<nnxnp\>";
					if ($matgen{$x}=~ /$snp_pattern_nuc[0]/)	{
						$snptype{$x}{n} = $snp_pattern_nuc[0];
						$snptype{$x}{p} = $snp_pattern_nuc[1];
					} else {
						$snptype{$x}{n} = $snp_pattern_nuc[1];
						$snptype{$x}{p} = $snp_pattern_nuc[0];
					}
	}
##########################################################################################################################
#	print out snp types for all snps and maternal and paternal genotypes in JoinMap format
	




