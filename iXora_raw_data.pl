#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: February 1, 2013
# Version		: ver. 1.0 based on joinmap_final_format_v2.pl
# History		: ver. 1.0
# Description	:  creates raw data file with SNPs in rows, DNA in columns for a population
# <input>		:  SNPchip data, DNA manifest		
# <output>		:   raw data file in iXora format
# 			
# Dependencies	: N/A
# Usage			: joinmap_rawdata.pl
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

#####################################################################################################################
#
my %hoh = ();
my %snp_hash = ();
my $snp = 0;
my $dna = 0;
my $geno1 = 0;
my $geno2 = 0;
my %WG = ();
my %pop1 = ();
my $popname = ();
my %mapsnp_hash = ();
my $mapsnp = 0;
my $seg_type = 0;
my $parent1 = 0;
my $parent2 = 0;
my %p1mapsnp = ();
my %p2mapsnp = ();
my %pop1_raw = ();
my $pop1_raw_header = "SNPname\t";

#####################################################################################################################
#	main	Read in SNP chip data.
open (SNPDATA,$ARGV[0]) || die "Can't open avo SNP chip data $!\n";
	while (<SNPDATA>)	{
		chomp;
		if (m/^SHRSPa/)	{	#if line has SNP name in first column, avocado SHRSPa  cacao contig
			($snp, $dna, $geno1, $geno2) = split /,/;	#split on tab for cacao, comma for avocado
			my @genotype = ($geno1, $geno2);	#put genotypes in an array
			my @genosort = sort (@genotype);	#sort genotypes ASCIIbetically A,C,G,T
			$hoh{$snp}{$dna} = "$genosort[0]$genosort[1]";	#hoh with SNP key and dna key and genotype value
			$snp_hash{$snp} = $snp;	
		}
}	
				
######################################################################################################################
# Read in manifest data and make hashes of parents  
open (DNADATA, $ARGV[1]) || die "Can't open avo manifest data $!\n";
	while (<DNADATA>)	{
		chomp;
		if (m/^\d/)	{	#if line has DNA name data  (WG cacao)  (\d avocado)
		my @dnadata = split /\t/;
		$WG{$dnadata[3]} = $dnadata[6];	#make hash with well numbers as keys and common names as values ([0] and [1] for cacao) ([3] [6] for avocado)
		}
	}	
#######################################################################################################################
#  Get all members of a population and make hashes of WG names
	$popname = $ARGV[2];
	foreach my $key (sort keys %WG)		{
		if ($WG{$key}=~ /$popname/)	{	#put in selector for population name (e.g. PNG_8)
			$pop1{$key} = $key;
				}
	}
######################################################################################################################

#########################################################################################################################
#	Go through SNPs for each genotype in a mapping population and make a file with genotype data.  snps are rows, genotypes are columns
	foreach my $key (sort keys %snp_hash)	{	#go through mappable SNPs one at a time
		foreach my $x (sort keys %pop1)	{  # go through all pop1 genotypes
			$pop1_raw{$key} .= "$hoh{$key}{$x}\t"; #creates a raw data file for pop1  $key is snp name, $x is DNA name
			}
		}

##########################################################################################################################
#	Print outfile with all genotype data for each genotype in FLAMP_SxT
	my $date = `date`;
	$date =~ s/\s//g;
open (RAW, ">>iXora_$popname\_raw\_$date\_.txt") || die "Can't open pop1_raw outfile $!\n";
		foreach my $x (sort keys %pop1)	{
			$pop1_raw_header .= "$x\t";
		}
		print RAW "$pop1_raw_header\n";
			foreach my $x (sort keys %pop1_raw)	{
			print RAW  "$x\t$pop1_raw{$x}\n";
			}
close RAW;
				
######################################################################################################################
#   

#######################################################################################################################
#   
 
#########################################################################################################################
#

##########################################################################################################################
#
 
