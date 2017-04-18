#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: October 4, 2012
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	:  creates raw data file in JoinMap format
# <input>		:  SNPchip data, DNA manifest, mappable SNPs for population, segregation type for each SNP		
# <output>		:   raw data file in JoinMap format
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
my %tonnage_parent = ();
my $tonnage_parent_header = "SNPname\t";
my %simmonds_parent = ();
my $simmonds_parent_header = "SNPname\t";
my %hass_parent = ();
my $hass_parent_header = "SNPname\t";
my %bacon_parent = ();
my $bacon_parent_header = "SNPname\t";
my %tonnage_parent_consensus = ();
my %simmonds_parent_consensus = ();
my %hass_parent_consensus = ();
my %bacon_parent_consensus = ();
my %FLAMP_SxT = ();
my %FLAMP_TxS = ();
my %mapsnp_hash = ();
my $mapsnp = 0;
my $seg_type = 0;
my $simmonds = 0;
my $tonnage = 0;
my %simparent = ();
my %tonparent = ();
my %FLAMP_SxT_raw = ();
my $FLAMP_SxT_raw_header = "SNPname\t";
my %FLAMP_TxS_raw = ();
my $FLAMP_TxS_raw_header = "SNPname\t";
#####################################################################################################################
#	main	Read in SNP chip data.
open (SNPDATA,$ARGV[0]) || die "Can't open avo SNP chip data $!\n";
	while (<SNPDATA>)	{
		chomp;
		if (m/^SHRSPa/)	{	#if line has SNP name in first column
			($snp, $dna, $geno1, $geno2) = split /,/;	#split on comma
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
		if (m/^\d/)	{
		my @dnadata = split /\t/;
		$WG{$dnadata[3]} = $dnadata[6];	#make hash with well numbers as keys and common names as values
		}
	}	
#######################################################################################################################
#  Get all parents and make hashes of WG names
	foreach my $key (sort keys %WG)		{
		if ($WG{$key}=~ /FLAMP_SxT/)	{
			$FLAMP_SxT{$key} = $key;
		} elsif ($WG{$key}=~ /FLAMP_TxS/)	{
			$FLAMP_TxS{$key} = $key;
		}
	}
######################################################################################################################
#	Read in mappable snps, segregation patterns and Simmonds and Tonnage genotypes at snp.  Create hashes with seg_type and simmonds parent genotype and tonnage parent genotype
open (MAPSNP, $ARGV[2]) || die "Can't open mappable snps file $!\n";
	while (<MAPSNP>)		{
		if (m/^SHRSPa/)	{
#			($mapsnp, $seg_type, $simmonds, $tonnage) = split /\t/;	#simmonds x tonnage cross
			($mapsnp, $seg_type, $tonnage, $simmonds) = split /\t/;	#tonnage x simmonds cross
			$mapsnp_hash{$mapsnp}= $seg_type;
			$simparent{$mapsnp} = $simmonds;
			$tonparent{$mapsnp} = $tonnage;
		}
	}
		
#########################################################################################################################
#	Go through mappable SNPs for each genotype in a mapping population and make a file with genotype data.  snps are rows, genotypes are columns
	foreach my $key (sort keys %mapsnp_hash)	{	#go through mappable SNPs one at a time
#		foreach my $x (sort keys %FLAMP_SxT)	{  # go through all FLAMP_SxT genotypes
#			$FLAMP_SxT_raw{$key} .= "$hoh{$key}{$x}\t"; #creates a raw data file for SxT
		foreach my $x (sort keys %FLAMP_TxS)	{  # go through all FLAMP_TxS genotypes
			$FLAMP_TxS_raw{$key} .= "$hoh{$key}{$x}\t"; #creates a raw data file for TxS
			}
		}

##########################################################################################################################
#	Print outfile with all genotype data for each genotype in FLAMP_SxT
open (RAW, ">>FLAMP_TxS_raw.txt") || die "Can't open FLAMP_TxS_raw outfile $!\n";
		foreach my $x (sort keys %FLAMP_TxS)	{
			$FLAMP_TxS_raw_header .= "$x\t";
		}
		print RAW "$FLAMP_TxS_raw_header\n";
			foreach my $x (sort keys %FLAMP_TxS_raw)	{
			print RAW  "$x\t$FLAMP_TxS_raw{$x}\n";
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
 
