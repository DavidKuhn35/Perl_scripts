#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: May 29, 2013
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: 	Read in SNP chip data.  Read in avo parents.  Make consensus genotypes for avo parents.  Determine segregation patterns for each SNP locus
# <input>		:  SNP chip data, avo parents WG#		
# <output>		:   consensus genotypes for parents, segregation patterns for each locus.
# 			
# Dependencies	: N/A
# Usage			: avo_parents_consensus_genotype_Oct3_2012.pl
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
		if ($WG{$key}=~ /Control_Simmonds/)	{
			$simmonds_parent{$key} = $key;
		} elsif ($WG{$key}=~ /Control_Tonnage/)	{
			$tonnage_parent{$key} = $key;
		} elsif ($WG{$key} =~ /Control_Bacon/)	{
			$bacon_parent{$key} = $key;
		} elsif ($WG{$key} =~ /Control_Hass/)	{
			$hass_parent{$key} = $key;
		}
	}

		
#########################################################################################################################
#	Go through SNPs for each parent and make a file with all genotype data for each parent
	foreach my $key (sort keys %snp_hash)	{	#go through SNPs one at a time
		foreach my $x (sort keys %simmonds_parent)	{  # go through all Simmonds parents
			$simmonds_parent_consensus{$key} .= "$hoh{$key}{$x}\t";
			}
		foreach my $x (sort keys %tonnage_parent)	{ # go through all Tonnage parents
			$tonnage_parent_consensus{$key} .= "$hoh{$key}{$x}\t";
			}
		foreach my $x (sort keys %bacon_parent)	{
			$bacon_parent_consensus{$key} .= "$hoh{$key}{$x}\t";
			}
		foreach my $x (sort keys %hass_parent)	{
			$hass_parent_consensus{$key} .= "$hoh{$key}{$x}\t";
			}
		}

##########################################################################################################################
#	Print outfile with all genotype data for each parent
open (BACON_OUT, ">>bacon_consensus_parent_genotypes.txt") || die "Can't open bacon consensus parent genotypes outfile $!\n";
		foreach my $x (sort keys %bacon_parent)	{
			$bacon_parent_header .= "$x\t";
		}
		print BACON_OUT "$bacon_parent_header\n";
			foreach my $x (sort keys %bacon_parent_consensus)	{
			print BACON_OUT "$x\t$bacon_parent_consensus{$x}\n";
			}
close BACON_OUT;
open (HASS_OUT, ">>hass_consensus_parent_genotypes.txt") || die "Can't open hass consensus parent genotypes outfile $!\n";
		foreach my $x (sort keys %hass_parent)	{
			$hass_parent_header .= "$x\t";
		}
		print HASS_OUT "$hass_parent_header\n";
			foreach my $x (sort keys %hass_parent_consensus)	{
				print HASS_OUT "$x\t$hass_parent_consensus{$x}\n";
			}
close HASS_OUT;
open (SIMMONDS_OUT, ">>simmonds_consensus_parent_genotypes.txt") || die "Can't open simmonds consensus parent genotypes outfile $!\n";
		foreach my $x (sort keys %simmonds_parent)	{
			$simmonds_parent_header .= "$x\t";
		}
		print SIMMONDS_OUT "$simmonds_parent_header\n";
			foreach my $x (sort keys %simmonds_parent_consensus)	{
				print SIMMONDS_OUT "$x\t$simmonds_parent_consensus{$x}\n";
			}
close SIMMONDS_OUT;				
open (TONNAGE_OUT, ">>tonnage_consensus_parent_genotypes.txt") || die "Can't open tonnage consensus parent genotypes outfile $!\n";
		foreach my $x (sort keys %tonnage_parent)	{
			$tonnage_parent_header .= "$x\t";
		}
		print TONNAGE_OUT "$tonnage_parent_header\n";
			foreach my $x (sort keys %tonnage_parent_consensus)	{
				print TONNAGE_OUT "$x\t$tonnage_parent_consensus{$x}\n";
			}
close TONNAGE_OUT;		 

