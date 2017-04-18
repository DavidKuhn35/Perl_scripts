#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: May 29, 2013
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: 	Read in SNP data.  Read in manifest data and parents' names.  Make genotypes for parents.  
# <input>		:  SNP data, parents names	
# <output>		:   consensus genotypes for parents, segregation patterns for each locus.
# 			
# Dependencies	: N/A
# Usage			: EP1_identify_parent_genotypes.pl
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
my %hoh2 = ();
my %snp_hash = ();
my %dna_hash = ();
my %Irwin = ();
my %KP = ();
my %map_snp = ();
my %IxKP = ();
my $snp = 0;
my $dna = 0;
my $geno1 = 0;
my $geno2 = 0;
my $print_header = "SNP name\tSegregation type\tIrwin genotype\tKP genotype\n";
my %segregation_type = ();
my $print_header2 = "SNP name\t";
my %snp_geno = ();

#####################################################################################################################
#	main	Read in SNP  data from EP-1 csv file.
open (SNPDATA,$ARGV[0]) || die "Can't open avo SNP chip data $!\n";
	while (<SNPDATA>)	{
		chomp;
		if (m/^S/)	{	#if line has unique SNP id in first column
			my @ep1 = split /,/;	#split on comma
			if (($ep1[9] ne "No Call")&&($ep1[9] ne "Invalid") && ($ep1[9] ne "NTC"))	{
			my @allele = split /:/, $ep1[9];	#get nucleotide genotypes and put in an array
			my @allesort = sort (@allele);	#sort genotypes ASCIIbetically A,C,G,T
			$hoh{$ep1[1]}{$ep1[4]} = "$allesort[0]$allesort[1]";	#hoh with SNP as first key and cultivar name key and genotype value
			$hoh2{$ep1[4]}{$ep1[1]} = "$allesort[0]$allesort[1]";	#hoh with cultivar name as first key, snp as key and genotype value
			$snp_hash{$ep1[1]} = $ep1[1];	#hash with snp name key and snp name value
			$dna_hash{$ep1[4]} = $ep1[4];	#hash with dna name key and dna name value
			}
		}
}	
				
close SNPDATA;
open (DATACHK, ">>datachk_$ARGV[1].txt") || die "Can't open data check $!\n";
		foreach my $key (sort keys %hoh2)	{
			foreach my $x (sort keys %{$hoh2{$key}})	{
			print DATACHK "$key\t$hoh2{$key}\t$hoh2{$key}{$x}\n";
		}
}
close DATACHK;	
#########################################################################################################################
#	Go through SNPs for each parent and make a file with all genotype data for each parent
	foreach my $key (sort keys %hoh2)	{	#go through dnas one at a time
		if ($key eq "Irwin")	{	# if cultivar is Irwin, make hash with all genotypes
			foreach my $x (sort keys %{$hoh2{$key}})	{
				$Irwin{$x} = $hoh2{$key}{$x};
			}
		}
	}
open (IRWINCHK, ">>irwinchk_May30_2013.txt") || die "Can't open data check $!\n";
		foreach my $key (sort keys %Irwin)	{
			print IRWINCHK "$key\t$Irwin{$key}\n";
		}
close IRWINCHK;	
##########################################################################################################################
#	Go through SNPs for KP and make a hash with all genotypes
	foreach my $key (sort keys %hoh2)	{
		if ($key eq "KP")	{	# if cultivar is KP, make hash with all genotypes
			foreach my $x (sort keys %{$hoh2{$key}})	{
				$KP{$x} = $hoh2{$key}{$x};
		}
	}
}

#####################################################################################################################
#	Go through all SNPs for both parents (Irwin and KP) and get mappable snps
		foreach my $key (sort keys %snp_hash)	{
			if ((exists $Irwin{$key}) && (exists $KP{$key}))	{
			unless ((($Irwin{$key} eq "AA") || ($Irwin{$key} eq "CC") || ($Irwin{$key} eq "GG") || ($Irwin{$key} eq "TT") || ($Irwin{$key} eq "--")) && (($KP{$key} eq "AA")|| ($KP{$key} eq "CC") || ($KP{$key} eq "GG") || ($KP{$key} eq "TT") || ($KP{$key} eq "--")))	{
		$map_snp{$key} = $key;
				}
	}	
}
						
######################################################################################################################
#  create list of mappable snps with segregation types
	foreach my $key (sort keys %map_snp)	{	
 		if ((($Irwin{$key} eq "AA") || ($Irwin{$key} eq "CC") || ($Irwin{$key} eq "GG") || ($Irwin{$key} eq "TT")) && (($KP{$key} ne "AA") && ($KP{$key} ne "CC") && ($KP{$key} ne "GG") && ($KP{$key} ne "TT")))	{
			$segregation_type{$key} = "\<nnxnp\>";  #parent1 x parent2 segregation type
#			$segregation_type{$key} = "\<lmxll\>";	#parent2 x parent1 segregation type
		}elsif ((($Irwin{$key} ne "AA") && ($Irwin{$key} ne "CC") && ($Irwin{$key} ne "GG") && ($Irwin{$key} ne "TT")) && (($KP{$key} eq "AA") || ($KP{$key} eq "CC") || ($KP{$key} eq "GG") || ($KP{$key} eq "TT")))	{
			$segregation_type{$key} = "\<lmxll\>";	#parent1 x parent2 segregation type
#			$segregation_type{$key} = "\<nnxnp\>";	#parent2 x parent1 segregation type
		}elsif ((($Irwin{$key} ne "AA") && ($Irwin{$key} ne "CC") && ($Irwin{$key} ne "GG") && ($Irwin{$key} ne "TT")) && (($KP{$key} ne "AA") && ($KP{$key} ne "CC") && ($KP{$key} ne "GG") && ($KP{$key} ne "TT")))	{
			$segregation_type{$key} = "\<hkxhk\>";
		}
	}


##########################################################################################################################
#	Print outfile with all genotype data for each parent and snp name and segregation type
open (SNP_OUT, ">>irwin_KP_parent_genotypes_May29_2013.txt") || die "Can't open irwin parent genotypes outfile $!\n";
		print SNP_OUT "$print_header";
		foreach my $x (sort keys %map_snp)	{
			print SNP_OUT "$x\t$segregation_type{$x}\t$Irwin{$x}\t$KP{$x}\n";
			}
close SNP_OUT;

###############################################################################################################################
#  Print outfile with all mappable genotype data for all cultivars on EP-1 chip
	foreach my $key (sort keys %hoh2)	{
		$print_header2 .= "$key\t";	#put name of all cultivars in header
		foreach my $x (sort keys %map_snp)		{
				$snp_geno{$x} .= "$hoh2{$key}{$x}\t";
			}
		}

open (GENO_OUT, ">>all_geno_out_May30_2013.txt") || die "Can't open all geno out $!\n";
	print GENO_OUT "$print_header2";
	foreach my $key (sort keys %snp_geno)		{
		print GENO_OUT	"$key\t$snp_geno{$key}\n";
	}
close GENO_OUT;
