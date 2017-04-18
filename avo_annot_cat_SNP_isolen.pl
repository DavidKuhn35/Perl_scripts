#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: January 29, 2013
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: Calculates average isotig length and average SNP number for different categories of goslim annotation
# <input>		: May25_avo_all_isolen_snpcount_snpfreq_agi_goslim_out and files with the unique goslim terms for C, F, and P (or APVO isotig and annotation file)		
# <output>		: file with average isotig length and average SNP number for different categories of goslim annotation  
# 			
# Dependencies	: N/A
# Usage			: avo_annot_cat_SNP_isolen.pl
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
my %isolen = ();
my %isosnp = ();
my %isoC = ();
my %isoF = ();
my %isoP = ();
my %Cterm = ();
my %Fterm = ();
my %Pterm = ();
my %Choh = ();
my %Fhoh = ();
my %Phoh = ();


#####################################################################################################################
#	main	
open (DATA, "$ARGV[0]") || die "Can't open input data file $! \n";
	while (<DATA>)
{
	chomp;
	if (m/^g/)	{	#if line starts with an isotig name
	my @x = split (/\t/);
	$isolen{$x[0]} = $x[1];	# hash with isotig name as key and isotig length as value
	$isosnp{$x[0]} = $x[2]; # hash with isotig name as key and number of SNPs in isotig as value
	$isoC{$x[0]} = $x[7]; # hash with isotig name as key and C annotation as value
	$isoF{$x[0]} = $x[13]; #hash with isotig name as key and F annotation as value
	$isoP{$x[0]} = $x[19]; #hash with isotig name as key and P annotation as value
	
}
}				
close DATA;				
######################################################################################################################
# read in files of C, F and P annotation terms  
open (CTERMS, "$ARGV[1]") || die "Can't open C terms file $! \n";
	while (<CTERMS>)
{
	chomp;
	if (m/^.\w/)	{
	my $x = "$_";
	$Cterm{$x} = $x;	# make a hash of unique Cterm annota
}
}
	my $Ckeys = (keys %Cterm);
	print "$Ckeys\n ";
close CTERMS;
open (FTERMS, "$ARGV[2]") || die "Can't open F terms file $! \n";
 	while (<FTERMS>)
{
	chomp;
	if (m/^.\w/)	{
	$Fterm{"$_"} = "$_"; #hash with unique  F term annotation 
}
}
	my $Fkeys = (keys %Fterm);
	print " $Fkeys\n";
close FTERMS;
open (PTERMS, "$ARGV[3]") || die "Can't open P terms file $! \n";
 	while (<PTERMS>)
{
	chomp;
	if (m/^.\w/)	{
	$Pterm{"$_"} = "$_"; #hash with unique  P terms 
}
}
	my $Pkeys = (keys %Pterm);	
	print " $Pkeys \n";
close PTERMS;
#######################################################################################################################
# Count number of isotigs in each category, calculate average isotig length and average number of SNPs
	foreach my $key (sort keys %isoC)		{
		foreach my $x (sort keys %Cterm)	{
			if ($isoC{$key} eq $Cterm{$x})	{
				$Choh{$x}{$key} = $key;	# if Cterms agree store isotig name in Choh with annotation term as first key and isotig name as second key 
			}
		}
	}
	foreach my $key (sort keys %isoF)		{
		foreach my $x (sort keys %Fterm)	{
			if ($isoF{$key} eq $Fterm{$x})	{
				$Fhoh{$x}{$key} = $key;	# if Fterms agree store isotig name in Fhoh with annotation term as first key and isotig name as second key 
			}
		}
	}

	foreach my $key (sort keys %isoP)		{
		foreach my $x (sort keys %Pterm)	{
			if ($isoP{$key} eq $Pterm{$x})	{
				$Phoh{$x}{$key} = $key;	# if Pterms agree store isotig name in Phoh with annotation term as first key and isotig name as second key 
			}
		}
	}
#########################################################################################################################
#Get average isotig length and average number of SNPs for each annotation term

	foreach my $key (sort keys %Choh)	{	#keys are annotation terms
		my $isotig_count = keys %{$Choh{$key}};	#get number of isotigs for annotation term
		my $iso_len_sum = 0;
		my $SNP_sum = 0;
		foreach my $x (sort keys %{$Choh{$key}})		{	#keys (x) are isotig names
			$iso_len_sum += $isolen{$x};	#get sum of isotig lengths for annotation term
			$SNP_sum += $isosnp{$x};	#get sum of SNPs for annotation term
		}
		$Choh{$key}{isotig_count} = $isotig_count;
		$Choh{$key}{ave_isotig} = $iso_len_sum/$isotig_count;	# get ave isotig length and put in Choh
		$Choh{$key}{ave_SNP} = $SNP_sum/$isotig_count;	#get ave isotig length and put in Choh
	} 	

	foreach my $key (sort keys %Fhoh)	{	#keys are annotation terms
		my $isotig_count = keys %{$Fhoh{$key}};	#get number of isotigs for annotation term
		my $iso_len_sum = 0;
		my $SNP_sum = 0;
		foreach my $x (sort keys %{$Fhoh{$key}})		{	#keys are isotig names
			$iso_len_sum += $isolen{$x};	#get sum of isotig lengths for annotation term
			$SNP_sum += $isosnp{$x};	#get sum of SNPs for annotation term
			}
		$Fhoh{$key}{isotig_count} = $isotig_count;
		$Fhoh{$key}{ave_isotig} = $iso_len_sum/$isotig_count;	# get ave isotig length and put in Fhoh
		$Fhoh{$key}{ave_SNP} = $SNP_sum/$isotig_count;	#get ave SNP per isotig and put in Fhoh
	} 	
	
	foreach my $key (sort keys %Phoh)	{	#keys are annotation terms
		my $isotig_count = keys %{$Phoh{$key}};	#get number of isotigs for annotation term
		my $iso_len_sum = 0;
		my $SNP_sum = 0;
		foreach my $x (sort keys %{$Phoh{$key}})		{	#keys are isotig names
			$iso_len_sum += $isolen{$x};	#get sum of isotig lengths for annotation term
			$SNP_sum += $isosnp{$x};	#get sum of SNPs for annotation term
			}
		$Phoh{$key}{isotig_count} = $isotig_count;
		$Phoh{$key}{ave_isotig} = $iso_len_sum/$isotig_count;	# get ave isotig length and put in Phoh
		$Phoh{$key}{ave_SNP} = $SNP_sum/$isotig_count;	#get ave isotig length and put in Phoh
	} 			
##########################################################################################################################
#Write out report with average isotig length and average number of SNPs for each annotation term, use printf to format averages

	open (OUTFILE, ">>avo_$ARGV[4]\_annot_cat_SNP_isolen_out.txt") || die "Can't open outfile $! \n";
		my $C = "C";
		my $F = "F";
		my $P = "P";
		printf OUTFILE "Type of annotation\tAnnotation term\tIsotig count\tAverage isotig length\tAverage SNP number\n";
		foreach my $key (sort keys %Choh)	{ 
			printf OUTFILE "$C\t$key\t$Choh{$key}{isotig_count}\t%6.3f\t%6.3f\n", $Choh{$key}{ave_isotig},$Choh{$key}{ave_SNP};
			}
		foreach my $key (sort keys %Fhoh)	{
			printf OUTFILE "$F\t$key\t$Fhoh{$key}{isotig_count}\t%6.3f\t%6.3f\n", $Fhoh{$key}{ave_isotig},$Fhoh{$key}{ave_SNP};
			}
		foreach my $key (sort keys %Phoh)	{
			printf OUTFILE "$P\t$key\t$Phoh{$key}{isotig_count}\t%6.3f\t%6.3f\n", $Phoh{$key}{ave_isotig},$Phoh{$key}{ave_SNP};
			}
close OUTFILE;

