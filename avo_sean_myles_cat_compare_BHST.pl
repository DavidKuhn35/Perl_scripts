#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: July 28, 2014
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	:  read in sean myles cat files for BHST, get SNPs common to two parents
# <input>		: MON_bacon_cat.log, etc.
# <output>		: list of SNPs by cat common to at least two parents
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
#################################################################################################################################
# declare global variables
my %hoh = ();
my %allsnp = ();
my %snpcat = ();
my %bacon = ();
my $bacon = "bacon";
my $tonnage = "tonnage";
my %tonnage = ();
my $simmonds = "simmonds";
my %simmonds = ();
my $hass = "hass";
my %hass = ();
my %bad307 = ();
my %badcat = ();
#####################################################################################################################
#  read in snp names for 307 bad chip snps
open (BAD, "bad_307_snps.txt") || die "can't open bad 307 snps $!\n";
	while (<BAD>)	{
	chomp;
	my ($name, $snp) = split /\t/;
	$bad307{$name} = $snp;
}
#####################################################################################################################
#	main	
open (CATB, "MON_bacon_cat.log")|| die "can't open bacon cat log $!\n";
	while (<CATB>)	{
	chomp;
	my ($bsnp, $bcat) = split /\t/;
	$allsnp{$bsnp} = $bsnp;			
	$bacon{$bsnp} = $bcat; 
	print "$bsnp\t$bacon{$bsnp}\n";
}			
######################################################################################################################
#   
open (CATH, "MON_hass_cat.log") || die "can't open hass cat log $!\n";
	while (<CATH>)	{
	chomp;
	my @h = split /\t/;
	my $hsnp = $h[0];			
	my $hcat = $h[1];
	$hass{$hsnp} = $hcat;
}
#######################################################################################################################
#   
open (CATS, "MON_simmonds_cat.log" )|| die "can't open simmonds cat log $!\n";
	while (<CATS>)	{
	chomp;
	my @s = split /\t/;
	my $ssnp = $s[0];			
	my $scat = $s[1];
	$simmonds{$ssnp} = $scat;
}
#########################################################################################################################
#
open (CATT, "MON_tonnage_cat.log") || die "can't open bacon cat log $!\n";
	while (<CATT>)	{
	chomp;
	my @t = split /\t/;
	my $tsnp = $t[0];			
	my $tcat = $t[1];
	$tonnage{$tsnp} = $tcat;
}
##########################################################################################################################
#  find cats of bad snps for each parent
open (BADOUT, ">>badcat_BHST.txt") || die "can't open badcat out $!\n";
foreach my $key (sort keys %bad307)	{
	;
	print BADOUT "$key\t$bacon{$key}\t$hass{$key}\t$simmonds{$key}\t$tonnage{$key}\n";
}


	

