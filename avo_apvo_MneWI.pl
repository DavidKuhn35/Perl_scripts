#!/usr/bin/perl -w
##################################################################################
# Author		: David Kuhn
# Date			: 6/21/2012
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: 	reads in file with apvo959 isotig name and annotation.  Reads in file with isotigs that contain Mexican Bacon) homo snp for one allele and West Indian
#			 (Simmonds) homo snp for the other allele (MneWI).  Puts out a file with apvo959 isotigs that have MneWI snps along with number of MneWI snps and frequency
# <input>		:  ZZMay25_apvo959_avo_out.txt 	ARGV[0], DDD_summary_MneWI_isotigs_snps.out ARGV[1]	
# <output>		:  MMM_apvo959_mnewi_out 
# 			
# Dependencies	: N/A
# Usage			: apvo_MneWI.pl
#  Still under development (6/22/12).                    
##################################################################################
#use strict;
use warnings;
use diagnostics;
use Statistics::Descriptive;
my %hoh = ();
my %hoh2 = ();
my %hoh3 = ();
my %isotig = ();
my $goodsnp = 0;
my $badsnp = 0;
my @snpdistribution = ();
my $snpname = "snpname";
my @bacon = ();
my @hass = ();
my @simmonds = ();
my @tonnage = ();
my %bacongeno = ();
my %hassgeno = ();
my %simmondsgeno = ();
my %tonnagegeno = ();
my $bacon_nosnp = 0;
my $bacon_hetsnp = 0;
my $bacon_homosnp = 0;
my $hass_nosnp = 0;
my $hass_homosnp = 0;
my $hass_hetsnp = 0;
my $simmonds_nosnp = 0;
my $simmonds_homosnp = 0;
my $simmonds_hetsnp = 0;
my $tonnage_nosnp = 0;
my $tonnage_homosnp = 0;
my $tonnage_hetsnp = 0;
my $baconuniq_nosnp = 0;
my $baconuniq_hetsnp = 0;
my $baconuniq_homosnp = 0;
my $hassuniq_nosnp = 0;
my $hassuniq_homosnp = 0;
my $hassuniq_hetsnp = 0;
my $simmondsuniq_nosnp = 0;
my $simmondsuniq_homosnp = 0;
my $simmondsuniq_hetsnp = 0;
my $tonnageuniq_nosnp = 0;
my $tonnageuniq_homosnp = 0;
my $tonnageuniq_hetsnp = 0;
my %avosnp = ();
my $baconuniq = 0;
my $hassuniq = 0;
my $simmondsuniq = 0;
my $tonnageuniq = 0;
my $notuniq = 0;
my $totaluniq = 0;
my %baconuniq = ();
my %hassuniq = ();
my %simmondsuniq = ();
my %tonnageuniq = ();
my %baconhomo = ();
my %baconhet = ();
my %hasshomo = ();
my %hasshet = ();
my %simmondshomo = ();
my %simmondshet = ();
my %tonnagehomo = ();
my %tonnagehet = ();
my %MneWI = ();
my $MneWI = 0;


#Read in the apvo959 data file with annotation ZZMay25_apvo959_out.txt
open (SNP, "$ARGV[0]") || die "Can't open apvo959 infile $!\n";
	while (<SNP>)
{
	chomp;
	if (m/^g/)  	{	#If line starts with isotig name
	my @snp = split (/\t/);
	my $isoname = $snp[0];
	$data959 = join "\t", @snp;
	$hoh{$isoname} = $data959;
	}
}
####################################################################################################################################
#Read in the MneWI datafile with isotigs, number of MneWI snps, total snps per isotig and ratio (DDD_summary_MneWI_isotigs_snps.out) output from avo_M_G_WI.pl
open (MNEWI, "$ARGV[1]") || die "Can't open MneWI infile $!\n";
	while (<MNEWI>)
{
	chomp;
	if (m/^g/)	{
		my @mnewi = split (/\t/);
		my $mnewiname = $mnewi[0];
		$mnewidata = join "\t", @mnewi;
		$hoh2{$mnewiname} = $mnewidata;
	}
}
###################################################################################################################################		
open (OUTFILE, ">>MMM_apvo959_mnewi_out") || die "Can't open apvo959 mnewi outfile $!\n";
	foreach my $key (sort keys %hoh)	{
		if (exists $hoh2{$key})		{
			$hoh3{$key} = "$hoh2{$key}\t$hoh{$key}";
		}
}
		$apvo = keys %hoh;
		$apvomnewi = keys %hoh3;
		print OUTFILE "There were $apvomnewi apvo isotigs that had MneWI snps out of $apvo apvo959 isotigs.\n";
	foreach my $key (sort keys %hoh3)	{
		print OUTFILE "$hoh3{$key}\n\n";
		}

