#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: 5/25/2012
# Version		: ver. 1.0.1
# History		: ver. 1.0
# Description	: reads in file with apvo AGI names, reads in avo snp data file with AGI names, SNP frequency, GOslim annotation, makes new hash with isotigs that match the apvo AGI and #	prints out that as a file
# <input>		: apvo959_out.txt file, May25_avo_all_isolen_snpcount_snpfreq_agi_goslim_out file (redid GOSlim annotation file due to bug in goslim_reduced 
#	program			
# <output>		:    ZZMay25_apvo959_avo_out.txt
# 			
# Dependencies	: N/A
# Usage			: avo_apvo.pl	apvo959_out.txt May25_avo_all_isolen_snpcount_snpfreq_agi_goslim_out
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
my %apvo = ();
my %avo_all = ();
my %apvo_avo = ();
open (APVO, "$ARGV[0]") || die "Can't open apvo infile $!\n";	#read in apvo file apvo959_out.txt which is basically the apvo959 database printed out
	while (<APVO>)
{
	chomp;
	(my $apvo_agi, my $apvo_annot) = split (/\t/);
	$apvo{$apvo_agi} = $apvo_annot;	#make hash with apvo annot data agi=>apvo_annot
}

open (AVO_ALL, "$ARGV[1]") || die "Can't open avo_all infile $!\n"; #read in avo_all_May3_isolen_snpcount_snpfreq_agi_goslim_out which is cleaned up version of all annotation data
	while (<AVO_ALL>)	
{
	chomp;
	if (m/^g/)	{	#look for lines that start with isotig name
	my @x = split (/\t/);
	$avo_all{$x[0]}{AGI} = $x[4];
	$avo_all{$x[0]}{annot} = join "\t", @x;
}
}

	for	my $key (sort keys %avo_all)	{
			my $z = $avo_all{$key}{AGI};
			if (exists $apvo{$z})	{
			$apvo_avo{$key} = $avo_all{$key}{annot};	#make new hash with only isotigs that match the apvo AGI names
			}
	}

open (OUTFILE, ">>ZZMay25_apvo959_avo_out.txt") || die "Can't open apvo959 avo outfile $!\n";
				my $goslimt_default = join "\t", ("C","Cdate","Cgoslimterm","Cgoid","Crelation","Cgoterm","F","Fdate","Fgoslimterm","Fgoid","Frelation","Fgoterm","P","Pdate","Pgoslimterm","Pgoid","Prelation","Pgoterm");

	print OUTFILE "isotig_name\tisotig_length\tsnp_count\tsnp_frequency\tAGI\t$goslimt_default\n";

	for my $q (sort keys %apvo_avo)	{
		print OUTFILE "$apvo_avo{$q}\n";
}


