#!/usr/bin/perl -w
##################################################################################
# Author		: David Kuhn
# Date			: 4/19/2012
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: 
# <input>		:  unifiedSNPs_Avocado_BLASTX_TAIR10_ANNOT.table_cvg.txt			
# <output>		:    
# 			
# Dependencies	: N/A
# Usage			: avo_snp_category.pl
#   STILL UNDER DEVELOPMENT AND MAY NO LONGER BE NECESSARY.  DOES CORRECTLY PASS HASH TO SUBROUTINE.                 
##################################################################################

use strict;
use warnings;
use diagnostics;
use Statistics::Descriptive;
my @snp = ();
my %hoh = ();
# %isotig = ();
# $goodsnp = 0;
my $badsnp = 0;
# %isotiglength = ();
# @snpdistribution = ();
my %annot = ();
my %curator = ();
my %fullannot = ();
my %prot_mod = ();
# %inverse_annot = ();
# %inverse_fullannot = ();
my %nodupannot = ();
my %nodupfullannot = ();
my %involved_in = ();
#my $inv_in = "lala";
#my $inv_in = ();
#Read in the avo snp file and get rid of rows with missing data or Ns
open (SNP, "$ARGV[0]") || die "Can't open SNP infile\n";
	while (<SNP>)
{
	chomp;
	if (m/^SEQID/)  	{	#If line has SNP name with Tcm
	@snp = split (/\t/);
#	if (($snp[2] eq "N")||($snp[3] eq "N")||($snp[5]<=6)||($snp[13]<=6)||($snp[21]<=6)||($snp[29]<= 6))	{
	if (($snp[2] eq "N")||($snp[3] eq "N"))	{
		$badsnp++;
	}	else	{
		$hoh{$snp[36]}{$snp[1]}= @snp;
		my $goodsnp++;
#		$isotiglength{$snp[36]} = $snp[37];
#		$prot_mod{$snp[36]} = $snp[51];	#arabidopsis gene model name, should give minimum number of unique annotations
#		$prot_type{$snp[36]} = $snp[52];	#all the same protein type, this should give 1 key when reversed
#		$annot{$snp[36]} = $snp[53];	#short description
#		$curator{$snp[36]} = $snp[54];	#curators description
		$fullannot{$snp[36]} = $snp[55];	#full description GO terms FUNCTIONS IN:  LOCATED IN:  INVOLVED IN:
		my $inv_in = $fullannot{$snp[36]};
#			if ($inv_in =~ m/INVOLVED IN:(.*); LOCATED/)	{
#			print "$1\n";
			if ($inv_in =~ m/INVOLVED IN: (.*),/)	{
			$involved_in{$snp[36]} = $1;	
		}
	}
}
}
	open (TESTER, ">>XX$ARGV[0]\_tester.out") || die "Can't open tester outfile";

my $hoh_count =&keycount(\%hoh);

print TESTER "Less stringent hoh key count is $hoh_count\nBadsnp is $badsnp\nGoodsnp is $goodsnp";

#my %annot_nodup = &inverthash(\%annot);
#my $annot_nodup_count = &keycount(\%annot_nodup);
#	my $annot_count = keys (%annot);
#	print "keys for annot is $annot_count\nkeys for annot_nodup_count is $annot_nodup_count\n";
#&hashprint("snp53",\%annot_nodup);

#my %prot_mod_nodup = &inverthash(\%prot_mod);
#my $prot_mod_nodup_count = &keycount(\%prot_mod_nodup);
#	print "keys for prot_mod nodup is $prot_mod_nodup_count\n";
#&hashprint ("snp51",\%prot_mod_nodup); 

#my %curator_nodup = &inverthash(\%curator);
#&hashprint ("snp54", \%curator_nodup);

#my %fullannot_nodup = &inverthash(\%fullannot);
#&hashprint ("snp55", \%fullannot_nodup);

my %involved_in_nodup = &inverthash(\%involved_in);
&hashprint ("involved_in", \%involved_in_nodup);


sub keycount	{
#	print "name of hash sent to subroutine is $_[0]\n";
	my $keycount = keys (%{$_[0]});
#	print "subroutine key count is $keycount\n";
}













sub inverthash	{
	my %z = reverse %{$_[0]};
	my %x = reverse %z;
	}

sub hashprint {
	print TESTER "name of hash is $_[0]\n";
	my $x = keys (%{$_[1]});
	print TESTER "hash $_[0] has $x keys\n";
	for my $key (sort keys %{$_[1]})	{
	print TESTER "$key\t${$_[1]}{$key}\n";
	}
}

