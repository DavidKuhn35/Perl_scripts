#!/usr/bin/perl -w
##################################################################################
# Author		: David Kuhn
# Date			: 4/25/2012
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: Get isotigs related to APVO genes, print out files so statistical analysis can be done for SNPfrequency distribution, isotig length distribution, etc.
# <input>		:  			
# <output>		:    
# 			
# Dependencies	: N/A
# Usage			: 
#    PROGRAM STILL UNDER DEVELOPMENT.  GAVE DATA IN ANOTHER FORMAT TO GUILIANA MUSTIGA FOR SAS ANALYSIS.                
##################################################################################
use Statistics::Descriptive;
my %hoh = ();
my %isotig = ();
my $goodsnp = 0;
my $badsnp = 0;
my @snpdistribution = ();
#Read in the avo snp file and get rid of rows with missing data or Ns
open (SNP, "$ARGV[0]") || die "Can't open SNP infile";
	while (<SNP>)
{
	chomp;
	if (m/^SEQID/)  	{	#If line has SNP name with Tcm
	@snp = split (/\t/);
#	if (($snp[2] eq "N")||($snp[3] eq "N")||($snp[5]<=6)||($snp[13]<=6)||($snp[21]<=6)||($snp[29]<= 6))	{
	if (($snp[2] eq "N")||($snp[3] eq "N"))	{
		$badsnp++;
	}	else	{
		$hoh{$snp[36]}{$snp[1]} = @snp;
		$isotig{$snp[36]}{isotiglength} = $snp[37];
		$isotig{$snp[36]}{atg_name} = $snp[40];	
		}
	}
}


open (OUTFILE, ">>QQ$ARGV[0]\_out.txt") || die "Can't open outfile";

for $x (sort keys %hoh) {
		my $count = keys %{$hoh{$x}};
		$isotig{$x}{snpcount} = $count;	#adds snp count value to hash for key of isotig
#		push (@snpdistribution, $count);
#		$goodsnp = $goodsnp + $count;
		$isotig{$x}{snpfrequency} = $count/$isotig{$x}{isotiglength};	#adds snp frequency as value to isotig hash with isotigs as keys
#		push (@snpfreqdistribution, $snpfrequency);

#		push (@isolendis, $isotig{$x}{isotiglength});
#		printf OUTFILE "$x\t$count\t$isotig{$x}{isotiglength}\t%1.3f\n", $snpfrequency;
	}


tie (%apvo959, "DB_File", "apvo959.db", O_CREAT|O_RDWR, 0644, $DB_HASH) || die "Cannot create or open apvo.db: $!";

for $x (sort keys %isotig)	{
	for $y (sort keys %apvo959)	{
		if  ($y =~ m/$isotig{$x}{atg_name}/)	{
			$apvo{$x} = $y;	#makes a new hash with isotigs as keys that have apvo959 annotation and the AT.G name as the value.  If this hash is sorted by keys, can get all the data (snp count, snp frequency) on the isotigs that match apvo959 genes
			}
		}
	}

	for $q  (sort keys %isotig)	{
		if (exists $apvo{$q})	{
			push (@snpdistribution, $isotig{$q}{snpcount});
			push (@snpfreqdistribution, $isotig{$q}{snpfrequency});
			push (@isolendis, $isotig{$q}{isotiglength});
			}
		}


open (SNPDISTRIBUTION, ">>APVO$ARGV[0]\_snpdistribution_out.txt") || die "Can't open snpdistribution outfile";
$stat = Statistics::Descriptive::Full->new();
$stat->add_data(@snpdistribution);
$stat->get_data();
%fsnpdistribution = $stat->frequency_distribution(100);
for (sort {$a <=> $b} keys %fsnpdistribution) 	{
	print SNPDISTRIBUTION "$_\t$fsnpdistribution{$_}\n";
}
open (FREQDISTRIBUTION, ">>APVO$ARGV[0]\_freqdistribution_out.txt") || die "Can't open freqdistribution outfile";
$freq = Statistics::Descriptive::Full->new();
$freq->add_data(@snpfreqdistribution);
$freq->get_data();
%freqdistribution = $freq->frequency_distribution(100);
for (sort {$a <=> $b} keys %freqdistribution)	{
	printf FREQDISTRIBUTION "%1.3f\t$freqdistribution{$_}\n", $_;
}

open (ISOTIGLENGTHDIS, ">>APVO$ARGV[0]\_isotiglengthdis_out.txt") || die "Can't open lengthdistribution outfile";
$isolen = Statistics::Descriptive::Full->new();
$isolen->add_data(@isolendis);
%isolendis = $isolen->frequency_distribution(100);
for (sort {$a <=> $b} keys %isolendis)	{
	printf ISOTIGLENGTHDIS "%5d\t$isolendis{$_}\n", $_;
	}
	$totalsnp = $goodsnp + $badsnp;
	$totalcontigs = keys %hoh;

open	(SUMMARY, ">>APVO$ARGV[0]\_summary.txt") || die "Can't open summary file";
	print SUMMARY "Summary of data in $ARGV[0]\nTotal contigs with SNPs is $totalcontigs\nTotal SNPs equals $totalsnp\nGood SNPs equals $goodsnp\nBad SNPs equals $badsnp\n";
