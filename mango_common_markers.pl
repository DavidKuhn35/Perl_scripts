#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: November 6, 2012
# Version		: ver. 2.0
# History		: ver. 1.0
# Description	: takes joinmap rawdata file and converts with segregation type to joinmap final format
# <input>		:  joinmap rawdata file	and file of mappable snps, segregation types and parental genotypes
# <output>		:  joinmap final format file 
# 			
# Dependencies	: N/A
# Usage			: joinmap_final_format_v2.pl
#                    Should have used femparent and malparent instead of cultivar names!
##################################################################################
use strict;
use warnings;
use diagnostics;
use Bio::Seq; 
use Bio::SeqIO;
use Bio::SeqUtils;
use Bio::Tools::OddCodes;
use Bio::PrimarySeq;
my $mapsnp = 0;
my $seg_type = 0;
my $parent1 = 0;
my $parent2 = 0;
my %mapsnp_hash = ();
my %P1mapsnp = ();
my %P2mapsnp = ();
my $header = "SNP\t";
my %common = ();
my %joinmap_final_data = ();
my $joinsnp = 0;
my $reformat_data = ();
my $x = ();

#####################################################################################################################
#	main	
#	Read in mappable snps, segregation patterns and parent1 and parent2 genotypes at snp.  Create hashes with seg_type and parent1 genotype and parent2 genotype
open (MAPSNP, $ARGV[0]) || die "Can't open mappable snps file $!\n";
	while (<MAPSNP>)		{
		chomp;
		if (m/^SNP/) {
			my($dummy, $pname) = split/\t/;
			$header .= "$pname\t";
		} 
		if ((m/^Mi/) || (m/^SSK/) || (m/^mango/) || (m/^Contig/))	{	#identifier for snp name (contig for cacao) (SHRSPa for avocado)
			($mapsnp, $seg_type) = split /\t/;	#parent1(female) x parent2(male) cross
				if ($seg_type =~ m/nomap/)	{
				$common{$mapsnp} .= "0\t";
				}else {
					$common{$mapsnp} .= "1\t";}
				}
					}
			
				
######################################################################################################################
#   	
		 my $date = `date`;
		$date =~ s/\s+/_/g;
		$date =~ s/\:/_/g;
		
open (FINAL, ">common/common_markers_$date.txt") || die "Can't open common markers outfile $!\n";
		print FINAL	"$header\n";
		foreach my $key (sort keys %common)	{
			print FINAL "$key\t$common{$key}\n";
			}
				
							
#######################################################################################################################
#   

#########################################################################################################################
#

##########################################################################################################################
#

