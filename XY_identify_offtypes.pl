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
my $pop_final_header = ();
my %joinmap_final_data = ();
my $joinsnp = 0;
my $reformat_data = ();
my $x = ();
my %offtype = ();
my %ot = (
	1 => "XY|YY",
	2=> "YY",
	3=> "XX|YY",
	4=> "YY",
	5=> "QQ",
	6=> "XX",
	7=> "XX|YY",
	8=> "XX",
	9=> "XX|XY",
	10=> "WW"
	);

#####################################################################################################################
#	main	
#	Read in mappable snps, segregation patterns and parent1 and parent2 genotypes at snp.  Create hashes with seg_type and parent1 genotype and parent2 genotype
open (MAPSNP, $ARGV[0]) || die "Can't open mappable snps file $!\n";
	while (<MAPSNP>)		{
		chomp;
		if ((m/^Mi/) || (m/^SSK/) || (m/^mango/) || (m/^Contig/))	{	#identifier for snp name (contig for cacao) (SHRSPa for avocado)
		($mapsnp, $seg_type, $parent1, $parent2) = split /\t/;	#parent1(female) x parent2(male) cross
			$mapsnp_hash{$mapsnp}= $seg_type;
			$P1mapsnp{$mapsnp} = $parent1;
			$P2mapsnp{$mapsnp} = $parent2;
		}
	}
					
				
######################################################################################################################
#   	Read in joinmap_rawdata file
open (RAW, $ARGV[1]) || die "Can't open joinmap rawdata file $!\n";
		while (<RAW>)		{
			chomp;
				my $dnk = $_;
				($joinsnp) = split /\t/, $_;
					if  (m/^SNP/)	{
#					$dnk =~ s/^SNP\t//;
					$pop_final_header .= $dnk;			
			}elsif ((m/^Mi/) || (m/^SSK/) || (m/^mango/) || (m/^Contig/))	{		#Get snp name (contig for cacao) (SHRSPa for avocado) (Mi, SSK, mango,Contig for mango)
				 
				my $x = $_;	#make a copy of current line
				$offtype{$joinsnp} = &otsub ($mapsnp_hash{$joinsnp},$x);
		}
	}

		 my $date = `date`;
		$date =~ s/\s+/_/g;
		$date =~ s/\:/_/g;
		my $popname = $ARGV[2];	#read in popname
open (FINAL, ">>$popname/$popname\_offtype_final_$date.txt") || die "Can't open $popname offtype final format outfile $!\n";
		print FINAL	"SNPname\tSegregation type\tParent1 geno\tParent2 geno\t$pop_final_header\n";
		foreach my $key (sort keys %offtype)	{
			print FINAL "$key\t$mapsnp_hash{$key}\t$P1mapsnp{$key}\t$P2mapsnp{$key}\t$offtype{$key}\n";
			}
				
							
#######################################################################################################################
#   sub ot gets segtype and rawgeno data for snp and global substitutes OT for offtype genotypes
sub otsub {
		my ($seg,$l) = @_;
			$x = $ot{$seg};
			$l =~ s/$x/OT/g;
			print "$x\t$l\n";
			$l;	
		
			

		}
	
	



#########################################################################################################################
#

##########################################################################################################################
#

