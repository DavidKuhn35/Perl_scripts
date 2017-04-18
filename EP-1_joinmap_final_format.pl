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
my $joinmap_final_header = ();
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
				($joinsnp) = split /\t/, $_;
			if  (m/^SNP/)	{
#		my		$dnk =~ s/^$joinsnp\t//;
				$joinmap_final_header = $_;			
			}elsif ((m/^Mi/) || (m/^SSK/) || (m/^mango/) || (m/^Contig/))	{		#Get snp name (contig for cacao) (SHRSPa for avocado) (Mi, SSK, mango,Contig for mango)
				 
				 $x = $_;	#make a copy of current line
					$x =~ s/Z/-/g; #global replacement of Z with Joinmap missing data symbol
					$x =~ s//-/g;
#					$x =~ s/^$joinsnp\t//;
			}		if ($mapsnp_hash{$joinsnp} eq "\<nnxnp\>")	{
					
				
					my $n = $P1mapsnp{$joinsnp};	#parent1 is homozygous for n nnxnp, this gets nn
					$n =~ s/^\b(\w)\w/$1/;  #identify nucleotide that is n
					my $q = $x;
					$q =~ s/$n/n/g;  #global substitution in the line for all nucleotides that are n
					my $a = $q;	#make a copy of the modified line
					my $p = $P2mapsnp{$joinsnp};	#parent2 should be np, can identify p by removing n
					$p =~ s/$n//;  #remove n, $p = p
					$a =~ s/$p/p/g;	#global substitution of p for nucleotide that is p in modified line
					$a =~ s/pp/--/g; #global substitution to change impossible genotypes to missing data
					$joinmap_final_data{$joinsnp} = $a;  #save modified line, should have all nucleotides replaced with n and p
					}elsif ($mapsnp_hash{$joinsnp} eq "\<lmxll\>")	{
					my $y = $x;	#make a copy of line
#					$y =~ s/Z/-/g, $x; #global replacement of Z with Joinmap missing data symbol
					my $ell = $P2mapsnp{$joinsnp};	  # parent 2 should be homozygous for l
					$ell =~ s/^\b(\w)\w/$1/;	#identify nucleotide that is l
					$y =~ s/$ell/l/g;	#substitute l for nucleotide in the line
					my $m = $P1mapsnp{$joinsnp};	#parent 1 should be heterozygous lm
					$m =~ s/$ell//;	#get rid of l, only m should be left
					$y =~ s/$m/m/g;	#substitute in line for each nucleotide that is m
					$y =~ s/mm/--/g; #global substitution to change impossible genotypes to missing data
					$joinmap_final_data{$joinsnp} = $y;
					}elsif ($mapsnp_hash{$joinsnp} eq "\<hkxhk\>")	{
					my $z = $x;
#					$z =~ s/Z/-/g; #global replacement of Z with Joinmap missing data symbol
					my $h = $P1mapsnp{$joinsnp};	#P1 genotype
					$h =~ s/^\b(\w)\w/$1/;	#identify first nucleotide
					$z =~ s/$h/h/g;	#substitute h for first nucleotide
					my $k = $P1mapsnp{$joinsnp};	#get second nucleotide
					$k =~ s/$h//;	#identify second nucleotide by removing h
					$z =~ s/$k/k/g;	#global replacement for k in line
					$joinmap_final_data{$joinsnp}= $z;
					
					}
			


		}
	

		 my $date = `date`;
		$date =~ s/\s+/_/g;
		$date =~ s/\:/_/g;
		my $popname = $ARGV[2];	#read in popname
open (FINAL, ">>$popname/$popname\_joinmap_final_$date.txt") || die "Can't open $popname joinmap final format outfile $!\n";
		print FINAL	"SNPname\tSegregation type\tParent1 geno\tParent2 geno\t$joinmap_final_header\n";
		foreach my $key (sort keys %joinmap_final_data)	{
			print FINAL "$key\t$mapsnp_hash{$key}\t$P1mapsnp{$key}\t$P2mapsnp{$key}\t$joinmap_final_data{$key}\n";
			}
				
							
#######################################################################################################################
#   

#########################################################################################################################
#

##########################################################################################################################
#

