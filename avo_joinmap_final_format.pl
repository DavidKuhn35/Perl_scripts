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
my $simmonds = 0;
my $tonnage = 0;
my %mapsnp_hash = ();
my %simparent = ();
my %tonparent = ();
my $joinmap_final_header = "SNPname\tSegregation type\t";
my %joinmap_final_data = ();
my $joinsnp = 0;
my $reformat_data = ();

#####################################################################################################################
#	main	
#	Read in mappable snps, segregation patterns and Simmonds and Tonnage genotypes at snp.  Create hashes with seg_type and simmonds parent genotype and tonnage parent genotype
open (MAPSNP, $ARGV[0]) || die "Can't open mappable snps file $!\n";
	while (<MAPSNP>)		{
		chomp;
		if (m/^SHRSPa/)	{
#			($mapsnp, $seg_type, $simmonds, $tonnage) = split /\t/;	#Simmonds x Tonnage cross
			($mapsnp, $seg_type, $tonnage, $simmonds) = split /\t/;	#Tonnage x Simmonds cross
			$mapsnp_hash{$mapsnp}= $seg_type;
			$simparent{$mapsnp} = $simmonds;
			$tonparent{$mapsnp} = $tonnage;
#			print "$mapsnp\t$tonparent{$mapsnp}\n";
		}
	}
					
				
######################################################################################################################
#   	Read in joinmap_rawdata file
open (RAW, $ARGV[1]) || die "Can't open joinmap rawdata file $!\n";
		while (<RAW>)		{
			chomp;
			if  (m/^SNPname/)	{
				$joinmap_final_header .= $_;			
			}elsif (m/^SHRSPa/)	{		#Get snp name
				($joinsnp) = split /\t/;
					if ($mapsnp_hash{$joinsnp} eq "\<nnxnp\>")	{
					my $x = $_;
#					my $n = $simparent{$joinsnp};	#Simmonds x Tonnage cross
					my $n = $tonparent{$joinsnp};	#Tonnage x Simmonds cross
					$n =~ s/^\b(\w)\w/$1/;
					$x =~ s/$n/n/g;
					my $a = $x;
#					my $p = $tonparent{$joinsnp};	#Simmonds x Tonnage cross
					my $p = $simparent{$joinsnp};	#Tonnage x Simmonds cross
					$p =~ s/$n//;
					my $pee = $p;
					$pee =~s/\s+(\w)\s+/$1/;
#					print "This $p is p for snp $joinsnp\n";
					$a =~ s/$pee/p/g;
#					print "This is the genotypes with substitutions for p\n$a\n";
					$joinmap_final_data{$joinsnp} = $a;
					}elsif ($mapsnp_hash{$joinsnp} eq "\<lmxll\>")	{
					my $y = $_;
#					my $ell = $tonparent{$joinsnp};  #Simmonds x Tonnage cross
					my $ell = $simparent{$joinsnp};	  #Tonnage x Simmonds cross
					$ell =~ s/^\b(\w)\w/$1/;
#					print "This is the memory parameter for l $1 for SNP $joinsnp\n";
#					print "This is $l for l for SNP $joinsnp\n";
					$y =~ s/$ell/l/g;
#					print "Here are the genotypes after substitution for l\n$y\n";
#					my $m = $simparent{$joinsnp};	#Simmonds x Tonnage cross
					my $m = $tonparent{$joinsnp};	#Tonnage x Simmonds cross
					$m =~ s/$ell//;
					$y =~ s/$m/m/g;
#					print "Here are the genotypes after substitution for m \n$y\n";
					$joinmap_final_data{$joinsnp} = $y;
					}elsif ($mapsnp_hash{$joinsnp} eq "\<hkxhk\>")	{
					my $z = $_;
#					my $h = $simparent{$joinsnp};	#Simmonds x Tonnage cross
					my $h = $tonparent{$joinsnp};	#Tonnage x Simmonds cross
					$h =~ s/^\b(\w)\w/$1/;
					$z =~ s/$h/h/g;
#					my $k = $simparent{$joinsnp};	#Simmonds x Tonnage cross
					my $k = $tonparent{$joinsnp};	#Tonnage x Simmonds cross
					$k =~ s/$h//;
					$z =~ s/$k/k/g;
					$joinmap_final_data{$joinsnp}= $z;


		}
	}
}

open (FINAL, ">>FLAMP_TxS_joinmap_final.txt") || die "Can't open FLAMP_TxS joinmap final format outfile $!\n";
		print FINAL	"SNPname\tSegregation type\tTonnage geno\tSimmonds geno\t$joinmap_final_header\n";
		foreach my $key (sort keys %joinmap_final_data)	{
			print FINAL "$key\t$mapsnp_hash{$key}\t$tonparent{$key}\t$simparent{$key}\t$joinmap_final_data{$key}\n";
			}
				
							
#######################################################################################################################
#   

#########################################################################################################################
#

##########################################################################################################################
#

