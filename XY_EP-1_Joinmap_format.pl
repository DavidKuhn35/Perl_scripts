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
my %jmap = (
	1 => "nomap",
	2=> "<lmxll>",
	3=> "nomap",
	4=> "<nnxnp>",
	5=> "<hkxhk>",
	6=> "<nnxnp>",
	7=> "nomap",
	8=> "<lmxll>",
	9=> "nomap",
	10=> "nomap"
	);
my %hoh = ();
my %ref = ();
my $q = ();
my $jseg = ();
my $snpheader = ();
my  $snip = ();
my $segr = ();

#####################################################################################################################
#	main	
#	Read in mappable snps, segregation patterns and parent1 and parent2 genotypes at snp.  Create hashes with seg_type and parent1 genotype and parent2 genotype

				
######################################################################################################################
#   	read in XY offtype file for pop  split into seg types
open (RAW, $ARGV[0]) || die "Can't open offtype rawdata file $!\n";
		while (<RAW>)		{
			chomp;

				 my @jmp = split /\t/;
				 my	$snip = shift @jmp;
				 my	$segr = shift @jmp;
#				 $segr =~ s/nodata/10/;
				 my	$p1 = shift @jmp;
				 my	$p2 = shift @jmp;
				 my $waste = shift @jmp;


				 if  ($snip eq "SNPname")	{
#				 	foreach my $x (@jmp)	{
#				 		 $snpheader .= "$x\t";
#				 	}
						$snpheader = join "\t", @jmp;				 		
# could also use $snpheader = join,\t,@jmp to get tabs between snp names				 		
				 		}else{
				 			foreach my $g (@jmp)	{
				 				$hoh{$segr}{$snip} .= "$g\t";

				 			}
				 			
				 		}
				 	}
				 $joinmap_final_header = "$snip\t$segr\t$snpheader\n";	
################################################################################################################################################################
#  convert to joinmap format

	foreach my $key (sort keys %hoh)	{
		if (($key == 2) || ($key == 8))	{
			 $jseg = "<lmxll>";
			foreach my $d (sort keys %{$hoh{$key}})	{
				$q = $hoh{$key}{$d};
				print "$q\n";
				$q =~ s/XY/lm/g;
				$q =~ s/XX|YY/ll/g;
				$ref{$d} = "$d\t$jseg\t$q";	#create a hash with snp as key and snp\tseg\treformatted geno as value
				}
		}  elsif (($key ==4)||($key == 6))	{
			 $jseg = "<nnxnp>";
			foreach my $d (sort keys %{$hoh{$key}})	{
				$q = $hoh{$key}{$d};
				$q =~ s/XY/np/g;
				$q =~ s/XX|YY/nn/g;
				$ref{$d} = "$d\t$jseg\t$q";	#create a hash with snp as key and snp\tseg\treformatted geno as value
				} 
		}	elsif ($key == 5)	{
			 $jseg = "<hkxhk>";
			foreach my $d (sort keys %{$hoh{$key}})	{
				$q = $hoh{$key}{$d};
				$q =~ s/XY/hk/g;
				$q =~ s/XX/hh/g;
				$q =~ s/YY/kk/g;
				$ref{$d} = "$d\t$jseg\t$q";	#create a hash with snp as key and snp\tseg\treformatted geno as value
				}
		}

	}



		 my $date = `date`;
		$date =~ s/\s+/_/g;
		$date =~ s/\:/_/g;
		my $popname = $ARGV[1];	#read in popname
open (FINAL, ">>$popname/$popname\_XYjoinmap_final_$date.txt") || die "Can't open $popname joinmap final format outfile $!\n";
		print FINAL	"$joinmap_final_header\n";
		foreach my $key (sort keys %ref)	{
			print FINAL "$ref{$key}\n";
			}
				
							
#######################################################################################################################
  

#########################################################################################################################
#

##########################################################################################################################
#

