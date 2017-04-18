#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: February 4, 2014
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: Find synonymous mango SNPs for Fluidigm design
# <input>		: Mangifera_Run2_VarCountsWithTAIR10blast.txt, FirstHalfMangifera_ALL_VarRun2_withTAIR10blast.txt (SecondHalf...), renamed_with_position_Mangifera_All_Run2_SNPflanks_201.fasta.txt (...401..) 		
# <output>		: File with IUPAC 201 or 401 mers to BLAST against mango CDS  
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

#####################################################################################################################
#	main	Make hashes of all contigs with number of variants and blast results
my @snp = ();
my $header = ();
my $header1 = ();
my %goodsnp = ();
my $snpreads = ();
my %refreads = ();
my %altreads = ();
my %refperc = ();
my %altperc = ();
my %reads = ();
my @dp4 = ();
my %score = ();
my %genotype = ();


####################################################################################################################
#  Read in list of good snps
open (GOODSNP, "$ARGV[0]") || die "Can't read in good snp file$!\n";
	while (<GOODSNP>)	{
		chomp;
		my $snpname = $_;
		$snpname =~ s/>//; #remove fasta header start
		$goodsnp{$snpname}= $snpname;
	}
close GOODSNP;				
######################################################################################################################
#  Look up reference and alternate allele read info from variant report and make hashes
open (SNP_ID, "$ARGV[1]") || die "Can't open Mango SNP variant report $!\n";
	while (<SNP_ID>)
{
	chomp;
	if (m/^i/)	{
	@snp = split /\t/, $_;
			if (exists $goodsnp{$snp[0]})	{
			$snp[15] =~ s/\"//g;
			@dp4 = split /,/, $snp[15];
				if (($snp[9] >= 200)&&($dp4[2]+$dp4[3]>=10))	{ 
			$reads{$snp[0]} = $snp[9]; 
			$refreads{$snp[0]} = $dp4[0] + $dp4[1];
			$altreads{$snp[0]} = $dp4[2] + $dp4[3];
			$refperc{$snp[0]} = ($dp4[0] + $dp4[1])/$snp[9];
			$altperc{$snp[0]} = ($dp4[2] + $dp4[3])/$snp[9];
				for (my $i = 20; $i <= 29; $i++)	{
				print "this is the value of i $i and this is the genotype in the main $snp[$i]\n";
				my @z = split /:/, $snp[$i];
				$z[0] =~ s/"//g;
				print "this is z[0] $z[0] and this is snp number $snp[0]\n";
				my $s = "$z[0]\t";
#			my $intgeno = &genotype ($snp[$i]);
				$genotype{$snp[0]} .= $s;
				print "this is the genotype $genotype{$snp[0]} for snp number $snp[0]\n";
			}
		}
		}
		
		}
}
###################################################################################################################
#  Sort hashes by value to find snps with greatest number of reads and greatest minor allele frequency
#
	my @readsort = sort by_reads keys %reads;
	my @altreadsort = sort by_altreads keys %altreads;
	my @altpercsort = sort by_altperc keys %altperc;
#	foreach my $y (@altpercsort)	{
#		print "$y\t$reads{$y}\t$altreads{$y}\t$altperc{$y}\n";
#		}
	
#######################################################################################################################
#print out results in descending order of minor allele frequency  
#
open (MAF,">>mango_maf_out.txt") || die "Can't print to mango MAF file$!\n";
	print MAF "snpname\treads\taltreads\taltpercent\n";

	foreach my $x (@altpercsort)	{
			if (($altperc{$x} <= 0.55) && ($altperc{$x} >= 0.25))	{
			
			printf MAF "$x\t$reads{$x}\t$altreads{$x}\t%.2f%%\t$genotype{$x}\n",$altperc{$x};
		}
	}	

#############################################################################################################
#subroutines to sort hashes by value
#

sub by_altperc { $altperc{$b} <=> $altperc{$a} }
sub by_reads	{$reads{$b} <=> $reads{$a}	}
sub by_altreads {$altreads{$b} <=> $altreads{$a} }
sub genotype {
	my @q = split /:/, $_;
	print "this is the genotype in the sub $q[0]\n";
	my $return = "$q[0]\t"; 
	}

