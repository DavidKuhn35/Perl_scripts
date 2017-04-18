#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: June 27, 2012
# Version		: ver. 1.1 
# History		: ver. 1.0
# Description	: Compare genotypes of Bacon, Hass, Simmonds, Tonnage from sequencing and Illumina chip
# <input>		:  	genotype files for B,H,S,T		
# <output>		:    table with comparison stats
# 			
# Dependencies	: N/A
# Usage			: avo_seq_chip_compare_table.pl
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
my %seq = ();
my %chip = ();
my %snp = ();
my @wg = ();

#####################################################################################################################
#	Read in files for comparison and make separate hashes for each wg in control

open (INFILE, "$ARGV[0]") || die "Can't open sequence genotype infile $!\n"; #read in sequence genotype data for a parent
	while (<INFILE>)	{
		my @x = split /\t/;
		$seq{$x[0]} = "$x[1]\t$x[2]";
		}
	
close INFILE;

open (INFILE, "$ARGV[1]") || die "Can't open chip genotype infile $!\n"; #read in chip genotype data for all controls for one parent
	while (<INFILE>)	{
		my @x = split /\t/;
		push (@wg, $x[3]);
		$chip{$x[3]}{$x[0]} = "$x[1]\t$x[2]";	#hoh with wg as first keys, snps as second keys and genotypes as values
		$snp{$x[0]}{$x[3]} = "$x[1]\t$x[2]";	#hoh with snps as first keys, wg as second keys and genotypes as values
		}	
close INFILE;
#######################################################################################################################
#	Compare chip snps to sequence snps, print out results with genotypes for differences chip first sequence second
open (OUTFILE, ">>1CCC_$ARGV[2]\_chip_seq_diff") || die "Can't open compare $ARGV[2] controls outfile $!\n";
		foreach my $y ( sort keys %snp)	{
			foreach my $z (sort keys %{$snp{$y}})	{
				if (exists $seq{$y})	{
					if ($snp{$y}{$z} ne $seq{$y})	{
					print OUTFILE "$y\t$z\t$snp{$y}{$z}\t$seq{$y}\tdiff_seq\n";
				}
			}
		}
	}
close OUTFILE;

#########################################################################################################################
#	Compare chip to other chip replicates, print out results with genotypes for differences, reference chip (first wg) second and other replicates
#	first
open (	OUTFILE, ">>1GGG_$ARGV[2]\_chip_chip_diff" ) || die "Can't open compare $ARGV[2] chip controls outfile $!\n";				
		foreach my $snp (sort keys %snp)	{
			foreach my $wg  (sort keys %{$snp{$snp}})	{
				if ($wg ne $wg[0])	{
					if ($snp{$snp}{$wg} ne $snp{$snp}{$wg[0]})	{
					print OUTFILE "$snp\t$wg\t$snp{$snp}{$wg}\t$snp{$snp}{$wg[0]}\tdiff_seq\n";
					}
				}
			}
		}
###########################################################################################################################


