#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: February 15, 2014
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: Filter mango SNPs for Fluidigm design
# <input>		: [0] All_Mangifera... variant file, [1] parsed blastx report from filter4, [2] 201mer fasta file
# <output>		: Fasta file with 201 mers through filter 6.  
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
my %filter1 = ();
my %cultreads = ();
my %snpname = ();
my %contigname = ();
my %snppos = ();
my %snpgeno = ();
my %filter2 = ();
my %orglength = ();
my %hoh = ();
my %snparray = ();
my %snparraynames = ();
my %filter3 = ();
my $IUPAC = ();
my %snpid = ();
my @blastx = ();
my %filter5 = ();
my %filter6 = ();

####################################################################################################################
#  		
######################################################################################################################
#  Look up reference and alternate allele read info from variant report and make hashes for filtering
open (SNP_ID, "$ARGV[0]") || die "Can't open Mango SNP variant report $!\n";
	while (<SNP_ID>)
{
	chomp;
	if (m/^i/)	{
	@snp = split /\t/, $_;
			$snpname{$snp[0]} = $snp[0];
			$contigname{$snp[0]} = $snp[1];
			$snppos{$snp[0]} = $snp[2];
			$snpgeno{$snp[0]} = "$snp[4]$snp[5]";
			$snp[15] =~ s/\"//g;
			@dp4 = split /,/, $snp[15];
#				if (($snp[9] >= 400)&&($dp4[2]+$dp4[3]>=10))	{ 
			$reads{$snp[0]} = $snp[9]; 
			$refreads{$snp[0]} = $dp4[0] + $dp4[1];
			$altreads{$snp[0]} = $dp4[2] + $dp4[3];
			$refperc{$snp[0]} = ($dp4[0] + $dp4[1])/$snp[9];
			$altperc{$snp[0]} = ($dp4[2] + $dp4[3])/$snp[9];
				for (my $i = 20; $i <= 29; $i++)	{	#for statement uses ; to separate initial, final and increment 
#				print "this is the value of i $i and this is the genotype in the main $snp[$i]\n";
				my @z = split /:/, $snp[$i];
				$z[0] =~ s/"//g;
#				print "this is z[0] $z[0] and this is snp number $snp[0]\n";
#				my $s = "$z[0]\t";	#gets genotype of 10 bins
#				my $t = "$z[2]\t";	#gets read number of 10 bins
				$genotype{$snp[0]} .= "$z[0]:";
				$cultreads{$snp[0]} .= "$z[2]:";
#			print "while i is $i the value of genotype is $genotype{$snp[0]} and the value of cultreads is $cultreads{$snp[0]}\n";
#				print "this is the genotype $genotype{$snp[0]} for snp number $snp[0]\n";
			}
		}
		}
		
#		}


#######################################################################################################################
#filter 5 MAF between 0.4-0.6
#
open (BLASTX, "$ARGV[1]") || die "Can't open parsed blastx report $!\n";
	while (<BLASTX>)	{	
	chomp;
	if (m/^i/)	{
	@blastx = split /\t/, $_;
			my $x = $blastx[0];
				if (($altperc{$x} >= 0.30) && ($altperc{$x} <= 0.70))	{
					$filter5{$x} = "$x\t$contigname{$x}\t$snppos{$x}\t$snpgeno{$x}\t$reads{$x}\t$altreads{$x}\t$altperc{$x}\t$genotype{$x}\t$cultreads{$x}";
#					print "the genotype for snp $x is $genotype{$x}\n";
#  filter for MAF 0.40-0.60 					

		}
	}	
}
close BLASTX;
##########################################################################################################################
# filter 6 A,B,C and TA with 0/1 genotype
#
	for my $d (sort keys %filter5)	{
		my @fil5 = split /\t/, $filter5{$d};
		my @snpgen = split /:/, $fil5[7]; #split into individual bin genotypes
		print "Snpgen array for snp $d is @snpgen\n";
		my $A = $snpgen[4];
		my $B = $snpgen[5];
		my $C = $snpgen[6];
		my $TA = $snpgen[8];
		my $indiv = "$snpgen[0]\t$snpgen[1]\t$snpgen[2]\t$snpgen[3]\t$snpgen[7]\t$snpgen[9]";
		
#		print "Set A genotype is $A\nSet B genotype is $B\nSet C genotype is $C\nTA genotype is $TA\n";
			if (((($A eq "0\/1")&&($B eq "0\/1"))|| (($A eq "0\/1")&&($C eq "0\/1")) || (($C eq "0\/1")&&($B eq "0\/1"))) && ($TA eq "0\/1"))	{
				if (($indiv=~ "1\/1")&&($indiv=~ "0\/0"))	{
#		if ($TA eq "0\/1")	{
			$filter6{$d} = $filter5{$d};
		}
	}
}
		my $filter5keys = keys %filter5;
		my $filter6keys = keys %filter6;

#open (REPORT, >>filter_report.txt)|| die "Can't open filter report$!\n";
	print  "The number of snps at filter 5 is $filter5keys\nThe number of snps at filter 6 is $filter6keys\n";
#######################################################################################################################
#make fasta file with 201 mers through filter 6
#
	my $xmerobj = Bio::SeqIO->new (-file =>"$ARGV[2]", -format => "fasta"); #get fasta objects
	while (my $xmer = $xmerobj->next_seq){
		 $header = $xmer->display_id;
		my $xmerseq = $xmer->seq;
		my $xmerdesc = $xmer->desc;
			if (exists $filter6{$header})	{
				my $mer_out = Bio::Seq->new (-seq => $xmerseq, -display_id => $header, -desc => $xmerdesc); 
				my $mer_outIO = Bio::SeqIO->new (-file => ">>IUPAC_filter6_201_out_less_stringent_30_70_MAF_indiv_00_11.fasta", -format => "fasta"); #creates new fasta file with IUPAC substitutions for Blast against mango CDS
				 $mer_outIO->write_seq($mer_out);
			}
		}	
