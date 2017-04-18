#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: April 2, 2015
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: 	Read in consensus genotypes for female then male parent.  Determine segregation type. Print out  
# <input>		:  SNP data, parents names	
# <output>		:   consensus genotypes for parents, segregation patterns for each locus.
# 			
# Dependencies	: N/A
# Usage			: EP1_identify_segregation_types.pl
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
my %hoh = ();
my %hoh2 = ();
my %snp_hash = ();
my %dna_hash = ();
my %fem = ();
my %mal = ();
my %snpname = ();
my %IxKP = ();

my $dna = 0;
my $geno1 = 0;
my $geno2 = 0;
my $print_header = "SNP name\tSegregation type\tfemale genotype\tmale genotype\n";
my %segregation_type = ();
my $print_header2 = "SNP name\t";
my %snp_geno = ();

#####################################################################################################################
	
#########################################################################################################################
#	


#####################################################################################################################
#	Read in female parent and male parent consensus genotype
		open (FEM, "$ARGV[0]") || die "can't open fem file $!\n";
		while (<FEM>)	{
			chomp;
			unless (m/SNP/) {
			my ($snp, $femgen) = split /\t/, $_;
			$fem{$snp} = $femgen;
			$snpname{$snp} = $snp;
			}
		} 
		open (MAL, "$ARGV[1]") || die "can't open mal file $!\n";
		while (<MAL>)	{
			chomp;
			unless (m/SNP/)	{
			my ($snp, $malgen) = split /\t/, $_;
			$mal{$snp} = $malgen;
			$snpname{$snp} = $snp;
			}
		} 
######################################################################################################################
#  create list of segregation types
	foreach my $key (sort keys %snpname)	{	
		if ((($fem{$key} eq "AA") || ($fem{$key} eq "CC") || ($fem{$key} eq "GG") || ($fem{$key} eq "TT" )) && (($mal{$key} ne "AA") && ($mal{$key} ne "CC") && ($mal{$key} ne "GG") && ($mal{$key} ne "TT") && ($mal{$key} ne "") && ($mal{$key} ne "ZZ")))	{
			$segregation_type{$key} = "\<nnxnp\>";  #parent1 x parent2 segregation type
#			$segregation_type{$key} = "\<lmxll\>";	#parent2 x parent1 segregation type
		}elsif ((($fem{$key} ne "AA") && ($fem{$key} ne "CC") && ($fem{$key} ne "GG") && ($fem{$key} ne "TT") && ($fem{$key} ne "") && ($fem{$key} ne "ZZ")) && (($mal{$key} eq "AA") || ($mal{$key} eq "CC") || ($mal{$key} eq "GG") || ($mal{$key} eq "TT")))	{
			$segregation_type{$key} = "\<lmxll\>";	#parent1 x parent2 segregation type
#			$segregation_type{$key} = "\<nnxnp\>";	#parent2 x parent1 segregation type
		}elsif ((($fem{$key} ne "AA") && ($fem{$key} ne "CC") && ($fem{$key} ne "GG") && ($fem{$key} ne "TT") && ($fem{$key} ne "") && ($fem{$key} ne "ZZ")) && (($mal{$key} ne "AA") && ($mal{$key} ne "CC") && ($mal{$key} ne "GG") && ($mal{$key} ne "TT") && ($mal{$key} ne "") && ($mal{$key} ne "ZZ")))	{
			$segregation_type{$key} = "\<hkxhk\>";
			} elsif ((($fem{$key} eq "AA") || ($fem{$key} eq "CC") || ($fem{$key} eq "GG") || ($fem{$key} eq "TT")|| ($fem{$key} eq "") || ($fem{$key} eq "ZZ")) && (($mal{$key} eq "AA") || ($mal{$key} eq "CC") || ($mal{$key} eq "GG") || ($mal{$key} eq "TT") ||($mal{$key} eq "") || ($mal{$key} eq "ZZ")))	{
			$segregation_type{$key} = "nomap";
		}else {$segregation_type{$key}= "nomap";
	}
	}


##########################################################################################################################
#	Print outfile with all genotype data for each parent and snp name and segregation type
my $date = `date`;
	$date =~ s/\s+/_/g;
	$date =~ s/\:/_/g;
open (SNP_OUT, ">>$ARGV[2]/$ARGV[2]\_parent_genotypes_segregation_types_$date.txt") || die "Can't open parent genotypes outfile $!\n";
		print SNP_OUT "$print_header";
		foreach my $x (sort keys %snpname)	{
			print SNP_OUT "$x\t$segregation_type{$x}\t$fem{$x}\t$mal{$x}\n";
			}
close SNP_OUT;

###############################################################################################################################
