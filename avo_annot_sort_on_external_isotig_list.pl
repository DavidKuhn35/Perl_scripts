#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: March 5, 2013
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: takes annotated isotigs and sorts into apvo, non-apvo or singlecopy, non-singlecopy depending on file loaded in to program
# <input>		:  edit_Mar4_avo_all_CPFannot.txt, files with names of isotigs that are apvo, singlecopy, whatever		
# <output>		:   files that have apvo, nonapvo isotigs with annotation to be run through avo_annot_cat_SNP_isolen.pl
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
my $header = ();
my %annot = ();
my %sort = ();
my %yes = ();
my %no = ();
#####################################################################################################################
#	main	read in annotation data and make a hash with isotigs as keys and annotation data as values
open (DATA, "$ARGV[0]") || die "Can't open annotation data for all avo isotigs $!\n";
	while (<DATA>)	
{
	chomp;	
		if (m/^isotig/)		{
		$header = $_;
		}elsif (m/^g/){
			my $x = $_;
			my @y = split /\t/;
			$annot{$y[0]} = $x;
		}
	}	
				
				
######################################################################################################################
#  read in file with isotig names to sort annotation data
open (SORT, "$ARGV[1]") || die "Can't open sort file with isotig names $!\n";
	while (<SORT>)
{
	chomp;
	if (m/^g/)	{
		my @x = split /\t/;
		$sort{$x[0]}= $x[0];
	}
}
		 

#######################################################################################################################
#   go through annotation hash and sort into two files based on isotig names in sort hash
foreach my $key (sort keys %annot)	{
	if (exists $sort{$key})	{
		$yes{$key} = $annot{$key};
	}else{
	$no{$key} = $annot{$key};
	}
}

#########################################################################################################################
#  print out files with annotation
open (YES, ">>YES_$ARGV[1]") || die "Can't open yes output file $!\n";
	print YES "$header\n";
	foreach my $key (sort keys %yes)	{
		print YES "$yes{$key}\n";
	}
close YES;

open (NO, ">>NO_$ARGV[1]") || die "Can't open no output file $!\n";
	print NO "$header\n";
	foreach my $key (sort keys %no)	{
		print NO "$no{$key}\n";
	}

##########################################################################################################################
#

