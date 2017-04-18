#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: March 7, 2013
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: Reformat isotig dataset to make parametric statistical analysis possible
# <input>		: edit_Mar4_avo_all_CPFannot.txt 		
# <output>		:   avo_all_param_reformat.txt
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
my %isolen = ();
my %isosnp = ();
my %isoC = ();
my %isoF = ();
my %isoP = ();
my %Cterm = ();
my %Fterm = ();
my %Pterm = ();
my %Choh = ();
my %Fhoh = ();
my %Phoh = ();
my %reformat = ();
my %apvo = ();
my %singlecopy = ();
my $header = "isotig\tlength\tSNPcount\tCterm\tFterm\tPterm\tSNPcat\tAPVOcat\tSinglecopycat";
#####################################################################################################################
#	main	read in datafile
open (APVO, "$ARGV[1]") || die "Can't open APVO input file $! \n";
	while (<APVO>)
{
	chomp;
	if (m/^g/)	{
		my @z = split /\t/;
		$apvo{$z[0]} = $z[0];
		}
}
close APVO;
open (SINGLECOPY, "$ARGV[2]") || die "Can't open singlecopy input file $! \n";
	while (<SINGLECOPY>)	
{
	chomp;
	if (m/^g/)	{
		my @z = split /\t/;
		$singlecopy{$z[0]} = $z[0];
		}
}
close SINGLECOPY;
	
open (DATA, "$ARGV[0]") || die "Can't open input data file $! \n";
	while (<DATA>)
{
	chomp;
	if (m/^g/)	{	#if line starts with an isotig name
	my @x = split (/\t/);
	$isolen{$x[0]} = $x[1];	# hash with isotig name as key and isotig length as value
	$isosnp{$x[0]} = $x[2]; # hash with isotig name as key and number of SNPs in isotig as value
	$isoC{$x[0]} = $x[7]; # hash with isotig name as key and C annotation as value
	$isoF{$x[0]} = $x[13]; #hash with isotig name as key and F annotation as value
	$isoP{$x[0]} = $x[19]; #hash with isotig name as key and P annotation as value
	my @y = ( "$x[0]", "$x[1]", "$x[2]", "$x[7]", "$x[13]", "$x[19]");  #make new array with isotig name, length, SNP number, C, F, P annotation
		if ($y[2] == 0)		{
			push (@y, "nosnp");
			}else{
			push (@y, "yessnp");
		}
		if (exists $singlecopy{$y[0]})		{
			push (@y, "yessc");
			}else{
			push (@y, "nosc");
		}
		if (exists $apvo{$y[0]})	{
			push (@y, "yesapvo");
			}else{
			push (@y, "noapvo");
		}
	$reformat{$y[0]} = join "\t", @y;
	}
}				
close DATA;		
open (OUT, ">>Mar7_avo_reformat_out.txt") || die "Can't open outfile $! \n";
		print OUT "$header\n";
	foreach my $key (sort keys %reformat)	{
		print OUT "$reformat{$key}\n";
		}
close OUT;			
				
######################################################################################################################
#   

#######################################################################################################################
#   

#########################################################################################################################
#

##########################################################################################################################
#

