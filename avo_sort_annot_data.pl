#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: March 5, 2013
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: 
# <input>		:  avo annotation dataset (edit_Mar4_avo_all_CPFannot.txt), file of APVO isotigs, file of single copy isotigs, 		
# <output>		:  single copy, non-single copy, APVO, non-APVO files with annotation to go into annotation counter program avo_annot_cat_SNP_isolen.pl 
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
my %nosnp = ();
my %snp = ();
my $header = ();

#####################################################################################################################
#	main  split data set into snp and nosnp	
open (DATA, "$ARGV[0]") || die "Can't open $ARGV[0] $!\n";
	while (<DATA>)
	{
		chomp;
		if (m/^isotig/) {
		$header = $_;
		} elsif (m/^g/)	{
			my @x = split /\t/;
				my $z = $x[0];
				my $y = join "\t", @x;	
				if ($x[2]== 0)	{
					$nosnp{$z} = $y;
				} else {
					$snp{$z} = $y;
				}
			}
		}
	
							
					
######################################################################################################################
#   print out snp and nosnp files to be entered into counting program
open (SNP,">>SNP_$ARGV[0]") || die "Can't open snp out $!\n";
	print SNP "$header \n";
	foreach my $key (sort keys %snp)	{
		print SNP "$snp{$key} \n";
	}
close SNP;
open (NOSNP, ">>NOSNP_$ARGV[0]") || die "Can't open nosnp out $!\n";
	print NOSNP "$header\n";
	foreach my $key (sort keys %nosnp)	{
		print NOSNP "$nosnp{$key}\n";
		}


#######################################################################################################################
#   

#########################################################################################################################
#

##########################################################################################################################
#

