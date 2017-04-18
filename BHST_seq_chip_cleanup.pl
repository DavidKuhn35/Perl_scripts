#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: March 8, 2013 (Happy Birthday to me!!!)
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: 	Cleans up BHST_6kSNP_chip_seq.txt file and adds actual genotype calls from sequence data and chip data
# <input>		:  Sequence calls from BHST_6kSNP_chip_seq.txt file, sequence calls from chip for all 6000 SNPs.		
# <output>		:  Data file with all sequence and chip data for 6000 SNPs 
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
#	main	Read in BHST_6kSNP_chip_seq.txt file and make a hash based on SNP name (SHRSPa) not contig name and position
open (BHSTIN, "$ARGV[0]") || die "Can't open bhst in file $!\n";
	while (<BHSTIN>)		
{
		chomp;
		my @x = $_;
		unshift (@x, $				
######################################################################################################################
#   

#######################################################################################################################
#   

#########################################################################################################################
#

##########################################################################################################################
#

