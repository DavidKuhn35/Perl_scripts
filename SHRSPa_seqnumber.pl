#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: July 24, 2014
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: read in Mar7_BHST_6kSNP_chip_seq_data file get seq SNP name and SHRSPa name
# <input>		:  Mar7_BHST_6kSNP_chip_seq_data		
# <output>		:  SHRSPachipname_seqname 
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
my $SHRS = 0;
my $seqno = 0;

#####################################################################################################################
#	main	
open (CHIP, "$ARGV[0]") || die "can't open $ARGV[0] $!\n";
	while (<CHIP>)	{
		chomp;
		my @x = split (/\t/);
		 $SHRS = $x[40];
	#	print "$SHRS\n";
		 $seqno = $x[1];	
	#	print "$seqno\n";
open (SCOUT, (">>SHRS_seqno_out_July24_2014.txt") || die "can't open outfile $! \n");
	print SCOUT "$SHRS\t$seqno\n";
}
		
				
######################################################################################################################
#   

#######################################################################################################################
#   

#########################################################################################################################
#

##########################################################################################################################
#

