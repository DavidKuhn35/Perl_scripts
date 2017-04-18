#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: March 4, 2013
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: 	does global substitution to change annotation names in avocado annotation files
# <input>		:  Mar4_avo_all_CPFannot file with all isotigs and their annotation		
# <output>		:   same data but Cgoslim_term, Fgoslim_term, Pgoslim_term have been renamed Cgoslimterm, Fgoslimterm, Pgoslimterm
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
#	main	
	open (DATA,"$ARGV[0]") || die "Can't open $ARGV[0] $!\n";
		while (<DATA>)
	{
	chomp;
	my $x = $_;
	$x =~ s/goslim_term/goslimterm/g;
	open (OUT, ">>edit_$ARGV[0]") || die "Can't open $ARGV[0]\_edit.txt \n";
		print OUT "$x\n";
	}		
				
######################################################################################################################
#   

#######################################################################################################################
#   

#########################################################################################################################
#

##########################################################################################################################
#

