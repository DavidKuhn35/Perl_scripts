#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: February 24,2017
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: take LG and psMol data and reformat parsed blast report for plotting LG vs psMol
# <input>		: parsed blast report		
# <output>		: file suitable for plotting LG cM vs position on psmol
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
my @parse = ();

#####################################################################################################################
#	main	

	open (INFILE, "$ARGV[0]") || die "Can't open infile $!\n";
				while (<INFILE>)	{
					chomp;
					if (m/^Marker/)	{		#don't take header, only markers
					@parse = split /\s+/;		# split on whitespace
					(my $psmol) = split /\_/ , $parse[5];	#split psmol name off length
					my $psnew = $psmol * (10**8);
					my $pshit = $psnew + $parse[18];	#put psmol onto hit start
					print "$parse[18]\n";
					print "$pshit\n";
					my $cmnew = $parse[2] + ($parse[1] * (10**3));  #put lg onto cm
						open (OUTFILE, ">>$ARGV[1]") || die "Can't open outfile $!\n";
						print OUTFILE "$cmnew,$pshit\n";
					}
					
				}
######################################################################################################################
#######################################################################################################################
#  subroutine 

#########################################################################################################################
#
 
##########################################################################################################################
#
  
