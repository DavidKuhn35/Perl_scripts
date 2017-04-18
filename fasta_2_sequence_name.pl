#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: March 6, 2013
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: 	takes a fasta file and makes a hash of sequence names
# <input>		:  fasta file of sequences		
# <output>		:  file of sequence names 
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
my %seqname = ();
my $seqio_obj = ();
my $seq_obj = ();
my $fasta_name = ();
#####################################################################################################################
#	main	get sequence names from fasta file

 	$seqio_obj = Bio::SeqIO->new(-file => $ARGV[0], -format => "fasta" );
	while ($seq_obj = $seqio_obj->next_seq)  {
 	$fasta_name = $seq_obj->display_id;
	$seqname{$fasta_name} = $fasta_name;
}			
				
######################################################################################################################
# print out file with sequence names
open (SEQOUT, ">>seqname_out_$ARGV[0]") || die "Can't open seqname out file $!\n";
	foreach my $key (sort keys %seqname)	{
	  print SEQOUT "$key\n";
	}

#######################################################################################################################
#   

#########################################################################################################################
#

##########################################################################################################################
#

