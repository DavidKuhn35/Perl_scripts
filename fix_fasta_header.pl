#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: November 7, 2016
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: 	takes a fasta file and extracts length info and puts in header
# <input>		:  fasta file of sequences		
# <output>		:  fasta file of sequences with length in header
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
#####################################################################################################################
#	main	get sequence names from fasta file
	my $newhead = Bio::SeqIO->new(-file => '>mapsnp.fasta', -format => 'fasta'); #create a new fasta file for output
 	my $seqio_obj = Bio::SeqIO->new(-file => $ARGV[0], -format => 'fasta' ); #read a fasta file for input
	while (my $seq_obj = $seqio_obj->next_seq)  { #read the seq files one at a time from the input fasta file
		my $fasta_name = $seq_obj->display_id;	
		
		unless ($fasta_name =~ /Tcm/) {
		
		
	$fasta_name = $seq_obj->display_id;
#	$seqname{$fasta_name} = $fasta_name;
#	my $len = $seq_obj->length;	#get the length of the sequence
#	my $di = $seq_obj->display_id; #get the display id of the sequence
#	my $newdi = "$di\_$len\_nts"; #put the length into the display id
	#$seq_obj->display_id("$newdi"); #set the new display id in the original seq file
#	$seq_obj->desc("$newdi"); #set the description to include the length, does not change display_id
	$newhead->write_seq($seq_obj); #write the new header and sequence to the output fasta file
	}
}
			
				
######################################################################################################################
# print out file with sequence names
#open (SEQOUT, ">>seqname_out_$ARGV[0]" || die "Can't open seqname out file $!\n";
#	foreach my $key (sort keys %seqname)	{
#	  print SEQOUT "$key\n";
#	}

#######################################################################################################################
#   

#########################################################################################################################
#

##########################################################################################################################
#

