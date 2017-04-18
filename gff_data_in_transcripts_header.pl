#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: January 26, 2017
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: 	takes data from a csv file and puts it in the description of the mango pseudomolecules fasta file header
# <input>		:  fasta file of sequences		
# <output>		:  fasta file of sequences with length in header
# 			
# Dependencies	: N/A
# Usage			: ARGV[0] = file with scaffold info for pseudomolecules, ARGV[1] is pseudomolecule fasta file
#                    
##################################################################################
use strict;
use warnings;
use diagnostics;
use Bio::Seq; 
use Bio::SeqIO;
my %scafhead = ();
my $x = ();
my $psmol = ();
my $y = ();
use Bio::SeqUtils;
#use Bio::Tools::OddCodes;
#use Bio::PrimarySeq;
######################################################################################################################
# get scaffold data from csv file

#get scaffold data from csv file
open (SCAFDATA, $ARGV[0]) || die "Can't open $ARGV[0]:$!\n";
#	my $scafdata = $ARGV[0];
	while (<SCAFDATA>) {
		chomp;
		my @scaf = split /\s/;
		$psmol = $scaf[0];
		my $scafstart = $scaf[3];
		my $scafend = $scaf[4];
		my $range = join "-", $scafstart, $scafend;
		(my $id, my $name) = split ";", $scaf[8];
#		$x = join ":", $psmol, $range;
#		$y = $x . "\|";
		$scafhead{$id} = $psmol . ":" . $range;
		} 		

#open (SCAFDATA, $ARGV[0]) || die "Can't open $ARGV[0]:$!\n";
#	my $scafdata = $ARGV[0];
#	while (<SCAFDATA>) {
#		chomp;
#		my @scaf = split /\t/;
#		$psmol = $scaf[0];
#		my $scafno = $scaf[1];
#		my $scafstart = $scaf[2];
#		my $scafend = $scaf[3];
#		my $range = join "-", $scafstart, $scafend;
#		
#		$x = join ":", $scafno, $range;
#		$y = $x . "\|";
#		$scafhead{$psmol}.= $y;
#		} 		



#####################################################################################################################
#	main	get sequence names from fasta file
	my $newhead = Bio::SeqIO->new(-file => ">$ARGV[1]\_psmol.start.stop.fasta", -format => 'fasta'); #create a new fasta file for output
 	my $seqio_obj = Bio::SeqIO->new(-file => $ARGV[1], -format => 'fasta' ); #read a fasta file for input
	while (my $seq_obj = $seqio_obj->next_seq)  { #read the seq files one at a time from the input fasta file
 #	$fasta_name = $seq_obj->display_id;
#	$seqname{$fasta_name} = $fasta_name;
#	my $len = $seq_obj->length;	#get the length of the sequence
	my $di = $seq_obj->display_id; #get the display id of the sequence
	if (exists $scafhead{$di}) {
		my $olddesc = $seq_obj->desc;
		my $newdesc = $scafhead{$di} . "\|" . $olddesc; #get the new header info from file read in above
#	my $newdi = "$di\_$len\_nts"; #put the length into the display id
#	$seq_obj->display_id("$newdi"); #set the new display id in the original seq file
		$seq_obj->desc("$newdesc"); #set the description to include the length, does not change display_id
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

