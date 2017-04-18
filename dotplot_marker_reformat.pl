#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: February 22, 2017
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: 	takes mango markers from fasta file, sorts by lg and cm, concatenates into single sequence with 10Ns between markers and 100Ns between lgs
# <input>		:  fasta file of sequences		
# <output>		:  single sequence file
# 			
# Dependencies	: N/A
# Usage			: ARGV[0] = file with markers, ARGV[1] name of output sequence
#                    
##################################################################################
use strict;
use warnings;
use diagnostics;
use Bio::Seq; 
use Bio::SeqIO;
my %seq = ();	#key is marker name, value is sequence
my %lg = ();		#key is marker name, value is linkage group
my %cm = ();		#key is marker name, value is cM position in linkage group
my @mn = ();		#array of all marker names to be sorted
my @concat = ();
my @lgseq = ();
my $N10 = "N" x 10;
my $N490 = "N" x 490;

use Bio::SeqUtils;

######################################################################################################################
# read in markers from fasta file, make hashes with marker name as key
open (SCAFDATA, $ARGV[0]) || die "Can't open $ARGV[0]:$!\n";

 	my $seqio_obj = Bio::SeqIO->new(-file => $ARGV[0], -format => 'fasta' ); #read a fasta file for input
	while (my $seq_obj = $seqio_obj->next_seq)  { #read the seq files one at a time from the input fasta file
	my $head = $seq_obj->display_id;
	
	print "head is $head\n";	
#		(my $mn, my $desc) = split /\|/, $head;		#split marker name from description for mangomap markers

#	print "desc is $desc\n";
#		(my $lg, my $cm) = split /:/, $desc;		#split lg and cm from description, discard rest of description
		


		(my $mn, my $lg, my $cm) = split /\_/, $head;  #split marker name, lg and cm for Luomap markers 

		print " marker is $mn lg is $lg cm is $cm\n";
		my $seq = $seq_obj->seq;				#get the sequence from the fasta file
		$seq{$mn}=$seq;
		$lg{$mn}=$lg;
		$cm{$mn}=$cm;
		unless ($lg{$mn} eq "N/A")	{			#if lg is N/A, don't add to marker name array
		push (@mn, $mn);						#add next mn to array
		}
}
			
close SCAFDATA;				
######################################################################################################################
# sort the marker name array by lg, cm and marker name

my @lgcmsort = sort 	{
	$lg{$a} <=> $lg{$b}  or		#sort by lg
	$cm{$a} <=> $cm{$b}  or  	#or sort by cm
	$a cmp $b				#or sort alphabetically by marker name
}	@mn;


#######################################################################################################################
#   concatenate marker sequences based on markers sorted by lg and cm

		foreach my $x (@lgcmsort)	{  #put sequence of marker into an array in lg cm order
			my $y = $lg{$x} - 1 ;	#make index by subtracting 1 from LG number
			$concat[$y] .= $seq{$x} . $N10;	#append sorted marker sequence to array of lgs with spacer of 10 Ns 
		}
			my $result = join $N490, @concat;  #create concatenated marker sequences with 500 Ns between LGs

		open (OUTPUT, ">>$ARGV[1]") || die "Can't open output $ARGV[1] :$!\n";
			print OUTPUT "$result";
			
	close OUTPUT;	


#########################################################################################################################
#

##########################################################################################################################
#

