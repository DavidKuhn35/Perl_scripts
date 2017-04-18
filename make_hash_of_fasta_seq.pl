#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: November 21, 2014
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: identifies a subset of sequences from a fasta file based on fasta id list
# <input>		:  fasta file and list of fasta names		
# <output>		:  file with snpname and sequence flanking SNP
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
my %fhash = ();
my %snpname =();


#####################################################################################################################
#	main	
#Takes a fasta file with multiple sequences from command line ARGV[0] and makes a hash of a subset of sequences based on a text list ARGV[1]


 	my $seqio_obj = Bio::SeqIO->new(-file => $ARGV[0], -format => "fasta" );
	while (my $seq_obj = $seqio_obj->next_seq)  {
 	my $fasta_name = $seq_obj->display_id;
	my $fasta_seq = $seq_obj-> seq;
	$fhash{$fasta_name} = $fasta_seq;
	}
	
open (SNP, "$ARGV[1]") || die "Can't open SNP name infile";
	while (<SNP>)
{
	chomp;
	if (m/^J/)  	{	#If line has J at start
	my @snp = split (/\t/);
	$snpname{$snp[0]} = $snp[0];	#make hash of snpnames	
	
}
}			

open (SNPOUT, ">>$ARGV[2].txt") || die "can't open snp sequence outfile $!\n";
	foreach my $key (sort keys %fhash)	{
		if (exists $snpname{$key})	{	
			print SNPOUT "$key\t$fhash{$key}\n";
}
}				
######################################################################################################################
#   

#######################################################################################################################
#   

#########################################################################################################################
#

##########################################################################################################################
#

