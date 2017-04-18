#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: July 3, 2013
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: Reads in parsed BLASTx report $ARGV[0] to get query frame and strand.  Uses this to correctly translate the full length query sequence and generate a peptide fasta #	#	file.
# <input>		: parsed BLASTx report of avo isotigs vs cacao gene models $ARGV[0].  avo isotig fasta file $ARGV[1] name of peptide outfile $ARGV[2]		
# <output>		:   avo peptide fasta file (avo_transcriptome_gene_models)
# 			
# Dependencies	: N/A
# Usage			: avo_transcriptome_gene_models.pl
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
my %isostrand = ();
my %isoframe = ();

#####################################################################################################################
#	main	
open (INFILE, "$ARGV[0]") || die "Can't open infile.\n"; #read in parsed BLASTx report of avo isotigs vs cacao gene models
	while (<INFILE>) #read in from file, line by line
{
	chomp;	#remove end of line character
	if (m/^g/) {	#skip any header information, start at genotype data, avo isotigs start with g
		my @iso = split (/\t/); #assign variables from genotype data
		my $isoname = shift(@iso);  #take the isotig name from the left end of the array
		$isostrand{$isoname} = pop(@iso);  #take the query strand from the right end of the array
		$isoframe{$isoname} = pop(@iso);  #take the query frame from the right end of the array
		}
}			
close INFILE;				
######################################################################################################################
# Read in avocado isotig fasta file
my $avo_iso = Bio::SeqIO->new(-file => $ARGV[1], -format => "fasta" );  #read in avocado isotig file and make sequence objects
my $avo_pep_out = Bio::SeqIO->new (-file => ">$ARGV[2].fasta", -format => 'fasta' ); #prepare output peptide fasta file
	while (my $avo_seq = $avo_iso->next_seq)  {
 	my $fasta_name = $avo_seq->display_id;
		my @iso_name = split /\s/, $fasta_name; #split the fasta_name on white space to just get the isotig name as the first array element
		my $g_name = $iso_name[0];
	if (exists $isostrand{$g_name})	{
		if ($isostrand{$g_name}>0)	{
		my @seqs = Bio::SeqUtils->translate_3frames($avo_seq); #translate three forward frames
		my $avo_pep = $seqs[$isoframe{$g_name}]; #gets the translation in the correct frame
		$avo_pep_out->write_seq ($avo_pep);
		}else{
		my $rev_avo_seq = $avo_seq->revcom;	#reverse and complement the sequence as it is on the negative strand
		my @rev_seqs = Bio::SeqUtils->translate_3frames($rev_avo_seq); #translate in the three forward frames the reversed and complemented seq
		my $rev_avo_pep = $rev_seqs[$isoframe{$g_name}];	#gets the translation in the correct frame
		$avo_pep_out->write_seq ($rev_avo_pep);
		}
	}	
}
#######################################################################################################################
#   

#########################################################################################################################
#

##########################################################################################################################
#

