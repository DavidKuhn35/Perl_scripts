#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: July 25, 2012
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: reads in TAIR peptide fasta file, reads in APVO ATxG names, makes file with APVO protein sequences
# <input>		:  	
# <output>		:   
# 			
# Dependencies	: N/A
# Usage			: get_APVO_peptides.pl
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
my %length = ();
my %snp = ();
my %cov_length = ();
my %percent_cov = ();
my %med_cov = ();
my %mean_cov = ();
my %std_dev = ();
my %min_cov = ();
my %max_cov = ();
my %hoh = ();
my %hoh2 = ();
my $fastaname = ();
my $seqio_obj = ();
my $name = ();
my $seq_obj = ();
my $seq = ();
my $out_obj = ();
#####################################################################################################################
#	main	
	open (INFILE, "$ARGV[0]") || die "Can't open infile with APVO names$!\n";
		while (<INFILE>)	{
		chomp;
		if (m/^AT/)	{
			$fastaname = $_;
			$hoh{$fastaname} = $fastaname;
		}
}
	$seqio_obj = Bio::SeqIO->new(-file => $ARGV[1], -format => "fasta" );
	while ($seq_obj = $seqio_obj->next_seq)  {
	 $name = $seq_obj->display_id;
	 $seq = $seq_obj->seq;
	foreach my $key (sort keys %hoh)	{
		if ($name =~ m/$key/)	{
		$out_obj = Bio::SeqIO-> new(-file => '>>APVO_959.fasta', -format => 'fasta');
		$out_obj->write_seq($seq_obj);
 	}
}
}


