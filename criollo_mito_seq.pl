#!/usr/bin/perl
#Takes a fasta file from command line ARGV[0] and a start ARGV[1] and stop ARGV[2] position and writes out a subsequence file.  ARGV[3] lets you add another name to subsequence.
use strict;
use warnings;
use diagnostics;

use Bio::Seq; 
use Bio::SeqIO;
 	my $seqio_obj = Bio::SeqIO->new(-file => "$ARGV[0]", -format => "fasta" );	#Read in Criollo linkage group 5 which has most of mito genome!
	my $Cri05 = $seqio_obj->next_seq;
	my $fasta_name = $Cri05->display_id;
#read in start and stop positions from parsed BLAST file between Matina mito genome and Criollo genome	

open (INFILE, "$ARGV[1]") || die "Can't open infile with parsed BLAST data $!\n";	#Mat_mito_vs_Criollo_blastn_top_106.txt
	while (<INFILE>)	
{
		chomp;
		if (m/^M16/)	{
			my @blast = split (/\t/);
			if ($blast[8]eq "Plus/Plus")	{
			my $start = $blast[11];
			my $stop = $blast[12];
			my $Mstart = $blast[9];
			my $Mstop = $blast[10];
			my $align_length = $blast[7];
			my $subseq = $Cri05->subseq($start,$stop);
			my $seq2_obj = Bio::Seq->new(-seq => $subseq, -display_id => "Criollo_subseq_from_$start\_to_$stop\_matches_Matina_mito_from_$Mstart\_to_$Mstop\_alignment_length_is\_$align_length");
			$seqio_obj = Bio::SeqIO->new(-file => ">Matina_$Mstart\_$Mstop.fasta", -format => 'fasta' );
			$seqio_obj->write_seq($seq2_obj);

			} elsif ($blast[8] eq "Plus/Miuns")	{
			my $start = $blast[12];
			my $stop  = $blast[11]; 
			my $Mstart = $blast[9];
			my $Mstop = $blast[10];
			my $align_length = $blast[7];
			my $subseq = $Cri05->subseq($start,$stop);
			my $seq2_obj = Bio::Seq->new(-seq => $subseq, -display_id => "Criollo_subseq_from_$start\_to_$stop\_matches_Matina_mito_from_$Mstart\_to_$Mstop\_alignment_length_is\_$align_length");
			my $seq2rev_obj = $seq2_obj->revcom;
			$seqio_obj = Bio::SeqIO->new(-file => ">Matina_$Mstart\_$Mstop.fasta", -format => 'fasta' );
			$seqio_obj->write_seq($seq2rev_obj);
}
}
}

