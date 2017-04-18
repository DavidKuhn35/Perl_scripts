#!/usr/bin/perl
#Takes a fasta file with multiple sequences from command line ARGV[0] and makes it into separate files for each sequence.  Use it to split genome fasta files up into pseudomolecules.  
#ARGV[1] can add specificity to names of each individual fasta file.
use Bio::Seq; 
use Bio::SeqIO;
 	$seqio_obj = Bio::SeqIO->new(-file => $ARGV[0], -format => "fasta" );
	while ($seq_obj = $seqio_obj->next_seq)  {
 	$fasta_name = $seq_obj->display_id;
#	print "I made it this far.The fasta id is $fasta_name. \n";
	$seqioa_obj = Bio::SeqIO->new(-file => ">$ARGV[1]\_$fasta_name.fasta", -format => 'fasta' );
	$seqioa_obj->write_seq($seq_obj);
}
