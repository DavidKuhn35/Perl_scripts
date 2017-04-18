#!/usr/bin/perl
#Takes a fasta file from command line ARGV[0] and a start ARGV[1] and stop ARGV[2] position and writes out a subsequence file.  ARGV[3] lets you add another name to subsequence.
use Bio::Seq; 
use Bio::SeqIO;
 	$seqio_obj = Bio::SeqIO->new(-file => "$ARGV[0]", -format => "fasta" );
	$seq_obj = $seqio_obj->next_seq;
	$fasta_name = $seq_obj->display_id;	
	$subseq = $seq_obj->subseq($ARGV[1],$ARGV[2]);
#	print "$ARGV[0]\n$ARGV[1]\n$ARGV[2]\n $subseq\n";
#	print $seq_obj->display_id ."\n". $seq_obj->subseq($ARGV[1], $ARGV[2]) . "\n";
#	print "$subseq is the subsequence of $fasta_name from $ARGV[1] to $ARGV[2]\n";
	$seq2_obj = Bio::Seq->new(-seq => $subseq, -display_id => "$fasta_name subseq from $ARGV[1] to $ARGV[2]");
	$seqio_obj = Bio::SeqIO->new(-file => ">$ARGV[3]\_subseq\_$fasta_name\_$ARGV[1]\_$ARGV[2].fasta", -format => 'fasta' );
	$seqio_obj->write_seq($seq2_obj);

