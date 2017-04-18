#!/usr/bin/perl
#Takes a fasta file from command line ARGV[0] and a start ARGV[1] and stop ARGV[2] position and writes out a reverse complemented subsequence file.  ARGV[3] gives cultivar name.
use Bio::Seq; 
use Bio::SeqIO;
 	$seqio_obj = Bio::SeqIO->new(-file => "$ARGV[0]", -format => "fasta" );
	$seq_obj = $seqio_obj->next_seq;
	$fasta_name = $seq_obj->display_id;	
	$subseq = $seq_obj->subseq($ARGV[1],$ARGV[2]);
#	print "$ARGV[0]\n$ARGV[1]\n$ARGV[2]\n $subseq\n";
#	print $seq_obj->display_id ."\n". $seq_obj->subseq($ARGV[1], $ARGV[2]) . "\n";
#	print "$subseq is the subsequence of $fasta_name from $ARGV[1] to $ARGV[2]\n";
	$seq2_obj = Bio::Seq->new(-seq => $subseq, -display_id => "$ARGV[0]\_$fasta_name\_revcom\_subseq\_from\_$ARGV[1]\_$ARGV[2]");
	$fasta2_name = $seq2_obj->display_id;
	$seq3_obj = $seq2_obj->revcom;
	$seqio_obj = Bio::SeqIO->new(-file => ">$ARGV[3]\_$fasta2_name.fasta", -format => 'fasta' );
	$seqio_obj->write_seq($seq3_obj);

