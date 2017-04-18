#!/usr/bin/perl
#Takes a fasta file with multiple sequences from command line ARGV[0] and takes the bracket format for SNPs into two separate sequences for translation and makes a hoh with the new sequences linked to a single SNP name.  To be used as a subroutine.
%hoh = ();


use Bio::Seq; 
use Bio::SeqIO;
 	$seqio_obj = Bio::SeqIO->new(-file => $ARGV[0], -format => "fasta" );
	$fastaoutput_obj = Bio::SeqIO->new (-file => '>separateseq.fasta', -format => 'fasta' );
	$pepoutput_obj = Bio::SeqIO->new (-file => '>separatepep.fasta', -format => 'fasta' );
	open (OUTFILE, ">>translate_bracket_$ARGV[0]") || die "Can't open outfile\n";
	while ($seq_obj = $seqio_obj->next_seq)  {
 	$fasta_name = $seq_obj->display_id;
	$seq = $seq_obj->seq;
	if ($seq =~ m!\[A/G\]!) {
	&translate_both ($seq,A,G,$fasta_name);
		}elsif ($seq =~ m!\[C/G\]!)  {
		&translate_both ($seq, C, G, $fasta_name);
	}elsif ($seq =~ m!\[T/G\]!)  {
		&translate_both ($seq, T, G, $fasta_name);
	}elsif ($seq =~ m!\[A/C\]!)  {
		&translate_both ($seq, A, C, $fasta_name);		
	}elsif ($seq =~ m!\[T/C\]!)  {
		&translate_both ($seq, T, C, $fasta_name);
	}elsif ($seq =~ m!\[A/T\]!) {
		&translate_both ($seq, A, T, $fasta_name);
	}	
}
	for $key (sort keys %hoh)	{
		print OUTFILE "$key\tOriginal bracket sequence is $hoh{$key}{bracketseq}\nSequence 1 is $hoh{$key}{seq1}\nSequence 2 is $hoh{$key}{seq2}\n";
	}
	
sub translate_both {
	my ($a,$b,$c,$d)=@_;
	$hoh{$d}{bracketseq}= $a;
	$seq1 = $a;
	$seq1 =~ s!\[$b/$c\]!$b!;
	$hoh{$d}{seq1} = $seq1;
	$seq2 = $a;
	$seq2 =~ s!\[$b/$c\]!$c!;
	$hoh{$d}{seq2} = $seq2;
}
