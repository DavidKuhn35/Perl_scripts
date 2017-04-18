#!/usr/bin/perl
#Takes a fasta file with multiple sequences from command line ARGV[0] and takes the bracket format for SNPs into two separate sequences for translation and makes a new fasta file with both sequences.
use Bio::Seq; 
use Bio::SeqIO;
 	$seqio_obj = Bio::SeqIO->new(-file => $ARGV[0], -format => "fasta" );
	$fastaoutput_obj = Bio::SeqIO->new (-file => '>separateseq.fasta', -format => 'fasta' );
	$pepoutput_obj = Bio::SeqIO->new (-file => '>separatepep.fasta', -format => 'fasta' );
	open (OUTFILE, ">>translate_outfile_$ARGV[0]") || die "Can't open outfile\n";
	while ($seq_obj = $seqio_obj->next_seq)  {
 	$fasta_name = $seq_obj->display_id;
	$seq = $seq_obj->seq;
	print "I made it this far.The fasta id is $fasta_name. \n";
	print "The sequence is $seq.\n";
	if ($seq =~ m!\[A/G\]!) {
		$seq1 = $seq;
		$seq2 = $seq;
		$seq1 =~ s!\[A/G\]!A!;
		$seq1_obj = Bio::Seq-> new(-seq => $seq1, display_id => "$fasta_name\_1", -alphabet => "dna");
		$seq1pep_obj = $seq1_obj->translate;
		$fastaoutput_obj->write_seq($seq1_obj);
		$pepoutput_obj->write_seq($seq1pep_obj);
		$seq2 =~ s!\[A/G\]!G!;
		$seq2_obj = Bio::Seq-> new(-seq => $seq2, display_id => "$fasta_name\_2", -alphabet => "dna");
		$seq2pep_obj = $seq2_obj->translate;
		$fastaoutput_obj->write_seq($seq2_obj);
		$pepoutput_obj->write_seq($seq2pep_obj);
	print "The first sequence is $seq1.\n  The second sequence is $seq2.\n";
	print OUTFILE ">$fasta_name\_1\n$seq1\n>$fasta_name\_2\n$seq2\n";	
	}elsif ($seq =~ m!\[C/G\]!)  {
		$seq1 = $seq;
		$seq2 = $seq;
		$seq1 =~ s!\[C/G\]!C!;
		$seq1_obj = Bio::Seq-> new(-seq => $seq1, display_id => "$fasta_name\_1", -alphabet => "dna");
		$seq1pep_obj = $seq1_obj->translate;
		$fastaoutput_obj->write_seq($seq1_obj);
		$pepoutput_obj->write_seq($seq1pep_obj);
		$seq2 =~ s!\[C/G\]!G!;
		$seq2_obj = Bio::Seq-> new (-seq => $seq2, display_id => "$fasta_name\_2", -alphabet => "dna");
		$seq2pep_obj = $seq2_obj->translate;
		$fastaoutput_obj->write_seq($seq2_obj);
		$pepoutput_obj->write_seq($seq2pep_obj);
	print "The first sequence is $seq1.\n  The second sequence is $seq2.\n";
	print OUTFILE ">$fasta_name\_1\n$seq1\n>$fasta_name\_2\n$seq2\n";	
	}elsif ($seq =~ m!\[T/G\]!)  {
		$seq1 = $seq;
		$seq2 = $seq;
		$seq1 =~ s!\[T/G\]!T!;
		$seq1_obj = Bio::Seq-> new(-seq => $seq1, display_id => "$fasta_name\_1", -alphabet => "dna");
		$seq1pep_obj = $seq1_obj->translate;
		$fastaoutput_obj->write_seq($seq1_obj);
		$pepoutput_obj->write_seq($seq1pep_obj);
		$seq2 =~ s!\[T/G\]!G!;
		$seq2_obj = Bio::Seq-> new (-seq => $seq2, display_id => "$fasta_name\_2", -alphabet => "dna");
		$seq2pep_obj = $seq2_obj->translate;
		$fastaoutput_obj->write_seq($seq2_obj);
		$pepoutput_obj->write_seq($seq2pep_obj);
	print "The first sequence is $seq1.\n  The second sequence is $seq2.\n";
	print OUTFILE ">$fasta_name\_1\n$seq1\n>$fasta_name\_2\n$seq2\n";	
	}elsif ($seq =~ m!\[A/C\]!)  {
		$seq1 = $seq;
		$seq2 = $seq;
		$seq1 =~ s!\[A/C\]!A!;
		$seq1_obj = Bio::Seq-> new(-seq => $seq1, display_id => "$fasta_name\_1", -alphabet => "dna");
		$seq1pep_obj = $seq1_obj->translate;
		$fastaoutput_obj->write_seq($seq1_obj);
		$pepoutput_obj->write_seq($seq1pep_obj);
		$seq2 =~ s!\[A/C\]!C!;
		$seq2_obj = Bio::Seq-> new (-seq => $seq2, display_id => "$fasta_name\_2", -alphabet => "dna");
		$seq2pep_obj = $seq2_obj->translate;
		$fastaoutput_obj->write_seq($seq2_obj);
		$pepoutput_obj->write_seq($seq2pep_obj);
	print "The first sequence is $seq1.\n  The second sequence is $seq2.\n";
	print OUTFILE ">$fasta_name\_1\n$seq1\n>$fasta_name\_2\n$seq2\n";	
	}elsif ($seq =~ m!\[T/C\]!)  {
		$seq1 = $seq;
		$seq2 = $seq;
		$seq1 =~ s!\[T/C\]!T!;
		$seq1_obj = Bio::Seq-> new(-seq => $seq1, display_id => "$fasta_name\_1", -alphabet => "dna");
		$seq1pep_obj = $seq1_obj->translate;
		$fastaoutput_obj->write_seq($seq1_obj);
		$pepoutput_obj->write_seq($seq1pep_obj);
		$seq2 =~ s!\[T/C\]!C!;
		$seq2_obj = Bio::Seq-> new (-seq => $seq2, display_id => "$fasta_name\_2", -alphabet => "dna");
		$seq2pep_obj = $seq2_obj->translate;
		$fastaoutput_obj->write_seq($seq2_obj);
		$pepoutput_obj->write_seq($seq2pep_obj);
	print "The first sequence is $seq1.\n  The second sequence is $seq2.\n";
	print OUTFILE ">$fasta_name\_1\n$seq1\n>$fasta_name\_2\n$seq2\n";	
	}elsif ($seq =~ m!\[A/T\]!) {
		$seq1 = $seq;
		$seq2 = $seq;
		$seq1 =~ s!\[A/T\]!A!;
		$seq1_obj = Bio::Seq-> new(-seq => $seq1, display_id => "$fasta_name\_1", -alphabet => "dna");
		$seq1pep_obj = $seq1_obj->translate;
		$fastaoutput_obj->write_seq($seq1_obj);
		$pepoutput_obj->write_seq($seq1pep_obj);
		$seq2 =~ s!\[A/T\]!T!;
		$seq2_obj = Bio::Seq-> new (-seq => $seq2, display_id => "$fasta_name\_2", -alphabet => "dna");
		$seq2pep_obj = $seq2_obj->translate;
		$fastaoutput_obj->write_seq($seq2_obj);
		$pepoutput_obj->write_seq($seq2pep_obj);
	print "The first sequence is $seq1.\n  The second sequence is $seq2.\n";
	print OUTFILE ">$fasta_name\_1\n$seq1\n>$fasta_name\_2\n$seq2\n";	
	}	
	else	{
	print OUTFILE ">$fasta_name\nATCG\n";
	}
}

