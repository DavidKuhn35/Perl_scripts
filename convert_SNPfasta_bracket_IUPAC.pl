#!/usr/bin/perl
#Takes a fasta file with multiple sequences from command line ARGV[0] converts the bracket format for SNPs into the IUPAC format and prints out the reformatted SNP file.
use Bio::Seq; 
use Bio::SeqIO;
 	$seqio_obj = Bio::SeqIO->new(-file => $ARGV[0], -format => "fasta" );
	open (OUTFILE, ">>IUPAC_outfile_$ARGV[0]") || die "Can't open outfile\n";
	while ($seq_obj = $seqio_obj->next_seq)  {
 	$fasta_name = $seq_obj->display_id;
	$seq = $seq_obj->seq;
	print "I made it this far.The fasta id is $fasta_name. \n";
	print "The sequence is $seq.\n";
	if (($seq =~ m!\[A/G\]!) || ($seq =~ m!\[G/A\]!)) {
		$seq =~	s!\[\w/\w\]!R! ;
	print "The converted sequence is $seq.\n";

	print OUTFILE ">$fasta_name\n$seq\n";	
	}elsif (($seq =~ m!\[C/G\]!)|| ($seq =~ m!\[G/C\]!)) {
		$seq =~ s!\[\w/\w\]!S!;
#	open (OUTFILE, ">>IUPAC_outfile$ARGV[0].fasta") || die "Can't open outfile\n";
	print OUTFILE ">$fasta_name\n$seq\n";
	print "The converted sequence is $seq.\n";
	}elsif (($seq =~ m!\[T/G\]!) || ($seq =~ m!\[G/T\]!)) {
		$seq =~ s!\[\w/\w\]!K!;
#	open (OUTFILE, ">>IUPAC_outfile$ARGV[0].fasta") || die "Can't open outfile\n";
	print OUTFILE ">$fasta_name\n$seq\n";
	print "The converted sequence is $seq.\n";
	}elsif (($seq =~ m!\[A/C\]!) || ($seq =~ m!\[C/A\]!)) {
		$seq =~ s!\[\w/\w\]!M!;
#	open (OUTFILE, ">>IUPAC_outfile$ARGV[0].fasta") || die "Can't open outfile\n";
	print OUTFILE ">$fasta_name\n$seq\n";
	print "The converted sequence is $seq.\n";
	}elsif (($seq =~ m!\[T/C\]!)|| ($seq =~ m!\[C/T\]!))  {
		$seq =~ s!\[\w/\w\]!Y!;
#	open (OUTFILE, ">>IUPAC_outfile$ARGV[0].fasta") || die "Can't open outfile\n";
	print OUTFILE ">$fasta_name\n$seq\n";
	print "The converted sequence is $seq.\n";
	}elsif (($seq =~ m!\[A/T\]!)|| ($seq =~ m!\[T/A\]!)) {
		$seq =~ s!\[\w/\w\]!W!;
	print "The converted sequence is $seq.\n";
#	open (OUTFILE, ">>IUPAC_outfile_$ARGV[0]") || die "Can't open outfile\n";
	print OUTFILE ">$fasta_name\n$seq\n";
	}	
#	else	{
#	open (OUTFILE, ">>IUPAC_outfile_".$ARGV[0].".fasta") || die "Can't open outfile\n";
#	print OUTFILE ">$fasta_name\nATCG\n";
#	$IUPAC_obj = Bio::Seq->new(-seq => $seq, -display_id => $fasta_name, -alphabet => "dna");
#	$IUPACseq = $IUPAC_obj->seq;
#	print "The IUPAC sequence is $IUPACseq.\n";
#	$seqioa_obj = Bio::SeqIO->new(-file => ">>$ARGV[0]\_IUPAC.fasta", -format => 'fasta' );
#	$seqioa_obj->write_seq($IUPAC_obj);
#	}
}

