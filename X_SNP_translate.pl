#!/usr/bin/perl
#Takes a parsed blastx file with the name of a 121mer that had an X in its translation and the translation frame and a fasta file with the bracket sequences for each SNP 121mer and finds the SNPs that cause major amino acid changes.
use Bio::Seq; 
use Bio::SeqIO;
use Bio::SeqUtils;
use Bio::Tools::OddCodes;
use Bio::PrimarySeq;
my $charge = ();
my $chemical = ();
my %hoh = ();
my @SNP = ();
##############################################################################################################################################################################
#
#		This section reads in the blastx file and makes an hoh based on SNP name with the original blastx translated match, a fragment of the original match and the SNP name.
#
#
open (INFILE, "$ARGV[0]") || die "Can't open infile.\n";
	while (<INFILE>) #read in from file, line by line
{
	chomp;	#remove end of line character
	if (m/^Tcm/) {	#skip any header information, start at genotype data
		@SNP = split (/\t/); #assign variables from genotype data
		$SNPname = shift(@SNP);  #take SNP name from the left end of the array
		$frame = pop(@SNP);  #take the frame from the right end of the array
		$SNPpep = pop(@SNP);  #take the original query translation with the X in it from the blastx report
		$hoh{$SNPname}{frame} = $frame;
		$hoh{$SNPname}{peptide_orig} = $SNPpep;
		$SNPfrag = $SNPpep;
		$SNPfrag =~ s/^\w{3}//;
		$SNPfrag =~ s/\w{27}$//;
		$hoh{$SNPname}{snpfrag} = $SNPfrag;
		}
}
###############################################################################################################################################################################
#
#		This section reads in the SNP 121mers in bracket format and makes an hoh based on SNP name.  It takes the brackets and makes two sequences with the different SNP #		nucleotides.  Uses the translate_both subroutine.
#
#
	$seqio_obj = Bio::SeqIO->new(-file => $ARGV[1], -format => "fasta" );
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

################################################################################################################################################################################	
#
#		This section sorts the hoh and translates the frame that matches the fragment from the blast report.  Different aa translations in the correct frame are added to
#		the hoh.  chemical, charge, functional.  Uses both the SNP_aa and framefind subroutines.  The chemical, charge and functional translations are compared by 		#		SNP_aa_compare and then a report is printed out.
#
#
#
	for $key  (sort keys %hoh)	{
	if (exists $hoh{$key}{snpfrag})	{
		&SNP_aa ($hoh{$key}{seq1}, seq1, $key);
		&SNP_aa ($hoh{$key}{seq2}, seq2, $key);
		&SNP_aa_compare ($key);
	}
}	
open (OUTFILE, ">>$ARGV[2]\.fasta") || die "Can't open outfile";  
	for $key (sort keys %hoh)	{
	if (exists $hoh{$key}{snpfrag})	{
#		if (exists $hoh{$key}{charge_diff} || $hoh{$key}{chemical_diff} || $hoh{$key}{functional_diff})	{
#		if (exists $hoh{$key}{charge_diff})	{				
#	print OUTFILE "$key\t$hoh{$key}{peptide_orig}\t$hoh{$key}{seq1_blastframe}\t$hoh{$key}{seq2_blastframe}\t$hoh{$key}{bracketseq}\t$hoh{$key}{seq1_diff_aa}\t$hoh{$key}{seq2_diff_aa}\n";
		if (exists $hoh{$key}{functional_diff})	{
	print OUTFILE "\>$key\n$hoh{$key}{bracketseq}\n";
		
        }
	}
}
####################################################################################################################################################################
sub 	SNP_aa	{
		my ($a,$b,$c) = @_;
		$seq_obj = Bio::Seq->new(-seq=> $a, -display_id=> $c, alphabet=> "dna");
		@seqs = Bio::SeqUtils->translate_6frames($seq_obj);
		$frame0 = $seqs[0]->seq;
		$hoh{$c}{"$b\_frame0"} = $frame0;
		&framefind ($frame0, $c, $seqs[0], $b);
		$frame1 = $seqs[1]->seq;
		$hoh{$c}{"$b\_frame1"} = $frame1;
		&framefind ($frame1, $c, $seqs[1]), $b;
		$frame2 = $seqs[2]->seq;
		$hoh{$c}{"$b\_frame2"} = $frame2;
		&framefind ($frame2, $c, $seqs[2], $b);
		$framerev0 = $seqs[3]->seq;
		$hoh{$c}{"$b\_revframe0"} = $framerev0;
		&framefind ($framerev0, $c, $seqs[3], $b);
		$framerev1 = $seqs[4]->seq;
		$hoh{$c}{"$b\_revframe1"} = $framerev1;
		&framefind ($framerev1, $c, $seqs[4], $b);
		$framerev2 = $seqs[5]->seq;
		$hoh{$c}{"$b\_revframe2"} = $framerev2;
		&framefind ($framerev2, $c, $seqs[5], $b);
		}

##############################################################################################################################################################################

sub framefind	{	#takes a frame from an input sequence, compares to the BLAST subfrag, if it matches, does a charge translation and a chemical translation of the amino acid sequence and stores it as part of the hash
	if ($_[0] =~ m/$hoh{$_[1]}{snpfrag}/)	{
		$hoh{$_[1]}{"$_[3]\_blastframe"} = $_[0];
		my $oddcode_obj = Bio::Tools::OddCodes->new(-seq => $_[2]);
		my $output = $oddcode_obj->charge();
		my $charge = $$output;
		$hoh{$_[1]}{"$_[3]\_charge"} = $charge;
		my $chemoutput = $oddcode_obj->chemical();
		my $chemical = $$chemoutput;
		$hoh{$_[1]}{"$_[3]\_chemical"} = $chemical;
		my $funcoutput = $oddcode_obj->functional();
		my $functional = $$funcoutput;
		$hoh{$_[1]}{"$_[3]\_functional"} = $functional;

	}
}

################################################################################################################################################################################

sub translate_both {
	my ($a,$b,$c,$d)= @_ ;
	$hoh{$d}{bracketseq}= $a;
	$seq1 = $a;
	$seq1 =~ s!\[$b/$c\]!$b!;
	$hoh{$d}{seq1} = $seq1;
	$seq2 = $a;
	$seq2 =~ s!\[$b/$c\]!$c!;
	$hoh{$d}{seq2} = $seq2;
}

#################################################################################################################################################################################
#
#		This subroutine compares the different amino acid translations and returns those that cause a significant change in charge, chemical or functional.
#
#
sub SNP_aa_compare	{
	my ($a) = @_ ;
	if ($hoh{$a}{seq1_charge} ne $hoh{$a}{seq2_charge})	{
		$hoh{$a}{charge_diff} = "Charge is different";
	&aa_id ($a);
	}
	if ($hoh{$a}{seq1_chemical} ne $hoh{$a}{seq2_chemical})	{
		$hoh{$a}{chemical_diff} = "Chemical is different";
	}
	if ($hoh{$a}{seq1_functional} ne $hoh{$a}{seq2_functional})	{
		$hoh{$a}{functional_diff} = "Functional is different"
	}
}	

##################################################################################################################################################################################


sub aa_id	{	#This subroutine identifies the amino acid difference between seq1 and seq2 for charge differences.
	my ($a) = @_;
	$seq1 = $hoh{$a}{seq1_blastframe};
	$seq2 = $hoh{$a}{seq2_blastframe};
	$seq_obj1 = Bio::Seq->new(-seq=>$hoh{$a}{seq1_blastframe}, -alphabet => "protein", -display_id => "seq1_blastframe" );
	$seq_obj2 = Bio::Seq->new(-seq=>$hoh{$a}{seq2_blastframe}, -alphabet => "protein",  -display_id => "seq2_blastframe");
	for ($i = 1; $i <= length $seq1; $i++) {
	$sq1 = $seq_obj1->subseq($i,$i);
	$sq2 = $seq_obj2->subseq($i,$i);
		if ($sq1 ne $sq2)	{
		$hoh{$a}{seq1_diff_aa} = $sq1;
		$hoh{$a}{seq2_diff_aa} = $sq2;
	}
}
}

	
