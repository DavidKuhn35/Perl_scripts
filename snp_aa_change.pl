#!/usr/bin/perl
#Takes a parsed blastx file with the name of a 121mer that had an X in its translation and the translation frame and a fasta file with the separated possible sequences for each SNP 121mer and outputs the two translated 121mers.
use Bio::Seq; 
use Bio::SeqIO;
use Bio::SeqUtils;
use Bio::Tools::OddCodes;
use Bio::PrimarySeq;
my $charge = ();
my $chemical = ();
my %hoh = ();
my @SNP = ();

open (INFILE, "$ARGV[0]") || die "Can't open infile.\n";
	while (<INFILE>) #read in from file, line by line
{
	chomp;	#remove end of line character
	if (m/^Tcm/) {	#skip any header information, start at genotype data
		@SNP = split (/\t/); #assign variables from genotype data
		$SNPframe = shift(@SNP);  #take SNP name from the left end of the array
#		print "$SNPframe\n";
		$frame = pop(@SNP);  #take the frame from the right end of the array
#		print "$frame\n";
		$SNPpep = pop(@SNP);  #take the original query translation with the X in it from the blastx report
		$SNP1 = "$SNPframe"."_1";
		$hoh{$SNP1}{frame} = $frame;
		$hoh{$SNP1}{peptide_orig} = $SNPpep;
#		my $seqobj = Bio::PrimarySeq->new(-seq=> $SNPpep, -alphabet => 'protein', -id => $SNPframe);
#		my $oddcode_obj = Bio::Tools::OddCodes->new(-seq => $seqobj);
#		my $output = $oddcode_obj->charge();
#		$charge = $$output;
#		my $output2 = $oddcode_obj->chemical();
#		$chemical = $$output2;
#		print "charge is $charge\tchemical is $chemical\n";
		$SNPfrag = $SNPpep;
#		print "SNPfrag is $SNPfrag\n";
		$SNPfrag =~ s/^\w{3}//;
#		print "SNPfrag is now three aas shorter: $SNPfrag\n";
		$SNPfrag =~ s/\w{27}$//;
#		print "SNPfrag is now 27 aas shorter from the carboxy terminus $SNPfrag\n";
		$hoh{$SNP1}{snpfrag} = $SNPfrag;
		$hoh{$SNP1}{origname} = $SNPframe;
		$SNP2 = "$SNPframe"."_2";
		$hoh{$SNP2}{frame} = $frame;
		$hoh{$SNP2}{peptide_orig} = $SNPpep;
		$hoh{$SNP2}{snpfrag} = $SNPfrag;
		$hoh{$SNP2}{origname} = $SNPframe;
#		print "The frame for $SNPframe is $hoh{$SNPframe}{frame}\nThe original query peptide sequence is $hoh{$SNPframe}{peptide_orig}\n";
#		print "The frame for $SNP2 is $hoh{$SNP2}{frame}\nThe original query peptide sequence is $hoh{$SNP2}{peptide_orig}\n";
		}
}
	
	$seqio_obj = Bio::SeqIO->new(-file => $ARGV[1], -format => "fasta" );  #get separate 121mer nucleotide sequences that represent both nucleotides at the SNP position
	$pepoutput_obj = Bio::SeqIO->new (-file => ">$ARGV[2]_separatepep.fasta", -format => 'fasta' );
	while ($seq_obj = $seqio_obj->next_seq)  {
 	$fasta_name = $seq_obj->display_id;
	if (exists $hoh{$fasta_name})	{
		$pepframe = $hoh{$fasta_name}{frame};
		@seqs = Bio::SeqUtils->translate_6frames($seq_obj);
		$frame0 = $seqs[0]->seq;
		$hoh{$fasta_name}{frame0} = $frame0;
		&framefind ($frame0, $fasta_name, $seqs[0]);
		$frame1 = $seqs[1]->seq;
		$hoh{$fasta_name}{frame1} = $frame1;
		&framefind ($frame1, $fasta_name, $seqs[1]);
		$frame2 = $seqs[2]->seq;
		$hoh{$fasta_name}{frame2} = $frame2;
		&framefind ($frame2, $fasta_name, $seqs[2]);
		$framerev0 = $seqs[3]->seq;
		$hoh{$fasta_name}{revframe0} = $framerev0;
		&framefind ($framerev0, $fasta_name, $seqs[3]);
		$framerev1 = $seqs[4]->seq;
		$hoh{$fasta_name}{revframe1} = $framerev1;
		&framefind ($framerev1, $fasta_name, $seqs[4]);
		$framerev2 = $seqs[5]->seq;
		$hoh{$fasta_name}{revframe2} = $framerev2;
		&framefind ($framerev2, $fasta_name, $seqs[5]);
	}
}	

	for $key (sort keys %hoh)	{
		

open (OUTFILE, ">>$ARGV[2]_.txt") || die "Can't open outfile";  
#Step through the new hoh, make two new hashes, one for the SNP names and one for the appended genotypes.
for $key (sort keys %hoh) {
	print OUTFILE "$key\t Frame 0 \t$hoh{$key}{frame0}\n$key\t Frame 1 \t$hoh{$key}{frame1}\n$key\t Frame 2 \t$hoh{$key}{frame2}\n$key\tRev Frame 0\t$hoh{$key}{revframe0}\n$key\tRev Frame 1\t$hoh{$key}{revframe1}\n$key\tRev Frame 2\t$hoh{$key}{revframe2}\nOriginal frame is $hoh{$key}{frame}\t$hoh{$key}{peptide_orig}\n";	
        }

sub framefind	{	#takes a frame from an input sequence, compares to the BLAST subfrag, if it matches, does a charge translation and a chemical translation of the amino acid sequence and stores it as part of the hash
	if ($_[0] =~ m/$hoh{$_[1]}{snpfrag}/)	{
		print "SNP name is $_[1] and SNP sequence is $_[0]\n";
		my $oddcode_obj = Bio::Tools::OddCodes->new(-seq => $_[2]);
		my $output = $oddcode_obj->charge();
		my $charge = $$output;
#		print "The charge aa sequence is $charge\n";
		$hoh{$_[1]}{charge} = $charge;
		my $chemoutput = $oddcode_obj->chemical();
		my $chemical = $$chemoutput;
#		print "The chemical aa sequence is $chemical\n";
		$hoh{$_[1]}{chemical} = $chemical;
		my $funcoutput = $oddcode_obj->functional();
		my $functional = $$funcoutput;
		$hoh{$_[1]}{functional} = $functional;

	}
}
