#!/usr/bin/perl
#Takes a parsed blastx file with the name of a 121mer that had an X in its translation and the translation frame and a fasta file with the separated possible sequences for each SNP 121mer and outputs the two translated 121mers.
use Bio::Seq; 
use Bio::SeqIO;

my %hoh = ();
my @SNP = ();

open (INFILE, "$ARGV[0]") || die "Can't open infile.\n";
	while (<INFILE>) #read in from file, line by line
{
	chomp;	#remove end of line character
	if (m/^Tcm/) {	#skip any header information, start at genotype data
		@SNP = split (/\t/); #assign variables from genotype data
		$SNPframe = shift(@SNP);  #take SNP name from the left end of the array
		print "$SNPframe\n";
		$frame = pop(@SNP);  #take the frame from the right end of the array
		print "$frame\n";
		$SNPpep = pop(@SNP);  #take the original query translation with the X in it from the blastx report
		$SNP1 = "$SNPframe"."_1";
		$hoh{$SNP1}{frame} = $frame;
		$hoh{$SNP1}{peptide_orig} = $SNPpep;
		$SNP2 = "$SNPframe"."_2";
		$hoh{$SNP2}{frame} = $frame;
		$hoh{$SNP2}{peptide_orig} = $SNPpep;
		print "The frame for $SNP1 is $hoh{$SNP1}{frame}\nThe original query peptide sequence is $hoh{$SNP1}{peptide_orig}\n";
		print "The frame for $SNP2 is $hoh{$SNP2}{frame}\nThe original query peptide sequence is $hoh{$SNP2}{peptide_orig}\n";
		}
}
	


	
	$seqio_obj = Bio::SeqIO->new(-file => $ARGV[1], -format => "fasta" );  #get separate 121mer nucleotide sequences that represent both nucleotides at the SNP position
	$pepoutput_obj = Bio::SeqIO->new (-file => ">$ARGV[2]_separatepep.fasta", -format => 'fasta' );
	while ($seq_obj = $seqio_obj->next_seq)  {
 	$fasta_name = $seq_obj->display_id;
	if (exists $hoh{$fasta_name})	{
		$pepframe = $hoh{$fasta_name}{frame};
		if ($pepframe = 0)	{
		print "The frame is $pepframe\n";
		$seqpep_obj = $seq_obj->translate(-frame => 0);
		$newpeptide = $seqpep_obj->seq;
		print "The translated sequence for $fasta_name is $newpeptide\n.  The original translated sequence is $hoh{$fasta_name}{peptide_orig}\n.";
		$pepoutput_obj->write_seq($seqpep_obj);
#		$peptide = $seq_obj->translate(-frame => "$pepframe");
		$hoh{$fasta_name}{peptide} = $newpeptide;
		}elsif ($pepframe = 1) {
		print "The frame is $pepframe\n";
		$seqpep_obj = $seq_obj->translate(-frame => 1);
		$newpeptide = $seqpep_obj->seq;
		print "The translated sequence for $fasta_name is $newpeptide\n.  The original translated sequence is $hoh{$fasta_name}{peptide_orig}\n.";
		$pepoutput_obj->write_seq($seqpep_obj);
#		$peptide = $seq_obj->translate(-frame => "$pepframe");
		$hoh{$fasta_name}{peptide} = $newpeptide;
		}elsif ($pepframe = 2) {
		print "The frame is $pepframe\n";
		$seqpep_obj = $seq_obj->translate(-frame => 2);
		$newpeptide = $seqpep_obj->seq;
		print "The translated sequence for $fasta_name is $newpeptide\n.  The original translated sequence is $hoh{$fasta_name}{peptide_orig}\n.";
		$pepoutput_obj->write_seq($seqpep_obj);
#		$peptide = $seq_obj->translate(-frame => "$pepframe");
		$hoh{$fasta_name}{peptide} = $newpeptide;
		}
	}
}
open (OUTFILE, ">>$ARGV[2]_name_peptides.txt") || die "Can't open outfile";  
#Step through the new hoh, make two new hashes, one for the SNP names and one for the appended genotypes.
for $key (sort keys %hoh) {
	print OUTFILE "$key\t$hoh{$key}{peptide}\n$key\t$hoh{$key}{peptide_orig}\n";	
        }
