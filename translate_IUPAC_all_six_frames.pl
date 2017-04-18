#!/usr/bin/perl
#Takes a parsed blastx file with the name of a 121mer that had an X in its translation and the translation frame and a fasta file with the separated possible sequences for each SNP 121mer and outputs the two translated 121mers.
use Bio::Seq; 
use Bio::SeqIO;
use Bio::SeqUtils;

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
#		$SNP1 = "$SNPframe"."_1";
		$hoh{$SNPframe}{frame} = $frame;
		$hoh{$SNPframe}{peptide_orig} = $SNPpep;
#		$SNP2 = "$SNPframe"."_2";
#		$hoh{$SNP2}{frame} = $frame;
#		$hoh{$SNP2}{peptide_orig} = $SNPpep;
		print "The frame for $SNPframe is $hoh{$SNPframe}{frame}\nThe original query peptide sequence is $hoh{$SNPframe}{peptide_orig}\n";
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
		$frame1 = $seqs[1]->seq;
		$hoh{$fasta_name}{frame1} = $frame1;
		$frame2 = $seqs[2]->seq;
		$hoh{$fasta_name}{frame2} = $frame2;
		$framerev0 = $seqs[3]->seq;
		$hoh{$fasta_name}{revframe0} = $framerev0;
		$framerev1 = $seqs[4]->seq;
		$hoh{$fasta_name}{revframe1} = $framerev1;
		$framerev2 = $seqs[5]->seq;
		$hoh{$fasta_name}{revframe2} = $framerev2;
	}
}	
#		print "The frame is $pepframe\n";
#		$seqpep_obj_frame0 = $seq_obj->translate(-frame => 0);
#		$newpeptide_frame0 = $seqpep_obj_frame0->seq;
#		$seqpep_obj_revframe0 = $seq_obj->translate(-frame =>)->revcom->seq;
#		$newpeptide_rev_frame0 = $seqpep_obj_revframe0->seq;
#		print "The forward translated sequence for frame 0 for $fasta_name is $newpeptide_frame0\n.  The reverse sequence for frame 0 for $fasta_name is $newpeptide_rev_frame0 .\nThe original frame is $pepframe and the originaltranslated sequence is $hoh{$fasta_name}{peptide_orig}\n.";
#		$pepoutput_obj->write_seq($seqpep_obj_frame0);
#		$pepoutput_obj->write_seq($seqpep_obj_revframe0);
#		$peptide = $seq_obj->translate(-frame => "$pepframe");
#		$hoh{$fasta_name}{peptideframe0} = $newpeptide_frame0;
#		$hoh{$fasta_name}{peptiderevframe0} = $newpeptide_revframe0;

#		$seqpep_obj_frame1 = $seq_obj->translate(-frame => 1);
#		$newpeptide_frame1 = $seqpep_obj_frame1->seq;
#		$seqpep_obj_revframe1 = $seq_obj->translate(-frame => 1)->revcom->seq;
#		$newpeptide_revframe1 = $seqpep_obj_revframe1->seq;
#		print "The forward translated sequence for frame 1 for $fasta_name is $newpeptide_frame1\n. The reverse sequence for frame 1 for $fasta_name is $newpeptide_rev_frame1. \n The original frame is $pepframe and the originaltranslated sequence is $hoh{$fasta_name}{peptide_orig}\n.";
#		$pepoutput_obj->write_seq($seqpep_obj_frame1);
#		$pepoutput_obj->write_seq($seqpep_obj_revframe1);
#		$peptide = $seq_obj->translate(-frame => "$pepframe");
#		$hoh{$fasta_name}{peptideframe1} = $newpeptide_frame1;
#		$hoh{$fasta_name}{peptiderevframe1} = $newpeptide_revframe1;

#		$seqpep_obj_frame2 = $seq_obj->translate(-frame => 2);
#		$newpeptide_frame2 = $seqpep_obj_frame2->seq;
#		$seqpep_obj_revframe2 = $seq_obj->translate(-frame => 1)->revcom->seq;
#		$newpeptide_revframe2 = $seqpep_obj_revframe2->seq;
#		print "The translated sequence for frame 2 for $fasta_name is $newpeptide_frame2\n. The reverse sequence for frame 1 for $fasta_name is $newpeptide_rev_frame1.\n The original frame is $pepframe and the originaltranslated sequence is $hoh{$fasta_name}{peptide_orig}\n.";
#		$pepoutput_obj->write_seq($seqpep_obj_frame2);
#		$pepoutput_obj->write_seq($seqpep_obj_revframe2);		
#		$peptide = $seq_obj->translate(-frame => "$pepframe");
#		$hoh{$fasta_name}{peptideframe2} = $newpeptide_frame2;
#		$hoh{$fasta_name}{peptiderevframe2} = $newpeptide_revframe2;		
#		}
#	}

open (OUTFILE, ">>$ARGV[2]_name_peptides.txt") || die "Can't open outfile";  
#Step through the new hoh, make two new hashes, one for the SNP names and one for the appended genotypes.
for $key (sort keys %hoh) {
	print OUTFILE "$key\t Frame 0 \t$hoh{$key}{frame0}\n$key\t Frame 1 \t$hoh{$key}{frame1}\n$key\t Frame 2 \t$hoh{$key}{frame2}\n$key\tRev Frame 0\t$hoh{$key}{revframe0}\n$key\tRev Frame 1\t$hoh{$key}{revframe1}\n$key\tRev Frame 2\t$hoh{$key}{revframe2}\nOriginal frame is $hoh{$key}{frame}\t$hoh{$key}{peptide_orig}\n";	
        }
