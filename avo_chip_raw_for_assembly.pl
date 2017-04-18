#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: 2/4/2013
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: reads in fasta file of all longest avo isotigs (/data/common/data/FASTAfiles/Avocado_longestIsotigs_110105.fna, makes them into sequence objects, reads in Avo1_14f_6kSNP_list_for_Illumina.txt which has both SNP names, 121 mer sequences and isotig names from which the 121mers came.  Add whole isotig sequence data to file with SNP name.
# <input>		:  fasta file of avocado isotigs, 
# 			
# Dependencies	: N/A
# Usage			: avo_isotig_length_from_fasta.pl /data/common/data/FASTAfiles/Avocado_longestIsotigs_110105.fna Avo1_14f_6kSNP_list_for_Illumina.txt
#	
#                 
##################################################################################
use strict;
use warnings;
use diagnostics;
use Bio::Seq; 
use Bio::SeqIO;
use Bio::SeqUtils;
use Bio::Tools::OddCodes;
use Bio::PrimarySeq;
my $charge = ();
my $chemical = ();
my %hoh = ();
my %hoh2 = ();
my @SNP = ();
my $goslimt_default = "foobar";
my $new_y = "ATfoobar";
my %goslim = ();


	my $seqio_obj = Bio::SeqIO->new(-file => $ARGV[0], -format => "fasta" );	#ARGV[0] is fasta file of all longest avo isotigs
	while (my $seq_obj = $seqio_obj->next_seq)  {
 	my $fasta_name = $seq_obj->display_id;
	my $seqlen = $seq_obj->length;
	my $seq = $seq_obj->seq;
	$hoh{$fasta_name}{isotig} = $fasta_name;
	$hoh{$fasta_name}{isolen} = $seqlen;
	$hoh{$fasta_name}{seq} = $seq;
	$hoh{$fasta_name}{snpcount} = 0;
	$hoh{$fasta_name}{snpfreq} = 0;
	
}

	open (INFILE, "$ARGV[1]")	|| die "Can't open infile $!\n";	#ARGV[1] is Avo1_14f_6kSNP_list_for_Illumina.txt 
	while (<INFILE>)	
{
	chomp;
	if (m/^SHRS/)	{
	my @chip = split (/\t/);
#	my $snpname = $chip[0];	#get snpname
#	my $seq121mer = $chip[1];	#get 121mer
#	my $isotigname = $chip[7];	#get isotigname
	$hoh2{$chip[0]}{seq121mer} = $chip[1];	#put 121mer sequence into hoh with snpname
	$hoh2{$chip[0]}{isotigname} = $chip[7]; # put isotigname with SNP name to show which isotig 121mer came from
	$hoh2{$chip[0]}{isotigseq} = $hoh{$chip[7]}{seq}; #put isotig sequence into hoh with snpname 
# Now should have hoh2 with all the information needed to print out SNP name, 121mer and original isotig sequence		
	}
}


open (OUTFILE, ">>avo_chip_data_for_assembly.txt" ) || die "Can't open outfile $!\n";#name of outfile for chip SNPs with 121mers and original isotig sequences
	print OUTFILE "SNPname\t121mer\tisotig_name\tisotig_sequence\n";

	for my $key (sort keys %hoh2)	{
		print OUTFILE "$key\t$hoh2{$key}{seq121mer}\t$hoh2{$key}{isotigname}\t$hoh2{$key}{isotigseq}\n";
	}
