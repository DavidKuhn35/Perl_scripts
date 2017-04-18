#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: 5/3/2012
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: reads in fasta file of all longest avo isotigs (/data/common/data/FASTAfiles/Avocado_longestIsotigs_110105.fna, makes them into sequence objects, reads in parsed
#	BLASTX report BLASTX_longestIsotigs_110105_TAIR10_with_ANNOT.table.refined.txt to put AGI name to isotig name, reads in file with SNP frequency linked to isotig name, reads in
#	 GOSlim annotation linked to AGI name, prints out table with all isotigs and GOSlim annotation plus SNP frequency for statistical analysis
# <input>		:  fasta file of avocado isotigs, parsed BLASTX report for annotation of avocado isotigs, file with SNP frequency linked to isotig name, file 
#			with GOSlim annotation linked to AGI name			
# <output>		:    
# 			
# Dependencies	: N/A
# Usage			: avo_isotig_length_from_fasta.pl /data/common/data/FASTAfiles/Avocado_longestIsotigs_110105.fna BLASTX_longestIsotigs_110105_TAIR10_with_ANNOT.table.refined.txt
#	QQunifiedSNPs_Avocado_BLASTX_TAIR10_ANNOT.table_cvg.txt_out.txt XXMay25_avo_all_cfp_goslim_out.txt avo_all_May3_isolen_snpcount_snpfreq_agi_goslim_out avo_all 
#	(or some identifier)
#     NAME OF PROGRAM DOES NOT REPRESENT CURRENT FUNCTION!               
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

	open (INFILE, "$ARGV[1]")	|| die "Can't open infile $!\n";	#ARGV[1] is BLASTX_TAIR10_...refined.txt which is the parsed BLASTX report for all the longest avo isotigs
	while (<INFILE>)	
{
	chomp;
	if (m/^g/)	{
	(my $isotig, my $AGI) = split (/\t/);
	my $new_AGI = $AGI;
	$new_AGI =~ s/(\w+).\d+/$1/;	#strip any numbers and . off AGI name to make it compatible with other AGI names
		if (exists $hoh{$isotig})	{
		$hoh{$isotig}{AGI} = $new_AGI;
	}	else 	{
		$hoh{$isotig}{AGI} = "ATXG00000";	#if no AGI name, give it default as a placeholder in the final output file
	}
}
}
open (SNPFREQ, "$ARGV[2]") || die "Can't open snpfreq infile$!\n";#ARGV[2] is the output from the variant report (see avo_snp.pl) giving the isotig name, snpcount, isotig length and snp frequency
	while (<SNPFREQ>)	
{
	chomp;
	(my $isofreq, my $snpcount, my $isofreqlen, my $snpfreq) = split (/\t/);
	if (exists $hoh{$isofreq})	{
		$hoh{$isofreq}{snpcount} = $snpcount;
		$hoh{$isofreq}{snpfreq} = $snpfreq;
	}
			
}

open (GOSLIM,"$ARGV[3]") || die "Can't open goslim infile\t$!\n";#goslim annotation for all longest isotigs
	while (<GOSLIM>)	
{
	chomp;
	my @goslim = split (/\|/);
	my $goslimagi = $goslim[0];
	shift @goslim;
	my $goslimt = join "\t", @goslim;
	$goslim{$goslimagi} = $goslimt;
}	
		for my $x (sort keys %hoh)	{
			my $z = $hoh{$x}{AGI};
			if (exists $goslim{$z})	{
				$hoh{$x}{goslim} = $goslim{$z};
				}  else   {			
			$goslimt_default = join "\t", ("C","Cdate","Cgoslimterm","Cgoid","Crelation","Cgoterm","F","Fdate","Fgoslimterm","Fgoid","Frelation","Fgoterm","P","Pdate","Pgoslimterm","Pgoid","Prelation","Pgoterm");
			$hoh{$x}{goslim} = $goslimt_default;
	}
}

open (OUTFILE, ">>$ARGV[4]\_isolen_snpcount_snpfreq_agi_goslim_out" ) || die "Can't open outfile $!\n";#name of outfile for most complete dataset
	print OUTFILE "isotig\tlength\tsnp_count\tsnp_frequency\tAGI\t$goslimt_default\n";

	for my $key (sort keys %hoh)	{
		print OUTFILE "$key\t$hoh{$key}{isolen}\t$hoh{$key}{snpcount}\t$hoh{$key}{snpfreq}\t$hoh{$key}{AGI}\t$hoh{$key}{goslim}\n";
	}
