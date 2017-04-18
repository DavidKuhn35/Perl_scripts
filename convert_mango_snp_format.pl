#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: February 4, 2014
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: Find synonymous mango SNPs for Fluidigm design
# <input>		: Mangifera_Run2_VarCountsWithTAIR10blast.txt, FirstHalfMangifera_ALL_VarRun2_withTAIR10blast.txt (SecondHalf...), renamed_with_position_Mangifera_All_Run2_SNPflanks_201.fasta.txt (...401..) 		
# <output>		: File with IUPAC 201 or 401 mers to BLAST against mango CDS  
# 			
# Dependencies	: N/A
# Usage			: 
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

#####################################################################################################################
#	main	Make hashes of all contigs with number of variants and blast results
my %var1 = ();
my %var3up = ();
my %varblast = ();
my %var1blast = ();
my %var2blast = ();
my %snpid = ();
my @snp = ();
my $IUPAC = ();
my @mer = ();
my @header = ();
my @varno = ();
my $xmerobj = ();
my $xmer = ();
my $header = ();
my $header1 = ();
my $desc = ();
my @xdesc = ();
my @id = ();
my $newpos = ();
my $contig = ();
my $merseq = ();
my $mer_out = ();
my $mer_outIO = ();
my $IUseq = ();
my $new1 = ();
my $newlen = ();
my %var2 = ();
my $var1count = ();
#Read in the avo snp file and get rid of rows with missing data or Ns
open (VAR_NO, "$ARGV[0]") || die "Can't open mango variant count infile $!\n";
	while (<VAR_NO>)
{
	s/#//;
open (OUT, ">>mango_snp_201_reformat.fasta") || die "Can't open var summary$!\n";
	print OUT;	
}

	
				
######################################################################################################################
#  Look up SNP for all 1 var contigs store as IUPAC code

