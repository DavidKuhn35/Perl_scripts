#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: February 4, 2014
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: Find synonymous mango SNPs for Fluidigm design from contigs with only one SNP
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
	chomp;
	if (m/^i/)  	{	#If line has contig name for mango but not unique identifier
	@varno = split /\s+/;
		if ($varno[1] == 1)	{ #if contig has only one variant, store in var1
		$var1{$varno[0]} = $varno[1];
			}elsif ($varno[1] == 2)	{
			$var2{$varno[0]} = $varno[1];
			}	else	{
				$var3up{$varno[0]} = $varno[1]; #if contig has more than two variant store in var2up
			}
	if (exists $varno[2] )	{ #if contig has blast result store in varblast
	$varblast{$varno[0]} = $varno[3];
	}
	if ((exists $var1{$varno[0]})&&(exists $varblast{$varno[0]}))	{
		$var1blast{$varno[0]} = $varno[3];
		} elsif ((exists $var2{$varno[0]})&&(exists $varblast{$varno[0]}))	{
		$var2blast{$varno[0]} = $varno[3];
		}			
	}
}
open (OUT, ">>var_summary.txt") || die "Can't open var summary$!\n";
	$var1count = keys %var1;	
	print OUT "The number of contigs with 1 variant is $var1count\n";
	my $var2count = keys %var2;
	print OUT "The number of contigs with 2 variants is $var2count\n";
	my $var3upcount = keys %var3up;
	print OUT "The number of contigs with 2 or more variants is $var3upcount\n";			
	my $varblastcount = keys %varblast;
	print OUT "The number of contigs with a BLAST report is $varblastcount\n";
	my $var1blastcount = keys %var1blast;
	print OUT "The number of contigs with one variant and a BLAST report is $var1blastcount\n";
	my $var2blastcount = keys %var2blast;
	print OUT "The number of contigs with two variants and a BLAST report is $var2blastcount\n"; 
				
######################################################################################################################
#  Look up SNP for all 1 var contigs store as IUPAC code
open (SNP_ID, "$ARGV[1]") || die "Can't open Mango SNP variant report $!\n";
	while (<SNP_ID>)
{
	chomp;
	if (m/^i/)	{
	@snp = split /\t/, $_;
	$IUPAC = "$snp[4]$snp[5]";	#get reference and alt nucleotides and join
	unless ($IUPAC eq "\.\.")	{
		if ((exists $var1{$snp[1]})||(exists $var2{$snp[1]}))	{
			if (($IUPAC eq "AC") || ($IUPAC eq "CA")) 	{	#replace dinucleotide SNP with IUPAC code
				$IUPAC = "M";
				$snpid{$snp[0]} = $IUPAC;
			}elsif (($IUPAC eq "AG") || ($IUPAC eq "GA"))	{
				$IUPAC = "R";
				$snpid{$snp[0]} = $IUPAC;
			}elsif (($IUPAC eq "AT") || ($IUPAC eq "TA")) 	{
				$IUPAC = "W";
				$snpid{$snp[0]} = $IUPAC;
			}elsif (($IUPAC eq "CG") || ($IUPAC eq "GC"))	{
				$IUPAC = "S";
				$snpid{$snp[0]} = $IUPAC;
			}elsif (($IUPAC eq "CT") || ($IUPAC eq "TC"))	{
				$IUPAC = "Y";
				$snpid{$snp[0]} = $IUPAC;
			}elsif (($IUPAC eq "GT") || ($IUPAC eq "TG"))	{
				$IUPAC = "K"; 
				$snpid{$snp[0]} = $IUPAC;
			}
		}
	}
}
}
#######################################################################################################################
#   Get 201 or 401 mer with SNP, substitute IUPAC for SNP
#open (XMER, "$ARGV[2]") || die "Can't open 201 or 401mer file $!\n";
	$xmerobj = Bio::SeqIO->new (-file =>"$ARGV[2]", -format => "fasta"); #get fasta objects
	while ($xmer = $xmerobj->next_seq){
		$header = $xmer->display_id;
		my $contigvar = $header;
		my $contigname = $header;
		$contigname =~ s/_#\d+//;
		if (exists $var1{$contigname})	{
#		print "$contig\n";
		$contigvar =~ s/#//;  #strip # from unique SNP identifier to make it compatible with variant report
			if (exists $snpid{$contigvar})	{
#		print "$contig\n";
#		print "$header\t $desc\n";
				@xdesc = split /\s+/, $xmer->desc; #get rest of header info
				$newpos = $xdesc[1];
				$newlen = $xdesc[0];
				$newlen =~ s/new_len=(\d+)/$1/;
				$newpos =~ s/snp_pos=(\d+)/$1/;	
					if (($newpos == 101)&&($newlen >= 201))	{  #get SNP position, only take those at least 101 from end of contig and where there is 201 nt of sequence
						$merseq = $xmer->seq;
						$IUseq = $merseq;
							if (exists $snpid{$contigvar})	{
#								print "$IUseq\n";
#								print "$snpid{$contig}\n";
								$new1 = $newpos - 1;
								substr($IUseq,$new1,1)= $snpid{$contigvar};
#								print "$IUseq\n";
								$mer_out = Bio::Seq->new (-seq => $IUseq, -display_id => $contigvar, -desc => join "\t", @xdesc); #may also use header for contigs with more than one SNP
								$mer_outIO = Bio::SeqIO->new (-file => ">>var1_ALL_201_filter_goodSNPs_out.fasta", -format => "fasta"); #creates new fasta file with IUPAC substitutions for Blast against mango CDS
								$mer_outIO->write_seq($mer_out);
					}
				}
			}
		}
}

