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
my %goodsnp = ();
####################################################################################################################
#  Read in list of good snps
open (GOODSNP, "$ARGV[0]") || die "Can't read in good snp file$!\n";
	while (<GOODSNP>)	{
		chomp;
		my $snpname = $_;
		$snpname =~ s/>//;
		$goodsnp{$snpname}= $snpname;
#	print "this is a goodsnp name $goodsnp{$snpname}\n";
	}				
######################################################################################################################
#  Look up SNP for all 1 var contigs store as IUPAC code
open (SNP_ID, "$ARGV[1]") || die "Can't open Mango SNP variant report $!\n";
	while (<SNP_ID>)
{
	chomp;
	if (m/^i/)	{
	@snp = split /\t/, $_;
			if (exists $goodsnp{$snp[0]})	{
			$IUPAC = "$snp[4]$snp[5]";	#get reference and alt nucleotides and join
			unless ($IUPAC eq "\.\.")	{
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
#
	$xmerobj = Bio::SeqIO->new (-file =>"$ARGV[2]", -format => "fasta"); #get fasta objects
	while ($xmer = $xmerobj->next_seq){
		$header = $xmer->display_id;
			if (exists $goodsnp{$header})	{
					@xdesc = split /\s+/, $xmer->desc; #get rest of header info
					$newpos = $xdesc[1];
					$newlen = $xdesc[0];
					$newlen =~ s/new_len=(\d+)/$1/;
					$newpos =~ s/snp_pos=(\d+)/$1/;	
						if (($newpos == 101)&&($newlen >= 201))	{  #get SNP position, only take those at least 101 from end of contig and where there is 201 nt of sequence
						$merseq = $xmer->seq;
						$IUseq = $merseq;
							if (exists $snpid{$header})	{
								$new1 = $newpos - 1;
								substr($IUseq,$new1,1)= $snpid{$header};
								$mer_out = Bio::Seq->new (-seq => $IUseq, -display_id => $header, -desc => join "\t", @xdesc); #may also use header for contigs with more than nne SNP
								$mer_outIO = Bio::SeqIO->new (-file => ">>IUPAC_goodsnps_201_out.fasta", -format => "fasta"); #creates new fasta file with IUPAC substitutions for Blast against mango CDS
								$mer_outIO->write_seq($mer_out);
				}
			}
		}
}
