#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: February 23, 2014
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: Read in fasta file with IUPAC codes, output txt file with SNP name and sequence with brackets around snp.
# <input>		: filter6 fasta file 		
# <output>		: txt file with snpname and sequence with brackets around snp
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

#######################################################################################################################
#   Get 201 or 401 mer with SNP, substitute IUPAC for SNP
#
	open (BRACKET, ">>$ARGV[1]") || die "Can't open output file for bracket snps $!\n";

	$xmerobj = Bio::SeqIO->new (-file =>"$ARGV[0]", -format => "fasta"); #get fasta objects
	while ($xmer = $xmerobj->next_seq){
		$header = $xmer->display_id;
		my $orgseq = $xmer->seq;
		$merseq = $xmer->seq;
			if ($merseq =~ m/R/)	{
				$merseq =~ s/R/\[A\/G\]/;
			} elsif ($merseq =~ m/Y/) {
				$merseq =~ s/Y/\[C\/T\]/;
			}  elsif ($merseq =~m/K/) {
				$merseq =~ s/K/\[G\/T\]/;
			}  elsif ($merseq =~ m/M/) {
				$merseq =~ s/M/\[A\/C\]/;
			}  elsif ($merseq =~ m/S/) {
				$merseq =~ s/S/\[G\/C\]/;
			}  elsif ($merseq =~ m/W/) {
				$merseq =~ s/W/\[A\/T\]/;
			}
			print BRACKET "$header\t$merseq\n";
									
		}
