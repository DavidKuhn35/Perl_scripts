#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: February 15, 2014
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: Filter mango SNPs for Fluidigm design
# <input>		: [0] All_Mangifera... variant file, [1] 201revised fasta 
# <output>		: File with IUPAC 201mers to BLAST against mango CDS  
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
my @snp = ();
my $header = ();
my $header1 = ();
my %goodsnp = ();
my $snpreads = ();
my %refreads = ();
my %altreads = ();
my %refperc = ();
my %altperc = ();
my %reads = ();
my @dp4 = ();
my %score = ();
my %genotype = ();
my %filter1 = ();
my %cultreads = ();
my %snpname = ();
my %contigname = ();
my %snppos = ();
my %snpgeno = ();
my %filter2 = ();
my %orglength = ();
my %hoh = ();
my %snparray = ();
my %snparraynames = ();
my %filter3 = ();
my $IUPAC = ();
my %snpid = ();

####################################################################################################################
#  		
######################################################################################################################
#  Look up reference and alternate allele read info from variant report and make hashes for filtering
open (SNP_ID, "$ARGV[0]") || die "Can't open Mango SNP variant report $!\n";
	while (<SNP_ID>)
{
	chomp;
	if (m/^i/)	{
	@snp = split /\t/, $_;
			$snpname{$snp[0]} = $snp[0];
			$contigname{$snp[0]} = $snp[1];
			$snppos{$snp[0]} = $snp[2];
			$snpgeno{$snp[0]} = "$snp[4]$snp[5]";
			$snp[15] =~ s/\"//g;
			@dp4 = split /,/, $snp[15];
#				if (($snp[9] >= 400)&&($dp4[2]+$dp4[3]>=10))	{ 
			$reads{$snp[0]} = $snp[9]; 
			$refreads{$snp[0]} = $dp4[0] + $dp4[1];
			$altreads{$snp[0]} = $dp4[2] + $dp4[3];
			$refperc{$snp[0]} = ($dp4[0] + $dp4[1])/$snp[9];
			$altperc{$snp[0]} = ($dp4[2] + $dp4[3])/$snp[9];
				for (my $i = 20; $i <= 29; $i++)	{	#for statement uses ; to separate initial, final and increment 
#				print "this is the value of i $i and this is the genotype in the main $snp[$i]\n";
				my @z = split /:/, $snp[$i];
				$z[0] =~ s/"//g;
#				print "this is z[0] $z[0] and this is snp number $snp[0]\n";
				my $s = "$z[0]\t";	#gets genotype of 10 bins
				my $t = "$z[2]\t";	#gets read number of 10 bins
				$genotype{$snp[0]} .= $s;
				$cultreads{$snp[0]} .= $t;
#				print "this is the genotype $genotype{$snp[0]} for snp number $snp[0]\n";
			}
		}
		}
		
#		}

###################################################################################################################
#  Sort hashes by value to find snps with greatest number of reads and greatest minor allele frequency
#
	my @readsort = sort by_reads keys %reads;
	my @altreadsort = sort by_altreads keys %altreads;
	my @altpercsort = sort by_altperc keys %altperc;
#	foreach my $y (@altpercsort)	{
#		print "$y\t$reads{$y}\t$altreads{$y}\t$altperc{$y}\n";
#		}
	
#######################################################################################################################
#do first filter for MAF between 0.05-0.95 and altreads >= 10
#

	foreach my $x (@altpercsort)	{
			if (($altperc{$x} >= 0.05) && ($altperc{$x} <= 0.95))	{
				if ($altreads{$x} >= 10)	{
					$filter1{$x} = "$x\t$contigname{$x}\t$snppos{$x}\t$snpgeno{$x}\t$reads{$x}\t$altreads{$x}\t$altperc{$x}\t$genotype{$x}\t$cultreads{$x}";
#  filter for MAF 0.05-0.95 and altreads >= 10 and create new hash filter1					
#			printf MAF "$x\t$reads{$x}\t$altreads{$x}\t%.2f%%\t$genotype{$x}\n",$altperc{$x};
		}
	}	
}
########################################################################################################################
#filter 2 total reads >= 400 (or whatever eventually seems reasonable
#
	for my $q (sort keys %filter1)	{
		my @f = split /\t/, $filter1{$q};
			if ($f[4]>= 400)	{
			$filter2{$q} = $filter1{$q};
		}
	}

#########################################################################################################################
# filter 3 100nt flanks open 201mer revised fasta file
#
#	open (SNP, "$ARGV[1]") || die "Can't open mango SNP 201mer fasta infile$!\n";
#	while (<SNP>)
#{#
#	chomp;
#	if (m/^>i/)  	{	#If line has fasta header
#		my @g = split (/\s+/); #split on white space
#		my $h = $g[0];		#Make copy of snp name to strip # off to make it compatible with variant report
#		my $h =~ s/#//;		#Strips out # so $h is same snp name as in variant report

	my $xmerobj = Bio::SeqIO->new (-file =>"$ARGV[1]", -format => "fasta"); #get fasta objects
	while (my $xmer = $xmerobj->next_seq){
		my $h = $xmer->display_id;
			if (exists $filter2{$h})	{	#reduce to snps in filter2
					my @xdesc = split /\s+/, $xmer->desc; #get rest of header info
					my $newpos = $xdesc[1];
					my $newlen = $xdesc[0];
					$newlen =~ s/new_len=(\d+)/$1/;
					$newpos =~ s/snp_pos=(\d+)/$1/;	
					my $orglen = $xdesc[3];
					$orglen =~ s/org_len=(\d+)/$1/;
					$orglength{$contigname{$h}}= $orglen;	#makes hash with contigname as key and contiglength as value
					$hoh{$contigname{$h}}{$snpname{$h}}= $snppos{$h};	#makes hash of hash contigname=>snpname=>snp position
						#value of orglen for contig, use to pass to subroutine 
						#hash of hashes contigname=>snpname=>$orgpos\t$orglen should be able to sort by contigname to get hash for each contig name of snps
				}
			}


#close SNP;	

	for my $key1 (sort keys %hoh)	{
		for my $key2 (sort keys %{$hoh{$key1}})	{
			$snparraynames{$key1}{$hoh{$key1}{$key2}} = $key2;  #contig=>snp_pos=>snp_name 
			$snparray{$key1} .= "$hoh{$key1}{$key2}\t";# for each contig name, make a single hash entry which is an array of snp pos for that contig
			}
				my $orglen = $orglength{$key1};		
				my @snpflank1 = split /\t/, $snparray{$key1};
				unshift (@snpflank1, 0);
				push (@snpflank1, $orglen);
				my @snpflank = sort {$a <=> $b} @snpflank1;
				my @returnsnp = &flank100 (@snpflank);
			foreach my $ret (@returnsnp)	{
						if (exists $snparraynames{$key1}{$ret})	{ #if snp pos is same as in returned array put in with goodsnps
						$filter3{$snparraynames{$key1}{$ret}} = $filter2{$snparraynames{$key1}{$ret}};	#identifies snp name of snp with 100nt flank region contig=>snpname=>snpname

# filter3 hash has snpnames for snps that have 100nt flanks as keys, filter2 {snpname} as value. Keys for filter3 give number of goodsnps so far
					}
				}
			}
#######################################################################################################################
#   filter 4 Get 201 revised fasta file (# out of snpname), substitute IUPAC for SNP
#
	$xmerobj = Bio::SeqIO->new (-file =>"$ARGV[1]", -format => "fasta"); #get fasta objects
	while (my $xmer = $xmerobj->next_seq){
		$header = $xmer->display_id;
			if (exists $filter3{$header})	{
					my @xdesc = split /\s+/, $xmer->desc; #get rest of header info
					my $newpos = $xdesc[1];
					my $newlen = $xdesc[0];
					$newlen =~ s/new_len=(\d+)/$1/;
					$newpos =~ s/snp_pos=(\d+)/$1/;	
						if (($newpos == 101)&&($newlen >= 201))	{  #get SNP position, only take those at least 101 from end of contig and where there is 201 nt of sequence
						my $merseq = $xmer->seq;
						my $IUseq = $merseq;

					$IUPAC = $snpgeno{$header};	#get reference and alt nucleotides and join
						unless ($IUPAC eq "\.\.")	{
							if (($IUPAC eq "AC") || ($IUPAC eq "CA")) 	{	#replace dinucleotide SNP with IUPAC code
							$IUPAC = "M";
							$snpid{$header} = $IUPAC;
							}elsif (($IUPAC eq "AG") || ($IUPAC eq "GA"))	{
							$IUPAC = "R";
							$snpid{$header} = $IUPAC;
							}elsif (($IUPAC eq "AT") || ($IUPAC eq "TA")) 	{
							$IUPAC = "W";
							$snpid{$header} = $IUPAC;
							}elsif (($IUPAC eq "CG") || ($IUPAC eq "GC"))	{
							$IUPAC = "S";
							$snpid{$header} = $IUPAC;
							}elsif (($IUPAC eq "CT") || ($IUPAC eq "TC"))	{
							$IUPAC = "Y";
							$snpid{$header} = $IUPAC;
							}elsif (($IUPAC eq "GT") || ($IUPAC eq "TG"))	{
							$IUPAC = "K"; 
							$snpid{$header} = $IUPAC;
					}
				}

								if (exists $snpid{$header})	{
									my $new1 = $newpos - 1;
									substr($IUseq,$new1,1)= $snpid{$header};
									my $mer_out = Bio::Seq->new (-seq => $IUseq, -display_id => $header, -desc => join "\t", @xdesc); 
									my $mer_outIO = Bio::SeqIO->new (-file => ">>IUPAC_filter4_201_out.fasta", -format => "fasta"); #creates new fasta file with IUPAC substitutions for Blast against mango CDS
									 $mer_outIO->write_seq($mer_out);
								}	
							
			}
		}
}
##########################################################################################################
#  Report for mango variant filter
#	
	open (MANGO, ">>mango_variant_filter_report.txt") || die "Can't open mango variant report$!\n";
		my $filter1count = keys %filter1;
		my $filter2count = keys %filter2;
		my $filter3count = keys %filter3;
		my $snpidcount = keys %snpid;
		print MANGO "The number of SNPs after filter 1 is $filter1count\nThe number of SNPs after filter 2 is $filter2count\nThe number of SNPs after filter 3 is $filter3count\nThe number of SNPs after IUPAC substitution is $snpidcount\n";
#############################################################################################################
#subroutines to sort hashes by value
#

sub by_altperc { $altperc{$b} <=> $altperc{$a} }
sub by_reads	{$reads{$b} <=> $reads{$a}	}
sub by_altreads {$altreads{$b} <=> $altreads{$a} }

######################################################################################################################
#  sub flank100 identifies snps with 100nt snp free flank regions
sub flank100 {
	my @b = ();
	my $i = ();
	my @c = ();
	my (@a) = @_;
		for ($i = 0; $i <= $#a-1; $i++)	{ #go through the array, calculating whether snp flanks are 100 nt
		  $b[$i] = ($a[$i+1] - $a[$i]);		#snp a0 doesn't have flank of 100 on left side of snp
		} 
		my $bno = @b;
		for ( $i = 0; $i <= $bno-2 ; $i++)	{
				if (($b[$i] >=100)&&($b[$i+1] >= 100))	{
					push (@c, $a[$i+1]);
					}else{
					push (@c, 0);
		}
	}
		my @return = @c;		#return array c
	}
