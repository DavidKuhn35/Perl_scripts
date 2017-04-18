#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: Sunday, February 9, 2014
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: Reads mango variant 201mer report, makes hoh of contigname=>snpname=>org_pos\torg_len, then filters for SNPs that have 100nt SNP free flanks to create a file for analysis using IUPAC_mango_201_401.pl
# <input>		:  renamed_with_position_All_Run2_SNPflanks_201_fasta.txt		
# <output>		:  mango_100nt_flankSNPs.fasta
# 			
# Dependencies	: N/A
# Usage:	 mango_100nt_flankSNPs.pl  Create fasta file with mango SNPs with 100nt SNP free flanking regions to be analyzed for synonymous SNPs
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
#####################################################################################
#  Define global variables
my @snpname = ();
my $contigname = ();
my %hoh = ();
my %cntgname = ();
my $varname = ();
my %poslen = ();
my %variantname = ();
my %goodsnp = ();
my $orgpos = ();
my $orglen = ();
my @returnsnp = ();
my @snpflank = ();
my %orglength = ();
my %snparray = ();
my @return = ();
my $snpno = ();
my %snparraynames = ();

#####################################################################################################################
#	Read in mango variant report with 201 mers for all SNPs
open (SNP, "$ARGV[0]") || die "Can't open mango SNP infile$!\n";
	while (<SNP>)
{
	chomp;
	if (m/^>i/)  	{	#If line has fasta header
		(@snpname) = split (/\s+/); #split on white space
		$contigname = $snpname[0];	#Make copy of snp name to strip off ending to get contig name
		$varname = $snpname[0];		#Make copy of snp name to strip # off to make it compatible with variant report
		$varname =~ s/#//;		#Strips out #
		$contigname =~ s/_#\d+//;	#Strips snp identifier off contig name	
		$cntgname{$varname} = $contigname;	#Load cntgname hash with contigname for every unique SNP identifier in variant report format
		$variantname{$varname} = $varname;
		$orgpos = $snpname[3];
		$orgpos =~ s/org_snp_pos=(\d+)/$1/;
		$orglen = $snpname[4];
		$orglen =~ s/org_len=(\d+)/$1/;
		$orglength{$contigname}= $orglen;
		$hoh{$cntgname{$varname}}{$variantname{$varname}}= $orgpos;
						#value of orglen for contig, use to pass to subroutine 
						#hash of hashes contigname=>snpname=>$orgpos\t$orglen should be able to sort by contigname to get hash for each contig name of snps
		
			}

}
close SNP;
#print out hash to see what has been stored
#	for my $f (sort keys %hoh)	{
#		for my $g (sort keys %{$hoh{$f}})	{
#		print "The first key is $f, the second key is $g and the value is $hoh{$f}{$g}\n";
#		push @snpflank, $hoh{$f}{$g};
#			foreach my $test (@snpflank)	{
#			print "$test is the $snpflank[$test] element in \@snpflank\n";
#		}
#		}
#	}
#	Get info from hash into format for subroutine
open (GOODSNP, ">>all_mango_goodsnp_name_Feb10_2014.txt") || die "Couldn't open mango goodsnp out file $!\n";
	for my $key1 (sort keys %hoh)	{
		for my $key2 (sort keys %{$hoh{$key1}})	{
			$snparraynames{$key1}{$hoh{$key1}{$key2}} = $key2;  #contig=>snp_pos=>snp_name 
			$snparray{$key1} .= "$hoh{$key1}{$key2}\t";# for each contig name, make a single hash entry which is an array of snp pos for that contig
#			print "the following is unsorted snpflank line 78\@snpflank\n";
			}
#		$snparray{$key1} .= "$orglength{$key1}\t";
				$orglen = $orglength{$key1};		
				my @snpflank1 = split /\t/, $snparray{$key1};
#				print "the following is unsorted snpflank line 88\n @snpflank1\n";
				unshift (@snpflank1, 0);
				push (@snpflank1, $orglen);
				@snpflank = sort {$a <=> $b} @snpflank1;
				print "the following is sorted snpflank line 90\n @snpflank\n";
				@returnsnp = &flank100 (@snpflank);
				print "the following is returnsnp at line 93\n@returnsnp\n";
			foreach my $ret (@returnsnp)	{
						if (exists $snparraynames{$key1}{$ret})	{ #if snp pos is same as in returned array put in with goodsnps
						$goodsnp{$key1}{$snparraynames{$key1}{$ret}} = $snparraynames{$key1}{$ret};	#identifies snp name of snp with 100nt flank region
								print "this is the name of a goodsnp $goodsnp{$key1}{$snparraynames{$key1}{$ret}}\n";
								print GOODSNP "$goodsnp{$key1}{$snparraynames{$key1}{$ret}}\n";
					}
				}
			}					
close GOODSNP;
######################################################################################################################
#  sub flank100 identifies snps with 100nt snp free flank regions
sub flank100 {
	my @b = ();
	my $i = ();
	my @c = ();
	my (@a) = @_;
#	print "this is what made it to sub: array a @a\n";
		for ($i = 0; $i <= $#a-1; $i++)	{ #go through the array, calculating whether snp flanks are 100 nt
		  $b[$i] = ($a[$i+1] - $a[$i]);		#snp a0 doesn't have flank of 100 on left side of snp
		} 
#		print "in sub, this is the value of \@b\n@b\n";
		my $bno = @b;
		for ( $i = 0; $i <= $bno-2 ; $i++)	{
				if (($b[$i] >=100)&&($b[$i+1] >= 100))	{
					push (@c, $a[$i+1]);
#					print "in sub this is the  value of $c[$i] at line 130\n";
					}else{
					push (@c, 0);
#					print "in sub, this is a value of c $c[$i] at line 133\n";
		}
	}
#		print "in sub, this is \@c: @c\n";
		my @return = @c;		#return array c
	}
#######################################################################################################################
#   

#########################################################################################################################
#

##########################################################################################################################
#

