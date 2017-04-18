#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: June 27, 2012
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: 	Analyzes diff_seq files between B,H,S,T chip marker calls and predicted sequence calls
# <input>		:  seq_diff file (output from avo_seq_chip_compare_table.pl		
# <output>		:   stats summary on types of differences (missing data, homo/homo and homo/het differences 
# 			
# Dependencies	: N/A
# Usage			: avo_diff_seq_chip.pl
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
my @wg = ();
my $seqnodata = 0;
my $chipnodata = 0;
my $homohomo = 0;
my $chiphomoseqhet = 0;
my $chiphetseqhomo = 0;
my @seqnodata = ();
my @chipnodata = ();
my @homohomo = ();
my @chiphomoseqhet = ();
my @chiphetseqhomo = ();
my %seqnodatacount = ();
my %chipnodatacount = ();
my %homohomocount = ();
my %chiphomoseqhetcount = ();
my %chiphetseqhomocount = ();
my %alldiff = ();

#####################################################################################################################
#	Read in diff file and score for types of differences
open	(DIFF, "$ARGV[0]") || die "Can't read in diff file $!\n";
	while (<DIFF>)	
{
	chomp;
		if (m/^SHRSPa/)	{	#if line starts with a snp name
		my @diff = split (/\t/)	;
		if ($diff[3] eq "X")	{
			push (@seqnodata, $diff[0]);
			$seqnodata++;
			$seqnodatacount{$diff[0]} = $seqnodatacount{$diff[0]}++;
			$alldiff{$diff[0]}++;
		}	elsif ($diff[1] eq "-")	{
			push (@chipnodata, $diff[0]);
			$chipnodata++;
			$chipnodatacount{$diff[0]} = $chipnodatacount{$diff[0]}++;
			$alldiff{$diff[0]}++;
		} elsif (($diff[1] eq $diff[2]) && ($diff[3] eq $diff[4]))	{
			push (@homohomo, $diff[0]);
			$homohomo++;
			$homohomocount{$diff[0]} = $homohomocount{$diff[0]}++;
			$alldiff{$diff[0]}++;
		} elsif (($diff[1] eq $diff[2]) && ($diff[3] ne $diff[4]))	{
			push (@chiphomoseqhet, $diff[0]);
			$chiphomoseqhet++;
			$chiphomoseqhetcount{$diff[0]} = $chiphomoseqhetcount{$diff[0]}++;
			$alldiff{$diff[0]}++;
		} elsif (($diff[1] ne $diff[2]) && ($diff[3] eq $diff[4]))	{
			push (@chiphetseqhomo, $diff[0]);
			$chiphetseqhomo++;
			$chiphetseqhomocount{$diff[0]} = $chiphetseqhomocount{$diff[0]}++;
			$alldiff{$diff[0]}++;
		}
	}
}
close DIFF;


#########################################################################################################################
#	Print out report on types of differences

open (OUT, ">>2all_diff_types_stats") || die "Can't open diff outfile $!\n";
		my $totaldiff = $chipnodata + $homohomo + $chiphomoseqhet + $chiphetseqhomo;
		my $seqnodatapercent = ($seqnodata/$totaldiff)* 100;
		my $chipnodatapercent = ($chipnodata/$totaldiff)* 100;
		my $homohomopercent = ($homohomo/$totaldiff) * 100;
		my $chiphomoseqhetpercent = ($chiphomoseqhet/$totaldiff) * 100;
		my $chiphetseqhomopercent = ($chiphetseqhomo/$totaldiff) * 100;
		printf OUT "Total differences analyzed\t $totaldiff\nSeq no data\t$seqnodata\t%3.1f%%\nChip no data\t $chipnodata\t%3.1f%%\nHomo homo diffs\t $homohomo\t%3.1f%%\nChip homo seq het\t $chiphomoseqhet\t%3.1f%%\nChip het seq homo\t $chiphetseqhomo\t%3.1f%%\n", $seqnodatapercent,$chipnodatapercent, $homohomopercent, $chiphomoseqhetpercent,$chiphetseqhomopercent;
		
		my @alldiffproblems = sort {	$alldiff{$b} <=> $alldiff{$a}
						or $a cmp $b 	}	%alldiff;
		my $k = keys %alldiff;
		for (my $i=0; $i < $k; $i++)	{
		print OUT "$alldiffproblems[$i]\t$alldiff{$alldiffproblems[$i]}\tall_diff\n";
		}
close OUT;
#############################################################################################################################	

open (OUT, ">>2indiv_diff")	|| die "Can't open individual diff file$!\n";	
		for	(my $i = 0; $i < $seqnodata; $i++)	{
			print OUT "$seqnodata[$i]\tseqnodata\n";
			}
		for 	(my $i = 0; $i < $chipnodata; $i++)		{
			print OUT "$chipnodata[$i]\tchipnodata\n";
			}
		for 	(my $i = 0; $i < $homohomo; $i++)		{
			print OUT "$homohomo[$i]\thomohomo\n";
			}
		for 	(my $i = 0; $i < $chiphomoseqhet; $i++)	{
			print OUT "$chiphomoseqhet[$i]\tchiphomoseqhet\n";
			}
		for 	(my $i = 0; $i < $chiphetseqhomo; $i++)	{
			print OUT "$chiphetseqhomo[$i]\tchiphetseqhomo\n";
			}
		my @seqnodataproblems = sort {	$seqnodatacount{$b} <=> $seqnodatacount{$a}
						or $a cmp $b 	}	keys %seqnodatacount;
		my $h = @seqnodataproblems;
		for (my $i=0; $i < $h; $i++)	{
		print OUT "$seqnodataproblems[$i]\t$seqnodatacount{$seqnodataproblems[$i]}\tseqnodata\n";
		}
		
		my @chipnodataproblems = sort { $chipnodatacount{$b} <=> $chipnodatacount{$a}
						or $a cmp $b	} keys %chipnodatacount;
		my $f = @chipnodataproblems;
		for (my $i=0; $i < $f; $i++)	{
		print OUT "$chipnodataproblems[$i]\t$chipnodatacount{$chipnodataproblems[$i]}\tchipnodata\n";
		}
		
		my @homohomoproblems = sort {	$homohomocount{$b} <=> $homohomocount{$a}
						or $a cmp $b	} keys %homohomocount;
		my $c = @homohomoproblems;
		for (my $i=0; $i < $c; $i++)	{
		print OUT "$homohomoproblems[$i]\t$homohomocount{$homohomoproblems[$i]}\thomohomo_diff\n";
		}
		
		my @chiphomoseqhetproblems = sort {$chiphomoseqhetcount{$b} <=> $chiphomoseqhetcount{$a}
						or $a cmp $b	}	keys %chiphomoseqhetcount;
		my $d = @chiphomoseqhetproblems;
		for (my $i=0; $i < $d; $i++)	{
		print OUT "$chiphomoseqhetproblems[$i]\t$chiphomoseqhetcount{$chiphomoseqhetproblems[$i]}\tchiphomo_seqhet\n";
		}
		
		my @chiphetseqhomoproblems = sort {	$chiphetseqhomocount{$b} <=> $chiphetseqhomocount{$a}
						or $a cmp $b	}	keys %chiphetseqhomocount;
		my $e = @chiphetseqhomoproblems;
		for (my $i=0; $i < $e; $i++)	{
		print OUT "$chiphetseqhomoproblems[$i]\t$chiphetseqhomocount{$chiphetseqhomoproblems[$i]}\tchiphet_seqhomo\n";	
		}
close OUT;


