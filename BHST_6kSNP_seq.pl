#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: March 7,2013
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: takes unifiedSNPs_Avocado_BLASTX_TAIR10_ANNOT.table_cvg.txt, strips off annotation, adds a unique identifier and pulls out the 6000 SNPs used for the chip.
# <input>		:  	unifiedSNPs_Avocado_BLASTX_TAIR10_ANNOT.table_cvg.txt, avo1_14f_6kSNP_list_for_Illumina.txt (to get SNP #s)	
# <output>		:   BHST_6kSNP_seq.txt
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
my %seqdata = ();
my %chipdata = ();
my %BHST_seq_chip = ();
my $newname = 0;
#####################################################################################################################
#	main	read in datafile
open (DATA,"$ARGV[0]") || die "Can't open input file $!\n";
	while (<DATA>)		{
		chomp;
		if (m/^SEQID/)	{
			my @x = split /\t/;
			my @y = ();
			unless ($x[2] =~ m/^N/)	{
			for (my $i = 0; $i <=39; $i+=1)	{
				push (@y, "$x[$i]");
				}
			$y[0] = "$y[36]\_$y[1]";
			$seqdata{$y[0]} = join "\t", @y;
		}
	}
}			
open (CHIP, "$ARGV[1]") || die "Can't open chip data file $!\n";
	while (<CHIP>)	
{
		chomp;
		if (m/^SHRSPa/)	{
			my @x = split /\t/;
				if ($x[10] =~ m/Bacon/)	{
					$x[11] = "Het Bacon";
				}else{
					$x[11] = "No Het Bacon";
			}
				if ($x[10] =~ m/Hass/)	{
					$x[12] = "Het Hass";
				}else{
					$x[12] = "No Het Hass";
			}
				if ($x[10] =~ m/Simmonds/)	{
					$x[13] = "Het Simmonds";
				}else{
					$x[13] = "No Het Simmonds";
			}
				if ($x[10] =~ m/Tonnage/)	{
					$x[14] = "Het Tonnage";
				}else{
					$x[14] = "No Het Tonnage";
			}
		$chipdata{$x[3]} = join "\t", @x;
	}
}
							
######################################################################################################################
#   sort seqdata and get file with 6k SNPs, seq data and chip data

		foreach my $key (sort keys %chipdata)	{
					if (exists $seqdata{$key})	{
					$BHST_seq_chip{$key} = "$seqdata{$key}\t$chipdata{$key}";
					}
		}	

#######################################################################################################################
#   print out file with BHST chip and seq data
open (OUT, ">>Mar7_BHST_6kSNP_chip_seq_data.txt") || die "Can't open outfile for chip seq data $!\n";
		foreach my $key (sort keys %BHST_seq_chip)	{
			print OUT "$BHST_seq_chip{$key}\n";
		}

#########################################################################################################################
#

##########################################################################################################################
#

