#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: 6/25/2012
# Version		: ver. 1.1  6/28/2012
# History		: ver. 1.0
# Description	: 	Compare genotypes for mapping population parents predicted from sequencing data and from Illumina SNP chip
# <input>		:  DNA manifest for SNP chip, Illumina SNP data w/ SNPs with leading zeros, SNP data from sequencing			
# <output>		:  Table comparing %identity between SNP chip and sequencing genotypes for each parent 
# 			
# Dependencies	: N/A
# Usage			: avo_seq_chip_compare.pl
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
my %man = ();
my %chip = ();
my %lookup = ();
my %chiplookup = ();
my %revlookup = ();
my %revchiplookup = ();
my %bacon = ();
my %hass = ();
my %simmonds = ();
my %tonnage = ();
my %baconreads = ();
my %hassreads = ();
my %simmondsreads = ();
my %tonnagereads = ();
my %baconchip = ();
my %hasschip = ();
my %simmondschip = ();
my %tonnagechip = ();
my %hasswg = ();
my %baconwg = ();
my %simmondswg = ();
my %tonnagewg = ();
#my %hasschip = ();
#my %baconchip = ();
#my %simmondschip = ();
#my %tonnagechip = ();
my %hasschipseqident = ();
my %baconchipseqident = ();
my %simmondschipseqident = ();
my %tonnagechipseqident = ();
my %hasschipseqdiff = ();
my %baconchipseqdiff = ();
my %simmondschipseqdiff = ();
my %tonnagechipseqdiff = ();
my $hasschipseqident = 0;
my $baconchipseqident = 0;
my $simmondschipseqident = 0;
my $tonnagechipseqident = 0;
my $hasschipseqdiff = 0;
my $baconchipseqdiff = 0;
my $simmondschipseqdiff = 0;
my $tonnagechipseqdiff = 0;
my @b = ();
my @h = ();
my @s = ();
my @t = ();
my @baconwg = ();
my @hasswg = ();
my @simmondswg = ();
my @tonnagewg = ();
my %snp = ();
my %wg = ();
my $snpname = "Pafubar";
my %bseq = ();
my %hseq = ();
my %sseq = ();
my %tseq = ();

############################################################################################
#Load in data from avocado DNA manifest file DNA_Manifest_USDA_ARS_Production_2012Apr.txt
open (AVO_MANIFEST, "$ARGV[0]") || die "Can't open DNA manifest file for avocado SNP chip $!\n"; #load in DNA manifest for avocado SNP chip, make %man hash with well number as key 
#			and name of accession as value
	while (<AVO_MANIFEST>)	
{
	chomp;
		my @x = split (/\t/);
			if ($x[3] =~ m/^WG/)	{
			$man{$x[3]} = $x[6];	#$x[3] is well number, $x[6] is comments or common name
		}
}
###################################################################################################
#Load in SNP chip data from USDA_Avocado_Apr2012_gtReport_PASS.csv				 
open 	(CHIP, "$ARGV[1]") || die "Can't open Illumina chip data infile $!\n";	
	while (<CHIP>)	
{
	chomp;
		if (m/^SHRS/)	{
			my @y = split (/,/);
			my @geno = ($y[2], $y[3]);
			my @sorted = sort (@geno);
			my $genotype = "$sorted[0]\t$sorted[1]";
			$chip{$y[1]}{$y[0]} = $genotype; #$y[1] is WG number, $y[0] is SHRS snp number
			$snp{$y[0]} = $y[0];	#Makes a hash of new snp names from chip
			$wg{$y[1]} = $y[1];	#Makes a hash of wg numbers
			}
}

#######################################################################################################
#  Load in SNP number converter Avo1_14f_6kSNP_list_for_ILLUMINA.txt, add leading zeros to old SNP names, make hash with keys as old snp names and values as new snp names
open 	(LOOKUP, "$ARGV[2]") || die "Can't open SNP name converter file $!\n";
	while (<LOOKUP>)
{
	chomp;
		if (m/^SHRS/)	{
			my @z = split (/\t/);
			$z[4] = $z[4] + 100000000;
			$lookup{$z[4]} = $z[0];	#$z[0] is new SNP name SHRSPa and $z[4] is the SNP position which had 100000000 added to it to generate leading zeros
			$revlookup{$z[0]} = $z[4];
		}
}
##############################################################################################################
	foreach my $x (sort keys %snp)	{	#make a lookup table only for snps on the chip!
		if (exists $revlookup{$x})	{
			$chiplookup{$revlookup{$x}} = $x; #chip hash lookup lz# and get back new snp name lz name => snpname
			$revchiplookup{$x} = $revlookup{$x}; #snpname => lz name
		}
	}
###############################################################################################################
#	Load in snp calls from sequencing with snps with leading zeros and create a hash for each parent with sorted genotypes by old snp name
#	PPP_avo_genotypes_out
open	(SEQ, "$ARGV[3]") ||	die "Can't open sequencing snp genotype data infile $!\n";
	while (<SEQ>)
{
	chomp;
		if (m/^g/)	{
			my @q = split (/\t/);
			my $x = $q[1];	#$q[1] is snp name with leading zeros
			if (exists $chiplookup{$x} )	{	#to get rid of run error perhaps "unitialized value for $snpname used in hash
			$snpname = $chiplookup{$x};
				if ($q[2] ne "X")	{ #puts new SNP name with bacon sequencing genotype if genotype is NOT a nosnp
				$bacon{$snpname} = "$q[2]\t$q[3]";
				$baconreads{$snpname} = $q[4];
				}  
				if ($q[5] ne "X") 	{	
				$hass{$snpname} = "$q[5]\t$q[6]";
				$hassreads{$snpname} = $q[7];
				}
				if ($q[8] ne "X")	{
				$simmonds{$snpname} = "$q[8]\t$q[9]";
				$simmondsreads{$snpname} = $q[10];
				}
				if ($q[11] ne "X")	{
				$tonnage{$snpname} = "$q[11]\t$q[12]";
				$tonnagereads{$snpname} = $q[13];
			}
		}
	}
}
##################################################################################################################
#  Print out sequencing genotype files
open (B_OUT, ">>2bacon_seq_geno_out") || die "Can't open sequencing genotype outfile $!\n";
		foreach my $key (sort keys %bacon)	{
			print B_OUT "$key\t$bacon{$key}\t$baconreads{$key}\tbacon_seq\n";
		}

open (H_OUT, ">>2hass_seq_geno_out") || die "Can't open sequencing genotype outfile $!\n";
		foreach my $key (sort keys %hass)	{
			print H_OUT "$key\t$hass{$key}\t$hassreads{$key}\thass_seq\n";
		}

open (S_OUT, ">>2simmonds_seq_geno_out") || die "Can't open sequencing genotype outfile $!\n";
		foreach my $key (sort keys %simmonds)	{
			print S_OUT "$key\t$simmonds{$key}\tsimmonds_seq\n";
		}
open (T_OUT, ">>2tonnage_seq_geno_out") || die "Can't open sequencing genotype outfile $!\n";
		foreach my $key (sort keys %tonnage)	{
			print T_OUT "$key\t$tonnage{$key}\ttonnage_seq\n";
		}
###########################################################################################################################
#  Make hash for bacon, hass, simmonds,and tonnage genotype data from snp chip, compare to seq genotype

	foreach my $key (sort keys %man)	{
		if ($man{$key} =~ /Control_Hass/)	{	#get all wg that are Hass controls and put them into arrays
		push @hasswg, $key;
	} elsif ($man{$key} =~ /Control_Bacon/)		{
		push @baconwg, $key;
	} elsif ($man{$key} =~ /Control_Simmonds/)	{
		push @simmondswg, $key;
	} elsif ($man{$key} =~ /Control_Tonnage/)	{
		push @tonnagewg, $key;
	}
}

#############################################################33
		foreach my $snp (sort keys %{$chip{$hasswg[0]}})	{
			$hasschip{$snp} = $chip{$hasswg[0]}{$snp};	#fill in genotypes from chip for each Hass control, make a hasschip hash
		}
		foreach my $snp (sort keys %{$chip{$baconwg[0]}})	{
			$baconchip{$snp} = $chip{$baconwg[0]}{$snp};
		}
		foreach my $snp (sort keys %{$chip{$simmondswg[0]}})	{
			$simmondschip{$snp} = $chip{$simmondswg[0]}{$snp};
		}
		foreach my $snp (sort keys %{$chip{$tonnagewg[0]}})	{
			$tonnagechip{$snp} = $chip{$tonnagewg[0]}{$snp};
		}
###########################################################################################################
#  Print out control genotype files for all the parents
open (BCONT, ">>2BBB_b_chip_geno_out") || die "Can't open bacon chip geno out file $!\n";
		foreach my $key (sort keys %baconchip)	{
			print BCONT "$key\t$baconchip{$key}\tbacon_chip\n";
		}
close BCONT;
open (HCONT, ">>2HHH_h_chip_geno_out") || die "Can't open hass chip geno out file $!\n";
		foreach my $key (sort keys %hasschip)	{
			print HCONT "$key\t$hasschip{$key}\thass_chip\n";
			}
close HCONT;

open (SCONT, ">>2SSS_s_chip_geno_out") || die "Can't open simmonds chip geno out file $!\n";
		foreach my $key (sort keys %simmondschip)	{
			print SCONT "$key\t$simmondschip{$key}\tsimmonds_chip\n";
			}
close SCONT;

open (TCONT, ">>2TTT_t_chip_geno_out") || die "Can't open tonnage chip geno out file $!\n";
		foreach my $key (sort keys %tonnagechip)	{
			print TCONT "$key\t$tonnagechip{$key}\ttonnage_chip\n";
			}
close TCONT;

####################################################################
#	compare chip control genotypes to sequencing genotype
	foreach my $snp (sort keys %baconchip)	{
				if (exists $bacon{$snp})	{
					if ($bacon{$snp} eq $baconchip{$snp})	{
					$baconchipseqident++;
					$baconchipseqident{$snp}= $baconchip{$snp}; 
					} else {
					$baconchipseqdiff++;
					$baconchipseqdiff{$snp} = $baconchip{$snp};
				}
			}
		}

	foreach my $snp (sort keys %hasschip)	{
				if (exists $hass{$snp})	{
					if ($hass{$snp} eq $hasschip{$snp})	{
					$hasschipseqident++;
					$hasschipseqident{$snp}= $hasschip{$snp}; 
					} else {
					$hasschipseqdiff++;
					$hasschipseqdiff{$snp} = $hasschip{$snp};
				}
			}
		}

	foreach my $snp (sort keys %simmondschip)	{
				if (exists $simmonds{$snp})	{
					if ($simmonds{$snp} eq $simmondschip{$snp})	{
					$simmondschipseqident++;
					$simmondschipseqident{$snp}= $simmondschip{$snp}; 
					} else {
					$simmondschipseqdiff++;
					$simmondschipseqdiff{$snp} = $simmondschip{$snp};
				}
			}
		}

	foreach my $snp (sort keys %tonnagechip)	{
				if (exists $tonnage{$snp})	{
					if ($tonnage{$snp} eq $tonnagechip{$snp})	{
					$tonnagechipseqident++;
					$tonnagechipseqident{$snp}= $tonnagechip{$snp}; 
					} else {
					$tonnagechipseqdiff++;
					$tonnagechipseqdiff{$snp} = $tonnagechip{$snp};
				}
			}
		}	

###################################################################################################################
#  Print out report on identity between genotypes of parents predicted by sequencing and scored by chip

open (OUTFILE, ">>2RRR_chip_seq_compare_out") || die "Can't open chip seq compare outfile $!\n";
	my $hasstotal = $hasschipseqident + $hasschipseqdiff;
	my $bacontotal = $baconchipseqident + $baconchipseqdiff;
	my $simmondstotal = $simmondschipseqident + $simmondschipseqdiff;
	my $tonnagetotal = $tonnagechipseqident + $tonnagechipseqdiff;

	print OUTFILE "Bacon total genotypes $bacontotal,seq and chip identical $baconchipseqident, seq and chip different $baconchipseqdiff\nHass total genotypes $hasstotal, seq and chip identical $hasschipseqident, seq and chip different $hasschipseqdiff\nSimmonds total genotypes $simmondstotal, seq and chip identical $simmondschipseqident, seq and chip different $simmondschipseqdiff\nTonnage total genotypes $tonnagetotal, seq and chip identical $tonnagechipseqident, seq and chip different $tonnagechipseqdiff\nSNP name\tchip genotype\tseq genotype\treads for seq genotype\tparent\n";

	foreach my $key (sort keys %baconchipseqdiff)	{
			print OUTFILE "$key\t$baconchipseqdiff{$key}\t$bacon{$key}\t$baconreads{$key}\tbacon\n";
		}

	foreach my $key (sort keys %hasschipseqdiff)	{
			print OUTFILE "$key\t$hasschipseqdiff{$key}\t$hass{$key}\t$hassreads{$key}\thass\n";
		}

	foreach my $key (sort keys %simmondschipseqdiff)	{
			print OUTFILE "$key\t$simmondschipseqdiff{$key}\t$simmonds{$key}\t$simmondsreads{$key}\tsimmonds\n";
		}

	foreach my $key (sort keys %tonnagechipseqdiff)	{
			print OUTFILE "$key\t$tonnagechipseqdiff{$key}\t$tonnage{$key}\t$tonnagereads{$key}\ttonnage\n";
		}


##############################################################################################################
#

