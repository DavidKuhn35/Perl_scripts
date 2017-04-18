#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: May 3, 2013
# Version		: ver. 1.1
# History		: ver. 1.0, April 2013 avo_SNPtype_counter.pl
# Description	: Counts transitions and transversions for SNPs, prints out file with all SNP types for 660k SNPs
# <input>		:  avocado variant report		
# <output>		:  list of all SNPs and a summary of number of transitions and transversions across all B H S T.  Plus, # of hets of B, H, S, T for a particular isotig
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
#use Statistics::Descriptive;
my %hoh = ();
my %hoh2 = ();
my %isotig = ();
my @snpdistribution = ();
my $snpname = ();
my @bacon = ();
my @hass = ();
my @simmonds = ();
my @tonnage = ();
my %bacongeno = ();
my %hassgeno = ();
my %simmondsgeno = ();
my %tonnagegeno = ();
my %refnuc = ();
my $bacon_nosnp = 0;
my $bacon_hetsnp = 0;
my $bacon_homosnp = 0;
my $hass_nosnp = 0;
my $hass_homosnp = 0;
my $hass_hetsnp = 0;
my $simmonds_nosnp = 0;
my $simmonds_homosnp = 0;
my $simmonds_hetsnp = 0;
my $tonnage_nosnp = 0;
my $tonnage_homosnp = 0;
my $tonnage_hetsnp = 0;
my %allhet = ();
my %AC = ();
my %AG = ();
my %AT = ();
my %CG = ();
my %CT = ();
my %GT = ();
my %transition = ();
my %transversion = ();
my %bhomo = ();
my %hhomo = ();
my %shomo = ();
my %thomo = ();
my %allsnp = ();
my %allhomo = ();
my %bindel = ();
my %hindel = ();
my %sindel = ();
my %tindel = ();
my %all_no_snp = ();
my %homosnpAC = ();
my %homosnpAG = ();
my %homosnpAT = ();
my %homosnpCG = ();
my %homosnpCT = ();
my %homosnpGT = ();
my %homotransition = ();
my %homotransversion = ();
my %N_snp = ();
my %isolength = ();
my %allN = ();
my %allindel = ();
my %bhet = ();
my %hhet = ();
my %shet = ();
my %thet = ();
my %het_transition = ();
my %het_transversion = ();
my %het_AC = ();
my %het_AG = ();
my %het_AT = ();
my %het_CG = ();
my %het_CT = ();
my %het_GT = ();
my %bcov = ();
my %hcov = ();
my %scov = ();
my %tcov = ();
my %annot = ();
my $isotig = ();
my %baconreads = ();
my %hassreads = ();
my %simmondsreads = ();
my %tonnagereads = ();
my %bnosnp = ();
my %hnosnp = ();
my %snosnp = ();
my %tnosnp = ();
########################################################################################################################################
#Read in the avo snp file and get rid of rows with Ns and indels
open (SNP, "$ARGV[0]") || die "Can't open SNP infile";
	while (<SNP>)
{
	chomp;
	if (m/^SEQID/)  	{	#If line has SNP name for avocado
	my @snp = split (/\t/);
	$snpname = $snp[1];	#snpname is snp name
	if (($snp[2] eq "N")||($snp[3] eq "N"))	{	# if most frequent nucleotide or reference nucleotide is N remove from further filtering
		$allN{$snpname} = $snpname;
	}	elsif (($snp[6] eq "-")||($snp[9] eq "-")||($snp[14] eq "-")|| ($snp[17] eq "-") || ($snp[22] eq "-") || ($snp[25] eq "-") || ($snp[30] eq "-") || ($snp[33] eq "-")) 	{
			$allindel{$snpname} = $snpname; #if any nucleotide is an indel "-" remove the variant from further filtering
			$isotig = $snp[36];	#isotig name
			$hoh{$isotig}{$snpname} = $snpname; #make hash of hash with isotig as first key and snpname as second key
	} else {
		 	 @bacon = ($snp[1], $snp[2], $snp[3], $snp[4], $snp[5], $snp[6], $snp[7], $snp[8], $snp[9], $snp[10], $snp[11]);
# 1 is SNP#, 2 most common nuc, 3 hass ref nuc, 4 bacon, 5 # of total reads for bacon, 6 bacon most common nuc, 7 # of reads supporting most common nuc
# 8 major allele freq, 9 minor allele nuc, 10 # of reads supporting minor allele, 11 minor allele frequency

			 @hass = ($snp[1], $snp[2], $snp[3], $snp[12], $snp[13], $snp[14], $snp[15], $snp[16], $snp[17], $snp[18], $snp[19]);
			 @simmonds = ($snp[1], $snp[2], $snp[3], $snp[20], $snp[21], $snp[22], $snp[23], $snp[24], $snp[25], $snp[26], $snp[27]);
			 @tonnage = ($snp[1], $snp[2], $snp[3], $snp[28], $snp[29], $snp[30], $snp[31], $snp[32], $snp[33], $snp[34], $snp[35]);
			$allsnp{$snpname} = $snpname;	#hash of snps without Ns or indels
			$isotig = $snp[36];	#isotig name
			$refnuc{$snpname} = $snp[3];
			$baconreads{$snpname} = $snp[5];
			$bacongeno{$snpname} = &genotype (\@bacon);
			$hassgeno{$snpname} = &genotype (\@hass);
			$hassreads{$snpname} = $snp[13];
	 		$simmondsgeno{$snpname} = &genotype (\@simmonds);
			$simmondsreads{$snpname} = $snp[21];
	 		$tonnagegeno{$snpname} = &genotype (\@tonnage);
			$tonnagereads{$snpname} = $snp[29];
			$hoh2{$isotig}{$snpname} = "$refnuc{$snpname}\t$bacongeno{$snpname}\t$baconreads{$snpname}\t$hassgeno{$snpname}\t$hassreads{$snpname}\t$simmondsgeno{$snpname}\t$simmondsreads{$snpname}\t$tonnagegeno{$snpname}\t$tonnagereads{$snpname}";	#hash hoh2 has isotig name key and snpname key for #
#       sorting on both sets of keys, values are all genotypes for bacon, hass, simmonds, tonnage and total reads for B,H,S,T
#	should make a database of the genotypes for each isotig => snpname BUT no Ns or indels in this database.

			}
	}
}

#########################################################################################################################################
open	(AVOGENO, ">>fixit_avo_May9_2013_refnuc_seq_genotypes_out_July21_2014.txt") || die "Can't open avo genotypes out file $!\n";
	print AVOGENO "Isotig\tSNPname\tRefNuc\tBacon genotype\tBacon reads\tHass genotype\tHass reads\tSimmonds genotype\tSimmonds reads\tTonnage genotype\tTonnage reads\n";
	foreach my $key (sort keys %hoh2)	{
		foreach my $snp	(sort keys %{$hoh2{$key}})	{
			print AVOGENO "$key\t$snp\t$hoh2{$key}{$snp}\n";
		}
	}
###########################################################################################################################################
############################################################################################################################################
##########################################################################################################################################3
#  subroutine to call and sort genotypes of SNPs
sub genotype	{	#determines genotype of SNP (XX if not enough reads) (homozygous major nucleotide) or heterozygous
	my @x = @{$_[0]};
#	print "This is the subroutine array @x\n";
#	if ($x[4] < 6)	{
#		my $genotype = "X\tX";
	if ($x[5]eq "\.")	{
		my $genotype = "$x[1]\t$x[1]";
	} elsif ($x[7]>= 0.80)	{
		my $genotype = "$x[5]\t$x[5]";
	} elsif	(($x[7] < 0.80) && ($x[10] > 0.20))	{
		my @z = ($x[5],$x[8]);
		my @zsort = sort @z;
		my $genotype = join "\t", @zsort;
	}
}



