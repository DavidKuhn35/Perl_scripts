#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: July 24, 2014
# Version		: ver. 1.1 copied from avo_fixit_May9_2013.pl and modified to get 5050 files of bacon, simmonds, hass, tonnage for sean myles table without removing <6 reads calls
#  edited to calculate SNP type from variant report, determine # of transitions and transversions, # of AC, AG, etc, #of total reads for parent, # average # of reads per SNP per parent, genotype coded for ref nuc, 
# History		: ver. 1.0, April 2013 avo_SNPtype_counter.pl, 
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
my $allhet = 0;
my $AC = 0;
my $AG = 0;
my $AT = 0;
my $CG = 0;
my $CT = 0;
my $GT = 0;
my $tsit = 0;
my $tver = 0;
my %bhomo = ();
my %hhomo = ();
my %shomo = ();
my %thomo = ();
my %bhom11 = ();
my %hhom11 = ();
my %shom11 = ();
my %thom11 = ();
my %bhom22 = ();
my %hhom22 = ();
my %shom22 = ();
my %thom22 = ();
my %allsnp = ();
my $allhom = 0;
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
my $totalbreads = 0;
my %hassreads = ();
my $totalhreads = 0;
my %simmondsreads = ();
my $totalsreads = 0;
my %tonnagereads = ();
my $totaltreads = 00000;
my %bnosnp = ();
my %hnosnp = ();
my %snosnp = ();
my %tnosnp = ();
my %SHRS = ();
my %seqno = ();
my %bseqsm = ();
my %hseqsm = ();
my %sseqsm = ();
my %tseqsm = ();
my %snptype = ();
my %chipsnptype = ();
my %seqsnptype = ();
my %chiprefnuc = ();
my %seqrefnuc = ();

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
			$snptype{$snpname} = &snptype (\@snp); #hash of calculated SNP type
		 	 
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
			$totalbreads = $totalbreads + $snp[5];
			$bacongeno{$snpname} = &genotype (\@bacon);
			$hassgeno{$snpname} = &genotype (\@hass);
			$hassreads{$snpname} = $snp[13];
			$totalhreads = $totalhreads + $snp[13];
	 		$simmondsgeno{$snpname} = &genotype (\@simmonds);
			$simmondsreads{$snpname} = $snp[21];
			$totalsreads = $totalsreads + $snp[21];
	 		$tonnagegeno{$snpname} = &genotype (\@tonnage);
			$tonnagereads{$snpname} = $snp[29];
			$totaltreads = $totaltreads + $snp[29];
			$hoh2{$isotig}{$snpname} = "$snptype{$snpname}\t$refnuc{$snpname}\t$bacongeno{$snpname}\t$baconreads{$snpname}\t$hassgeno{$snpname}\t$hassreads{$snpname}\t$simmondsgeno{$snpname}\t$simmondsreads{$snpname}\t$tonnagegeno{$snpname}\t$tonnagereads{$snpname}";	#hash hoh2 has isotig name key and snpname key for #
#       sorting on both sets of keys, values are all genotypes for bacon, hass, simmonds, tonnage and total reads for B,H,S,T
#	should make a database of the genotypes for each isotig => snpname BUT no Ns or indels in this database.

			}
	}
}

##########################################################################################################################################
#  read in SHRS and snp id for 6000 snps from SHRS_seqno_out_July24_2014.txt
open (CHIP6000, "$ARGV[1]") || die "can't open SHRS snp id 6000 file $!\n";
	while (<CHIP6000>) 	{
		chomp;
		my @x = split (/\t/) ;
		$SHRS{$x[0]} = $x[1]; #hash with shrs name as key and seq id as value
		$seqno{$x[1]} = $x[0];	#hash with seq id as key and shrs name as value
}


#########################################################################################################################################
open	(AVOGENO, ">>$ARGV[2]\_SNPtype_refnuc_seq_genotypes_from_variant_July31_2014.txt") || die "Can't open avo genotypes out file $!\n";
	print AVOGENO "Isotig\tSNPname\tSNPtype\tRefNuc\tBacon genotype\tBacon reads\tHass genotype\tHass reads\tSimmonds genotype\tSimmonds reads\tTonnage genotype\tTonnage reads\n";
	foreach my $key (sort keys %hoh2)	{
		foreach my $snp	(sort keys %{$hoh2{$key}})	{
			print AVOGENO "$key\t$snp\t$seqno{$snp}\t$hoh2{$key}{$snp}\n";
		}
	}

###########################################################################################################################################
#get 6000 bacon seq calls for sean myles table
open (BSEQ, ">>$ARGV[2]\_bacon_seq_6000_sean_myles_v2.txt") || die "can't open bacon 6000 seq outfile $!\n";
	foreach my $key (sort keys %seqno)	{
		if (exists $bacongeno{$key})	{	
			$bseqsm{$seqno{$key}} = $bacongeno{$key}; #sorted by seq no not SHRS no, can't compare to SHRS easily
			print BSEQ "$seqno{$key}\t$snptype{$key}\t$refnuc{$key}\t$bacongeno{$key}\t$baconreads{$key}\n";
}
}
###########################################################################################################################################
#get 6000 hass seq calls for sean myles table
open (HSEQ, ">>$ARGV[2]\_hass_seq_6000_sean_myles_v2.txt") || die "can't open hass 6000 seq outfile $!\n";
	foreach my $key (sort keys %seqno)	{
		if (exists $hassgeno{$key})	{	
			$hseqsm{$seqno{$key}} = $hassgeno{$key};
	print HSEQ "$seqno{$key}\t$snptype{$key}\t$refnuc{$key}\t$hassgeno{$key}\t$hassreads{$key}\n";
}
}
###########################################################################################################################################
#get 6000 simmonds seq calls for sean myles table
open (SSEQ, ">>$ARGV[2]\_simmonds_seq_6000_sean_myles_v2.txt") || die "can't open simmonds 6000 seq outfile $!\n";
	foreach my $key (sort keys %seqno)	{
		if (exists $simmondsgeno{$key})	{	
			$sseqsm{$seqno{$key}} = $simmondsgeno{$key};
	print SSEQ "$seqno{$key}\t$snptype{$key}\t$refnuc{$key}\t$simmondsgeno{$key}\t$simmondsreads{$key}\n";
}
}
###########################################################################################################################################
#get 6000 bacon seq calls for sean myles table
open (TSEQ, ">>$ARGV[2]\_tonnage_seq_6000_sean_myles_v2.txt") || die "can't open tonnage 6000 seq outfile $!\n";
	foreach my $key (sort keys %seqno)	{
		if (exists $tonnagegeno{$key})	{	
			$tseqsm{$seqno{$key}} = $tonnagegeno{$key};
	print TSEQ "$seqno{$key}\t$snptype{$key}\t$refnuc{$key}\t$tonnagegeno{$key}\t$tonnagereads{$key}\n";
}
}

	
############################################################################################################################################
#calculate # of each SNP type from %snptype (AC, AG, AT, CG, CT, GT) and transitions and transversions
foreach my $key (sort keys %snptype)	{
	my $x = $snptype{$key};
	if ($x eq "A\tC")	{
	$AC = $AC + 1;
	$tver = $tver + 1;
		}elsif ($x eq "A\tG")	{
		$AG = $AG + 1;
		$tsit = $tsit + 1;
			} elsif ($x eq "A\tT")	{
			$AT = $AT + 1;
			$tver = $tver + 1;
				} elsif ($x eq "C\tG")	{
				$CG = $CG + 1;
				$tver = $tver + 1;
					} elsif ($x eq "C\tT")	{
					$CT = $CT + 1;
					$tsit = $tsit + 1;
						}elsif ($x eq "G\tT")	{
						$GT = $GT + 1;
						$tver = $tver + 1;
					}
	}
##########################################################################################################################################
# Calculate homozygosity, coded homozygosity and heterozygosity for each parent
foreach my $key (sort keys %bacongeno)	{
	my @x = split /\t/, $bacongeno{$key};
	if (($x[0] eq $x[1]) && ($x[0] eq $refnuc{$key}))	{
	$bhomo{$key} = $key;
	$bhom11{$key} = $key;
		}elsif (($x[0] eq $x[1]) && ($x[0] ne $refnuc{$key}))	{
		$bhomo{$key} = $key;
		$bhom22{$key} = $key;
		}else{
		$bhet{$key} = $key;
	}
}
#######################################################################################################################################
#  	
foreach my $key (sort keys %hassgeno)	{
	my @x = split /\t/, $hassgeno{$key};
	if (($x[0] eq $x[1]) && ($x[0] eq $refnuc{$key}))	{
	$hhomo{$key} = $key;
	$hhom11{$key} = $key;
		}elsif (($x[0] eq $x[1]) && ($x[0] ne $refnuc{$key}))	{
		$hhomo{$key} = $key;
		$hhom22{$key} = $key;
		}else{
		$hhet{$key} = $key;
	}
}
############################################################################################################################
foreach my $key (sort keys %simmondsgeno)	{
	my @x = split /\t/, $simmondsgeno{$key};
	if (($x[0] eq $x[1]) && ($x[0] eq $refnuc{$key}))	{
	$shomo{$key} = $key;
	$shom11{$key} = $key;
		}elsif (($x[0] eq $x[1]) && ($x[0] ne $refnuc{$key}))	{
		$shomo{$key} = $key;
		$shom22{$key} = $key;
		}else{
		$shet{$key} = $key;
	}
}
##########################################################################################################################
foreach my $key (sort keys %tonnagegeno)	{
	my @x = split /\t/, $tonnagegeno{$key};
	if (($x[0] eq $x[1]) && ($x[0] eq $refnuc{$key}))	{
	$thomo{$key} = $key;
	$thom11{$key} = $key;
		}elsif (($x[0] eq $x[1]) && ($x[0] ne $refnuc{$key}))	{
		$thomo{$key} = $key;
		$thom22{$key} = $key;
		}else{
		$thet{$key} = $key;
	}
}
###############################################################################################################################
# calculate number of all hets and all homos for SNPs
foreach my $key (sort keys %allsnp)	{
	if (($bhomo{$key} eq $hhomo{$key}) && ($bhomo{$key} eq $shomo{$key}) && ($bhomo{$key} eq $thomo{$key}))	{
		$allhom = $allhom + 1;
	} elsif ((exists $bhet{$key}) && (exists $hhet{$key}) && (exists $shet{$key}) && (exists $thet{$key}))	{
		$allhet = $allhet + 1;
	}
}
##########################################################################################################################################3
#  subroutine to call and sort genotypes of SNPs
sub genotype	{	#determines genotype of SNP (XX if not enough reads) (homozygous major nucleotide) or heterozygous
	my @x = @{$_[0]};
#	print "This is the subroutine array @x\n";
#	if ($x[4] < 6)	{
#		my $genotype = "X\tX";
	if ($x[5]eq "\.")	{
		my $genotype = "$x[1]\t$x[1]"; #  if ., means that genotype is homozygous for most common allele $x[1]
	} elsif ($x[7]>= 0.80)	{ #if allele freq is >0.8, genotype is homozygous for less comm[on allele
		my $genotype = "$x[5]\t$x[5]"; #only have allele frequency if it is the less common allele
	} elsif	(($x[7] < 0.80) && ($x[10] > 0.20))	{ #conditions for heterozygous, both alleles are stated
		my @z = ($x[5],$x[8]);
		my @zsort = sort @z;
		my $genotype = join "\t", @zsort; #sort heterozygotes ascibetically
	}elsif (($x[7]<= 0.2)&& ($x[10] >= 0.8))	{	 #major allele freq low, minor allele freq high, this may not be possible
#								 but it shouldn't matter
		my $genotype = "$x[8]\t$x[8]" 
	}
}

#############################################################################################################################################
# subroutine to call snptype from variant report for all snps
sub snptype 	{
	my @s = @{$_[0]};
	my $x = 0;
	if ($s[2] ne $s[3])	{
	 $x = "$s[2]\t$s[3]";
		} elsif (($s[2] ne $s[6]) && ($s[6] ne "\."))	{
		$x = "$s[2]\t$s[6]";
			} elsif (($s[2] ne $s[9]) && ($s[9] ne "\.")) {
			$x = "$s[2]\t$s[9]";
				} elsif (($s[2] ne $s[14]) && ($s[14] ne "\.")) {
				$x = "$s[2]\t$s[14]";
					} elsif (($s[2] ne $s[17]) && ($s[17] ne "\.")) {
					$x = "$s[2]\t$s[17]";
						} elsif (($s[2] ne $s[22]) && ($s[22] ne "\.")) {
						$x = "$s[2]\t$s[22]";
							} elsif (($s[2] ne $s[25]) && ($s[25] ne "\.")) {
							$x = "$s[2]\t$s[25]";
								} elsif (($s[2] ne $s[30]) && ($s[30] ne "\.")) {
								$x = "$s[2]\t$s[30]";
									} elsif (($s[2] ne $s[33]) && ($s[33] ne "\.")) {
									$x = "$s[2]\t$s[33]";
										} else {
										$x = "X\tX";
										}
		}
##########################################################################################
# subroutine to give ref allele 1 or alt allele 2 coding
sub code {
    my ( $s, $ref ) = @_;
	if (($s eq "X") || ($s eq "-"))	{ # codes for missing data seq X chip -
	$s = 9;
    } elsif ( $s eq $refnuc{$ref} ) {
        $s = 1;
    }
    else {
        $s = 2;
    }
}	
