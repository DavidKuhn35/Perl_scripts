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
use Statistics::Descriptive;
my %hoh = ();
my %isotig = ();
my @snpdistribution = ();
my $isoname = ();
my @bacon = ();
my @hass = ();
my @simmonds = ();
my @tonnage = ();
my %bacongeno = ();
my %hassgeno = ();
my %simmondsgeno = ();
my %tonnagegeno = ();
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
my %breads = ();
my %hreads = ();
my %sreads = ();
my %treads = ();
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
	$isoname = $snp[1];	#isoname is snp name
	if (($snp[2] eq "N")||($snp[3] eq "N"))	{	# if most frequent nucleotide or reference nucleotide is N remove from further filtering
		$allN{$isoname} = $isoname;
	}	elsif (($snp[6] eq "-")||($snp[9] eq "-")||($snp[14] eq "-")|| ($snp[17] eq "-") || ($snp[22] eq "-") || ($snp[25] eq "-") || ($snp[30] eq "-") || ($snp[33] eq "-")) 	{
			$allindel{$isoname} = $isoname; #if any nucleotide is an indel "-" remove the variant from further filtering
			$isotig = $snp[36];	#isotig name
			$hoh{$isotig}{$isoname} = $isoname; #make hash of hash with isotig as first key and snpname as second key
	} else {
		 @bacon = ($snp[1], $snp[2], $snp[3], $snp[4], $snp[5], $snp[6], $snp[7], $snp[8], $snp[9], $snp[10], $snp[11]);
		 @hass = ($snp[1], $snp[2], $snp[3], $snp[12], $snp[13], $snp[14], $snp[15], $snp[16], $snp[17], $snp[18], $snp[19]);
		 @simmonds = ($snp[1], $snp[2], $snp[3], $snp[20], $snp[21], $snp[22], $snp[23], $snp[24], $snp[25], $snp[26], $snp[27]);
		 @tonnage = ($snp[1], $snp[2], $snp[3], $snp[28], $snp[29], $snp[30], $snp[31], $snp[32], $snp[33], $snp[34], $snp[35]);
	 		$bacongeno{$isoname} = &genotype (\@bacon);
			$hassgeno{$isoname} = &genotype (\@hass);
	 		$simmondsgeno{$isoname} = &genotype (\@simmonds);
	 		$tonnagegeno{$isoname} = &genotype (\@tonnage);
			$allsnp{$isoname} = $isoname;	#hash of snps without Ns or indels
			$isotig = $snp[36];	#isotig name
			$hoh{$isotig}{$isoname} = $isoname; #make hash of hash with isotig as first key and snpname as second key
			}
	}
}
################################################################################################################################################3
#   Get genotypes for each cultivar and determine if they are homozygous or heterozygous for the SNP
	foreach my $key (sort keys %bacongeno)	{
			(my $z1,my $z2) = split /\t/, $bacongeno{$key};
			if ($z1 eq "X") {
				$bnosnp{$key} = $bacongeno{$key};
			} elsif ($z1 eq $z2)	{
				$bhomo{$key} = $bacongeno{$key}; #genotype for snp is homozygous
			}else	{
				$bhet{$key} = $bacongeno{$key};	#genotype for snp is heterozygous
				$allhet{$key} = $bacongeno{$key};	#store genotype in all heterozygotes hash	
			}
} 

	foreach my $key (sort keys %hassgeno)	{
			(my $z1,my $z2) = split /\t/, $hassgeno{$key};
			if ($z1 eq "X") {
				$hnosnp{$key} = $hassgeno{$key};
			} elsif ($z1 eq $z2)	{
				$hhomo{$key} = $hassgeno{$key};
			}else	{
				$hhet{$key} = $hassgeno{$key};
				$allhet{$key} = $hassgeno{$key};
			}
} 
 
	foreach my $key (sort keys %simmondsgeno)	{
			(my $z1,my $z2) = split /\t/, $simmondsgeno{$key};
			if ($z1 eq "X") {
				$snosnp{$key} = $simmondsgeno{$key};
			} elsif ($z1 eq $z2)	{
				$shomo{$key} = $simmondsgeno{$key};
			}else	{
				$shet{$key} = $simmondsgeno{$key};
				$allhet{$key} = $simmondsgeno{$key};
			}
} 

foreach my $key (sort keys %tonnagegeno)	{
			(my $z1,my $z2) = split /\t/, $tonnagegeno{$key};
			if ($z1 eq "X") {
				$tnosnp{$key} = $tonnagegeno{$key};
			 }elsif ($z1 eq $z2)	{
				$thomo{$key} = $tonnagegeno{$key};
			}else	{
				$thet{$key} = $tonnagegeno{$key};
				$allhet{$key} = $tonnagegeno{$key};
			}
} 
			
#################################################################################################################################################
#  Determine type of SNPs (transition, transversion and actual type for all SNPs) using %allhet

		foreach my $key (sort keys %allhet)	{
			if ($allhet{$key} eq "A\tC")	{
				$het_AC{$key} = $key;
				$het_transversion{$key} = $key;
			}elsif ( $allhet{$key} eq "A\tG")	{
				$het_AG{$key} = $key;
				$het_transition{$key} = $key;
			}elsif ( $allhet{$key} eq "A\tT")	{
				$het_AT{$key} = $key;
				$het_transversion{$key} = $key;	
			}elsif ( $allhet{$key} eq "C\tG")	{
				$het_CG{$key} = $key;
				$het_transversion{$key} = $key;
			}elsif ( $allhet{$key} eq "C\tT")	{
				$het_CT{$key} = $key;
				$het_transition{$key} = $key;
			}elsif ( $allhet{$key} eq "G\tT")	{
				$het_GT{$key} = $key;
				$het_transversion{$key} = $key;	
	}
}

#####################################################################################################################		
#   Determine the number of homozygous SNPs (all parents homozygous) and make a hash
		foreach my $key (sort keys %allsnp)		{
			if ((exists $bhomo{$key}) && (exists $hhomo{$key}) && (exists $shomo{$key}) && (exists $thomo{$key}))	{
				$allhomo{$key} = "$bhomo{$key}\t$hhomo{$key}\t$shomo{$key}\t$thomo{$key}";
			}
		}
		foreach my $key (sort keys %allhomo)	{
				if (($allhomo{$key} =~ m/A/) && ($allhomo{$key} =~ m/C/))	{
				$homosnpAC{$key} = $key;
				$homotransversion{$key} = $key;
			}elsif (($allhomo{$key} =~ m/A/) && ( $allhomo{$key} =~ m/G/))	{
				$homosnpAG{$key} = $key;
				$homotransition{$key} = $key;
			}elsif  (($allhomo{$key} =~ m/A/) && ($allhomo{$key} =~ m/T/))	{
				$homosnpAT{$key} = $key;
				$homotransversion{$key} = $key;
			}elsif (($allhomo{$key} =~ m/C/) && ( $allhomo{$key} =~ m/G/))	{
				$homosnpCG{$key} = $key;
				$homotransversion{$key} = $key;
			}elsif (($allhomo{$key} =~ m/C/) && ($allhomo{$key} =~ m/T/))	{
				$homosnpCT{$key} = $key;
				$homotransition{$key} = $key;
			}elsif (($allhomo{$key} =~ m/G/) && ($allhomo{$key} =~ m/T/))	{
				$homosnpGT{$key} = $key;
				$homotransversion{$key} = $key;
			}
}	

#			my $totalhomo = keys %homotransition + keys %homotransversion;

################################################################################################################################
#  Read in coverage files for B H S T
open (B, "/home/david/avocado_Illumina_data/avo5_isotig_cvg_snp.stats") || die "Can't open Bacon coverage file $!\n";
	while (<B>)		
{
	chomp;
	if (m/^g/)	{
	my @b = split /\t/;
	$isolength{$b[0]} = $b[1];
	$bcov{$b[0]} = $b[2];
	$breads{$b[0]} = $b[5];
	}
}
close B;

open (H, "/home/david/avocado_Illumina_data/avo6_isotig_cvg_snp.stats") || die "Can't open Hass coverage file $!\n";
	while (<H>)		
{
	chomp;
	if (m/^g/)	{
	my @h = split /\t/;
	$hcov{$h[0]} = $h[2];
	$hreads{$h[0]} = $h[5];
	}
}
close H;
open (S, "/home/david/avocado_Illumina_data/avo7_isotig_cvg_snp.stats") || die "Can't open Simmonds coverage file $!\n";
	while (<S>)		
{
	chomp;
	if (m/^g/)	{
	my @s = split /\t/;
	$scov{$s[0]}= $s[2];
	$sreads{$s[0]} = $s[5];
	}
}
close S;

open (T, "/home/david/avocado_Illumina_data/avo8_isotig_cvg_snp.stats") || die "Can't open Tonnage coverage file $!\n";
	while (<T>)		
{
	chomp;
	if (m/^g/)	{
	my @t = split /\t/;
	$tcov{$t[0]} = $t[2];
	$treads{$t[0]} = $t[5];
	}
}
close T;

############################################################################################################################################
#	Read in annotation for isotigs from Mar7_avo_reformat_CFPcodes.txt
open (ANNOT, "/home/david/avocado/SNPs/April2_2013_avo_reformat_CFPcodes.txt") || die "Can't open annotation file $!\n";
	while	(<ANNOT>)	
{
	chomp;
	if (m/^g/)	{
		my @y = split /\t/;
		my $annot_isotig = $y[0];
		my $APVO = pop @y;
		my $singlecopy = pop @y;
		pop @y;
		shift @y;
		shift @y;
		shift @y;
		my $annotation = join "\t", @y, $singlecopy, $APVO;
		$annot{$annot_isotig} = $annotation;
	}
}
close ANNOT;
###########################################################################################################################
#   Go back through each isotig and count up each type of SNP in each isotig
	
	foreach my $key (sort keys %annot)	{  #Should have all 33,957 isotigs
		unless (exists $hoh{$key})	{
		$hoh{$key}{IND} = 0;
		$hoh{$key}{het_AC} = 0;
		$hoh{$key}{het_AG} = 0;
		$hoh{$key}{het_AT} = 0;
		$hoh{$key}{het_CG} = 0;
		$hoh{$key}{het_CT} = 0;
		$hoh{$key}{het_GT} = 0;
		$hoh{$key}{homoAC} = 0;
		$hoh{$key}{homoAG} = 0;
		$hoh{$key}{homoAT} = 0;
		$hoh{$key}{homoCG} = 0;
		$hoh{$key}{homoCT} = 0;
		$hoh{$key}{homoGT} = 0;
		}else{
		$hoh{$key}{IND} = 0;
		$hoh{$key}{het_AC} = 0;
		$hoh{$key}{het_AG} = 0;
		$hoh{$key}{het_AT} = 0;
		$hoh{$key}{het_CG} = 0;
		$hoh{$key}{het_CT} = 0;
		$hoh{$key}{het_GT} = 0;
		$hoh{$key}{homoAC} = 0;
		$hoh{$key}{homoAG} = 0;
		$hoh{$key}{homoAT} = 0;
		$hoh{$key}{homoCG} = 0;
		$hoh{$key}{homoCT} = 0;
		$hoh{$key}{homoGT} = 0;
		foreach my $x (sort keys %{$hoh{$key}})	{
			if (exists $allindel{$x})	{
				$hoh{$key}{IND} = $hoh{$key}{IND} + 1;
			}elsif (exists $homosnpAC{$x})	{
				$hoh{$key}{homoAC} = $hoh{$key}{homoAC} + 1;
			}elsif (exists $homosnpAG{$x})	{
				$hoh{$key}{homoAG} = $hoh{$key}{homoAG} + 1;
			}elsif (exists $homosnpAT{$x})  {
				$hoh{$key}{homoAT} = $hoh{$key}{homoAT} + 1;
			}elsif (exists $homosnpCG{$x}) 	{
				$hoh{$key}{homoCG} = $hoh{$key}{homoCG} + 1;
			}elsif (exists $homosnpCT{$x}) {
				$hoh{$key}{homoCT} = $hoh{$key}{homoCT} + 1;
			}elsif (exists $homosnpGT{$x}) {
				$hoh{$key}{homoGT} = $hoh{$key}{homoGT} + 1;
			}elsif (exists $het_AC{$x})	{
				$hoh{$key}{het_AC} = $hoh{$key}{het_AC} + 1;
			}elsif (exists $het_AG{$x})	{
				$hoh{$key}{het_AG} = $hoh{$key}{het_AG} + 1;
			}elsif (exists $het_AT{$x})  {
				$hoh{$key}{het_AT} = $hoh{$key}{het_AT} + 1;
			}elsif (exists $het_CG{$x}) 	{
				$hoh{$key}{het_CG} = $hoh{$key}{het_CG} + 1;
			}elsif (exists $het_CT{$x}) {
				$hoh{$key}{het_CT} = $hoh{$key}{het_CT} + 1;
			}elsif (exists $het_GT{$x}) {
				$hoh{$key}{het_GT} = $hoh{$key}{het_GT} + 1;
				}
			}
		}
}

			my $totalN = keys %allN;
			my $totalindel = keys %allindel;
			my $totalhets = keys %allhet;
			my $totalhomo = keys %allhomo;
			my $totalsnp = keys %allsnp;
##############################################################################################################################################
#  Calculate het SNPs homo SNPs and no SNPs for B H S T
		my $bhettotal = keys %bhet;
		my $bhomototal = keys %bhomo;
		my $bnosnptotal = keys %bnosnp;

		my $hhettotal = keys %hhet;
		my $hhomototal = keys %hhomo;
		my $hnosnptotal = keys %hnosnp;

		my $shettotal = keys %shet;
		my $shomototal = keys %shomo;
		my $snosnptotal = keys %snosnp;

		my $thettotal = keys %thet;
		my $thomototal = keys %thomo;
		my $tnosnptotal = keys %tnosnp;

		my $bsnptotal = $bhettotal + $bhomototal + $bnosnptotal ;

		my $hsnptotal = $hhettotal + $hhomototal + $hnosnptotal;

		my $ssnptotal = $shettotal + $shomototal + $snosnptotal;

		my $tsnptotal = $thettotal + $thomototal + $tnosnptotal;

		my $bperhet = $bhettotal/$bsnptotal;
		my $bperhomo = $bhomototal/$bsnptotal;
		my $bpernosnp = $bnosnptotal/$bsnptotal;

		my $hperhet = $hhettotal/$hsnptotal;
		my $hperhomo = $hhomototal/$hsnptotal;
		my $hpernosnp = $hnosnptotal/$hsnptotal;

		my $sperhet = $shettotal/$ssnptotal;
		my $sperhomo = $shomototal/$ssnptotal;
		my $spernosnp = $snosnptotal/$ssnptotal;

		my $tperhet = $thettotal/$tsnptotal;
		my $tperhomo = $thomototal/$tsnptotal;
		my $tpernosnp = $tnosnptotal/$tsnptotal;

open (HETOUT, ">>avo_SNPtype_het_BHST_counter_May3_out_check.txt") || die "Can't open HETOUT $!\n";
		printf HETOUT "%20s\t" x 8 . "\n", "Cultivar","SNP total", "Het total","Het percent", "Homo total", "Homo percent", "No SNP total", "No SNP percent";
		printf HETOUT "%20s\t" . "%20d\t" . ("%20d\t%20.3f" x 3). "\n", "Bacon", "$bsnptotal","$bhettotal","$bperhet","$bhomototal","$bperhomo","$bnosnptotal","$bpernosnp";   
		printf HETOUT "%20s\t" . "%20d\t" . ("%20d\t%20.3f" x 3). "\n", "Hass", "$hsnptotal","$hhettotal","$hperhet","$hhomototal","$hperhomo","$hnosnptotal","$hpernosnp"; 
		printf HETOUT "%20s\t" . "%20d\t" . ("%20d\t%20.3f" x 3). "\n", "Simmonds", "$ssnptotal","$shettotal","$sperhet","$shomototal","$sperhomo","$snosnptotal","$spernosnp"; 
		printf HETOUT "%20s\t" . "%20d\t" . ("%20d\t%20.3f" x 3). "\n", "Tonnage", "$tsnptotal","$thettotal","$tperhet","$thomototal","$tperhomo","$tnosnptotal","$tpernosnp"; 
close HETOUT;
#############################################################################################################################################
#  Print out datafile with all columns
	
open (BIGOUT, ">>avo_SNPtype_counter_May7_out_noSNP.txt") || die "Can't open data outfile $!\n";
	print BIGOUT "Total Ns\tTotal Indels\tTotal SNPs\tTotal het SNPs\tTotal homo SNPs\n$totalN\t$totalindel\t$totalsnp\t$totalhets\t$totalhomo\n";
	print BIGOUT "isotig name\tisotig length\tbacon cov length\tbacon mean reads\thass covered length\thass mean reads\tsimmonds covered length\tsimmonds mean reads\ttonnage covered length\ttonnage mean reads\thet AC\thet AG\thet AT\thet CG\thet CT\thet GT\thomo AC\thomo AG\thomo AT\thomo CG\thomo CT\thomo GT\tNumber of indels\tC annotation\tC code\tF annotation\tF code\tP annotation\tP code\tsinglecopy status\tAPVO status\n";
	foreach my $key (sort keys %hoh)		{
		print BIGOUT "$key\t$isolength{$key}\t$bcov{$key}\t$breads{$key}\t$hcov{$key}\t$hreads{$key}\t$scov{$key}\t$sreads{$key}\t$tcov{$key}\t$treads{$key}\t$hoh{$key}{het_AC}\t$hoh{$key}{het_AG}\t$hoh{$key}{het_AT}\t$hoh{$key}{het_CG}\t$hoh{$key}{het_CT}\t$hoh{$key}{het_GT}\t$hoh{$key}{homoAC}\t$hoh{$key}{homoAG}\t$hoh{$key}{homoAT}\t$hoh{$key}{homoCG}\t$hoh{$key}{homoCT}\t$hoh{$key}{homoGT}\t$hoh{$key}{IND}\t$annot{$key}\n";
	}
close BIGOUT;	

##########################################################################################################################################3
#  subroutine to call and sort genotypes of SNPs
sub genotype	{	#determines genotype of SNP (XX if not enough reads) (homozygous major nucleotide) or heterozygous
	my @x = @{$_[0]};
#	print "This is the subroutine array @x\n";
	if ($x[4] < 6)	{
		my $genotype = "X\tX";
	} elsif ($x[5]eq "\.")	{
		my $genotype = "$x[1]\t$x[1]";
	} elsif ($x[7]>= 0.80)	{
		my $genotype = "$x[5]\t$x[5]";
	} elsif	(($x[7] < 0.80) && ($x[10] > 0.20))	{
		my @z = ($x[5],$x[8]);
		my @zsort = sort @z;
		my $genotype = join "\t", @zsort;
	}
}



