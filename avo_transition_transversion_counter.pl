#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: March 26, 2013
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: Counts transitions and transversions for SNPs
# <input>		:  avocado variant report		
# <output>		:  list of all SNPs and a summary of number of transitions and transversions 
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
my $goodsnp = 0;
my $badsnp = 0;
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
my %bno_snp = ();
my %hno_snp = ();
my %sno_snp = ();
my %tno_snp = ();
my %all_no_snp = ();
my $homosnpAC = 0;
my $homosnpAG = 0;
my $homosnpAT = 0;
my $homosnpCG = 0;
my $homosnpCT = 0;
my $homosnpGT = 0;
my $homotransition = 0;
my $homotransversion = 0;

#Read in the avo snp file and get rid of rows with missing data or Ns
open (SNP, "$ARGV[0]") || die "Can't open SNP infile";
	while (<SNP>)
{
	chomp;
	if (m/^SEQID/)  	{	#If line has SNP name with Tcm
	my @snp = split (/\t/);
	if (($snp[2] eq "N")||($snp[3] eq "N"))	{
		my $badsnp++;
	}	else	{
		my $goodsnp++;
		 $isoname = $snp[1];
		 @bacon = ($snp[1], $snp[2], $snp[3], $snp[4], $snp[5], $snp[6], $snp[7], $snp[8], $snp[9], $snp[10], $snp[11]);
		 @hass = ($snp[1], $snp[2], $snp[3], $snp[12], $snp[13], $snp[14], $snp[15], $snp[16], $snp[17], $snp[18], $snp[19]);
		 @simmonds = ($snp[1], $snp[2], $snp[3], $snp[20], $snp[21], $snp[22], $snp[23], $snp[24], $snp[25], $snp[26], $snp[27]);
		 @tonnage = ($snp[1], $snp[2], $snp[3], $snp[28], $snp[29], $snp[30], $snp[31], $snp[32], $snp[33], $snp[34], $snp[35]);
	 $bacongeno{$isoname} = &genotype (\@bacon);
	$hassgeno{$isoname} = &genotype (\@hass);
	 $simmondsgeno{$isoname} = &genotype (\@simmonds);
	 $tonnagegeno{$isoname} = &genotype (\@tonnage);
	$allsnp{$isoname} = $isoname;
		}
	}
}

	foreach my $key (sort keys %bacongeno)	{
			(my $z1,my $z2) = split /\t/, $bacongeno{$key};
			if (($z1 eq "X")||($z1 eq "-"))	{
				 $bacon_nosnp++;
				$bacongeno{nosnp} = $bacon_nosnp;
				$bno_snp{$key} = $bacongeno{$key};
			}elsif ($z1 eq $z2)	{
				 $bacon_homosnp++;
				$bacongeno{homosnp} = $bacon_homosnp;
				$bhomo{$key} = $bacongeno{$key};
			}else	{
				 $bacon_hetsnp++;
				$bacongeno{hetsnp} = $bacon_hetsnp;
				$allhet{$key} = $bacongeno{$key};	
			}
} 

	foreach my $key (sort keys %hassgeno)	{
			(my $z1,my $z2) = split /\t/, $hassgeno{$key};
			if (($z1 eq "X")||($z1 eq "-"))	{
				 $hass_nosnp++;
				$hassgeno{nosnp} = $hass_nosnp;
				$hno_snp{$key} = $hassgeno{$key};
			}elsif ($z1 eq $z2)	{
				 $hass_homosnp++;
				$hassgeno{homosnp} = $hass_homosnp;
				$hhomo{$key} = $hassgeno{$key};
			}else	{
				 $hass_hetsnp++;
				$hassgeno{hetsnp} = $hass_hetsnp;
				$allhet{$key} = $hassgeno{$key};
			}
} 
 
	foreach my $key (sort keys %simmondsgeno)	{
			(my $z1,my $z2) = split /\t/, $simmondsgeno{$key};
			if (($z1 eq "X")||($z1 eq "-"))	{
				 $simmonds_nosnp++;
				$simmondsgeno{nosnp} = $simmonds_nosnp;
				$sno_snp{$key} = $simmondsgeno{$key};
			}elsif ($z1 eq $z2)	{
				 $simmonds_homosnp++;
				$simmondsgeno{homosnp} = $simmonds_homosnp;
				$shomo{$key} = $simmondsgeno{$key};
			}else	{
				$simmonds_hetsnp++;
				$simmondsgeno{hetsnp} = $simmonds_hetsnp;
				$allhet{$key} = $simmondsgeno{$key};
			}
} 

foreach my $key (sort keys %tonnagegeno)	{
			(my $z1,my $z2) = split /\t/, $tonnagegeno{$key};
			if (($z1 eq "X")||($z1 eq "-"))	{
				$tonnage_nosnp++;
				$tonnagegeno{nosnp} = $tonnage_nosnp;
				$tno_snp{$key} = $tonnagegeno{$key};
			}elsif ($z1 eq $z2)	{
				$tonnage_homosnp++;
				$tonnagegeno{homosnp} = $tonnage_homosnp;
				$thomo{$key} = $tonnagegeno{$key};
			}else	{
				$tonnage_hetsnp++;
				$tonnagegeno{hetsnp} = $tonnage_hetsnp;
				$allhet{$key} = $tonnagegeno{$key};
			}
} 
			my $bacon_total_snp = $bacongeno{hetsnp} + $bacongeno{homosnp} + $bacongeno{nosnp};
			my $bacon_percent_het = $bacongeno{hetsnp}/ ($bacongeno{hetsnp} + $bacongeno{homosnp} + $bacongeno{nosnp});
#			print "bacon nosnp $bacongeno{nosnp}\n";
#			print "Bacon percent het is $bacon_percent_het \n";
#			print "Total SNPs are $bacon_total_snp\n";
			my $hass_total_snp = $hassgeno{hetsnp} + $hassgeno{homosnp} + $hassgeno{nosnp};
			my $hass_percent_het = $hassgeno{hetsnp}/($hassgeno{hetsnp} + $hassgeno{homosnp} + $hassgeno{nosnp});
			my $simmonds_total_snp = $simmondsgeno{hetsnp} + $simmondsgeno{homosnp} + $simmondsgeno{nosnp};
			my $simmonds_percent_het = $simmondsgeno{hetsnp}/($simmondsgeno{hetsnp} + $simmondsgeno{homosnp} + $simmondsgeno{nosnp});
			my $tonnage_total_snp = $tonnagegeno{hetsnp} + $tonnagegeno{homosnp} + $tonnagegeno{nosnp};
			my $tonnage_percent_het = $tonnagegeno{hetsnp}/($tonnagegeno{hetsnp} + $tonnagegeno{homosnp} + $tonnagegeno{nosnp});

	open (OUTFILE, ">>Mar28_2013_percent_hetout_missing_data_equal.txt") || die "Can't open outfile";

		
		printf OUTFILE "The percent heterozygosity in Bacon is %1.3f\nThe number of Bacon nosnps is $bacongeno{nosnp}\nThe number of Bacon hetsnps is $bacongeno{hetsnp}\nThe number of Bacon homosnps is $bacongeno{homosnp}\nThe total number of Bacon snps is $bacon_total_snp\n. The percent heterozygosity in Hass is %1.3f\nThe number of Hass nosnps is $hassgeno{nosnp}\nThe number of Hass hetsnps is $hassgeno{hetsnp}\nThe number of Hass homosnps is $hassgeno{homosnp}\nThe total number of Hass snps is $hass_total_snp\n.The percent heterozygosity in Simmonds is %1.3f\nThe number of Simmonds nosnps is $simmondsgeno{nosnp}\nThe number of Simmonds hetsnps is $simmondsgeno{hetsnp}\nThe number of Simmonds homosnps is $simmondsgeno{homosnp}\nThe total number of Simmonds snps is $simmonds_total_snp\n.The percent heterozygosity in Tonnage is %1.3f\nThe number of Tonnage nosnps is $tonnagegeno{nosnp}\nThe number of Tonnage hetsnps is $tonnagegeno{hetsnp}\nThe number of Tonnage homosnps is $tonnagegeno{homosnp}\nThe total number of Tonnage snps is $tonnage_total_snp\n.", $bacon_percent_het, $hass_percent_het, $simmonds_percent_het, $tonnage_percent_het;

#################################################################################################################################################
#  Determine type of SNPs (transition, transversion and actual type for all SNPs) using %allhet

		foreach my $key (sort keys %allhet)	{
			if ($allhet{$key} eq "A\tC")	{
				$AC{$key} = $allhet{$key};
				$transversion{$key} = $allhet{$key};
			}elsif ( $allhet{$key} eq "A\tG")	{
				$AG{$key} = $allhet{$key};
				$transition{$key} = $allhet{$key};
			}elsif ( $allhet{$key} eq "A\tT")	{
				$AT{$key} = $allhet{$key};
				$transversion{$key} = $allhet{$key};	
			}elsif ( $allhet{$key} eq "C\tG")	{
				$CG{$key} = $allhet{$key};
				$transversion{$key} = $allhet{$key};
			}elsif ( $allhet{$key} eq "C\tT")	{
				$CT{$key} = $allhet{$key};
				$transition{$key} = $allhet{$key};
			}elsif ( $allhet{$key} eq "G\tT")	{
				$GT{$key} = $allhet{$key};
				$transversion{$key} = $allhet{$key};	
	}
}

		my $transition_count = keys %transition;
		print "The transition count is $transition_count\n";
		my $transversion_count = keys %transversion;
		print "The transversion count is $transversion_count\n";
		my $allSNP_count = keys %allhet;
		print "The number of SNPs is $allSNP_count\n";
		my $AC_count = keys %AC;
		my $AG_count = keys %AG;
		my $AT_count = keys %AT;
		my $CG_count = keys %CG;
		my $CT_count = keys %CT;
		my $GT_count = keys %GT;
		my $SNP_total = $transition_count + $transversion_count;
		print "The SNP total of transition + transversion SNPs is $SNP_total \n";

		foreach my $key (sort keys %allsnp)		{
			if (((exists $bhomo{$key})||(exists $bno_snp{$key})) && ((exists $hhomo{$key})|| (exists $hno_snp{$key})) && ((exists $shomo{$key})|| (exists $sno_snp{$key})) && ((exists $thomo{$key})||(exists $tno_snp{$key})))	{
				$allhomo{$key} = $key;
			}
		}
		my $goodsnp_count = keys %allsnp;
		my $allhomo_count = keys %allhomo;
		

		print "The all homozygous SNP total is $allhomo_count \n";
		print "The total number of SNPs is $goodsnp_count \n";
#########################################################################################################################################################################
#  Print out summary of transitions and transversions and types

open (OUT, ">>Mar29_2013_Transition_transversion_stats_missing_data_equal.txt") || die "Can't open transition_transversion_stats out file $!\n";
		printf OUT  "%-25s %-25s %-25s %-25s %-25s \n", "Type","Number", "Percent of transitions", "Percent of transversions", "Percent of total";
		printf OUT  "%-25s %-25.2f %-25.2f %-25.2f %-25.2f\n", "Het Transitions",$transition_count, 1, 0, ($transition_count/$SNP_total);
		printf OUT  "%-25s %-25.2f %-25.2f %-25.2f %-25.2f\n", "Het Transversions", $transversion_count, 0, 1, $transversion_count/$SNP_total;
		printf OUT  "%-25s %-25.2f %-25.2f %-25.2f %-25.2f\n",  "Het AC",$AC_count, 0, $AC_count/$transversion_count, $AC_count/$SNP_total;
		printf OUT  "%-25s %-25.2f %-25.2f %-25.2f %-25.2f\n",  "Het AG", $AG_count, $AG_count/$transition_count, 0, $AG_count/$SNP_total;
		printf OUT  "%-25s %-25.2f %-25.2f %-25.2f %-25.2f\n", "Het AT" , $AT_count, 0, $AT_count/$transversion_count, $AT_count/$SNP_total;
		printf OUT  "%-25s %-25.2f %-25.2f %-25.2f %-25.2f\n", "Het CG", $CG_count, 0, $CG_count/$transversion_count, $CG_count/$SNP_total;
		printf OUT  "%-25s %-25.2f %-25.2f %-25.2f %-25.2f\n", "Het CT", $CT_count, $CT_count/$transition_count, 0, $CT_count/$SNP_total;
		printf OUT  "%-25s %-25.2f %-25.2f %-25.2f %-25.2f\n", "Het GT", $GT_count, 0, $GT_count/$transversion_count, $GT_count/$SNP_total;	
		printf OUT  "The number of all homozygous SNPs is $allhomo_count \n";
close OUT;

open (ALLHET, ">>Mar29_2013_allhet_hash_out_missing_data_equal.txt") || die "Can't open all het hash out $!\n";
		printf ALLHET "%-25s %-25s \n", "SNP name", "Het genotype";
		foreach my $key (sort keys %allhet)		{
			printf ALLHET "%25d %25s \n",$key, $allhet{$key};
		}
close ALLHET;
open (ALLHOMO, ">>Mar29_2013_allhomo_hash_out_missing_data_removed_equal.txt") || die "Can't open all homo hash out $!\n";
		printf ALLHOMO "%-25s %-25s %-25s %-25s %-25s %-25s %-25s %-25s %-25s\n", "SNP name", "Bacon geno","Bacon no_snp","Hass geno","Hass no_snp", "Simmonds geno","Simmonds no snp", "Tonnage geno", "Tonnage no snp" ;
		foreach my $key (sort keys %allhomo)	{
			printf ALLHOMO "%-25s %-25s %-25s %-25s %-25s %-25s %-25s %-25s %-25s \n", "$key", "$bhomo{$key}","$bno_snp{$key}", "$hhomo{$key}","$hno_snp{$key}", "$shomo{$key}","$sno_snp{$key}", "$thomo{$key}", "$tno_snp{$key}";			
	}
close ALLHOMO;
open (HOMOSNP, "Mar29_2013_allhomo_hash_out_missing_data_removed_equal.txt") || die "Can't open allhomo snp file $!\n";
	while (<HOMOSNP>)	{
		chomp;
			if ((m/A/) && (m/C/))	{
				$homosnpAC++;
				$homotransversion++;
			}elsif ((m/A/) && (m/G/))	{
				$homosnpAG++;
				$homotransition++;
			}elsif  ((m/A/) && (m/T/))	{
				$homosnpAT++;
				$homotransversion++;
			}elsif ((m/C/) && (m/G/))	{
				$homosnpCG++;
				$homotransversion++;
			}elsif ((m/C/) && (m/T/))	{
				$homosnpCT++;
				$homotransition++;
			}elsif ((m/G/) && (m/T/))	{
				$homosnpGT++;
				$homotransversion++;
			}
	}
close HOMOSNP;
			my $totalhomo = $homotransition + $homotransversion;
open (SNPREPORT, ">>Mar29_2013_snpreport.txt") || die "Can't open snp report out $!\n";
		printf SNPREPORT "%-25s %-25s %-25s %-25s %-25s \n", "Type","Number", "Percent of transitions", "Percent of transversions", "Percent of total";
		printf SNPREPORT "%-25s %-25.2f %-25.2f %-25.2f %-25.2f\n", "Homo transition", $homotransition, $homotransition/$homotransition, 0, $homotransition/$totalhomo;
		printf SNPREPORT "%-25s %-25.2f %-25.2f %-25.2f %-25.2f \n", "Homo transversion", $homotransversion, 0, $homotransversion/$homotransversion, $homotransversion/$totalhomo;
		printf SNPREPORT "%-25s %-25.2f %-25.2f %-25.2f %-25.2f \n",	"Homo AC", $homosnpAC, 0, $homosnpAC/$homotransversion, $homosnpAC/$totalhomo;
		printf SNPREPORT "%-25s %-25.2f %-25.2f %-25.2f %-25.2f \n", "Homo AG", $homosnpAG, $homosnpAG/$homotransition, 0, $homosnpAG/$totalhomo;
		printf SNPREPORT "%-25s %-25.2f %-25.2f %-25.2f %-25.2f \n", "Homo AT", $homosnpAT, 0, $homosnpAT/$homotransversion, $homosnpAT/$totalhomo;
		printf SNPREPORT "%-25s %-25.2f %-25.2f %-25.2f %-25.2f \n", "Homo CG", $homosnpCG, 0, $homosnpCG/$homotransversion, $homosnpCG/$totalhomo;
		printf SNPREPORT "%-25s %-25.2f %-25.2f %-25.2f %-25.2f \n", "Homo CT", $homosnpCT, $homosnpCT/$homotransition, 0, $homosnpCT/$totalhomo;
		printf SNPREPORT "%-25s %-25.2f %-25.2f %-25.2f %-25.2f \n", "Homo GT", $homosnpGT, 0, $homosnpGT/$homotransversion, $homosnpGT/$totalhomo;
		printf SNPREPORT "%-25s %-25d \n", "Total homo SNPs", $totalhomo;
close SNPREPORT;

		
open (NOSNP, ">>Mar29_2013_nosnp_hash_out_equal.txt") || die "Can't open nosnp hash out $! \n";
		printf NOSNP "%-25s %-25s %-25s %-25s %-25s \n", "SNP name", "Bacon geno", "Hass geno", "Simmonds geno", "Tonnage geno" ;
		foreach my $key (sort keys %allsnp)	{
			if ((exists $bno_snp{$key}) && (exists $hno_snp{$key}) && (exists $sno_snp{$key}) && (exists $tno_snp{$key}))	{
			$all_no_snp{$key} = $key;
			printf NOSNP "%-25s %-25s %-25s %-25s %-25s \n", "$key", "$bno_snp{$key}", "$hno_snp{$key}", "$sno_snp{$key}", "$tno_snp{$key}";			
	}
}
			my $all_no_snp_count = keys (%all_no_snp);
			print 	"The total number of loci with no snp for B, H, S, T is $all_no_snp_count \n";
			printf NOSNP "The total number of loci with no snp for B, H, S, T is $all_no_snp_count \n";
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

#	open (OUTFILE, ">>avo_SNPtype_counter_BHSThet.txt") || die "Can't open outfile";
#
#		
#		printf OUTFILE "%-25s %-25s %-25s %-25s %-25s %-25s %-25s %-25s %-25s\n", "Cultivar name", "No of N SNPs", "No of indels", "No of homo SNPs", "percentage of homo snps", "No of het snps", "percentage of hetsnps", "total indels","total snps";
#		printf OUTFILE "%-25s, %-25d, %-25d, %-25d, %-25.3f, %-25d, %-25.3f, %-25d %-25d\n", "Bacon", keys %allN, keys %allindel, keys %bhomo, (keys %bhomo)/ ((keys %bhet) + (keys %bhomo)), keys %bhet, $bacon_percent_het, keys %allindel,$bacon_total_snp;
#		printf OUTFILE "%-25s, %-25d, %-25d, %-25d, %-25.3f, %-25d, %-25.3f, %-25d %-25d\n", "Hass", keys %allN, keys %allindel, keys %hhomo, (keys %hhomo)/ ((keys %hhet) + (keys %hhomo)), keys %hhet, $hass_percent_het, keys %allindel,keys %N_snp + keys %hindel + $hass_total_snp;
#		printf OUTFILE "%-25s, %-25d, %-25d, %-25d, %-25.3f, %-25d, %-25.3f, %-25d %-25d\n", "Simmonds", keys %N_snp, keys %sindel, keys %shomo, (keys %shomo)/ ((keys %shet) + (keys %shomo)), keys %shet, $simmonds_percent_het, keys %allindel, keys %N_snp + keys %sindel + $simmonds_total_snp;
#		printf OUTFILE "%-25s, %-25d, %-25d, %-25d, %-25.3f, %-25d, %-25.3f, %-25d %-25d\n", "Tonnage", keys %N_snp, keys %tindel, keys %thomo, (keys %thomo)/ ((keys %thet) + (keys %thomo)), keys %thet, $tonnage_percent_het, keys %allindel, keys %N_snp + keys %tindel + $tonnage_total_snp;

#  Print out summary of transitions and transversions and types

#open (OUT, ">>avo_SNPtype_counter_het_homo_SNP_types.txt") || die "Can't open transition_transversion_stats out file $!\n";
#		printf OUT  "%-25s %-25s %-25s %-25s %-25s \n", "Type","Number", "Percent of transitions", "Percent of transversions", "Percent of total";
#		printf OUT  "%-25s %-25.2f %-25.2f %-25.2f %-25.2f\n", "Het Transversions", $het_transversion_count, 0, 1, $het_transversion_count/$het_SNP_total;
#		printf OUT  "%-25s %-25.2f %-25.2f %-25.2f %-25.2f\n",  "Het AC",$het_AC_count, 0, $het_AC_count/$het_transversion_count, $AC_count/$het_SNP_total;
#		printf OUT  "%-25s %-25.2f %-25.2f %-25.2f %-25.2f\n",  "Het AG", $het_AG_count, $het_AG_count/$het_transition_count, 0, $AG_count/$het_SNP_total;
#		printf OUT  "%-25s %-25.2f %-25.2f %-25.2f %-25.2f\n", "Het AT" , $het_AT_count, 0, $het_AT_count/$het_transversion_count, $AT_count/$het_SNP_total;
#		printf OUT  "%-25s %-25.2f %-25.2f %-25.2f %-25.2f\n", "Het CG", $het_CG_count, 0, $het_CG_count/$het_transversion_count, $CG_count/$het_SNP_total;
#		printf OUT  "%-25s %-25.2f %-25.2f %-25.2f %-25.2f\n", "Het CT", $het_CT_count, $het_CT_count/$het_transition_count, 0, $CT_count/$het_SNP_total;
#		printf OUT  "%-25s %-25.2f %-25.2f %-25.2f %-25.2f\n", "Het GT", $het_GT_count, 0, $het_GT_count/$het_transversion_count, $GT_count/$het_SNP_total;	
#		printf OUT "%-25s %-25.2f %-25.2f %-25.2f %-25.2f\n", "Homo transition", $homotransition, $homotransition/$homotransition, 0, $homotransition/$totalhomo;
#		printf OUT "%-25s %-25.2f %-25.2f %-25.2f %-25.2f \n", "Homo transversion", $homotransversion, 0, $homotransversion/$homotransversion, $homotransversion/$totalhomo;
#		printf OUT "%-25s %-25.2f %-25.2f %-25.2f %-25.2f \n",	"Homo AC", $homosnpAC, 0, $homosnpAC/$homotransversion, $homosnpAC/$totalhomo;
#		printf OUT "%-25s %-25.2f %-25.2f %-25.2f %-25.2f \n", "Homo AG", $homosnpAG, $homosnpAG/$homotransition, 0, $homosnpAG/$totalhomo;
#		printf OUT "%-25s %-25.2f %-25.2f %-25.2f %-25.2f \n", "Homo AT", $homosnpAT, 0, $homosnpAT/$homotransversion, $homosnpAT/$totalhomo;
#		printf OUT "%-25s %-25.2f %-25.2f %-25.2f %-25.2f \n", "Homo CG", $homosnpCG, 0, $homosnpCG/$homotransversion, $homosnpCG/$totalhomo;
#		printf OUT "%-25s %-25.2f %-25.2f %-25.2f %-25.2f \n", "Homo CT", $homosnpCT, $homosnpCT/$homotransition, 0, $homosnpCT/$totalhomo;
#		printf OUT "%-25s %-25.2f %-25.2f %-25.2f %-25.2f \n", "Homo GT", $homosnpGT, 0, $homosnpGT/$homotransversion, $homosnpGT/$totalhomo;
#		printf OUT "%-25s %-25.2f %-25.2f %-25.2f %-25.2f\n", "Total transition", $homotransition + $het_transition_count, 1, 0, ($homotransition + $het_transition_count)/$totalhomo + $het_SNP_total;
#		printf OUT "%-25s %-25.2f %-25.2f %-25.2f %-25.2f \n", "Total transversion", $homotransversion + $het_transversion_count, 0, 1, $homotransversion + $het_transversion_count/$totalhomo + $het_SNP_total;
#		printf OUT "%-25s %-25.2f %-25.2f %-25.2f %-25.2f \n",	"Total AC", $homosnpAC + $het_AC_count, 0, $homosnpAC + $het_AC_count/$homotransversion + $het_transversion_count, $homosnpAC + $het_AC_count/$totalhomo + $het_SNP_total;
#		printf OUT "%-25s %-25.2f %-25.2f %-25.2f %-25.2f \n", "Total AG", $homosnpAG, $homosnpAG/$homotransition, 0, $homosnpAG/$totalhomo + $het_SNP_total;
#		printf OUT "%-25s %-25.2f %-25.2f %-25.2f %-25.2f \n", "Total AT", $homosnpAT, 0, $homosnpAT/$homotransversion + $het_transversion_count, $homosnpAT/$totalhomo + $het_SNP_total;
#		printf OUT "%-25s %-25.2f %-25.2f %-25.2f %-25.2f \n", "Total CG", $homosnpCG, 0, $homosnpCG/$homotransversion + $het_transversion_count, $homosnpCG/$totalhomo + $het_SNP_total;
#		printf OUT "%-25s %-25.2f %-25.2f %-25.2f %-25.2f \n", "Total CT", $homosnpCT, $homosnpCT/$homotransition, 0, $homosnpCT/$totalhomo + $het_SNP_total;
#		printf OUT "%-25s %-25.2f %-25.2f %-25.2f %-25.2f \n", "Total GT", $homosnpGT, 0, $homosnpGT/$homotransversion + $het_transversion_count, $homosnpGT/$totalhomo + $het_SNP_total;
#close OUT;

		
#open (NOSNP, ">>Mar29_2013_nosnp_hash_out_equal.txt") || die "Can't open nosnp hash out $! \n";
#		printf NOSNP "%-25s %-25s %-25s %-25s %-25s \n", "SNP name", "Bacon geno", "Hass geno", "Simmonds geno", "Tonnage geno" ;
#		foreach my $key (sort keys %allsnp)	{
#			if ((exists $bno_snp{$key}) && (exists $hno_snp{$key}) && (exists $sno_snp{$key}) && (exists $tno_snp{$key}))	{
#			$all_no_snp{$key} = $key;
#			printf NOSNP "%-25s %-25s %-25s %-25s %-25s \n", "$key", "$bno_snp{$key}", "$hno_snp{$key}", "$sno_snp{$key}", "$tno_snp{$key}";			
#	}
#}
#			my $all_no_snp_count = keys (%all_no_snp);
#			print 	"The total number of loci with no snp for B, H, S, T is $all_no_snp_count \n";
#			printf NOSNP "The total number of loci with no snp for B, H, S, T is $all_no_snp_count \n";




