#!/usr/bin/perl -w
##################################################################################
# Author		: David Kuhn
# Date			: 5/9/2012
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: 	reads in avocado variant report, counts number of heterozygous SNPs for each parent of avocado mapping populations and outputs counts
# <input>		:  unifiedSNPs_Avocado_BLASTX_TAIR10_ANNOT.table_cvg.txt			
# <output>		:   QQ_percent_hetout.txt 
# 			
# Dependencies	: N/A
# Usage			: 
#  Still under development (5/24/12).                    
##################################################################################
use strict;
use warnings;
use diagnostics;
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
#		print "Bacon array is @bacon\n";
		 @hass = ($snp[1], $snp[2], $snp[3], $snp[12], $snp[13], $snp[14], $snp[15], $snp[16], $snp[17], $snp[18], $snp[19]);
		 @simmonds = ($snp[1], $snp[2], $snp[3], $snp[20], $snp[21], $snp[22], $snp[23], $snp[24], $snp[25], $snp[26], $snp[27]);
		 @tonnage = ($snp[1], $snp[2], $snp[3], $snp[28], $snp[29], $snp[30], $snp[31], $snp[32], $snp[33], $snp[34], $snp[35]);
	 $bacongeno{$isoname} = &genotype (\@bacon);
	$hassgeno{$isoname} = &genotype (\@hass);
	 $simmondsgeno{$isoname} = &genotype (\@simmonds);
	 $tonnagegeno{$isoname} = &genotype (\@tonnage);
		}
	}
}
#	foreach my $key (sort keys %bacongeno)	{
#			print "bacon geno $key \t $bacongeno{$key}\n";
#		}

	foreach my $key (sort keys %bacongeno)	{
			(my $z1,my $z2) = split /\t/, $bacongeno{$key};
#			print "for $key, z1 is $z1, z2 is $z2\n";
			if ($z1 eq "X")	{
				 $bacon_nosnp++;
				$bacongeno{nosnp} = $bacon_nosnp;
			}elsif ($z1 eq $z2)	{
				 $bacon_homosnp++;
				$bacongeno{homosnp} = $bacon_homosnp;
			}else	{
				 $bacon_hetsnp++;
				$bacongeno{hetsnp} = $bacon_hetsnp;
			}
} 

	foreach my $key (sort keys %hassgeno)	{
			(my $z1,my $z2) = split /\t/, $hassgeno{$key};
			if ($z1 eq "X")	{
				 $hass_nosnp++;
				$hassgeno{nosnp} = $hass_nosnp;
			}elsif ($z1 eq $z2)	{
				 $hass_homosnp++;
				$hassgeno{homosnp} = $hass_homosnp;
			}else	{
				 $hass_hetsnp++;
				$hassgeno{hetsnp} = $hass_hetsnp;
			}
} 
 
	foreach my $key (sort keys %simmondsgeno)	{
			(my $z1,my $z2) = split /\t/, $simmondsgeno{$key};
			if ($z1 eq "X")	{
				 $simmonds_nosnp++;
				$simmondsgeno{nosnp} = $simmonds_nosnp;
			}elsif ($z1 eq $z2)	{
				 $simmonds_homosnp++;
				$simmondsgeno{homosnp} = $simmonds_homosnp;
			}else	{
				$simmonds_hetsnp++;
				$simmondsgeno{hetsnp} = $simmonds_hetsnp;
			}
} 

foreach my $key (sort keys %tonnagegeno)	{
			(my $z1,my $z2) = split /\t/, $tonnagegeno{$key};
			if ($z1 eq "X")	{
				$tonnage_nosnp++;
				$tonnagegeno{nosnp} = $tonnage_nosnp;
			}elsif ($z1 eq $z2)	{
				$tonnage_homosnp++;
				$tonnagegeno{homosnp} = $tonnage_homosnp;
			}else	{
				$tonnage_hetsnp++;
				$tonnagegeno{hetsnp} = $tonnage_hetsnp;
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

open (HETOUT, ">>avo_SNPtype_het_BHST_counter_May3_out.txt") || die "Can't open HETOUT $!\n";
		printf HETOUT ("%20s\t" x 8) . "\n", "Cultivar","SNP total", "Het total","Het percent", "Homo total", "Homo percent", "No SNP total", "No SNP percent";
		printf HETOUT "%20s\t" . "%20d\t" . ("%20d\t%20.3f" x 3). "\n", "Bacon", $bacon_total_snp,$bacongeno{hetsnp},$bacon_percent_het,$bacongeno{homosnp}, $bacongeno{homosnp}/$bacon_total_snp, $bacongeno{nosnp}, $bacongeno{nosnp}/$bacon_total_snp;   
		printf HETOUT "%20s\t" . "%20d\t" . ("%20d\t%20.3f" x 3). "\n", "Hass", $hass_total_snp,$hassgeno{hetsnp},$hass_percent_het,$hassgeno{homosnp}, $hassgeno{homosnp}/$hass_total_snp, $hassgeno{nosnp}, $hassgeno{nosnp}/$hass_total_snp;   

		printf HETOUT "%20s\t" . "%20d\t" . ("%20d\t%20.3f" x 3). "\n", "Simmonds", $simmonds_total_snp,$simmondsgeno{hetsnp},$simmonds_percent_het,$simmondsgeno{homosnp}, $simmondsgeno{homosnp}/$simmonds_total_snp, $simmondsgeno{nosnp}, $simmondsgeno{nosnp}/$simmonds_total_snp;   

		printf HETOUT "%20s\t" . "%20d\t" . ("%20d\t%20.3f" x 3). "\n", "Tonnage", $tonnage_total_snp,$tonnagegeno{hetsnp},$tonnage_percent_het,$tonnagegeno{homosnp}, $tonnagegeno{homosnp}/$tonnage_total_snp, $tonnagegeno{nosnp}, $tonnagegeno{nosnp}/$tonnage_total_snp;   
close HETOUT;

	open (OUTFILE, ">>avo_het_counter_May3_out.txt") || die "Can't open outfile";

		
		printf OUTFILE "The percent heterozygosity in Bacon is %1.3f\nThe number of Bacon nosnps is $bacongeno{nosnp}\nThe number of Bacon hetsnps is $bacongeno{hetsnp}\nThe number of Bacon homosnps is $bacongeno{homosnp}\nThe total number of Bacon snps is $bacon_total_snp\n. The percent heterozygosity in Hass is %1.3f\nThe number of Hass nosnps is $hassgeno{nosnp}\nThe number of Hass hetsnps is $hassgeno{hetsnp}\nThe number of Hass homosnps is $hassgeno{homosnp}\nThe total number of Hass snps is $hass_total_snp\n.The percent heterozygosity in Simmonds is %1.3f\nThe number of Simmonds nosnps is $simmondsgeno{nosnp}\nThe number of Simmonds hetsnps is $simmondsgeno{hetsnp}\nThe number of Simmonds homosnps is $simmondsgeno{homosnp}\nThe total number of Simmonds snps is $simmonds_total_snp\n.The percent heterozygosity in Tonnage is %1.3f\nThe number of Tonnage nosnps is $tonnagegeno{nosnp}\nThe number of Tonnage hetsnps is $tonnagegeno{hetsnp}\nThe number of Tonnage homosnps is $tonnagegeno{homosnp}\nThe total number of Tonnage snps is $tonnage_total_snp\n.", $bacon_percent_het, $hass_percent_het, $simmonds_percent_het, $tonnage_percent_het;
	


sub genotype	{	#determines genotype of SNP (XX if not enough reads) (homozygous major nucleotide) or heterozygous
	my @x = @{$_[0]};
#	print "This is the subroutine array @x\n";
	if ($x[4] < 6)	{
		my $genotype = "X\tX";
	} elsif ($x[5]eq "\.")	{
		my $genotype = "$x[1]\t$x[1]";
	} elsif ($x[7]> 0.80)	{
		my $genotype = "$x[5]\t$x[5]";
	} elsif	(($x[7] < 0.80) && ($x[10] > 0.20))	{
		my $genotype = "$x[5]\t$x[8]";
	}
}


