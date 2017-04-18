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
#use strict;
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
my $baconuniq_nosnp = 0;
my $baconuniq_hetsnp = 0;
my $baconuniq_homosnp = 0;
my $hassuniq_nosnp = 0;
my $hassuniq_homosnp = 0;
my $hassuniq_hetsnp = 0;
my $simmondsuniq_nosnp = 0;
my $simmondsuniq_homosnp = 0;
my $simmondsuniq_hetsnp = 0;
my $tonnageuniq_nosnp = 0;
my $tonnageuniq_homosnp = 0;
my $tonnageuniq_hetsnp = 0;
my %avosnp = ();
my $baconuniq = 0;
my $hassuniq = 0;
my $simmondsuniq = 0;
my $tonnageuniq = 0;
my $notuniq = 0;
my $totaluniq = 0;
my %baconuniq = ();
my %hassuniq = ();
my %simmondsuniq = ();
my %tonnageuniq = ();


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
		$avosnp{$isoname} = $isoname;
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

	foreach my $key (sort keys %avosnp)	{
		if (($bacongeno{$key} ne $hassgeno{$key}) && ($bacongeno{$key} ne $simmondsgeno{$key}) && ($bacongeno{$key} ne $tonnagegeno{$key}))	{
		$avosnp{$key}{unique} = "bacon";
		$baconuniq++;
		$avosnp{$key}{genotype} = $bacongeno{$key};
		$baconuniq{$key} = $bacongeno{$key};
			}elsif	(( $hassgeno{$key} ne $bacongeno{$key}) && ($hassgeno{$key} ne $simmondsgeno{$key}) && ($hassgeno{$key}	ne $tonnagegeno{$key}))	{
		$avosnp{$key}{unique} = "hass";
		$avosnp{$key}{genotype} = $hassgeno{$key};
		$hassuniq++;
		$hassuniq{$key} = $hassgeno{$key};
			}elsif	(( $simmondsgeno{$key} ne $bacongeno{$key}) && ($simmondsgeno{$key} ne $hassgeno{$key}) && ($simmondsgeno{$key}	ne $tonnagegeno{$key}))	{
		$avosnp{$key}{unique} = "simmonds";
		$simmondsuniq++;
		$avosnp{$key}{genotype} = $simmondsgeno{$key};
		$simmondsuniq{$key} = $simmondsgeno{$key};
			}elsif	(( $tonnagegeno{$key} ne $bacongeno{$key}) && ($tonnagegeno{$key} ne $simmondsgeno{$key}) && ($tonnagegeno{$key} ne $hassgeno{$key}))	{
		$avosnp{$key}{unique} = "tonnage";
		$tonnageuniq++;
		$avosnp{$key}{genotype} = $tonnagegeno{$key};
		$tonnageuniq{$key} = $tonnagegeno{$key};
			}else{
		$avosnp{$key}{unique} = "not unique";
		$notuniq++;
	}
}
		$totaluniq = $baconuniq + $hassuniq + $simmondsuniq + $tonnageuniq + $notuniq;
	open (UNIQSNP, ">>ZZuniqsnp.out") || die "Can't open uniqsnp outfile $!";
		print UNIQSNP "The total number of Bacon unique SNPs is $baconuniq\nThe total number of Hass unique SNPs is $hassuniq\nThe total number of Simmonds unique SNPs is $simmondsuniq\nThe total number of Tonnage unique SNPs is $tonnageuniq\nThe total number of not unique SNPs is $notuniq\nThe total number of SNPs analyzed is $totaluniq\n";

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

	foreach my $key (sort keys %baconuniq)	{
			(my $z1,my $z2) = split /\t/, $baconuniq{$key};
#			print "for $key, z1 is $z1, z2 is $z2\n";
			if ($z1 eq "X")	{
				 $baconuniq_nosnp++;
			}elsif ($z1 eq $z2)	{
				 $baconuniq_homosnp++;
			}else	{
				 $baconuniq_hetsnp++;
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
	foreach my $key (sort keys %hassuniq)	{
			(my $z1,my $z2) = split /\t/, $hassuniq{$key};
			if ($z1 eq "X")	{
				 $hassuniq_nosnp++;
			}elsif ($z1 eq $z2)	{
				 $hassuniq_homosnp++;
			}else	{
				 $hassuniq_hetsnp++;
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
	foreach my $key (sort keys %simmondsuniq)	{
			(my $z1,my $z2) = split /\t/, $simmondsuniq{$key};
			if ($z1 eq "X")	{
				 $simmondsuniq_nosnp++;
			}elsif ($z1 eq $z2)	{
				 $simmondsuniq_homosnp++;
			}else	{
				$simmondsuniq_hetsnp++;
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

	foreach my $key (sort keys %tonnageuniq)	{
			(my $z1,my $z2) = split /\t/, $tonnageuniq{$key};
			if ($z1 eq "X")	{
				$tonnageuniq_nosnp++;
			}elsif ($z1 eq $z2)	{
				$tonnageuniq_homosnp++;
			}else	{
				$tonnageuniq_hetsnp++;
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

	open (UNIQUE, ">>YY_type_of_unique_SNPs_out") || die "Can't open type of unique SNPs outfile $!\n";

	print UNIQUE  "Bacon unique nosnps is $baconuniq_nosnp\nBacon unique homosnps is $baconuniq_homosnp\nBacon unique heterosnps is $baconuniq_hetsnp\nHass unique nosnps is $hassuniq_nosnp\nHass unique homosnp is $hassuniq_homosnp\nHass unique hetsnp is $hassuniq_hetsnp\nSimmonds unique nosnp is $simmondsuniq_nosnp\nSimmonds unique homosnp is $simmondsuniq_homosnp\nSimmonds unique hetsnp is $simmondsuniq_hetsnp\nTonnage unique nosnp is $tonnageuniq_nosnp\nTonnage unique homosnp is $tonnageuniq_homosnp\nTonnage unique hetsnp is $tonnageuniq_hetsnp\n";

#	open (OUTFILE, ">>QQ_percent_hetout.txt") || die "Can't open outfile";

		
#		printf OUTFILE "The percent heterozygosity in Bacon is %1.3f\nThe number of Bacon nosnps is $bacongeno{nosnp}\nThe number of Bacon hetsnps is $bacongeno{hetsnp}\nThe number of Bacon homosnps is $bacongeno{homosnp}\nThe total number of Bacon snps is $bacon_total_snp\n. The percent heterozygosity in Hass is %1.3f\nThe number of Hass nosnps is $hassgeno{nosnp}\nThe number of Hass hetsnps is $hassgeno{hetsnp}\nThe number of Hass homosnps is $hassgeno{homosnp}\nThe total number of Hass snps is $hass_total_snp\n.The percent heterozygosity in Simmonds is %1.3f\nThe number of Simmonds nosnps is $simmondsgeno{nosnp}\nThe number of Simmonds hetsnps is $simmondsgeno{hetsnp}\nThe number of Simmonds homosnps is $simmondsgeno{homosnp}\nThe total number of Simmonds snps is $simmonds_total_snp\n.The percent heterozygosity in Tonnage is %1.3f\nThe number of Tonnage nosnps is $tonnagegeno{nosnp}\nThe number of Tonnage hetsnps is $tonnagegeno{hetsnp}\nThe number of Tonnage homosnps is $tonnagegeno{homosnp}\nThe total number of Tonnage snps is $tonnage_total_snp\n.", $bacon_percent_het, $hass_percent_het, $simmonds_percent_het, $tonnage_percent_het;
	


sub genotype	{	#determines genotype of SNP (XX if not enough reads) (homozygous major nucleotide) or heterozygous
	my @x = @{$_[0]};
#	print "This is the subroutine array @x\n";
	if ($x[4] < 6)	{
		my $genotype = "X\tX";
	} elsif ($x[5]eq "\.")	{
		my $genotype = "$x[1]\t$x[1]";	#genotype is most frequent nucleotide
	} elsif ($x[7]> 0.80)	{
		my $genotype = "$x[5]\t$x[5]";	#genotype is nucleotide referenced in $x[5]
	} elsif	(($x[7] < 0.80) && ($x[10] > 0.20))	{
		my $genotype = "$x[5]\t$x[8]";	#genotype is heterozygous
	}
}


