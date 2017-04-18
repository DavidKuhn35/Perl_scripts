#!/usr/bin/perl -w
##################################################################################
# Author		: David Kuhn
# Date			: January 25, 2013
# Version		: ver. 1.2
# History		: ver. 1.0
# Description	: 	reads in avocado variant report, counts number of heterozygous SNPs for each parent of avocado mapping populations and outputs counts, calculates number of SNPs
#			 unique to each cultivar, identifies SNPs that are homozygous in Bacon for one allele and homozygous in Simmonds for the other allele.
# <input>		:  KK_avo_snp_lzsnpname			
# <output>		:   fixit_avo_genotypes (B,H,S,T genotypes for all 642961 snps; fixit_uniq_snp (snps uniq to each of B,H,S,T); fixit_type_of_uniq_snp (gives stats on
#	type of snps for B,H,S,T nosnp, homosnp, hetsnp);fixit_percent_het_out (%heterozygosity for B,H,S,T); fixit_MneWI_out (snps that are homo and opposite in Bacon (Mexican) and #Simmonds (West Indian)); fixit_MneWI_isotigs_snps; fixit_MneWI_summary. Prints out unique SNPs for B,H,S,T under new definition.  E.g.  B,H,S homozygous and identical, T ne B => SNP #unique to T.
# 			
# Dependencies	: N/A
# Usage			: 
#  Still under development (6/29/2012) January 25,2013.                    
##################################################################################
use strict;
use warnings;
use diagnostics;
use Statistics::Descriptive;
my %hoh = ();
my %hoh2 = ();
my %isotig = ();
my $goodsnp = 0;
my $badsnp = 0;
my @snpdistribution = ();
my $snpname = "snpname";
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
my %baconhomo = ();
my %baconhet = ();
my %hasshomo = ();
my %hasshet = ();
my %simmondshomo = ();
my %simmondshet = ();
my %tonnagehomo = ();
my %tonnagehet = ();
my %MneWI = ();
my $MneWI = 0;
my %notuniq = ();
my @baconuniq = ();
my @hassuniq = ();
my @simmondsuniq = ();
my @tonnageuniq = ();
my @notuniq = ();
my %baconreads = ();
my %hassreads = ();
my %simmondsreads = ();
my %tonnagereads = ();
my %bacon_unique = ();
my %hass_unique = ();
my %simmonds_unique = ();
my %tonnage_unique = ();
my %bacon_unique_homo = ();
my %hass_unique_homo = ();
my %simmonds_unique_homo = ();
my %tonnage_unique_homo = ();
my %bacon_unique_het = ();
my %hass_unique_het = ();
my %simmonds_unique_het = ();
my %tonnage_unique_het = ();
my $b_homo_unik = 0;
my $b_het_unik = 0;
my $h_homo_unik = 0;
my $h_het_unik = 0;
my $s_homo_unik = 0;
my $s_het_unik = 0;
my $t_homo_unik = 0;
my $t_het_unik = 0;


#Read in the avo snp file with leading zeros for the snp names (KK_avo_snp_lzsnpname) and gets rid of rows with missing data or Ns
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
		 $snpname = $snp[1];
		$avosnp{$snpname} = $snpname;	#hash with snpname keys for all 642k snps for sorting on keys
		my $isotigname = $snp[36];
		 @bacon = ($snp[1], $snp[2], $snp[3], $snp[4], $snp[5], $snp[6], $snp[7], $snp[8], $snp[9], $snp[10], $snp[11]);
# 1 is SNP#, 2 most common nuc, 3 hass ref nuc, 4 bacon, 5 # of total reads for bacon, 6 bacon most common nuc, 7 # of reads supporting most common nuc
# 8 major allele freq, 9 minor allele nuc, 10 # of reads supporting minor allele, 11 minor allele frequency

		 @hass = ($snp[1], $snp[2], $snp[3], $snp[12], $snp[13], $snp[14], $snp[15], $snp[16], $snp[17], $snp[18], $snp[19]);
		 @simmonds = ($snp[1], $snp[2], $snp[3], $snp[20], $snp[21], $snp[22], $snp[23], $snp[24], $snp[25], $snp[26], $snp[27]);
		 @tonnage = ($snp[1], $snp[2], $snp[3], $snp[28], $snp[29], $snp[30], $snp[31], $snp[32], $snp[33], $snp[34], $snp[35]);
	$baconreads{$snpname} = $snp[5];
	 $bacongeno{$snpname} = &genotype (\@bacon);
	$hassgeno{$snpname} = &genotype (\@hass);
	$hassreads{$snpname} = $snp[13];
	 $simmondsgeno{$snpname} = &genotype (\@simmonds);
	$simmondsreads{$snpname} = $snp[21];
	 $tonnagegeno{$snpname} = &genotype (\@tonnage);
	$tonnagereads{$snpname} = $snp[29];
	$hoh{$isotigname}{$snpname} = "$bacongeno{$snpname}\t$baconreads{$snpname}\t$hassgeno{$snpname}\t$hassreads{$snpname}\t$simmondsgeno{$snpname}\t$simmondsreads{$snpname}\t$tonnagegeno{$snpname}\t$tonnagereads{$snpname}";	#hash with isotig name key and snpname key for #
#       sorting on both sets of keys, values are all genotypes for bacon, hass, simmonds, tonnage and total reads for B,H,S,T
#	should make a database of the genotypes for each isotig => snpname
			}
	}
}
#########################################################################################################################################
open	(AVOGENO, ">>fixit_avo_genotypes_out") || die "Can't open avo genotypes out file $!\n";
	print AVOGENO "Bacon genotype\tBacon reads\tHass genotype\tHass reads\tSimmonds genotype\tSimmonds reads\tTonnage genotype\tTonnage reads\n";
	foreach my $key (sort keys %hoh)	{
		foreach my $snp	(sort keys %{$hoh{$key}})	{
			print AVOGENO "$key\t$snp\t$hoh{$key}{$snp}\n";
		}
	}

####################################################################################################################################
#  Count as unique SNPs only those that are identical in the other three parents
	foreach my $key (sort keys %avosnp)	{
		if (($bacongeno{$key} eq $hassgeno{$key}) && ($hassgeno{$key} eq $simmondsgeno{$key}) && ($simmondsgeno{$key} ne $tonnagegeno{$key}))	{
		$tonnageuniq++;
		$tonnageuniq{$key} = "$bacongeno{$key}\t$hassgeno{$key}\t$simmondsgeno{$key}\t$tonnagegeno{$key}";
			}elsif	(( $bacongeno{$key} eq $hassgeno{$key}) && ($hassgeno{$key} eq $tonnagegeno{$key}) && ($tonnagegeno{$key}ne $simmondsgeno{$key}))	{
		$simmondsuniq++;
		$simmondsuniq{$key} = "$bacongeno{$key}\t$hassgeno{$key}\t$simmondsgeno{$key}\t$tonnagegeno{$key}";
			}elsif	(( $bacongeno{$key} eq $simmondsgeno{$key}) && ($simmondsgeno{$key} eq $tonnagegeno{$key}) && ($tonnagegeno{$key}ne $hassgeno{$key}))	{
		$hassuniq++;
		$hassuniq{$key} = "$bacongeno{$key}\t$hassgeno{$key}\t$simmondsgeno{$key}\t$tonnagegeno{$key}";
			}elsif	(( $hassgeno{$key} eq $simmondsgeno{$key}) && ($simmondsgeno{$key} eq $tonnagegeno{$key}) && ($tonnagegeno{$key} ne $bacongeno{$key}))	{
		$baconuniq++;
		$baconuniq{$key} = "$bacongeno{$key}\t$hassgeno{$key}\t$simmondsgeno{$key}\t$tonnagegeno{$key}";
			}else{
		$notuniq++;
		$notuniq{$key} = "$bacongeno{$key}\t$hassgeno{$key}\t$simmondsgeno{$key}\t$tonnagegeno{$key}";
	}
}
		$totaluniq = $baconuniq + $hassuniq + $simmondsuniq + $tonnageuniq + $notuniq;
	open (UNIQSNP, ">>fixit_uniqsnp.out") || die "Can't open uniqsnp outfile $!";
		print UNIQSNP "The total number of Bacon unique SNPs is $baconuniq\nThe total number of Hass unique SNPs is $hassuniq\nThe total number of Simmonds unique SNPs is $simmondsuniq\nThe total number of Tonnage unique SNPs is $tonnageuniq\nThe total number of not unique SNPs is $notuniq\nThe total number of SNPs analyzed is $totaluniq\n";

	foreach my $key (sort keys %baconuniq)	{
			print UNIQSNP "$key\t$baconuniq{$key}\tbacon\n";
		}
	foreach my $key (sort keys %hassuniq)	{
			print UNIQSNP "$key\t$hassuniq{$key}\thass\n";
		}
	foreach my $key (sort keys %simmondsuniq)	{
			print UNIQSNP "$key\t$simmondsuniq{$key}\tsimmonds\n";
		}
	foreach my $key (sort keys %tonnageuniq)	{
			print UNIQSNP "$key\t$tonnageuniq{$key}\ttonnage\n";
		}
	foreach my $key (sort keys %notuniq)	{
			print UNIQSNP "$key\t$notuniq{$key}\tnotuniq\n";
	}
#########################################################################################################################	
	foreach my $key (sort keys %bacongeno)	{
			(my $z1,my $z2) = split /\t/, $bacongeno{$key};
			if ($z1 eq "X")	{
				 $bacon_nosnp++;
			}elsif ($z1 eq $z2)	{
				 $bacon_homosnp++;
				$baconhomo{$key} = $bacongeno{$key};
			}elsif  ($z1 ne $z2)	{
				 $bacon_hetsnp++;
				$baconhet{$key} = $bacongeno{$key};
			}
} 

	foreach my $key (sort keys %baconuniq)	{
			my @baconuniq = split /\t/, $baconuniq{$key};
			my $z1 = $baconuniq[0];
			my $z2 = $baconuniq[1];
			if ($z1 eq "X")	{
				 $baconuniq_nosnp++;
			}elsif ($z1 eq $z2)	{
				 $baconuniq_homosnp++;
			}elsif ($z1 ne $z2)	{
				 $baconuniq_hetsnp++;
			}
} 
	foreach my $key (sort keys %hassgeno)	{
			(my $z1,my $z2) = split /\t/, $hassgeno{$key};
			if ($z1 eq "X")	{
				 $hass_nosnp++;
			}elsif ($z1 eq $z2)	{
				 $hass_homosnp++;
				$hasshomo{$key} = $hassgeno{$key};
			}elsif ($z1 ne $z2)	{
				 $hass_hetsnp++;
				$hasshet{$key} = $hassgeno{$key};
			}
} 
	foreach my $key (sort keys %hassuniq)	{
			my @hassuniq = split /\t/, $hassuniq{$key};
			my $z1 = $hassuniq[2];
			my $z2 = $hassuniq[3];
			if ($z1 eq "X")	{
				 $hassuniq_nosnp++;
			}elsif ($z1 eq $z2)	{
				 $hassuniq_homosnp++;
			}elsif ($z1 ne $z2)	{
				 $hassuniq_hetsnp++;
			}
}  

	foreach my $key (sort keys %simmondsgeno)	{
			(my $z1,my $z2) = split /\t/, $simmondsgeno{$key};
			if ($z1 eq "X")	{
				 $simmonds_nosnp++;
			}elsif ($z1 eq $z2)	{
				 $simmonds_homosnp++;
				$simmondshomo{$key} = $simmondsgeno{$key};
			}elsif ($z1 ne $z2)	{
				$simmonds_hetsnp++;
				$simmondshet{$key} = $simmondsgeno{$key};
			}
} 
	foreach my $key (sort keys %simmondsuniq)	{
			@simmondsuniq = split /\t/, $simmondsuniq{$key};
			my $z1 = $simmondsuniq[4];
			my $z2 = $simmondsuniq[5];
			if ($z1 eq "X")	{
				 $simmondsuniq_nosnp++;
			}elsif ($z1 eq $z2)	{
				 $simmondsuniq_homosnp++;
			}elsif ($z1 ne $z2)	{
				$simmondsuniq_hetsnp++;
			}
} 
	foreach my $key (sort keys %tonnagegeno)	{
			(my $z1,my $z2) = split /\t/, $tonnagegeno{$key};
			if ($z1 eq "X")	{
				$tonnage_nosnp++;
			}elsif ($z1 eq $z2)	{
				$tonnage_homosnp++;
				$tonnagehomo{$key} = $tonnagegeno{$key};
			}elsif ($z1 ne $z2)	{
				$tonnage_hetsnp++;
				$tonnagehet{$key} = $tonnagegeno{$key};
			}
} 

	foreach my $key (sort keys %tonnageuniq)	{
			@tonnageuniq = split /\t/, $tonnageuniq{$key};
			my $z1 = $tonnageuniq[6];
			my $z2 = $tonnageuniq[7];
			if ($z1 eq "X")	{
				$tonnageuniq_nosnp++;
			}elsif ($z1 eq $z2)	{
				$tonnageuniq_homosnp++;
			}elsif ($z1 ne $z2)	{
				$tonnageuniq_hetsnp++;
			}
} 
#####################################################################################################################################################

			

############################################################################################################################################
##########################################################################################################################################
#  Calculate unique SNPs by first calculating when three parents are homozygous and identical  January 25, 2013
	foreach my $key (sort keys %avosnp)	{
		if ((exists $baconhomo{$key}) && ( exists $hasshomo{$key}) && ( exists $simmondshomo{$key})&& ( exists $tonnagehomo{$key}) && ($baconhomo{$key} eq $hasshomo{$key}) && ($hasshomo{$key} eq $simmondshomo{$key}) && ($tonnagehomo{$key} ne $baconhomo{$key} ))	{ 
		$tonnage_unique_homo{$key} = "$baconhomo{$key}\t\t$tonnagehomo{$key}";	# creates a hash where only Tonnage unique SNPs.  Homozygous SNP of other three parents is shown first. 
		$t_homo_unik++;
		}elsif ((exists $baconhomo{$key}) && ( exists $hasshomo{$key}) && ( exists $simmondshomo{$key}) && ($baconhomo{$key} eq $hasshomo{$key}) && ($hasshomo{$key} eq $simmondshomo{$key}) && (exists $tonnagehet{$key} ))	{ 
		$tonnage_unique_het{$key} = "$baconhomo{$key}\t\t$tonnagehet{$key}";
		$t_het_unik++;
		}elsif ((exists $baconhomo{$key}) && ( exists $hasshomo{$key}) && (exists $simmondshomo{$key})&& ( exists $tonnagehomo{$key}) && ($baconhomo{$key} eq $hasshomo{$key}) && ($hasshomo{$key} eq $tonnagehomo{$key}) && ($simmondshomo{$key} ne $baconhomo{$key}))	{ 
		$simmonds_unique_homo{$key} = "$baconhomo{$key}\t\t$simmondshomo{$key}";
		$s_homo_unik++;
		}elsif ((exists $baconhomo{$key}) && ( exists $hasshomo{$key}) && ( exists $tonnagehomo{$key}) && ($baconhomo{$key} eq $hasshomo{$key}) && ($hasshomo{$key} eq $tonnagehomo{$key}) && (exists $simmondshet{$key}))	{ 
		$simmonds_unique_het{$key} = "$baconhomo{$key}\t\t$simmondshet{$key}";
		$s_het_unik++;
		}elsif ((exists $baconhomo{$key}) && (exists $hasshomo{$key}) && ( exists $simmondshomo{$key}) && ( exists $tonnagehomo{$key}) && ($baconhomo{$key} eq $simmondshomo{$key}) && ($simmondshomo{$key} eq $tonnagehomo{$key}) && ($hasshomo{$key} ne $baconhomo{$key}))	{ 
		$hass_unique_homo{$key} = "$baconhomo{$key}\t\t$hasshomo{$key}";
		$h_homo_unik++;
		}elsif ((exists $baconhomo{$key}) && ( exists $simmondshomo{$key}) && ( exists $tonnagehomo{$key}) && ($baconhomo{$key} eq $simmondshomo{$key}) && ($simmondshomo{$key} eq $tonnagehomo{$key}) && (exists $hasshet{$key}))	{ 
		$hass_unique_het{$key} = "$baconhomo{$key}\t\t$hasshet{$key}";
		$h_het_unik++;
		}elsif ((exists $baconhomo{$key}) && ( exists $hasshomo{$key}) && ( exists $simmondshomo{$key}) && (exists $tonnagehomo{$key}) &&($hasshomo{$key} eq $simmondshomo{$key}) && ($simmondshomo{$key} eq $tonnagehomo{$key}) && ($baconhomo{$key} ne $hasshomo{$key}))	{
		$bacon_unique_homo{$key} = "$hasshomo{$key}\t\t$baconhomo{$key}"; 
		$b_homo_unik++; 
		}elsif (( exists $hasshomo{$key}) && ( exists $simmondshomo{$key}) && (exists $tonnagehomo{$key}) &&($hasshomo{$key} eq $simmondshomo{$key}) && ($simmondshomo{$key} eq $tonnagehomo{$key}) && (exists $baconhet{$key}))	{
		$bacon_unique_het{$key} = "$hasshomo{$key}\t\t$baconhet{$key}"; 
		$b_het_unik++; 
		}
	}
##############################################################################################################################################
			my $bacon_total_snp = $bacon_nosnp + $bacon_homosnp + $bacon_hetsnp;
			my $bacon_percent_het = $bacon_hetsnp/$bacon_total_snp;
			my $hass_total_snp = $hass_nosnp + $hass_homosnp + $hass_hetsnp;
			my $hass_percent_het = $hass_hetsnp/$hass_total_snp;
			my $simmonds_total_snp = $simmonds_nosnp + $simmonds_homosnp + $simmonds_hetsnp;
			my $simmonds_percent_het = $simmonds_hetsnp/$simmonds_total_snp;
			my $tonnage_total_snp = $tonnage_nosnp + $tonnage_homosnp + $tonnage_hetsnp;
			my $tonnage_percent_het = $tonnage_hetsnp/$tonnage_total_snp;	

################################################################################################################################################
#  Print out files with unique SNPs for each parent January 25, 2013
	open (BACON, ">>fixit_bacon_unique_SNPs.out") || die "Can't open fixit_bacon_unique_SNPs.out $!\n"; 
		printf BACON	"The number of unique homo SNPs for Bacon is $b_homo_unik\n The percentage of homo SNPs that are unique is %6.3f%% .\n The number of unique het SNPs for Bacon is $b_het_unik and the percentage of het SNPs that are unique is %6.3f%%.\nThe overall percent heterozygosity in Bacon is %6.3f%%\nThe number of Bacon nosnps is $bacon_nosnp\t%6.3f%%\nThe number of Bacon homosnps is $bacon_homosnp\t%6.3f%%\nThe number of Bacon hetsnps is $bacon_hetsnp\t%6.3f%%\nThe total number of Bacon snps is $bacon_total_snp .\nSNP name\t Hass Simmonds Tonnage genotype \t Bacon homogenotype\n", 100 * $b_homo_unik/$bacon_homosnp, 100 * $b_het_unik/$bacon_hetsnp, 100 * $bacon_percent_het,100 * $bacon_nosnp/$bacon_total_snp, 100 * $bacon_homosnp/$bacon_total_snp, 100 * $bacon_hetsnp/$bacon_total_snp;
		foreach my $key (sort keys %bacon_unique_homo)		{
			print BACON "$key\t$bacon_unique_homo{$key}\n";
		}
		print BACON	"SNP name\t Hass Simmonds Tonnage genotype \t Bacon hetgenotype\n";
		foreach my $key (sort keys %bacon_unique_het)		{
			print BACON "$key\t$bacon_unique_het{$key}\n";
		}
	close BACON;

	open (HASS, ">>fixit_hass_unique_SNPs.out") || die "Can't open fixit_hass_unique_SNPs.out $!\n";
		printf HASS	"The number of unique homo SNPs for Hass is $h_homo_unik\nThe percentage homo SNPs that are unique is %6.3f%%.\n The number of unique het SNPs for Hass is $h_het_unik .\n The percentage het SNPs that are unique is %6.3f%%\nThe overall percent heterozygosity in Hass is %6.3f%%\nThe number of Hass nosnps is $hass_nosnp\t%6.3f%%\nThe number of Hass homosnps is $hass_homosnp\t%6.3f%%\nThe number of Hass hetsnps is $hass_hetsnp\t%6.3f%%\nThe total number of Hass snps is $hass_total_snp.\nSNP name\tBacon Simmonds Tonnage genotype \t Hass homogenotype\n", 100 * $h_homo_unik/$hass_homosnp, 100 * $h_het_unik/$hass_hetsnp, 100 * $hass_percent_het,100 * $hass_nosnp/$hass_total_snp,100 * $hass_homosnp/$hass_total_snp,100 * $hass_hetsnp/$hass_total_snp ;
		foreach my $key (sort keys %hass_unique_homo)		{
			print HASS "$key\t$hass_unique_homo{$key}\n";
		}
		print HASS "SNP name\t Bacon Simmonds Tonnage genotype \t Hass hetgenotype";
		foreach my $key (sort keys %hass_unique_het)		{
			print HASS "$key\t$hass_unique_het{$key}\n";
		}
	close HASS;

	open (SIMMONDS, ">>fixit_simmonds_unique_SNPs.out") || die "Can't open fixit_simmonds_unique_SNPs.out $!\n";
		printf SIMMONDS	"The number of unique homo SNPs for SIMMONDS is $s_homo_unik.\nThe percentage of unique homo SNPs for Simmonds is %6.3f%%.\nThe number of unique het SNPs for Simmonds is $s_het_unik.\n The percentage of unique het SNPs for Simmonds is %6.3f%%\nThe percent heterozygosity in Simmonds is %6.3f%%\nThe number of Simmonds nosnps is $simmonds_nosnp\t%6.3f%%\nThe number of Simmonds homosnps is $simmonds_homosnp\t%6.3f%%\nThe number of Simmonds hetsnps is $simmonds_hetsnp\t%6.3f%%\nThe total number of Simmonds snps is $simmonds_total_snp.\nSNP name\t Bacon Hass Tonnage genotype \t Simmonds homogenotype \n", ($s_homo_unik/$simmonds_homosnp)*100, 100 * $s_het_unik/$simmonds_hetsnp, 100 *$simmonds_percent_het,100 * $simmonds_nosnp/$simmonds_total_snp,100 * $simmonds_homosnp/$simmonds_total_snp,100 * $simmonds_hetsnp/$simmonds_total_snp ;
		foreach my $key (sort keys %simmonds_unique_homo)		{
			print SIMMONDS "$key\t$simmonds_unique_homo{$key}\n";
		}
		print SIMMONDS "SNP name\tBacon, Hass, Tonnage genotype\t Simmonds hetgenotype";
		foreach my $key (sort keys %simmonds_unique_het)		{
			print SIMMONDS "$key\t$simmonds_unique_het{$key}\n";
		}
	close SIMMONDS;

	open (TONNAGE, ">>fixit_tonnage_unique_SNPs.out") || die "Can't open fixit_tonnage_unique_SNPs.out $! \n";
		printf TONNAGE "The number of unique homo SNPs for Tonnage is $t_homo_unik .\nThe percentage of unique homo SNPs for Tonnage is %6.3f%%\nThe number of unique het SNPs for Tonnage is $t_het_unik.\n The percentage of unique het SNPs for Tonnage is %6.3f%%\nThe percent heterozygosity in Tonnage is %6.3f%%\nThe number of Tonnage nosnps is $tonnage_nosnp\t%6.3f%%\nThe number of Tonnage homosnps is $tonnage_homosnp\t%6.3f%%\nThe number of Tonnage hetsnps is $tonnage_hetsnp\t%6.3f%%\nThe total number of Tonnage snps is $tonnage_total_snp\nSNP name\tBacon Hass Simmonds genotype \t Tonnage homogenotype \n", ($t_homo_unik/$tonnage_homosnp)*100, ($t_het_unik/$tonnage_hetsnp)*100, $tonnage_percent_het * 100, ($tonnage_nosnp/$tonnage_total_snp) * 100, ($tonnage_homosnp/$tonnage_total_snp) * 100, ($tonnage_hetsnp/$tonnage_total_snp) * 100 ;
		foreach my $key (sort keys %tonnage_unique_homo)		{
			print TONNAGE " $key\t$tonnage_unique_homo{$key} \n";
		}
			print TONNAGE "SNP name\t Bacon, Hass, Simmonds genotype\t Tonnage hetgenotype";
		foreach my $key (sort keys %tonnage_unique_het)		{
			print TONNAGE " $key\t$tonnage_unique_het{$key} \n";
		}
	close TONNAGE;

#####################################################################################################################################################

############################################################################################################################################


	open (OUTFILE, ">>fixit_percent_hetout.txt") || die "Can't open outfile";

		
		printf OUTFILE "The percent heterozygosity in Bacon is %6.3f%%\nThe number of Bacon nosnps is $bacon_nosnp\t%6.3f%%\nThe number of Bacon homosnps is $bacon_homosnp\t%6.3f%%\nThe number of Bacon hetsnps is $bacon_hetsnp\t%6.3f%%\nThe total number of Bacon snps is $bacon_total_snp\n. The percent heterozygosity in Hass is %6.3f%%\nThe number of Hass nosnps is $hass_nosnp\t%6.3f%%\nThe number of Hass homosnps is $hass_homosnp\t%6.3f%%\nThe number of Hass hetsnps is $hass_hetsnp\t%6.3f%%\nThe total number of Hass snps is $hass_total_snp\n.The percent heterozygosity in Simmonds is %6.3f%%\nThe number of Simmonds nosnps is $simmonds_nosnp\t%6.3f%%\nThe number of Simmonds homosnps is $simmonds_homosnp\t%6.3f%%\nThe number of Simmonds hetsnps is $simmonds_hetsnp\t%6.3f%%\nThe total number of Simmonds snps is $simmonds_total_snp\n.The percent heterozygosity in Tonnage is %6.3f%%\nThe number of Tonnage nosnps is $tonnage_nosnp\t%6.3f%%\nThe number of Tonnage homosnps is $tonnage_homosnp\t%6.3f%%\nThe number of Tonnage hetsnps is $tonnage_hetsnp\t%6.3f%%\nThe total number of Tonnage snps is $tonnage_total_snp\n.", 100 * $bacon_percent_het,100 * $bacon_nosnp/$bacon_total_snp, 100 * $bacon_homosnp/$bacon_total_snp, 100 * $bacon_hetsnp/$bacon_total_snp, 100 * $hass_percent_het,100 * $hass_nosnp/$hass_total_snp,100 * $hass_homosnp/$hass_total_snp,100 * $hass_hetsnp/$hass_total_snp, 100 * $simmonds_percent_het,$simmonds_nosnp/$simmonds_total_snp,$simmonds_homosnp/$simmonds_total_snp,$simmonds_hetsnp/$simmonds_total_snp, $tonnage_percent_het, $tonnage_nosnp/$tonnage_total_snp, 100 * $tonnage_homosnp/$tonnage_total_snp, 100 * $tonnage_hetsnp/$tonnage_total_snp;
	
###############################################################################################################################################
			foreach my $key (sort keys %avosnp)	{
				if ((exists $baconhomo{$key}) && (exists $simmondshomo{$key}) && ($baconhomo{$key} ne $simmondshomo{$key}))	{
					$MneWI{$key} = 	"$baconhomo{$key}\t$simmondshomo{$key}";
					$MneWI++;
				}
			}			

	open (MneWI, ">>fixit_MneWI_out.txt") || die "Can't open MneWI outfile $!\n";
			my $totalMneWIsnps = keys %MneWI;
			print MneWI	"Total number of SNPs that are homozygous in Bacon for one allele and homozygous in Simmonds for the other allele is $totalMneWIsnps\n";
			print MneWI	"The following are the snpnames and genotypes for those snps, Bacon genotype first.\n";
			foreach my $key (sort keys %MneWI)	{
			print MneWI	"$key\t$MneWI{$key}\n";
			}
#################################################################################################################################################
	open (OUTFILE, ">>fixit_MneWI_isotigs_snps") || die "Can't open MneWI isotig snp outfile $!\n";

			foreach my $isotigname (sort keys %hoh)	{
				foreach my $snpname (sort keys %{$hoh{$isotigname}})	{
					if (exists $MneWI{$snpname})	{
						$hoh2{$isotigname}{$snpname} = $MneWI{$snpname};
						print OUTFILE "$isotigname\t$snpname\t$MneWI{$snpname}\n";
					}
				}
			}
	open (HOH2, ">>fixit_summary_MneWI_isotigs_snps.out") || die "Can't open MneWI summary file $!\n";
		my $hoh2count = keys %hoh2;
		print HOH2 "There were $hoh2count isotigs that had MneWI snps\n";
		print HOH2 "Isotig\tMneWIsnps\tTotal SNPs\tRatio\n";
			foreach my $key (sort keys %hoh2)		{
				my $keycount = keys (%{$hoh2{$key}});
				my $keycounttotal = keys (%{$hoh{$key}});
				my $ratio = $keycount/$keycounttotal;
				printf HOH2 "$key\t$keycount\t$keycounttotal\t%6.3f%%\n", 100 * $ratio;
		}

sub genotype	{	#determines genotype of SNP (XX if not enough reads) (homozygous major nucleotide) or heterozygous or QQ if there are blanks
	my @x = @{$_[0]};
	if ($x[4] < 6)	{	#if total reads are less than 6, score as no data
		my $genotype = "X\tX";
	}elsif (($x[5] eq "-") || ($x[8] eq "-"))	{
		my $genotype = "X\tX";
	} elsif ($x[5] eq "\.")	{	#if total reads more than 6 and a . in most common nucleotide, score as homozygous most common nucleotide
		my $genotype = "$x[1]\t$x[1]";	#genotype is most frequent nucleotide
	} elsif (($x[5] ne "\.") && ($x[8] eq ".")) 	{ #if $x[5] doesn't have a dot and $x[8] does, genotype is homozygous $x[5]
		my $genotype = "$x[5]\t$x[5]";
	} elsif (($x[5] ne "\.") && ($x[8] ne "\.")) 	{ #if maj allele is not dot and minor allele is not dot, genotype is het
		my @geno = ($x[5], $x[8]);
		my @sorted = sort @geno;
		my $genotype = "$sorted[0]\t$sorted[1]";		
	}
	
}
