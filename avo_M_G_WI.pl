#!/usr/bin/perl -w
##################################################################################
# Author		: David Kuhn
# Date			: 5/9/2012
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: 	reads in avocado variant report, counts number of heterozygous SNPs for each parent of avocado mapping populations and outputs counts, calculates number of SNPs
#			 unique to each cultivar, identifies SNPs that are homozygous in Bacon for one allele and homozygous in Simmonds for the other allele.
# <input>		:  KK_avo_snp_lzsnpname			
# <output>		:   QQ_percent_hetout.txt 
# 			
# Dependencies	: N/A
# Usage			: 
#  Still under development (6/22/12).                    
##################################################################################
#use strict;
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
		$isotigname = $snp[36];
		 @bacon = ($snp[1], $snp[2], $snp[3], $snp[4], $snp[5], $snp[6], $snp[7], $snp[8], $snp[9], $snp[10], $snp[11]);
#		print "Bacon array is @bacon\n";
		 @hass = ($snp[1], $snp[2], $snp[3], $snp[12], $snp[13], $snp[14], $snp[15], $snp[16], $snp[17], $snp[18], $snp[19]);
		 @simmonds = ($snp[1], $snp[2], $snp[3], $snp[20], $snp[21], $snp[22], $snp[23], $snp[24], $snp[25], $snp[26], $snp[27]);
		 @tonnage = ($snp[1], $snp[2], $snp[3], $snp[28], $snp[29], $snp[30], $snp[31], $snp[32], $snp[33], $snp[34], $snp[35]);
	 $bacongeno{$snpname} = &genotype (\@bacon);
	$hassgeno{$snpname} = &genotype (\@hass);
	 $simmondsgeno{$snpname} = &genotype (\@simmonds);
	 $tonnagegeno{$snpname} = &genotype (\@tonnage);
	$hoh{$isotigname}{$snpname} = "$bacongeno{$snpname}\t$hassgeno{$snpname}\t$simmondsgeno{$snpname}\t$tonnagegeno{$snpname}";	#hash with isotig name key and snpname key for sorting #	on both sets of keys, values are all genotypes for bacon, hass, simmonds, tonnage
#	should make a database of the genotypes for each isotig => snpname
			}
	}
}
#########################################################################################################################################
open	(AVOGENO, ">>PPP_avo_genotypes_out") || die "Can't open avo genotypes out file $!\n";
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
	open (UNIQSNP, ">>ZZZ_uniqsnp.out") || die "Can't open uniqsnp outfile $!";
		print UNIQSNP "The total number of Bacon unique SNPs is $baconuniq\nThe total number of Hass unique SNPs is $hassuniq\nThe total number of Simmonds unique SNPs is $simmondsuniq\nThe total number of Tonnage unique SNPs is $tonnageuniq\nThe total number of not unique SNPs is $notuniq\nThe total number of SNPs analyzed is $totaluniq\n";

	foreach my $key (sort keys %baconuniq)	{
			print UNIQSNP "bacon geno $key\t$baconuniq{$key}\n";
		}
	foreach my $key (sort keys %hassuniq)	{
			print UNIQSNP "hass geno $key\t$hassuniq{$key}\n";
		}
	foreach my $key (sort keys %simmondsuniq)	{
			print UNIQSNP "simmonds geno $key\t$simmondsuniq{$key}\n";
		}
	foreach my $key (sort keys %tonnageuniq)	{
			print UNIQSNP "tonnage geno $key\t$tonnageuniq{$key}\n";
		}
	foreach my $key (sort keys %notuniq)	{
			print UNIQSNP "notuniq $key\t$notuniq{$key}\n";
	}
#########################################################################################################################	
	foreach my $key (sort keys %bacongeno)	{
			(my $z1,my $z2) = split /\t/, $bacongeno{$key};
#			print "for $key, z1 is $z1, z2 is $z2\n";
			if ($z1 =~ /X|Q/)	{
				 $bacon_nosnp++;
				$bacongeno{nosnp} = $bacon_nosnp;
			}elsif ($z1 eq $z2)	{
				 $bacon_homosnp++;
				$bacongeno{homosnp} = $bacon_homosnp;
				$baconhomo{$key} = $bacongeno{$key};
			}else	{
				 $bacon_hetsnp++;
				$bacongeno{hetsnp} = $bacon_hetsnp;
				$baconhet{$key} = $bacongeno{$key};
			}
} 

	foreach my $key (sort keys %baconuniq)	{
			my @baconuniq = split /\t/, $baconuniq{$key};
			my $z1 = $baconuniq[0];
			my $z2 = $baconuniq[1];
			if ($z1 =~ /X|Q/)	{
				 $baconuniq_nosnp++;
			}elsif ($z1 eq $z2)	{
				 $baconuniq_homosnp++;
			}else	{
				 $baconuniq_hetsnp++;
			}
} 
	foreach my $key (sort keys %hassgeno)	{
			(my $z1,my $z2) = split /\t/, $hassgeno{$key};
			if ($z1 =~ /X|Q/)	{
				 $hass_nosnp++;
				$hassgeno{nosnp} = $hass_nosnp;
			}elsif ($z1 eq $z2)	{
				 $hass_homosnp++;
				$hassgeno{homosnp} = $hass_homosnp;
				$hasshomo{$key} = $hassgeno{$key};
			}else	{
				 $hass_hetsnp++;
				$hassgeno{hetsnp} = $hass_hetsnp;
				$hasshet{$key} = $hassgeno{$key};
			}
} 
	foreach my $key (sort keys %hassuniq)	{
			my @hassuniq = split /\t/, $hassuniq{$key};
			my $z1 = $hassuniq[2];
			my $z2 = $hassuniq[3];
			if ($z1 =~ /X|Q/)	{
				 $hassuniq_nosnp++;
			}elsif ($z1 eq $z2)	{
				 $hassuniq_homosnp++;
			}else	{
				 $hassuniq_hetsnp++;
			}
}  

	foreach my $key (sort keys %simmondsgeno)	{
			(my $z1,my $z2) = split /\t/, $simmondsgeno{$key};
			if ($z1 =~ /X|Q/)	{
				 $simmonds_nosnp++;
				$simmondsgeno{nosnp} = $simmonds_nosnp;
			}elsif ($z1 eq $z2)	{
				 $simmonds_homosnp++;
				$simmondsgeno{homosnp} = $simmonds_homosnp;
				$simmondshomo{$key} = $simmondsgeno{$key};
			}else	{
				$simmonds_hetsnp++;
				$simmondsgeno{hetsnp} = $simmonds_hetsnp;
				$simmondshet{$key} = $simmondsgeno{$key};
			}
} 
	foreach my $key (sort keys %simmondsuniq)	{
			@simmondsuniq = split /\t/, $simmondsuniq{$key};
			my $z1 = $simmondsuniq[4];
			my $z2 = $simmondsuniq[5];
			if ($z1 =~ /X|Q/)	{
				 $simmondsuniq_nosnp++;
			}elsif ($z1 eq $z2)	{
				 $simmondsuniq_homosnp++;
			}else	{
				$simmondsuniq_hetsnp++;
			}
} 
	foreach my $key (sort keys %tonnagegeno)	{
			(my $z1,my $z2) = split /\t/, $tonnagegeno{$key};
			if ($z1 =~ /X|Q/)	{
				$tonnage_nosnp++;
				$tonnagegeno{nosnp} = $tonnage_nosnp;
			}elsif ($z1 eq $z2)	{
				$tonnage_homosnp++;
				$tonnagegeno{homosnp} = $tonnage_homosnp;
				$tonnagehomo{$key} = $tonnagegeno{$key};
			}else	{
				$tonnage_hetsnp++;
				$tonnagegeno{hetsnp} = $tonnage_hetsnp;
				$tonnagehet{$key} = $tonnagegeno{$key};
			}
} 

	foreach my $key (sort keys %tonnageuniq)	{
			@tonnageuniq = split /\t/, $tonnageuniq{$key};
			my $z1 = $tonnageuniq[6];
			my $z2 = $tonnageuniq[7];
			if ($z1 =~ /X|Q/)	{
				$tonnageuniq_nosnp++;
			}elsif ($z1 eq $z2)	{
				$tonnageuniq_homosnp++;
			}else	{
				$tonnageuniq_hetsnp++;
			}
} 
##########################################################################################################################################
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

	open (UNIQUE, ">>YYY_type_of_unique_SNPs_out") || die "Can't open type of unique SNPs outfile $!\n";

	print UNIQUE  "Bacon unique nosnps is $baconuniq_nosnp\nBacon unique homosnps is $baconuniq_homosnp\nBacon unique heterosnps is $baconuniq_hetsnp\nHass unique nosnps is $hassuniq_nosnp\nHass unique homosnp is $hassuniq_homosnp\nHass unique hetsnp is $hassuniq_hetsnp\nSimmonds unique nosnp is $simmondsuniq_nosnp\nSimmonds unique homosnp is $simmondsuniq_homosnp\nSimmonds unique hetsnp is $simmondsuniq_hetsnp\nTonnage unique nosnp is $tonnageuniq_nosnp\nTonnage unique homosnp is $tonnageuniq_homosnp\nTonnage unique hetsnp is $tonnageuniq_hetsnp\n";
############################################################################################################################################


	open (OUTFILE, ">>QQQ_percent_hetout.txt") || die "Can't open outfile";

		
		printf OUTFILE "The percent heterozygosity in Bacon is %1.3f\nThe number of Bacon nosnps is $bacongeno{nosnp}\nThe number of Bacon hetsnps is $bacongeno{hetsnp}\nThe number of Bacon homosnps is $bacongeno{homosnp}\nThe total number of Bacon snps is $bacon_total_snp\n. The percent heterozygosity in Hass is %1.3f\nThe number of Hass nosnps is $hassgeno{nosnp}\nThe number of Hass hetsnps is $hassgeno{hetsnp}\nThe number of Hass homosnps is $hassgeno{homosnp}\nThe total number of Hass snps is $hass_total_snp\n.The percent heterozygosity in Simmonds is %1.3f\nThe number of Simmonds nosnps is $simmondsgeno{nosnp}\nThe number of Simmonds hetsnps is $simmondsgeno{hetsnp}\nThe number of Simmonds homosnps is $simmondsgeno{homosnp}\nThe total number of Simmonds snps is $simmonds_total_snp\n.The percent heterozygosity in Tonnage is %1.3f\nThe number of Tonnage nosnps is $tonnagegeno{nosnp}\nThe number of Tonnage hetsnps is $tonnagegeno{hetsnp}\nThe number of Tonnage homosnps is $tonnagegeno{homosnp}\nThe total number of Tonnage snps is $tonnage_total_snp\n.", $bacon_percent_het, $hass_percent_het, $simmonds_percent_het, $tonnage_percent_het;
	
###############################################################################################################################################
			foreach my $key (sort keys %avosnp)	{
				if ((exists $baconhomo{$key}) && (exists $simmondshomo{$key}) && ($baconhomo{$key} ne $simmondshomo{$key}))	{
					$MneWI{$key} = 	"$baconhomo{$key}\t$simmondshomo{$key}";
					$MneWI++;
				}
			}			

	open (MneWI, ">>KKK_MneWI_out.txt") || die "Can't open MneWI outfile $!\n";
			$totalMneWIsnps = keys %MneWI;
			print MneWI	"Total number of SNPs that are homozygous in Bacon for one allele and homozygous in Simmonds for the other allele is $totalMneWIsnps\n";
			print MneWI	"The following are the snpnames and genotypes for those snps, Bacon genotype first.\n";
			foreach my $key (sort keys %MneWI)	{
			print MneWI	"$key\t$MneWI{$key}\n";
			}
#################################################################################################################################################
	open (OUTFILE, ">>NNN_MneWI_isotigs_snps") || die "Can't open MneWI isotig snp outfile $!\n";

			foreach my $isotigname (sort keys %hoh)	{
				foreach my $snpname (sort keys %{$hoh{$isotigname}})	{
					if (exists $MneWI{$snpname})	{
						$hoh2{$isotigname}{$snpname} = $MneWI{$snpname};
						print OUTFILE "$isotigname\t$snpname\t$MneWI{$snpname}\n";
					}
				}
			}
	open (HOH2, ">>DDD_summary_MneWI_isotigs_snps.out") || die "Can't open MneWI summary file $!\n";
		$hoh2count = keys %hoh2;
		print HOH2 "There were $hoh2count isotigs that had MneWI snps\n";
		print HOH2 "Isotig\tMneWIsnps\tTotal SNPs\tRatio\n";
			foreach my $key (sort keys %hoh2)		{
				my $keycount = keys (%{$hoh2{$key}});
				my $keycounttotal = keys (%{$hoh{$key}});
				my $ratio = $keycount/$keycounttotal;
				printf HOH2 "$key\t$keycount\t$keycounttotal\t%1.3f\n", $ratio;
		}

sub genotype	{	#determines genotype of SNP (XX if not enough reads) (homozygous major nucleotide) or heterozygous or QQ if there are blanks
	my @x = @{$_[0]};
#	print "This is the subroutine array @x\n";
	if ($x[4] < 6)	{
		my $genotype = "X\tX";
	} elsif ($x[5]eq "\.")	{
		my $genotype = "$x[1]\t$x[1]";	#genotype is most frequent nucleotide
	} elsif ($x[7]> 0.80)	{
		my $genotype = "$x[5]\t$x[5]";	#genotype is nucleotide referenced in $x[5]
	} elsif	(($x[7] < 0.80) && ($x[10] > 0.20))	{
		my @geno = ($x[5], $x[8]);
		my @sorted = sort (@geno);
		my $genotype = "$sorted[0]\t$sorted[1]";	#genotype is heterozygous
	} else {
		my $genotype = "Q\tQ";
	}
}


