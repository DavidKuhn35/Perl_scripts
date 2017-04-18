	#!/usr/bin/perl
	##################################################################################
	# Author		: David Kuhn
	# Date			: March 18, 2015
	# Version		: ver. 1.0
	# History		: ver. 1.0
	# Description	: 	Read in SNP data.  Create hashes of genotypes for each mapping population.  Export genotypes in Excel format, loci as rows, individuals as columns
	# <input>		:  SNP data as csv file (DTR export from EP-1), population members, parents	
	# <output>		: data for each population and all parents in format with loci as rows and individuals as columns
	# 			
	# Dependencies	: N/A
	# Usage			: EP-1_recode_SNP-genotypes_for_map_pops.pl
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
	my %hoh1 = ();
	my %hoh2 = ();
	my %hoh3 = ();
	my %snp_hash = ();
	my %dna_hash = ();
	my %Brazil = ();
	my %TAxKP = ();
	my %map_snp = ();
	my %IxKP = ();
	my %CrxKP = ();
	my %1243xKP = ();
	my $snp = 0;
	my $dna = 0;
	my $geno1 = 0;
	my $geno2 = 0;
	my $date = 
	#my $print_header = "SNP name\tSegregation type\tIrwin genotype\tKP genotype\n";
	#my %segregation_type = ();
	#my $print_header2 = "SNP name\t";
	my %snp_geno = ();
	my %popo = ("Br","Br", "IxKP","IxKP", "TAxKP", "TAxKP", "KP1243xKP","KP1243xKP","CrxKP", "CrxKP", "TomAtSelfs", "TomAtSelfs", "HadenSelfs", "HadenSelfs", "SHRSGermplasm", "SHRSGermplasm", "TA_Parent", "TA_Parent", "Creeper_Parent", "Creeper_Parent" , "Haden_Parent", "Haden_Parent", "Irwin_Parent","Irwin_Parent" ,"KP_Parent", "KP_Parent", "KP1243_Parent", "KP1243_Parent",);
	###############################################################################################################################
	# make a hash of individuals in each population
	my $date = `date`;
	$date =~ s/\s+//g;
	
	open (POPS, $ARGV[0]) || die "Can't open pop file $!\n";
		while (<POPS>)	{
			chomp;
			($pop, $name) = split /\t/ ; #split on tab
			if (m/^Br) $pop 	{
				$Brpop{$name} = $pop;
				} elsif (m/^IxKP/ $pop)	{
					$IxKPpop{$name} = $pop;
					}  elsif (m/^TAxKP/ $pop)	{
						$TAxKPpop{$name} = $pop;
						} elsif (m/^KP1243xKP/ $pop)	{
							$KP1243xKPpop{$name} = $pop;
							} elsif(m/^CrxKP/ $pop)	{
							$CrxKPpop{$name} = $pop;
							} elsif(m/^TomAtSelfs/ $pop)	{
								$TomAtSelfspop{$name} = $pop;
								} elsif(m/^HadenSelfs/ $pop)	{
									$HadenSelfspop{$name} = $pop;
									} elsif(m/^SHRSGermplasm/ $pop) {
										$SHRSGermplasmpop{$name} = $pop;
										} elsif(m/^TA_Parent/ $pop) {
											$TA_Parentpop{$name} = $pop;
											}elsif(m/^Creeper_Parent/ $pop) {
												$Creeper_Parentpop{$name} = $pop;
												}elsif(m/^Haden_Parent/ $pop) {
													$Haden_Parentpop{$name} = $pop;
													}elsif(m/^Irwin_Parent/ $pop) {
														$Irwin_Parentpop{$name} = $pop;
														}elsif(m/^KP_Parent/ $pop) {
															$KP_Parentpop{$name} = $pop;
															}elsif(m/^KP1243_Parent/ $pop) {
																$KP1243_Parentpop{$name} = $pop;
															}
														}
													
												
	
	#####################################################################################################################
	#	main	Read in SNP  data from EP-1 csv file.
	open (SNPDATA,$ARGV[1]) || die "Can't open avo SNP chip data $!\n";
		while (<SNPDATA>)	{
			chomp;
			if (m/^S/)	{	#if line has unique SNP id in first column
				my @ep1 = split /,/;	#split on comma
				$snp = $ep1[3];
				$name = $ep1[6];
				my ($s, $dummy) = split /-/, $ep1[0]; #get S number
				if (($name ne "No Call")&&($name ne "Invalid") && ($name ne "NTC") && ($name ne ""))	{ 

				 #checks sample name for NTC or blank
				my @allele = split /:/, $ep1[11];	#get nucleotide genotypes and put in an array
				my @allesort = sort (@allele);	#sort genotypes ASCIIbetically A,C,G,T
				my $genotype = "$allesort[0]$allesort[1]";
				
				$hoh1{$snp}{$name} = 	#hoh with SNP as first key and cultivar name key and genotype value
				$hoh2{$name}{$snp} = $genotype;	#hoh with cultivar name as first key, snp as key and genotype value
				$hoh3{$id}{$snp} = $genotype;	#hoh with ID name as first key, snp as key and genotype value need for multiple genotypes of parents
				$snp_hash{$snp} = $snp;	#hash with snp name key and snp name value
				$dna_hash{$name} = $s;	#hash with dna name key and dna s value
				$snum_hash{$s} = $name;
				my @snpXY = sort ($ep1[4],$ep1[5]);
				my $snptype{$snp} = "$snpXY[0]$snpXY[1]";  #hash of snptypes for each locus
				if (exists $Brpop{$name}) {
					$Br{$snp}{$name} = $genotype;
					}elsif(exists $IxKPpop{$name}) {
						$IxKP{$snp}{$name} = $genotype;
						}elsif(exists $TAxKPpop{$name}) {
							$TAxKP{$snp}{$name} = $genotype;
							}elsif(exists $KP1243xKPpop{$name}) {
								$KP1243xKP{$snp}{$name} = $genotype;
								}elsif(exists $CrxKPpop{$name}) {
									$CrxKP{$snp}{$name} = $genotype;
									}elsif(exists TomAtSelfspop{$name}) {
										$TomAtSelfs{$snp}{$name} = $genotype;
										}elsif(exists $HadenSelfspop{$name}) {
											$HadenSelfs{$snp}{$name} = $genotype;
												}elsif(exists $SHRSGermplasmpop{$name}) {
													$SHRSGermplasm{$snp}{$name} = $genotype;
													}elsif(exists $Haden_Parentpop{$s}) {
														$Haden_Parent{$snp}{$s} = $genotype;
														}elsif(exists $TA_Parent{$s}) {
															$TA_Parent{$snp}{$s} = $genotype;
															}elsif(exists $Creeper_Parentpop{$s}) {
																$Creeper_Parent{$snp}{$s} = $genotype;
																}elsif(exists $Irwin_Parentpop{$s}) {
																	$Irwin_Parent{$snp}{$s} = $genotype;
																	}elsif(exists $KP_Parentpop{$s}) {
																		$KP_Parent{$snp}{$s} = $genotype;
																		}elsif(exists $KP1243_Parentpop{$s}) {
																			$KP1243_Parent{$snp}{$s} = $genotype;
					}
				
				
			
			
	}
	}
	}
					
	close SNPDATA;
################################################################################################################################################################ get pop genotypes for export
	&popprint (\%Br,\%Brpop,"Br");
	&popprint (\%IxKP,\%IxKPpop,"IxKP");
	&popprint	(\%TAxKP, \%TAxKPpop, "TAxKP");
		&popprint (\%KP1243xKP, \%KP1243xKPpop, "KP1243xKP");
			&popprint(\%CrxKP,\%CrxKPpop, "CrxKP");
				&popprint(\%TomAtSelfs,\%TomAtSelfspop,"TomAtSelfs");
					&popprint(\%HadenSelfs,\%HadenSelfspop,"HadenSelfs");
						&popprint(\%SHRSGermplasm,\%SHRSGermplasmpop, "SHRSGermplasm");
							&popprint(\%TA_Parent,\%TA_Parentpop, "TA_Parent");
								&popprint(\%Creeper_Parent,\%Creeper_Parentpop,"Creeper_Parent");
									&popprint(\%Haden_Parent,\%Haden_Parentpop,"Haden_Parent");
										&popprint(\%Irwin_Parent,\%Irwin_Parentpop,"Irwin_Parent");
											&popprint(\%KP_Parent,\%KP_Parentpop,"KP_Parent");
												&popprint(\%KP1243_Parent,\%KP1243_Parentpop,"KP1243_Parent");
	#####################################################################################################################
	#	Go through all SNPs for both parents (Irwin and KP) and get mappable snps
			foreach my $key (sort keys %snp_hash)	{
				if ((exists $Irwin{$key}) && (exists $KP{$key}))	{
				unless ((($Irwin{$key} eq "AA") || ($Irwin{$key} eq "CC") || ($Irwin{$key} eq "GG") || ($Irwin{$key} eq "TT") || ($Irwin{$key} eq "--")) && (($KP{$key} eq "AA")|| ($KP{$key} eq "CC") || ($KP{$key} eq "GG") || ($KP{$key} eq "TT") || ($KP{$key} eq "--")))	{
			$map_snp{$key} = $key;
					}
		}	
	}
							
	######################################################################################################################
	#  create list of mappable snps with segregation types
		foreach my $key (sort keys %map_snp)	{	
	 		if ((($Irwin{$key} eq "AA") || ($Irwin{$key} eq "CC") || ($Irwin{$key} eq "GG") || ($Irwin{$key} eq "TT")) && (($KP{$key} ne "AA") && ($KP{$key} ne "CC") && ($KP{$key} ne "GG") && ($KP{$key} ne "TT")))	{
				$segregation_type{$key} = "\<nnxnp\>";  #parent1 x parent2 segregation type
	#			$segregation_type{$key} = "\<lmxll\>";	#parent2 x parent1 segregation type
			}elsif ((($Irwin{$key} ne "AA") && ($Irwin{$key} ne "CC") && ($Irwin{$key} ne "GG") && ($Irwin{$key} ne "TT")) && (($KP{$key} eq "AA") || ($KP{$key} eq "CC") || ($KP{$key} eq "GG") || ($KP{$key} eq "TT")))	{
				$segregation_type{$key} = "\<lmxll\>";	#parent1 x parent2 segregation type
	#			$segregation_type{$key} = "\<nnxnp\>";	#parent2 x parent1 segregation type
			}elsif ((($Irwin{$key} ne "AA") && ($Irwin{$key} ne "CC") && ($Irwin{$key} ne "GG") && ($Irwin{$key} ne "TT")) && (($KP{$key} ne "AA") && ($KP{$key} ne "CC") && ($KP{$key} ne "GG") && ($KP{$key} ne "TT")))	{
				$segregation_type{$key} = "\<hkxhk\>";
			}
		}


	##########################################################################################################################
	#	Print outfile with all genotype data for each parent and snp name and segregation type
	open (SNP_OUT, ">>irwin_KP_parent_genotypes_May29_2013.txt") || die "Can't open irwin parent genotypes outfile $!\n";
			print SNP_OUT "$print_header";
			foreach my $x (sort keys %map_snp)	{
				print SNP_OUT "$x\t$segregation_type{$x}\t$Irwin{$x}\t$KP{$x}\n";
				}
	close SNP_OUT;

	###############################################################################################################################
	#  Print outfile with all mappable genotype data for all cultivars on EP-1 chip
		foreach my $key (sort keys %hoh2)	{
			$print_header2 .= "$key\t";	#put name of all cultivars in header
			foreach my $x (sort keys %map_snp)		{
					$snp_geno{$x} .= "$hoh2{$key}{$x}\t";
				}
			}

	open (GENO_OUT, ">>all_geno_out_May30_2013.txt") || die "Can't open all geno out $!\n";
		print GENO_OUT "$print_header2";
		foreach my $key (sort keys %snp_geno)		{
			print GENO_OUT	"$key\t$snp_geno{$key}\n";
		}
	close GENO_OUT;
############################################################################################################################  subroutine to print out pop genotypes in wide format
sub popprint	{
my ($href,$href2,$pname) = @_; 
 $header = "Snp name\t";
foreach my $b (sort keys %$href2)	{
	$header .= "$b\t";
}
open (POPEYE, ">>$pname\_$date\_rawgeno.txt") || die "Can't open $pname rawgeno out $!\n";
print POPEYE "$header\n";
foreach my $key (sort keys %snp_hash)	{
	foreach my $z (sort keys %href2)	{
		my $gen .= "$href{$key}{$z}";
	}
	print POPEYE "$key\t$gen\n";
}

}