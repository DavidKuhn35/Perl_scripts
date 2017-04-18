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
	my %hoh4 = ();
	my	%hoh5 = ();
	my %hoh6 = ();
	my %snp_hash = ();
	my %dna_hash = ();
	my %snum_hash = ();
	my %Br = ();
	my %TAxKP = ();
#	my %map_snp = ();
	my %IxKP = ();
	my %CrxKP = ();
	my %KP1243xKP = ();
	my %TomAtSelfs = ();
	my %HadenSelfs = ();
	my %SHRSGermplasm = ();
	my %TA_Parent = ();
	my %Creeper_Parent = ();
	my %Haden_Parent = ();
	my %Irwin_Parent = ();
	my %KP_Parent = ();
	my %KP1243_Parent = ();
		my %Brpop = ();
	my %TAxKPpop = ();
#	my %map_snp = ();
	my %IxKPpop = ();
	my %CrxKPpop= ();
	my %KP1243xKPpop = ();
	my %TomAtSelfspop = ();
	my %HadenSelfspop = ();
	my %SHRSGermplasmpop = ();
	my %TA_Parentpop = ();
	my %Creeper_Parentpop = ();
	my %Haden_Parentpop = ();
	my %Irwin_Parentpop = ();
	my %KP_Parentpop = ();
	my %KP1243_Parentpop = ();
	my %snptype = ();

	my $snp = 0;
	my $dna = 0;
	my $geno1 = 0;
	my $geno2 = 0;
	my %popdir = ();
	my %snp_geno = ();
	my %popo = ("Br","Br", "IxKP","IxKP", "TAxKP", "TAxKP", "KP1243xKP","KP1243xKP","CrxKP", "CrxKP", "TomAtSelfs", "TomAtSelfs", "HadenSelfs", "HadenSelfs", "SHRSGermplasm", "SHRSGermplasm", "TA_Parent", "TA_Parent", "Creeper_Parent", "Creeper_Parent" , "Haden_Parent", "Haden_Parent", "Irwin_Parent","Irwin_Parent" ,"KP_Parent", "KP_Parent", "KP1243_Parent", "KP1243_Parent",);
	###############################################################################################################################
	# make a hash of individuals in each population read in mapping_pops_individuals.txt from ~/mango/manmap
	my $date = `date`;
	$date =~ s/\s+/_/g;
	$date =~ s/\:/_/g;
	
	open (POPS, $ARGV[0]) || die "Can't open pop file $!\n";
		while (<POPS>)	{
			chomp;
			my ($pop, $cultname) = split /\t/ ; #split on tab
			
				$hoh5{$pop}{$cultname} = $pop;
			
			}
		foreach my $d (sort keys %hoh5) {
				`mkdir $d`;
			}						
	#####################################################################################################################################################
	#  Read in STA ordered SNP names
	open (STA, $ARGV[2]) || die "Can't open STA list$!\n";
		while(<STA>)	{
			chomp;
			my ($qsnp,$qsta,$qchip) = split /\t/;
			$hoh4{$qsta}{$qsnp} = $qchip;

		}
 	#####################################################################################################################
	#	main	Read in SNP  data from EP-1 csv file.
	open (SNPDATA,$ARGV[1]) || die "Can't open EP-1 SNP chip data $!\n";
		while (<SNPDATA>)	{
			chomp;
			if (m/^S\d+/)	{	#if line has unique SNP id in first column
				my @ep1 = split /,/;	#split on comma
				my $snp = $ep1[3];
				my $name = $ep1[6];
				if ($name =~ m/\d+\s\w+/ ){
					$name =~ s/(\d+)\s(\w+)/$1$2/;  #remove space from Australian names
				}
				my ($s, $dummy) = split /-/, $ep1[0]; #get S number
#				if (($name ne "No Call")&&($name ne "Invalid") && ($name ne "NTC") && ($name ne "") && ($name ne "BLANK") && ($name ne "blank"))	{ 

				 #checks sample name for NTC or blank
				 my $nc = $ep1[11];
				 unless ($nc =~ m/:/) {
				 		 $nc = "Z:Z"; #turn any non genotype values into missing data
				 }
				
				my @allele = split /:/, $nc;	#get nucleotide genotypes and put in an array

				my @allesort = sort (@allele);	#sort genotypes ASCIIbetically A,C,G,T

				my $genotype = "$allesort[0]$allesort[1]"; #adds tab to genotype, perhaps this should occur at printing time
				
				$hoh1{$snp}{$name} = $genotype;	#hoh with SNP as first key and cultivar name key and genotype value
				$hoh2{$name}{$snp} = $genotype;	#hoh with cultivar name as first key, snp as key and genotype value
				$hoh3{$snp}{$s} = $genotype;	#hoh with ID name as first key, snp as key and genotype value need for multiple genotypes of parents
				$snp_hash{$snp} = $snp;	#hash with snp name key and snp name value
				$dna_hash{$name} = $s;	#hash with dna name key and dna s value
				$snum_hash{$s} = $name;
				my @snpXY = sort ($ep1[4],$ep1[5]);
				$snptype{$snp} = "$snpXY[0]$snpXY[1]";  #hash of snptypes for each locus

												
												if(exists $hoh5{Haden_Parent}{$name}) {
														$hoh6{Haden_Parent}{$snp} .= "$genotype\t";
														}elsif(exists $hoh5{TA_Parent}{$name}) {
															$hoh6{TA_Parent}{$snp} .= "$genotype\t";
															}elsif(exists $hoh5{Creeper_Parent}{$name}) {
																$hoh6{Creeper_Parent}{$snp} .= "$genotype\t";
																}elsif(exists $hoh5{Irwin_Parent}{$name}) {
																	$hoh6{Irwin_Parent}{$snp} .= "$genotype\t";
																	}elsif(exists $hoh5{KP_Parent}{$name}) {
																		$hoh6{KP_Parent}{$snp}.= "$genotype\t";
																		}elsif(exists $hoh5{KP1243_Parent}{$name}) {
																			$hoh6{KP1243_Parent}{$snp}.= "$genotype\t";
																			}elsif(exists $hoh5{SHRSGermplasm}{$name}) {
																			$hoh6{SHRSGermplasm}{$snp}{$name} = "$genotype\t";
					}
				
				
			
			
	}
	}
	
					
	close SNPDATA;
################################################################################################################################################################ get pop genotypes for export
	&popprint (\%hoh1,\%{$hoh5{Br}},"Br");
	&popprint (\%hoh1,\%{$hoh5{IxKP}},"IxKP");
	&popprint	(\%hoh1, \%{$hoh5{TAxKP}}, "TAxKP");
		&popprint (\%hoh1, \%{$hoh5{KP1243xKP}}, "KP1243xKP");
			&popprint(\%hoh1,\%{$hoh5{CrxKP}}, "CrxKP");
			&popprint(\%hoh1,\%{$hoh5{TomAtSelfs}},"TomAtSelfs");
					&popprint(\%hoh1,\%{$hoh5{HadenSelfs}},"HadenSelfs");
						&popprint(\%hoh1,\%{$hoh5{SHRSGermplasm}}, "SHRSGermplasm");
							&parprint(\%{$hoh6{TA_Parent}}, "TA_Parent");
								&parprint(\%{$hoh6{Creeper_Parent}},"Creeper_Parent");
									&parprint(\%{$hoh6{Haden_Parent}},"Haden_Parent");
										&parprint(\%{$hoh6{Irwin_Parent}},"Irwin_Parent");
											&parprint(\%{$hoh6{KP_Parent}},"KP_Parent");
												&parprint(\%{$hoh6{KP1243_Parent}},"KP1243_Parent");
	#####################################################################################################################
############################################################################################################################  subroutine to print out pop genotypes in wide format
sub popprint	{
my ($href,$href2,$pname) = @_; 
 my $header = "SNP\t";
 my $gen = ();
 my $st = ();
foreach my $b (sort keys %$href2)	{
	$header .= "$b\t";
}
open (POPEYE, ">>$pname/$pname\_$date\_STA1-10\_rawgeno.txt") || die "Can't open $pname rawgeno out $!\n";
print POPEYE "$header\n";
foreach my $key (sort keys %hoh4)	{	#sort by STA
		foreach  $st (sort keys %{$hoh4{$key}})	{	#sort by SNP in STA
			foreach my $z (sort keys %$href2)	{	#sort by name in pop
				$gen .= "$$href{$st}{$z}\t";	#genotype of all members of pop for one snp
	}
	
	print POPEYE "$st\t$gen\n";
	$gen = ();
}
}
close POPEYE;
}
##############################################################################################################################################
#  print out parent genotypes
sub parprint	{
	my ($parref,$parname) = @_;
	my $header = "SNP\tgenotype";
	open (PAR, ">>$parname/$parname\_$date\_STA1-10\_rawgeno.txt") || die "Can't open parent printout file $!\n";
	print PAR "$header\n";
	 	foreach my $key (sort keys %hoh4)	{	#sort by STA
	 		foreach my $st (sort keys %{$hoh4{$key}})	{	#sort by SNP in STA
	 			print PAR "$st\t$$parref{$st}\n";	#
	 		}
	}
close PAR
}
#close PAR