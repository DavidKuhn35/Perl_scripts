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
	my %Criollopop = ();
	my %Criollo = ();

	my $snp = 0;
	my $dna = 0;
	my $geno1 = 0;
	my $geno2 = 0;
	
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
			my ($pop, $name) = split /\t/ ; #split on tab
			if (m/^Criollo/ )  	{
				$Criollopop{$name} = $pop;
				} 
#			open (LOG, ">>Criollolog.txt") || die "Can't open Criollolog $!\n";										
			}									
#			foreach my $key (sort keys %Criollopop)	{
#				print LOG "$key\n";

#			}
#			close LOG;
	#####################################################################################################################
	#	main	Read in SNP  data from EP-1 csv file.
	open (SNPDATA,$ARGV[1]) || die "Can't open EP-1 SNP chip data $!\n";
		while (<SNPDATA>)	{
			chomp;
			if (m/^S\d+/)	{	#if line has unique SNP id in first column
				my @ep1 = split /,/;	#split on comma
				my $snp = $ep1[1];
				my $name = $ep1[4];
#				if ($name =~ m/\d+\s\w+/ ){
#					$name =~ s/(\d+)\s(\w+)/$1$2/;  #remove space from Australian names
#				}
				my ($s, $dummy) = split /-/, $ep1[0]; #get S number
				if (($name ne "No Call")&&($name ne "Invalid") && ($name ne "NTC") && ($name ne "") && ($name ne "BLANK") && ($name ne "blank"))	{ 

				 #checks sample name for NTC or blank
				 my $nc = $ep1[9];
				 	unless ($nc =~ m/:/) {
				 		 $nc = "Z:Z"; #turn any non genotype values into missing data
				 	}
				foreach my $n (sort keys %Criollopop)	{
					if ($name =~ m/$n/)	{

				
				my @allele = split /:/, $nc;	#get nucleotide genotypes and put in an array

				my @allesort = sort (@allele);	#sort genotypes ASCIIbetically A,C,G,T

				my $genotype = "$allesort[0]$allesort[1]"; #adds tab to genotype, perhaps this should occur at printing time
				
				$hoh1{$snp}{$name} = $genotype;	#hoh with SNP as first key and cultivar name key and genotype value

				print "$name\t$snp\t$genotype\t$hoh1{$snp}{$name}\n";
				$hoh2{$name}{$snp} = $genotype;	#hoh with cultivar name as first key, snp as key and genotype value
				$hoh3{$snp}{$s} = $genotype;	#hoh with ID name as first key, snp as key and genotype value need for multiple genotypes of parents
				$snp_hash{$snp} = $snp;	#hash with snp name key and snp name value
				$dna_hash{$name} = $s;	#hash with dna name key and dna s value
				$snum_hash{$s} = $name;
				my @snpXY = sort ($ep1[4],$ep1[5]);
				$snptype{$snp} = "$snpXY[0]$snpXY[1]";  #hash of snptypes for each locus
#				if (exists $Criollopop{$name}) {
#					$Criollo{$snp}{$name} = $genotype;
#				}
					}
				
				
			
	}		
	}
	
	}
					
	close SNPDATA;
################################################################################################################################################################ get pop genotypes for export
	&popprint (\%hoh1,"Criollo");
#	&popprint (\%hoh1,\%IxKPpop,"IxKP");
#	&popprint	(\%hoh1, \%TAxKPpop, "TAxKP");
#		&popprint (\%hoh1, \%KP1243xKPpop, "KP1243xKP");
#			&popprint(\%hoh1,\%CrxKPpop, "CrxKP");
#			&popprint(\%hoh1,\%TomAtSelfspop,"TomAtSelfs");
#					&popprint(\%hoh1,\%HadenSelfspop,"HadenSelfs");
#						&popprint(\%hoh1,\%SHRSGermplasmpop, "SHRSGermplasm");
#							&parprint(\%TA_Parent, "TA_Parent");
#								&parprint(\%Creeper_Parent,"Creeper_Parent");
#									&parprint(\%Haden_Parent,"Haden_Parent");
#										&parprint(\%Irwin_Parent,"Irwin_Parent");
#											&parprint(\%KP_Parent,"KP_Parent");
#												&parprint(\%KP1243_Parent,"KP1243_Parent");
	#####################################################################################################################
############################################################################################################################  subroutine to print out pop genotypes in wide format
sub popprint	{
my ($href,$pname) = @_; 
 my $header = "Snp\t";
 my $gen = ();
foreach my $b (sort keys %$href2)	{
	$header .= "$b\t";
}
open (POPEYE, ">>$pname/$pname\_$date\_newDTR\_rawgeno.txt") || die "Can't open $pname rawgeno out $!\n";
print POPEYE "$header\n";
foreach my $key (sort keys %snp_hash)	{
	foreach my $z (sort keys %$href2)	{
		$gen .= "$$href{$key}{$z}\t";
	}
	print POPEYE "$key\t$gen\n";
	$gen = ();
}
close POPEYE;
}
##############################################################################################################################################
#  print out parent genotypes
sub parprint	{
	my ($parref,$parname) = @_;
	my $header = "SNP\tgenotype";
	open (PAR, ">>$parname/$parname\_$date\_newDTR\_rawgeno.txt") || die "Can't open parent printout file $!\n";
	print PAR "$header\n";
	 	foreach my $b (sort keys %$parref) {
	 		print PAR "$b\t$$parref{$b}\n";

	}
close PAR;
}