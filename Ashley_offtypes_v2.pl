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
	my %snp_hash = ();
	my %dna_hash = ();
	my %snptype = ();
	my $snp = 0;
	my $dna = 0;
	my $geno1 = 0;
	my $geno2 = 0;
	my %pops = ();
	my %popdir = ();
	my %popref = ();
	my %snp_geno = ();
	my $date = `date`;
	$date =~ s/\s+/_/g;
	$date =~ s/\:/_/g;
	###############################################################################################################################
	# make a hash of individuals in each population read in mapping_pops_individuals.txt from ~/mango/manmap
	open (POPS, $ARGV[0]) || die "Can't open pop file $!\n";
		while (<POPS>)	{
			chomp;
			my ($pop, $name) = split /\t/ ; #split on tab
			 $pops{$name} = $pop;
			 $popdir{$pop} = $pop;
			}
			foreach my $d (sort keys %popdir) {
				`mkdir $d`;
			}
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
#				my ($s, $dummy) = split /-/, $ep1[0]; #get S number
				if (($name ne "No Call")&&($name ne "Invalid") && ($name ne "NTC") && ($name ne "") && ($name ne "BLANK") && ($name ne "blank"))	{ 

				 #checks sample name for NTC or blank
				 my $nc = $ep1[9];
				 	unless ($nc =~ m/:/) {
				 		 $nc = "Z:Z"; #turn any non genotype values into missing data
				 	}
				my @allele = split /:/, $nc;	#get nucleotide genotypes and put in an array
				my @allesort = sort (@allele);	#sort genotypes ASCIIbetically A,C,G,T
				my $genotype = "$allesort[0]$allesort[1]"; #adds tab to genotype, perhaps this should occur at printing time
				
				$hoh1{$snp}{$name} .= "$genotype\t";	#hoh with SNP as first key and cultivar name key and genotype value, this hash contains the actual names
				$snp_hash{$snp} = $snp;	#hash with snp name key and snp name value
				my @snpXY = sort ($ep1[4],$ep1[5]);
				$snptype{$snp} = "$snpXY[0]$snpXY[1]";  #hash of snptypes for each locus
					}
		}		
	}
				close SNPDATA;
################################################################################################################################################################ get pop genotypes for export
# print out just the desired genotypes
	foreach my $x (sort keys %pops)	{	#get names of clones to look for
		foreach my $y (sort keys %hoh1)	{	#get snp from main data
			foreach my $z (sort keys %{$hoh1{$y}})	{ #go through individual names to get genotypes
				if ($z =~ m/$x/)	{	#if specific name matches the general name or TC#
					open (OUT, ">>$pops{$x}/$x") || die "Can't open output $!\n";	#need to mkdir for all pops above
						print OUT "$pops{$x}\t$x\t$y\t$z\t$hoh1{$y}{$z}\n";
						close OUT;
				}

			}

		}

	}

#####################################################################################################################
