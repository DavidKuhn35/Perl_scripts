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
			if (m/^Br/ )  	{
				$Brpop{$name} = $pop;
				} elsif (m/^IxKP/ )	{
					$IxKPpop{$name} = $pop;
					}  elsif (m/^TAxKP/)	{
						$TAxKPpop{$name} = $pop;
						} elsif (m/^KP1243xKP/ )	{
							$KP1243xKPpop{$name} = $pop;
							} elsif(m/^CrxKP/ )	{
							$CrxKPpop{$name} = $pop;
							} elsif(m/^TomAtSelfs/)	{
								$TomAtSelfspop{$name} = $pop;
								} elsif(m/^HadenSelfs/ )	{
									$HadenSelfspop{$name} = $pop;
									} elsif(m/^SHRSGermplasm/) {
										$SHRSGermplasmpop{$name} = $pop;
										} elsif(m/^TA_Parent/ ) {
											$TA_Parentpop{$name} = $pop;
											}elsif(m/^Creeper_Parent/ ) {
												$Creeper_Parentpop{$name} = $pop;
												}elsif(m/^Haden_Parent/ ) {
													$Haden_Parentpop{$name} = $pop;
													}elsif(m/^Irwin_Parent/ ) {
														$Irwin_Parentpop{$name} = $pop;
														}elsif(m/^KP_Parent/ ) {
															$KP_Parentpop{$name} = $pop;
															}elsif(m/^KP1243_Parent/ ) {
																$KP1243_Parentpop{$name} = $pop;
															}
														}
													
												
	
	#####################################################################################################################
	

#####################################################################################################################
#  to get count for number of times an individual has been genotyped
&grepprint (\%Brpop, "Br");
	&grepprint (\%IxKPpop,"IxKP");
	&grepprint	( \%TAxKPpop, "TAxKP");
		&grepprint (\%KP1243xKPpop, "KP1243xKP");
			&grepprint(\%CrxKPpop, "CrxKP");
			&grepprint(\%TomAtSelfspop,"TomAtSelfs");
					&grepprint(\%HadenSelfspop,"HadenSelfs");
						&grepprint(\%SHRSGermplasmpop, "SHRSGermplasm");
#############################################################################################################################################################
# greps name and prints out count
sub grepprint	{
	my ($href, $pname) =@_;
	open (GREP, ">>test/greplog\_$date.txt") || die "Can't open grep log $!\n";
	foreach my $g (sort keys %$href) {
		my $count = `grep -c "$g"  /home/david/mango/mangomap/DTR_files_for_mapping_pops/STA1-8_144.csv`;
		print GREP "$$href{$g}\t$g\t$count\n";

	}
}











