#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: 
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: 
# <input>		:  		
# <output>		:   
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
my $snp = "bufu";
my %seen = ();
my %parent = ();
my $q = ();
my @z = ();
my @y = ();
my $madonna = undef;
my @d = ();
########################################################################################################
#	main
my $date = `date`;
	$date =~ s/\s+/_/g;
	$date =~ s/\:/_/g;

	open (PARENT, "$ARGV[0]") || die "Can't open parent genotype file $! \n"		;
		while (<PARENT>)	{
			chomp;
			if ((m/^Mi/) || (m/^SSK/) || (m/^mango/) || (m/^Contig/))	{
				my @x = split /\t/, $_;
			 	$snp = shift @x;
			 	foreach my $d (@x)	{
			 		if (($d ne "ZZ")&&(defined $d)) {
			 			unshift ( @d, $d); 
			 		}
			 	}

			$seen{$_}++ foreach (@d);	
			@z = sort by_score keys %seen;
			%seen = ();
		
#			@y = ();
#			open (LOG, ">>Consensus_log") || die "no log $!\n";
#			print LOG "$snp\t@x\t\t@d\t\t@z\t\t$z[0]\n";
#			$q = $z[0];
			
				@d = ();
#			if (($z[0] ne "ZZ") && (defined $z[0] )) 	{
#			 	$parent{$snp} = shift @z ;
			
#			 } elsif 
#			 	(($z[1] ne "ZZ") && ( defined $z[1] ))	{
			 		$parent{$snp} = $z[0] ;
#			 		open (LOG, ">>Consensus_log") || die "no log $!\n";
#			 		print LOG "$parent{$snp} is the value of zzero\n";
#			 		}else {
#			 			$parent{$snp} = $z[2];
#			 		}
					}
			 	}		 
			 

			 @z = ();
			
			 
			open (CONGENO,">>$ARGV[1]/$ARGV[1]\_consensus_genotype_$date.txt") || die "Can't open consensus genotype outfile $!\n";
			print CONGENO "SNP \tGenotype\n";

			foreach my $key (sort keys %parent)	{
				print CONGENO "$key\t$parent{$key}\n";
			}
		
			
######################################################################################################################
#   sub by_score
sub by_score {$seen{$b} <=> $seen{$a} or $a cmp $b }

#######################################################################################################################
#   

#########################################################################################################################
#

##########################################################################################################################
#

