#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: July 17, 2012
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: Count the number of isotigs greater than 500 nt in length, count the number of isotigs with read coverage but no SNPs, count number of isotigs with covered length
#		greater than a threshhold set as an argument
# <input>		:  avo5_isotig_cvg_SNP.stats		
# <output>		:   summary file 
# 			
# Dependencies	: N/A
# Usage			: avo_isotig_SNP_counter.pl avo5_isotig_cvg.stat 1000
#                    
##################################################################################
#	Reads in avocado SNP stat file for a parent, then counts number of isotigs greater or less than a threshhold, counts number of those isotigs with no SNPs
#
#
use strict;
use warnings;
use diagnostics;
use Bio::Seq; 
use Bio::SeqIO;
use Bio::SeqUtils;
use Bio::Tools::OddCodes;
use Bio::PrimarySeq;
use Statistics::Descriptive;
my %length = ();
my %snp = ();
my %cov_length = ();
my %percent_cov = ();
my %med_cov = ();
my %mean_cov = ();
my %std_dev = ();
my %min_cov = ();
my %max_cov = ();
my %isoname = ();
my @snp_freq_dis = ();
my @snpdistribution = ();
#####################################################################################################################
#	main  contains calls to subroutines to keep the organization simple!	
			&readin ($ARGV[0]);
			my @isotig_under = &isotig_under_500 (\%isoname, $ARGV[2]);
			my @snp_isotig = &isotig_w_no_snp (\%isoname, $ARGV[2]);
			@snp_freq_dis = &snp_freq_dis (\%isoname, $ARGV[2]);
			&print_stats (\@snp_freq_dis, $ARGV[1]);
			printf "The number of isotigs greater than $ARGV[2] nts is $isotig_under[1]\nThe number of isotigs over $ARGV[2] with no SNPs is $snp_isotig[1]\nThe percent no SNPs is %4.2f%%\n.  Number of isotigs for $ARGV[1] less than $ARGV[2] nts long is $isotig_under[0]\nThe number of isotigs under $ARGV[2] with no SNPs is $snp_isotig[0]\nThe percent no SNPs is %4.2f%%\n", (($snp_isotig[1]/$isotig_under[1]) * 100), (($snp_isotig[0]/$isotig_under[0]) * 100);
				
######################################################################################################################
#  subroutine to read data into a hashes based on columns in 
sub readin	{
	my $infile = $_[0];
	open (IN, "$infile") || die "Can't open infile with isotig coverage data $!\n";
		while (<IN>)	{
			chomp;
			if (m/^g/)	{	#if line has isotig in first column
			my @x = split /\t/;
			my $isotig = $x[0];
			$isoname{$isotig} = $isotig;
			$length{$isotig} = $x[1];
			$snp{$isotig} = $x[9];
			$cov_length{$isotig} = $x[2];
			$percent_cov{$isotig} = $x[3];
			$med_cov{$isotig} = $x[4];
			$mean_cov{$isotig} = $x[5];
			$std_dev{$isotig} = $x[6];
			$min_cov{$isotig} = $x[7];
			$max_cov{$isotig} = $x[8];
		}
	}
}
#######################################################################################################################
#  subroutine which counts number of isotigs under 500 nts in length
sub isotig_under_500	{
		my  %x = %{$_[0]};
		my $z = $_[1];
		my $gt_eq_500 = 0;
		my $less_than_500 = 0;
		foreach my $key (sort keys %x)	{
			if (($length{$key} < $z) && ($percent_cov{$key} == 100))	{
			$less_than_500++;
			}  elsif (($length{$key}>= $z) && ($percent_cov{$key} == 100)) {
			$gt_eq_500++;
 		}
	}
		my @y = ($less_than_500, $gt_eq_500);
}
#########################################################################################################################
##########################################################################################################################
#
sub isotig_w_no_snp		{
 	my (%x) = (%{$_[0]});
	my $k = $_[1];
	my $iso_no_snp_under_500 = 0;
	my $iso_no_snp_over_500 = 0;
	my @nosnp = (0,0);
		foreach my $key (sort keys %x)	{
			if (($length{$key} < $k) && ($percent_cov{$key} == 100) && ($snp{$key} == 0))	{
			$iso_no_snp_under_500++;
			} elsif (($length{$key}>= $k) && ($percent_cov{$key} == 100) && ($snp{$key} == 0)) {
			$iso_no_snp_over_500++;
			}
		}
	@nosnp = ($iso_no_snp_under_500, $iso_no_snp_over_500);
}
###############################################################################################################################
#  make snp frequency distribution array from subset of isotigs
sub  snp_freq_dis		{
		my (%x) = (%{$_[0]});
		my $k = $_[1];
			foreach my $key (sort keys %x)	{
				if (($length{$key}>= $k) && ($percent_cov{$key} == 100)) {
				my $snp_freq = $snp{$key}/$length{$key};
				push (@snpdistribution, $snp_freq);
			}
		}
		my @s = @snpdistribution;
}			
###################################################################################################################################
#  makes distribution of stats in array and prints them out
sub 	print_stats					{
	my @x = $_[0];
	my $name = $_[1];
	my @bin = (0.0005,0.001,0.0015,0.002,0.0025,0.003,0.0035,0.004,0.0045,0.005,0.0055,0.006,0.0065,0.007,0.0075,0.008,0.0085,0.009,0.0095,0.010,0.011,0.012,0.013,0.014,0.015,0.016,0.017,0.018,0.019,0.02,0.0225,0.025,0.0275,0.03,0.035,0.04,0.045,0.05,0.06,0.07,0.08,0.09,0.10,0.20);
open (OUTFILE, ">>$ARGV[1]\_snp_freq_stats") || die "Can't open $ARGV[1] snp_freq_stat outfile $!\n";
my $stat = Statistics::Descriptive::Full->new();
$stat->add_data(@snpdistribution);
$stat->get_data();
my %fsnpdistribution = $stat->frequency_distribution(\@bin);
for (sort {$a <=> $b} keys %fsnpdistribution) 	{
	print OUTFILE "$_\t$fsnpdistribution{$_}\n";
	}
}
####################################################################################################################################
