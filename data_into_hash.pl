#!/usr/bin/perl -w
##################################################################################
# Author		: David Kuhn
# Date			: 4/25/2012
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: Loads data into hash with isotig names as keys.  Reverses and re-reverses the hash to identify unique isotigs
# <input>		:  tab delimited data file			
# <output>		:  output to screen  
# 			
# Dependencies	: N/A
# Usage			: 
#  PROGRAM STILL UNDER DEVELOPMENT.                  
##################################################################################
use strict;
use warnings;
use diagnostics;
my @data = ();
my %isotig = ();
my %y =  ();
my %x = ();
my %z = ();
my $isotig_count = 0;
my $uniq_isotig = 0;




#Read in the tab delimited data file and make a hash
open (DATA, "$ARGV[0]") || die "Can't open DATA infile\n";
	while (<DATA>)
{
	chomp;
	if (m/^g/)  	{	#enter identifier in first column for data
	 @data = split (/\t/);
	 $isotig{$data[0]} = $data[1];	# make a hash with isotig names as keys and AT.G names as values
	}
}

 $isotig_count = keys (%isotig);
print "the total number of isotigs is $isotig_count\n";

 %y = reverse (%isotig);
 %z = reverse (%y);
 $uniq_isotig = keys (%z);
print "the total number of uniq isotigs is $uniq_isotig\n";
