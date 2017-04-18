#!/usr/bin/perl -w
##################################################################################
# Author		: David Kuhn
# Date			: 5/25/2012
# Version		: ver. 1.0.1
# History		: ver. 1.0
# Description	: reads in the GOSlim data file "ATH_GO_GOSLIM_ref_date_sorted.txt" which has GOSlim annotation but with the date converted to a number, puts most recent GOSlim 
# annotation on a single line with place holders so it is a complete data array and prints it out
# <input>		:  ATH_GO_GOSLIM_ref_date_sorted.txt	name for output file		
# <output>		:   XX$ARGV[1]\_cfp_goslim_out.txt 
# 			
# Dependencies	: N/A
# Usage			: perl goslim_reduced.pl ATH_GO_GOSLIM_ref_date_sorted.txt	name for output file
#                    
##################################################################################
use strict;
use warnings;
use diagnostics;
my @data = ();
my %goslim = ();
my @yoyo = ();
my $clastdate = 0;
my $flastdate = 0;
my $plastdate = 0;
my $lastAGI = ();
#Read in the tab delimited data file and make a hash
open (DATA, "$ARGV[0]") || die "Can't open DATA infile\n";
	while (<DATA>)
{
	chomp;
	if (m/^AT/)  	{	#enter identifier in first column for data
	 @data = split (/\t/);
#	my $x = join ":", @data;
#	$goslim{$data[0]}.= $x;
	$goslim{$data[0]}{AGI} = $data[0];
			if (($data[7] eq "C") && ($data[14] > $clastdate)) 	{ #identify which type of goslim term (C,F,P) and the date of the annotation
			$goslim{$data[0]}{"Crelat"} = $data[3];
			$goslim{$data[0]}{"Cgoterm"} = $data[4];
			$goslim{$data[0]}{"Cgoid"} = $data[5];
			$goslim{$data[0]}{"Cgoslim"} = $data[8];
			$goslim{$data[0]}{"Cdate"} = $data[14];
			my $C = join "|", ($data[7],$data[14], $data[8], $data[5], $data[3], $data[4]); #saves the latest annotation for component as a joined array
			$goslim{$data[0]}{"Carray"} = $C;
	}		elsif	(($data[7] eq "F") && ($data[14] > $flastdate))	{	#identify which type of goslim term (C,F,P) and the date of the annotation
			$goslim{$data[0]}{"Frelat"} = $data[3];
			$goslim{$data[0]}{"Fgoterm"} = $data[4];
			$goslim{$data[0]}{"Fgoid"} = $data[5];
			$goslim{$data[0]}{"Fgoslim"} = $data[8];
			$goslim{$data[0]}{"Fdate"} = $data[14];
			my $F = join "|", ($data[7],$data[14], $data[8], $data[5], $data[3], $data[4]); #saves the latest annotation for function as a joined array
			$goslim{$data[0]}{"Farray"} = $F;

	}		elsif (($data[7] eq "P") && ($data[14] > $plastdate))	{    #identify which type of goslim term (C,F,P) and the date of the annotation
			$goslim{$data[0]}{"Prelat"} = $data[3];
			$goslim{$data[0]}{"Pgoterm"} = $data[4];
			$goslim{$data[0]}{"Pgoid"} = $data[5];
			$goslim{$data[0]}{"Pgoslim"} = $data[8];
			$goslim{$data[0]}{"Pdate"} = $data[14];
			my $P = join "|", ($data[7],$data[14], $data[8], $data[5], $data[3], $data[4]); #saves the latest annotation for process as a joined array
			$goslim{$data[0]}{"Parray"} = $P;
			
		}
	}
}

open	(OUTFILE, ">>XXgoslim_intermediate.txt") || die "Can't print outfile $!";
for my $key (sort keys %goslim)	{
#		for my $x (sort keys %{$goslim{$key}})	{
#		push ( @yoyo, $goslim{$key}{$x});
#		}
#	my $q = join "\t", @yoyo;
	
	print OUTFILE "$key\|$goslim{$key}{Carray}\|$goslim{$key}{Farray}\|$goslim{$key}{Parray}\n";	#prints out AGI name C array F array and P array
#	$q = 0;
#	@yoyo = ();
}

close OUTFILE;

 open (INFILE, "XXgoslim_intermediate.txt") || die "Can't open goslim_intermediate.txt $!\n";	#opens the file that was just printed out back up and looks for missing data, which it will replace with placeholders so that all the annotation is in the same place in the array
	while (<INFILE>)
{
	my @y = split (/\|/);
	if 	(($y[1] ne "C") && ($y[2] eq "F") && ($y[8] eq "P"))	{
			 $goslim{$y[0]}{Carray} = join "|", ("C","Cdate", "Cgoslim_term", "Cgo_id", "Crelation", "Cgoterm");
	}	
	if 	(($y[1] ne "C") && ($y[2] ne "F")  && ($y[3] eq "P") )	{
			 $goslim{$y[0]}{Carray} = join "|", ("C","Cdate", "Cgoslim_term", "Cgo_id", "Crelation", "Cgoterm");
			 $goslim{$y[0]}{Farray} = join "|", ("F","Fdate", "Fgoslim_term", "Fgo_id", "Frelation", "Fgoterm"); 
	}	
	if	( ($y[1] ne "C")&&($y[2] eq "F")&& ($y[8] ne "P")) 	{
 			 $goslim{$y[0]}{Carray} = join "|", ("C","Cdate", "Cgoslim_term", "Cgo_id", "Crelation", "Cgoterm");
			 $goslim{$y[0]}{Parray} = join "|", ("P","Pdate", "Pgoslim_term", "Pgo_id", "Prelation", "Pgoterm");
	}
	if 	( ($y[1] eq "C") && ($y[7] ne "F") && ($y[8] ne "P"))	{
			 $goslim{$y[0]}{Farray} = join "|", ("F","Fdate", "Fgoslim_term", "Fgo_id", "Frelation", "Fgoterm");
			 $goslim{$y[0]}{Parray} = join "|", ("P","Pdate", "Pgoslim_term", "Pgo_id", "Prelation", "Pgoterm");
	}
	if 	( ($y[1] eq "C") && ($y[7] ne "F") && ($y[8] eq "P"))	{
			 $goslim{$y[0]}{Farray} = join "|", ("F","Fdate", "Fgoslim_term", "Fgo_id", "Frelation", "Fgoterm");
	}	
	if 	(($y[1] eq "C") && ($y[7] eq "F") && ($y[13] ne "P"))	{
			 $goslim{$y[0]}{Parray} = join "|", ("P","Pdate", "Pgoslim_term", "Pgo_id", "Prelation", "Pgoterm");						
	}
}

open (OUTFILE, ">>XX$ARGV[1]\_cfp_goslim_out.txt") || die "Can't open goslim_cpf_complete out $!\n";	#

for my $z (sort keys %goslim)  {
	print OUTFILE "$z\|$goslim{$z}{Carray}\|$goslim{$z}{Farray}\|$goslim{$z}{Parray}\n";
}
	
