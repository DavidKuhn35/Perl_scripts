#!/usr/bin/perl -w
#This program makes a hash of hashes from asap.py allele specific primer data read in from a file.
my %hoh = ();

open (INFILE, "[ARGV[0]")|| die "Can't open infile";  #open the file to be read.
	while (<INFILE>) #read in from file, line by line
{
	chomp;	#remove end of line character
	if (/^<snp id="(\w+)">/) {	#match SNP name and store in memory variable
		$SNPname = $1;
		$hoh{$SNPname}{SNPname} = $1;  #In hash of hashes, first key is SNP name.
		} elsif (/^<template>/) {
			next;
		 }elsif (/^<sequence id="id1" length="(\d+)">/)  {
			print $_ ;
			$hoh{$SNPname}{seqlength} = $1;
		} elsif  (%^(\w+[\w+/\w+]\w+)</sequence>$%)  {
			$hoh{$SNPname}{sequence} = $1;
		}elsif (/^<properties | ^item value/)	{
			next;
		}elsif (/^oligo id="\w"\s+\w+="(.+)"\s+position="(.*)"\s+length="(.*)"\s+seq="(w+)">$/)  {
			$oligo_id = $1;
			$hoh{$SNPname}{$oligo_id\_primertype} = $2;
			$hoh{$SNPname}{$oligo_id\_position} = $3;
			$hoh{$SNPname}{$oligo_id\_primerlength} = $4;
			$hoh{$SNPname}{$oligo_id\_primerseq} = $5;
		}elsif (/^Tm=(\.+),\s+GC%=(\.+),\s+Self Any=(\.+),\s+Self End=(\.+)$/)	{
			$hoh{$SNPname}{$oligo_id\_Tm} =  $1;
			$hoh{$SNPname}{$oligo_id\_GC} =	 $2;
			$hoh{$SNPname}{$oligo_id\_Self_Any} = $3;
			$hoh{$SNPname}{$oligo_id\_Self_End} = $4;
