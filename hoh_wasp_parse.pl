#!/usr/bin/perl -w
#This program makes a hash of hashes from asap.py allele specific primer data read in from a file.
my %hoh = ();

open (INFILE, "[ARGV[0]")|| die "Can't open infile";  #open the file to be read.
	while (<INFILE>) #read in from file, line by line
{
	chomp;	#remove end of line character
	if (/^<snp id="(\w\d)">/) {	#match SNP name and store in memory variable
		$SNPname = $1;
		$hoh{$SNPname}{SNPname} = $1;  #In hash of hashes, first key is SNP name.
		} elsif (/^<template>/) {
			next;
		 }elsif (/^<sequence id="id1" length="(\d)">/)  {
			$hoh{$SNPname}{seqlength} = $1;
		  
		
	}
}
	open (LIST, "/home/dkuhn/CoreSNPs/$ARGV[1]")|| die "Can't open infile";  #open the file to be read.
	while (<LIST>)
{
	chomp;
	$cultivar1 = $_; 
#	print "The value of cultivar $ARGV[0] and SNP $SNP is $hoh{$ARGV[0]}{$SNP}\n";
	for $key (sort keys %{ $hoh{$ARGV[0]} }) {
		if ($hoh{$ARGV[0]}{$key} eq $hoh{$cultivar1}{$key}) {
		$same_genotype{$key} = $hoh{$ARGV[0]}{$key};
#	print "The value of same_genotype is $same_genotype{$key}\n";
		$same_SNP{$key} += 1; #make a counter to see how many times SNP is same genotype
		print "SNP $key is the same SNP $same_SNP{$key} times\n";
		$same_SNP_cultivar{$key} .= "$cultivar1\t"; #says which cultivars have the same genotype at the SNP
		print "SNP $key is the same SNP in these cultivars $same_SNP_cultivar{$key}\n";
		}else {
		$different_genotype{$key} = $hoh{$ARGV[0]}{$key};
		$different_SNP{$key} += 1; #make a counter to see how many times SNP is a different genotype
		print "SNP $key is a different SNP $different_SNP{$key} times\n";
		$different_SNP_cultivar{$key} .= "$cultivar1\t";  #says which cultivars have a different genotype at the SNP
#		print "The SNP where the genotypes differ is $key.\n";
		print "SNP $key is a different SNP in these cultivars $different_SNP_cultivar{$key}\n";
	}
}
	$same_count = keys %same_genotype;
	$diff_count = keys %different_genotype; 
open (OUTFILE, ">>/home/dkuhn/CoreSNPs/$ARGV[0]_vs_$ARGV[1]")|| die "Can't open outfile";
	printf OUTFILE "$ARGV[0] vs $cultivar1\tsame is ($same_count / 5388)\t %5.2f \t different is ($diff_count / 5388) %5.2f \n", $same_count/5388, $diff_count/5388;
	%same_genotype = ();
	%different_genotype = (); 
}
		
open (OUTFILE, ">>/home/dkuhn/CoreSNPs/$ARGV[0]_vs$ARGV[1]_SNP_out.txt")|| die "Can't open outfile";
	foreach $key (sort keys %same_SNP)
{
	print OUTFILE "$key\t$same_SNP{$key}\n";
	print OUTFILE "$key\t$same_SNP_cultivar{$key}\n";
	print OUTFILE "$key\t$different_SNP{$key}\n";
	print OUTFILE "$key\t$different_SNP_cultivar{$key}\n";
}	 
