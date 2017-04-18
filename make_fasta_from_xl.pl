#!/usr/bin/perl -w
#This program makes a fasta file from an xl format.

open (INFILE, "$ARGV[0]")|| die "Can't open infile";  #open the file to be read.
	while (<INFILE>) #read in from file, line by line
{
	chomp;	#remove end of line character
	if (m/^g/) {	#skip any header information, start at genotype data
		($name, $dum1, $dum2, $dum3, $dum4, $sequence) = split (/\t/); #assign variables from genotype data
	open (OUTFILE, ">>$ARGV[1]")|| die "Can't open outfile";
	print "\>$name\n$sequence\n";
	print OUTFILE "\>$name \n$sequence\n";
}
		
}	 
