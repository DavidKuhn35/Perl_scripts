#!/usr/bin/perl -w
#This program makes a fasta file from an xl format.
#January 15,2017
#Modify to create fasta file from modified excel file saved as csv
#Fills in spaces with _
#Converts bracket format to IUPAC
#Reads in format >header:sequence and outputs a fasta file 
#program ARGV[0]  ARGV[1] output fasta
my %map = ();
my %header = ();
open (MARKER, "$ARGV[0]")|| die "Can't open MARKER";  #open the file to be read.
	while (<MARKER>) #read in from file, line by line
{
	chomp;	#remove end of line character
	if (m/^Marker\d+/) {	#skip any header information, start at genotype data
		(my$marker, my$lg, my $cm) = split (/\t/), $_; #assign variables from marker, linkage group and cM position of Luo map
    	$map{$marker} = "$marker\_$lg\_$cm";	
}
}
close MARKER;		
###########################################################################3

open (HEADER, "$ARGV[1]")|| die "Can't open HEADER";  #open the file to be read.
	while (<HEADER>) #read in from file, line by line
{
	chomp;	#remove end of line character
	if (m/^Marker\d+/) {	#skip any header information, start at genotype data
		(my$header, my$seq) = split (/\t/), $_; #assign variables from header and sequence of Luo map markers	
        $header{$header} = $seq;
        }
}

###############################################################################
    

	open (OUTFILE, ">>$ARGV[2]")|| die "Can't open outfile";
    foreach my $key (sort keys %map) {
	    
	print OUTFILE "\>$map{$key}\n$header{$key}\n";
}
		

