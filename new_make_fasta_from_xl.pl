#!/usr/bin/perl -w
#This program makes a fasta file from an xl format.
#January 15,2017
#Modify to create fasta file from modified excel file saved as csv
#Fills in spaces with _
#Converts bracket format to IUPAC
#Reads in format >header:sequence and outputs a fasta file 
#program ARGV[0] input csv ARGV[1] output fasta

open (INFILE, "$ARGV[0]")|| die "Can't open infile";  #open the file to be read.
	while (<INFILE>) #read in from file, line by line
{
	chomp;	#remove end of line character
	if (m/^>/) {	#skip any header information, start at genotype data
		my $line= s/\s/_/g;
#		my $line=~ s/\s+/_/g ; #remove any spaces and replace with underscores so that info will be in header as id, not id and description
		(my$header, my$seq) = split (/\:/), $line; #assign variables from single line header:sequence format of csv
#  Convert bracket format in sequence to IUPAC		
	if (($seq =~ m!\[A/G\]!) || ($seq =~ m!\[G/A\]!)) {
		$seq =~	s!\[\w/\w\]!R! ;
#	print "The converted sequence is $seq.\n";

#	print OUTFILE ">$fasta_name\n$seq\n";	
	}elsif (($seq =~ m!\[C/G\]!)|| ($seq =~ m!\[G/C\]!)) {
		$seq =~ s!\[\w/\w\]!S!;
#	open (OUTFILE, ">>IUPAC_outfile$ARGV[0].fasta") || die "Can't open outfile\n";
#	print OUTFILE ">$fasta_name\n$seq\n";
#	print "The converted sequence is $seq.\n";
	}elsif (($seq =~ m!\[T/G\]!) || ($seq =~ m!\[G/T\]!)) {
		$seq =~ s!\[\w/\w\]!K!;
#	open (OUTFILE, ">>IUPAC_outfile$ARGV[0].fasta") || die "Can't open outfile\n";
#	print OUTFILE ">$fasta_name\n$seq\n";
#	print "The converted sequence is $seq.\n";
	}elsif (($seq =~ m!\[A/C\]!) || ($seq =~ m!\[C/A\]!)) {
		$seq =~ s!\[\w/\w\]!M!;
#	open (OUTFILE, ">>IUPAC_outfile$ARGV[0].fasta") || die "Can't open outfile\n";
#	print OUTFILE ">$fasta_name\n$seq\n";
#	print "The converted sequence is $seq.\n";
	}elsif (($seq =~ m!\[T/C\]!)|| ($seq =~ m!\[C/T\]!))  {
		$seq =~ s!\[\w/\w\]!Y!;
#	open (OUTFILE, ">>IUPAC_outfile$ARGV[0].fasta") || die "Can't open outfile\n";
#	print OUTFILE ">$fasta_name\n$seq\n";
#	print "The converted sequence is $seq.\n";
	}elsif (($seq =~ m!\[A/T\]!)|| ($seq =~ m!\[T/A\]!)) {
		$seq =~ s!\[\w/\w\]!W!;
#	print "The converted sequence is $seq.\n";
#	open (OUTFILE, ">>IUPAC_outfile_$ARGV[0]") || die "Can't open outfile\n";
#	print OUTFILE ">$fasta_name\n$seq\n";
	}	
		
		
		
	open (OUTFILE, ">>$ARGV[1]")|| die "Can't open outfile";
#	print "\>$name\n$seq\n";
	print OUTFILE "$header\n$seq\n";
}
		
}	 
