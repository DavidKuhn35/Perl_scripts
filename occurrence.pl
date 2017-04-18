#!/usr/bin/perl -w
my %hash1 = ();
my $line_count = ();
		open (INFILE, "/home/dkuhn/unique_library/xtotal_unique_count.txt"); 
			while (<INFILE>) 	
{
			chomp;
			($a,$b) = split(/\t/);
			$hash1{$b} += 1;
			$line_count += 1;
}
	foreach $key (sort keys %hash1) 
{
	open (OUTFILE, ">> /home/dkuhn/unique_library/xoccurrence.txt");
	print OUTFILE "$key \t $hash1{$key}\n";  #prints out each key value pair, can be sent to a file
}
