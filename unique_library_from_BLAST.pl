#!/usr/bin/perl -w

# Loading data into hash arrays and initializing the hashes
%hash1_query = ();
%hash2_CGD_counter = ();
open(fp, $ARGV[0]) || die "Cant open file\n";  #  Open the parsed BLAST file to read in the queries and hits
my $library_name = $ARGV[1];
while(<fp>)
{
	if  (/^[^qQ]/) #Gets rid of the header which starts with Query
{
	($query, $dummy1, $dummy2, $CGDnumber) = split(/\t/);  #The query is in column 1 and the CGD hit is in column 4, the other data is lost
	$hash1_query{$query} = $CGDnumber;  #Load the first hash with a query key and a CGD value
	$hash2_CGD_counter{$CGDnumber} += 1; #Load the second hash with a CGD key and a count number for a value
#	print "This is the CGD number $CGDnumber \n";
#	print "This is the value of the CGD counter $hash2_CGD_counter{$CGDnumber} for the CGD number $CGDnumber \n";
}	
}
close(fp);  #Just being tidy

	while ( ($key, $value) = each %hash1_query) #looks all the way through the hash
{
	open (OUTFILE, ">>$library_name\_queryhash.txt";
	print OUTFILE "$key \t $value \n";  #prints out each key value pair, can be sent to a file
}
	open (OUTFILE, ">>$library_name\_uniqueCGD_count");
	while ( ($key1, $value1) = each %hash2_CGD_counter) #prints out each CGD number and the number of times it was seen.  Can be sent to a file.
{
	print OUTFILE "$key1 \t $value1 \n";
}
# Using the data in hash arrays to create a report
	open (OUTFILE ">> $library_name\_unique_summary");
	my $count = keys %hash2_CGD_counter; #counts the number of keys in each hash.
	my $total = keys %hash1_query;	print OUTFILE "There were $total number of queries and $count number of unique hits in the $library_name library. \n";
#  Print out to a file the list of unique CGD numbers from the library.
	my @CGD = keys %hash2_CGD_counter;
	print "@CGD \n"; #should print out list of keys with new lines but doesn't!
#  Check to see that inverse hash is the same as CGD hash
	%inverse_query = reverse %hash1_query;
	my $inv_count = keys %inverse_query;
#	print "The unique hits from the CGD hash is $count.  The unique hits from the inverse hash is $inv_count\n"; 
open (OUTFILE, ">>$library_name\_unique.out") || die "Can't open exist.out\n"; #Name the outfile and set that data will be appended
	while ( ($key3, $value3) = each %inverse_query)
{	
	if (exists $hash2_CGD_counter{$key3}) #Compare keys from one hash to another and print out the keys that exist in both hashes, should be able to imbed more loops to compare more hashes with elsif
	{
	print OUTFILE "$key3 exists in both CGD counter and inverse hash \n";# may be able to print out the names of all the libraries that have this CGD in it at the end of the imbedded loops.  
	}
}
