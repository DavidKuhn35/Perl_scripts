#!/usr/bin/perl -w
my $library = ();  #initializing library variable as undef
chdir "/home/dkuhn/FrostyPod_2011_working/sub_fasta";
system "ls -a *.fasta >> fasta_list";
open(LIBRARY_LIST, "/home/dkuhn/FrostyPod_2011_working/sub_fasta/fasta_list" ) || die "Cant open file\n";  #open filehandle for library list
while (<LIBRARY_LIST>)
{
	chomp ($library = $_) ;  #read in library names and do blast search
	
system "blastall -p blastn -d /home/dkuhn/FrostyPod_2011_working/All_ESTtik_fasta/CGDv09 -i /home/dkuhn/FrostyPod_2011_working/sub_fasta/$library -e 1e-15 -o /home/dkuhn/FrostyPod_2011_working/sub_fasta/$library\_vs_CGDv09\.blastn -v 5 -b 5"; #parameters for blast search are set only library name and output file will vary
}
my $blast_report = ();
chdir "/home/dkuhn/sub_fasta" or die "Cannot chdir /home/dkuhn/sub_fasta";
system "ls -a *.blastn >>blast_list.txt";
open (BLAST_REPORT, "blast_list.txt")|| die "Cant open file\n";  #open filehandle for library list
while (<BLAST_REPORT>)
{
	chomp ($blast_report = $_) ;  #read in blast report names and parse blast search
system "/home/dkuhn/perl_training/parse_blastReport_TopHit.pl $blast_report";
}	#output files go to current directory

	system "cp *_TopHit*parsed.txt /home/dkuhn/FrostyPod_2011_working/sub_fasta/";
	chdir "/home/dkuhn/FrostyPod_2011_working/sub_fasta";
	system "ls -a *_TopHit*parsed.txt >> parse_list.txt";
	my $parse_list = ();
	
	my $library_name = ();
open (PARSE_LIST, "parse_list.txt")  || die "Cant open parse list file\n";
	while (<PARSE_LIST>)
{
	chomp ($parse_list = $_) ;

# Loading data into hash arrays and initializing the hashes
%hash1_query = ();
%hash2_CGD_counter = ();
		open(PARSE_FILE, "$parse_list") || die "Cant open $parse_list file\n";  #  Open the parsed BLAST file to read in the queries and hits
		
		if (m/^(\w+)\./) {
			$library_name = $1;
}
		while(<PARSE_FILE>)
{
			if  (/^[^qQ]/) #Gets rid of the header which starts with Query
{
			($query, $dummy1, $dummy2, $CGDnumber) = split(/\t/);  #The query is in column 1 and the CGD hit is in column 4, the other data is lost
			$hash1_query{$query} = $CGDnumber;  #Load the first hash with a query key and a CGD value
			$hash2_CGD_counter{$CGDnumber} += 1; #Load the second hash with a CGD key and a count number for a value
#	print "This is the CGD number $CGDnumber \n";
#	print "This is the value of the CGD counter $hash2_CGD_counter{$CGDnumber} for the CGD number $CGDnumber \n";
}
#close (PARSE_FILE);	
}
#close (PARSE_LIST);  #Just being tidy

	while ( ($key, $value) = each %hash1_query) #looks all the way through the hash
{
	open (OUTFILE, ">>$library_name\_queryhash.txt");
	print OUTFILE "$key \t $value \n";  #prints out each key value pair, can be sent to a file
}
	open (OUTFILE, ">>$library_name\_uniqueCGD_count");
	while ( ($key1, $value1) = each %hash2_CGD_counter) #prints out each CGD number and the number of times it was seen.  Can be sent to a file.
{
	print OUTFILE "$key1 \t $value1 \n";
}
# Using the data in hash arrays to create a report
	open (OUTFILE, ">>$library_name\_unique_summary");
	my $count = keys %hash2_CGD_counter; #counts the number of keys in each hash.
	my $total = keys %hash1_query;	print OUTFILE "There were $total number of queries and $count number of unique hits in the $library_name library. \n";
#  Print out to a file the list of unique CGD numbers from the library.
#	my @CGD = keys %hash2_CGD_counter;
#	print "@CGD \n"; #should print out list of keys with new lines but doesn't!
#  Check to see that inverse hash is the same as CGD hash
#	%inverse_query = reverse %hash1_query;
#	my $inv_count = keys %inverse_query;
#	print "The unique hits from the CGD hash is $count.  The unique hits from the inverse hash is $inv_count\n"; 
#open (OUTFILE, ">>$library_name\_unique.out") || die "Can't open exist.out\n"; #Name the outfile and set that data will be appended
#	while ( ($key3, $value3) = each %inverse_query)
#{	
#	if (exists $hash2_CGD_counter{$key3}) 
#{
#	print OUTFILE "$key3 exists in both CGD counter and inverse hash \n";
#}
#}
}
