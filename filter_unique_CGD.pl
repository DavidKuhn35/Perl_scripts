#!/usr/bin/perl -w
#Reads in unique CGDs from libraries and produces a list of all unique CGDs for all libraries analyzed.  In addition, counts number of libraries each CGD occurs in.  For single occurrence CGDs, will give identity of library from which it came.
my $library_name = ();  #initializing library variable as undef
my $library_file_name = ();
my $library_count = ();
my @library = ();
%unique_count = ();
%unique_library = ();
chdir "/home/dkuhn/unique_library/";
system "ls -a *_uniqueCGD_count >> xunique_library_list.txt";
open(UNIQUE_LIBRARY_LIST, "/home/dkuhn/unique_library/xunique_library_list.txt" ) || die "Cant open file\n";  #open filehandle for library list
while (<UNIQUE_LIBRARY_LIST>)
{
	chomp ($library_file_name = $_) ;
	if (m/^(\w+)\_uniqueCGD/) {
		$library_name = $1;
		}
		open (UNIQUE_LIBRARY_FILE, "$library_file_name"); 
			while (<UNIQUE_LIBRARY_FILE>) 	
{
			chomp;
			($CGD ) = split(/\t/);
			$unique_count{$CGD} += 1;
			$unique_library{$CGD} .= "$library_name\t";
}
		$library_count += 1;
#		push @library, $library_name;	
}
#	while (($key, $value) = each %unique_count) #looks all the way through the hash
	foreach $key (sort keys %unique_count)
{
	open (OUTFILE, ">> xtotal_unique_count.txt")|| die "Cant open xtotal_unique_count.txt file\n";  #open filehandle for library list;
	print OUTFILE "$key\t$unique_count{$key} \n";  #prints out each key value pair, can be sent to a file
}
	foreach $key (sort keys %unique_library)
{
	open (OUTFILE, ">> xtotal_unique_library_count.txt")|| die "Cant open xtotal_unique_library_count.txt file\n";  #open filehandle for library list;
	print OUTFILE "$key\t$unique_library{$key} \n";  #prints out each key value pair, can be sent to a file
}		
#		if ($unique_count{$key} == 1)  {
#		open (OUTFILE, ">> xunique_CGD_library_id.txt")|| die "Cant open xunique_CGD_library_id.txt file\n";  #open filehandle for library list;
#		print OUTFILE "$key was found only in $unique_library{$key} library \n"; 
#}
#		if ($unique_count{$key} >= ($library_count)) {
#		open (OUTFILE, ">> xcommon_CGD.txt");
#		print OUTFILE "$key was found in more than " . ($library_count) . " libraries searched. \n";
#}

#my $CGD_single_library = ();
#my $CGD_all_library = ();
#my $CGD_unique_count = keys %unique_count;
#chomp ($CGD_single_library = `grep -c "was" /home/dkuhn/unique_library/xunique_CGD_library_id.txt`);
#chomp  ($CGD_all_library = `grep -c "was" /home/dkuhn/unique_library/xcommon_CGD.txt`);
open (OUTFILE, ">> xunique_summary.txt")|| die "Cant open xunique_summary.txt file\n";  #open filehandle for library list;
	print OUTFILE "There were $library_count libraries searched for unique CGD hits.\n  There were $CGD_unique_count CGDs found in the $library_count libraries.\n";

#count how many unique sequences were in each library
%unique_CGD = ();
%unique_library_count = ();
open (INFILE, "/home/dkuhn/unique_library/xunique_CGD_library_id.txt")|| die "Cant open infile\n";  #open filehandle for library list;
	while (<INFILE>)
{
	if (m/^(\w+)\s+was\s+found\s+only\s+in\s+(\w+)\s+library/)	{
		$unique_CGD{$1} = $2;
		$unique_library_count{$2} += 1;
#	print "$1 \t $2\n";
		}
}
open (OUTFILE, ">> /home/dkuhn/unique_library/xunique_library_uniqueCGDs.txt")|| die "Cant open outfile\n";  
while (($key, $value) = each %unique_library_count)
{
	print	OUTFILE "$key\t\t$value \n";
}







