#!/usr/bin/perl -w
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
