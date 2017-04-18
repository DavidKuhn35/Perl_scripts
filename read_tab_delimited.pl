#!c:\perl\bin

$line_no=0;
open(fp1, "tab_delimited_data.txt");
while(<fp1>)
{
	$line_no = $line_no + 1;
	print "Reading line $line_no\n";
	
	@line_elements = split(/\t/, $_);
	# $line_elements[0] contains the first data item
	# $line_elements[1] contains the second data item
	# $line_elements[2] contains the third data item
	print "@line_elements[3]\n";
}
print "Done!\n";
close(fp1)