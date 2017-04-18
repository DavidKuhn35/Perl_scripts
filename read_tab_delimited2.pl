#!c:\perl\bin

$line_no=0;
open(fp1, "tab_delimited_data.txt");
open(fp2, ">tab_delimited_data2.txt");
while(<fp1>)
{
	$line_no = $line_no + 1;
	print "Reading line $line_no\n";
	$current_line = $_;
	chomp $current_line;
	@line_elements = split(/\t/, $current_line);
	# $line_elements[0] contains the first data item
	# $line_elements[1] contains the second data item
	# $line_elements[2] contains the third data item
	$bottom_element = pop (@line_elements);
	unshift(@line_elements, $bottom_element);
	foreach $name (@line_elements)
	{
		print "$name\t";
	}
	print "\n";
	print fp2 "$line_no\t$line_elements[2]\t$line_elements[3]\n";
}
print "Done!\n";
close(fp1);
close(fp2);