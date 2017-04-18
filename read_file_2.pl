#!c:\perl\bin

$line_no=0;
open(fp1, "gene_info.txt");
while(<fp1>)
{
	$line_no = $line_no + 1;
	print "Reading line $line_no\n";
	
	@info = $_;
	# $line_elements[0] contains the first data item
	# $line_elements[1] contains the second data item
	# $line_elements[2] contains the third data item
	print "@info";
}
print "Done!\n";
close(fp1);

open(fp2, "gene_data.txt");
while(<fp2>)
{
	@data = $_;
	print "@data";
}
print "Done and done!\n";
close(fp2);
