#!c:\perl\bin
$line=0;
$line_no=0;
open(fp1, "gene_info.txt");
while(<fp1>)
{
	
	print "Reading line $line_no\n";
	chomp;
	@info[$line_no] = $_;
	# $line_elements[0] contains the first data item
	# $line_elements[1] contains the second data item
	# $line_elements[2] contains the third data item
	$line_no = $line_no + 1;
}
print "@info\n";
print "Done!\n";
close(fp1);

open(fp,"gene_data.txt");
while(<fp>)
{
	
	print "Reading line $line\n";
	chomp;
	@data[$line] = $_;
	# $line_elements[0] contains the first data item
	# $line_elements[1] contains the second data item
	# $line_elements[2] contains the third data item
	$line = $line + 1;
}
print "@data\n";
print "Done!\n";
close(fp);

$i=0;
foreach $name (@info)
{
	$merge = join("\t", $i+1,$info[$i],$data[$i]);
	@merge[$i] = $merge;
	$i= $i+1;
	print "$merge\n";
}
print "@merge\n";