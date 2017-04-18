open (fp, "gene_info.txt");
while(<fp>)
{
	chomp($_);
	@info = $_;
}	
close(fp);

print "@info\n";
print "this is the info part\n";

#print "@info[2]\n";

open (fp1, "gene_data.txt");
while (<fp1>)
{
	$current_line = $_;
	chomp($current_line);	
	@data = <fp1>;
}
close(fp1);
print "@data\n";
$i=0;
foreach $name (@info)
{
	$merge = join("\t", $i+1,$info[$i],$data[$i]);
	$i= $i+1;
	print "$merge\n";
}
