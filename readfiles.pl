#!:\perl\bin\perl

#This program reads in variable names, and creates a report to file

open(fp1, "gene_info.txt");
print "<fp1>";
$i=0;
while(<fp1>)
{
	chomp;
	@info[$i]=$_;
	$i=$i+1;
	print "@info[$i]\n";
}
end
close(fp1)


#open(fp2, ">example6_new.txt");
#print fp2 "The variables are:\n";
#for $i=1:length(@var_names)
#	print fp2 "@var_names($i)\n";
#end
#close(fp2)