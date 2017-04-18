#!c:\perl\bin\perl


#######################################################
# Loading sequence data
#######################################################
print "Loading sequence data from file into a Perl variable...\n";
$dna = "";
open(fp1, "sequence2.txt") || die "Cant open file!\n";
while(<fp1>)
{
	chomp $_;
	$dna = $dna.$_;
}
close(fp1);


#######################################################
# Cutting DNA
#######################################################
print "Cutting with enzyme $enzyme...\n";
$enzyme = "GATT";
@dna_fragments = split(/$enzyme/g, $dna);


#######################################################
# Displaying
#######################################################
print "Printing the sequence of resulting DNA fragments...\n";
$i=0;
foreach $cut_dna (@dna_fragments)
{
	$i++;
	print "$i\t$cut_dna\n"
}

#######################################################
# Identifying fragments
#######################################################
print "\n\n\nFinding fragments that have >1 CG, and are 100-200 bp...\n";
foreach $cut_dna (@dna_fragments)
{
	if((/CG/,$cut_dna) && length($cut_dna)>100 && length($cut_dna)<200)
	{
		print "$cut_dna\n";
	}
}