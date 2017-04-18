#!c:\perl\bin\perl


$dna = "";
open(fp1, "sequence2.txt") || die "Cant open file!\n";
while(<fp1>)
{
	chomp $_;
	$dna = $dna.$_;
}
close(fp1);
print "The DNA sequence is:\n$dna\n";
print "The length of DNA is sequence is ".length($dna);

$count_a = 0;
$count_c = 0;
$count_t = 0;
$count_g = 0;

for($i=0; $i<length($dna); $i++)
{
	if(substr($dna,$i,1) eq 'A')
	{ $count_a = $count_a + 1; }
	
	if(substr($dna,$i,1) eq 'C')
	{ $count_c = $count_c + 1; }
	
	if(substr($dna,$i,1) eq 'T')
	{ $count_t = $count_t + 1; }
	
	if(substr($dna,$i,1) eq 'G')
	{ $count_g = $count_g + 1; }
}

print "\n";
print "Number of As are: $count_a\n";
print "Number of Cs are: $count_c\n";
print "Number of Gs are: $count_g\n";
print "Number of Ts are: $count_t\n";