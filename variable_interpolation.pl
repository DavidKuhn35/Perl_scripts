#!C:\perl\bin

$a = 2;
$b = 3;
print "The first number is $a, second number is $b\n";
print "The first number is \$a = $a, second number is \$b = $b\n";


@dna_fragments = ("atagaggactaa", "aaatttggc", "atatagagcagcagcagcatagcagatagcagac");
$i = 0;
foreach (@dna_fragments)	
{ 
	$i = $i + 1;
	print "Restriction fragment number $i is $_\n";
}