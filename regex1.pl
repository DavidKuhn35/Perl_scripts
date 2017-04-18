#!C:\perl\bin

$dna = "GCAGTACAGCAGAATATGACGACGATAGAGAGCAGCGATAGCAGAT";
$rest_enzyme = $ARGV[0];

if($dna =~ /$rest_enzyme/i)
{
	print "Enzyme found!\n";
}
else
{
	print "Enzyme not found!\n";
}