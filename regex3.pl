#!C:\perl\bin

$dna = "TAGCAgcgaatagcaTGcgatagcGCAtagcagatagcatTTTAATTATA";
print "Dna is $dna\n";

$dna =~ tr/[a-z]/[A-Z]/, $dna;
print "New dna is $dna\n";
