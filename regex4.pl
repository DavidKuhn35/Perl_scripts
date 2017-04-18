#!C:\perl\bin

$dna = "nnnTAGCAgcgaatagcaTGcgatagcGCAtagcagatagcatTTTAATTATA";
print "Old dna is $old_dna\n";

$_ = $dna;
tr/[a-z]/[A-Z]/;
print "DNA in uppercase letters is $_\n";

s/N/T/g;
print "After substituting N with T, DNA is $_\n";

tr/ATGC/UACG/;
print "After transcription, RNA is $_\n";