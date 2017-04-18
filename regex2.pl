#!C:\perl\bin

$dna = "GCAGTACAGCAGAATACAGCTGACGACGATAGAGACAGCGCAGCGATAGCAGAT";
$rest_enzyme = $ARGV[0];

@dna_fragments = split(/$rest_enzyme/i, $dna);

print "The number of hits is ".(scalar(@dna_fragments)-1)."\n";