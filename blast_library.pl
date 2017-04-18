#!/usr/bin/perl -w

my $library = $ARGV[0];
system "blastall -p blastn -d CGDv09 -i $library -e 1e-15 -o $library\_vs_CGDv09\.blastn -v 5 -b 5";
