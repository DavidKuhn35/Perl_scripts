#!/usr/bin/perl -w
my $blast_report = ();
chdir "/home/dkuhn/library_blast" or die "Cannot chdir /home/dkuhn/library_blast";
system "ls -a *.blastn >>blast_list.txt";
open (BLAST_REPORT, "blast_list.txt")|| die "Cant open file\n";  #open filehandle for library list
while (<BLAST_REPORT>)
{
	chomp ($blast_report = $_) ;  #read in library names and do blast search
system "parse_blastReport_TopHit.pl $blast_report";
}
