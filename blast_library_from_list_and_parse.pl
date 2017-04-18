#!/usr/bin/perl -w
my $library = ();  #initializing library variable as undef
open(LIBRARY_LIST, "/home/dkuhn/FrostyPod_2011_working/All_ESTtik_fasta/lib_list_test" ) || die "Cant open file\n";  #open filehandle for library list
while (<LIBRARY_LIST>)
{
	chomp ($library = $_) ;  #read in library names and do blast search
	
system "blastall -p blastn -d /home/dkuhn/FrostyPod_2011_working/All_ESTtik_fasta/CGDv09 -i /home/dkuhn/FrostyPod_2011_working/All_ESTtik_fasta/$library -e 1e-15 -o /home/dkuhn/library_blast/$library\_vs_CGDv09\.blastn -v 5 -b 5"; #parameters for blast search are set only library name and output file will vary
}
my $blast_report = ();
chdir "/home/dkuhn/library_blast" or die "Cannot chdir /home/dkuhn/library_blast";
system "ls -a *.blastn >>blast_list.txt";
open (BLAST_REPORT, "blast_list.txt")|| die "Cant open file\n";  #open filehandle for library list
while (<BLAST_REPORT>)
{
	chomp ($blast_report = $_) ;  #read in library names and do blast search
system "/home/dkuhn/perl_training/parse_blastReport_TopHit.pl $blast_report";
}
