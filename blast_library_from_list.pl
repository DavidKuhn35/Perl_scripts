#!/usr/bin/perl -w
my $library = ();
open(LIBRARY_LIST, "/home/dkuhn/FrostyPod_2011_working/All_ESTtik_fasta/lib_list_test" ) || die "Cant open file\n";
while (<LIBRARY_LIST>)
{
	chomp ($library = $_) ;
	
system "blastall -p blastn -d /home/dkuhn/FrostyPod_2011_working/All_ESTtik_fasta/CGDv09 -i /home/dkuhn/FrostyPod_2011_working/All_ESTtik_fasta/$library -e 1e-15 -o /home/dkuhn/library_blast/$library\_vs_CGDv09\.blastn -v 5 -b 5";
}
