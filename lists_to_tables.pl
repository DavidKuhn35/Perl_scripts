#!/usr/bin/perl -w
my $dummy1 = ();
my $SNP = ();
my $SNPpos1 = ();
#my $newSNPpos1 = ();
my $cultivar = ();
my $cultivar1 = ();
my $allele1 = ();
my $allele2 = ();
my %SNP_lookup = ();
my %hoh = ();
my %hoh2 = ();
my %SNPlist = ();
my %cultivarlist = ();
my $header = "Tree name\t";
#my $genotype = ();
#Read in the SNP conversion file and make a hash
open (SNP, "/home/dkuhn/6kSNPchip/6kSNPs_sortable_Aug11.txt") || die "Can't open SNP infile";
	while (<SNP>)
{
	chomp;
	if (m/^Tcm/)  {
		($dummy1,$SNPpos1, $SNP) = split (/\t/);
		$SNP_lookup{$SNP} = $SNPpos1;
	}
}
#	open (OUTFILE, ">>/home/dkuhn/6kSNPchip/SNPlookupout.txt") || die "Can't open SNP lookup outfile";
#	foreach $key (sort keys %SNP_lookup) 	{
#		print OUTFILE "$key\t$SNP_lookup{$key}\n";
#	}

#Read in the SNPchip data and convert SNP names to SNP position names using hash
open (INFILE, "/home/dkuhn/6kSNPchip/20tree_testset.txt")|| die "Can't open infile";  #open the file to be read.
	while (<INFILE>) #read in from file, line by line
{
	chomp;	#remove end of line character
	if (m/^contig/) {	#skip any header information, start at genotype data
		($SNP, $cultivar, $allele1,$allele2) = split (/\t/); #assign variables from genotype data
	$hoh{$cultivar}{$SNP_lookup{$SNP}} = "$allele1$allele2";
#	print "For $cultivar, the $SNP genotype is $hoh{$cultivar}{$SNP}\n";
	}
}
open (CULTIVAR_LIST, "/home/dkuhn/CoreSNPs/diversity_DNA_list.txt") || die "Can't open Cultivarlist_test";
	while (<CULTIVAR_LIST>)
{		
	($cultivarlist) = split /\t/ ;
	$cultivarlist{$cultivarlist} = $cultivarlist;
}
		open (SNP_LIST, "/home/dkuhn/6kSNPchip/SNPlist_test10.txt") || die "Can't open SNP list";
		while (<SNP_LIST>)
		{
			chomp;
		($SNPlist) = split /\t/ ;
		$SNPlist{$SNP_lookup{$SNPlist}} = $SNP_lookup{$SNPlist};
}
	open (OUTFILE, ">>/home/dkuhn/6kSNPchip/list_out.txt" || die "Can't open list_out.txt");
for $cultivar (sort keys %cultivarlist) {
		if (exists $hoh{$cultivar})  {
			for $SNPlist (sort keys %SNPlist) {
			if (exists $hoh{$cultivar}{$SNPlist}) {
			print OUTFILE "$cultivar\t$SNPlist\t$hoh{$cultivar}{$SNPlist}\n";
			$hoh2{$cultivar}{$SNPlist} = $hoh{$cultivar}{$SNPlist};	
			}
			}
		}
		}
	close OUTFILE;
#Call coreSNPs as a subroutine to get table printed out!!!!
	&list_to_table;

sub list_to_table {	
my $SNP = ();	#Initialize variables as undefined.
my $cultivar = ();
my $genotype = ();
my $header = "Cultivar_name\t";  #Initialize header variable with cultivar name
my %hash1 = ();
my %hash2 = ();
my %hash3 = ();

open (INFILE, "/home/dkuhn/6kSNPchip/list_out.txt")|| die "Can't open infile";  #open the file to be read.
	while (<INFILE>) #read in from file, line by line
{
	chomp;	#remove end of line character
	if (m/^WG/)   {
	($cultivar, $SNP, $genotype) = split (/\t/); #assign variables from genotype data
	$hash1{$cultivar} .= "$genotype\t";  #append genotype for first individual to first SNP
	$hash2{$cultivar} += 1; #set a counter to see when this cultivar is read in again
	$hash3{$hash2{$cultivar}+10000} = $SNP;  #assign the SNP name for the genotypes for this cultivar and use the counter to assign the next SNP
#	print "This is the value of hash2 $hash2{$SNP} and this is the value of hash3 $hash3{$hash2{$SNP}+10000}\n";
	print "I see $cultivar and $hash1{$cultivar}\n";
		}
	}
	open (OUTFILE, ">>/home/dkuhn/6kSNPchip/table_out.txt")|| die "Can't open outfile";
	foreach $key (sort keys %hash3)
{
	print "This is the hash3 key $key and this is the hash3 value $hash3{$key}\n";

	$header .= "$hash3{$key}\t";  #append all the SNP names to make the header row and print to file.
}
	print OUTFILE "$header\n";
	foreach $key (sort keys %hash1)
{	
	
	print OUTFILE "$key\t$hash1{$key}\n";  #Print to the same file the SNPs for each genotype with a tab between the genotype for each SNP.
}	
}
