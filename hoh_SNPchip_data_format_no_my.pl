#!/usr/bin/perl -w

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

for $cultivar (sort keys %cultivarlist) {
		if (exists $hoh{$cultivar})  {
			for $SNPlist (sort keys %SNPlist) {
			if (exists $hoh{$cultivar}{$SNPlist}) {
			print "$cultivar\t$SNPlist\t$hoh{$cultivar}{$SNPlist}\n";
			$hoh2{$cultivar}{$SNPlist} = $hoh{$cultivar}{$SNPlist};	
			}
			}
		}
		}


open (OUTFILE, ">>SNPchip_genotypes_out.txt") || die "Can't open outfile";
for $key (sort keys %{$hoh2{$cultivar}})  {
	$header .= "$key\t";
#	print OUTFILE "$key\n";
}
	print OUTFILE "$header\n";
my %hash1 = ();
for $cultivar (sort keys %hoh2) {
	print "$cultivar:  ";
		for $key (sort keys %{ $hoh2{$cultivar} }) {
#		print OUTFILE "$cultivar\t$key\t$hoh2{$cultivar}{$key}\n";
		$genotype = $hoh2{$cultivar}{$key};
		$hash1{$cultivar} .= "$genotype\t";
		print "$key=\t$hoh2{$cultivar}{$key}\n";
		print "$hash1{$cultivar}\n";
		}
#open (OUTFILE, ">>append_genotypes_out.txt") || die "Can't open genotypes outfile";
#		print OUTFILE "$cultivar\t  $genotype  \n";
	print "\n";
	}

	for $key (sort keys %hash1) {	
		print "$key\t$hash1{$key}\n";
		print OUTFILE "$key\t$hash1{$key}\n";
}
