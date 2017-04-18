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
my $snpheader = "SNP name\t";

#Read in the SNPchip data 
open (INFILE, "/home/david/coreSNPs/6k_valid_data.txt")|| die "Can't open infile";  
	while (<INFILE>) #read in from file, line by line
{
	chomp;	#remove end of line character
	if (m/^contig/) {	#skip any header information, start at genotype data
	($SNP, $wellnumber, $allele1,$allele2) = split (/\t/); #assign variables from genotype data
	$hoh{$wellnumber}{$SNP} = "$allele1$allele2";
	$hoh3{$SNP}{$wellnumber} = "$allele1$allele2";
			}
}
	

open (OUTFILE, ">>all_6k_SNP_out_WG_rows.txt") || die "Can't open outfile";  
#Step through the new hoh, make two new hashes, one for the SNP names and one for the appended genotypes.
for $cultivar (sort keys %hoh) {
		for $key (sort keys %{ $hoh{$cultivar} }) {
		$hash1{$cultivar} .= "$hoh{$cultivar}{$key}\t";
		$hash2{$key} = $key;		
        }
	}
	for $key (sort keys %hash2) {
		$header .= "$key\t";
		}
	print OUTFILE "$header\n";
	for $key (sort keys %hash1) {	
	print OUTFILE "$key\t$hash1{$key}\n";
}

open (OUTFILE, ">>all_6k_SNP_out_WG_columns.txt") || die "Can't open outfile";  
#Step through the new hoh, make two new hashes, one for the SNP names and one for the appended genotypes.
for $snp (sort keys %hoh3) {
		for $key (sort keys %{ $hoh3{$snp} }) {
		$hash_1{$snp} .= "$hoh3{$snp}{$key}\t";	#appends genotype to snp
		$hash_2{$key} = $key;			#makes hash of well numbers
        }
	}
	for $key (sort keys %hash_2) {
		$snpheader .= "$key\t";
		}
	print OUTFILE "$snpheader\n";
	for $key (sort keys %hash_1) {	
	print OUTFILE "$key\t$hash_1{$key}\n";
}


