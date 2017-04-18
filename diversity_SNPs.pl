#!/usr/bin/perl -w
#This program makes a table out of the genotype data from the SNP chip.
my $SNP = ();	#Initialize variables as undefined.
my $cultivar = ();
my $allele1 = ();
my $allele2 = ();
my $header = "SNPname\t";  #Initialize header variable with SNPname
#my $lastSNP = ARGV[0];
%hash1 = ();
%hash2 = ();
%hash3 = ();

open (INFILE, "/home/dkuhn/CoreSNPs/Cacoa_Validation_FinalReport.txt")|| die "Can't open infile";  #open the file to be read.
	while (<INFILE>) #read in from file, line by line
{
	chomp;	#remove end of line character
	if (m/^contig/) {	#skip any header information, start at genotype data
		($SNP, $cultivar, $allele1,$allele2) = split (/\t/); #assign variables from genotype data
		$hash1{$SNP} .= "$allele1$allele2\t";  #append genotype for first individual to first SNP
		$hash2{$SNP} += 1; #set a counter to see when this SNP is read in again
		$hash3{$hash2{$SNP}} = $cultivar;  #assign the cultivar name for this set of SNPs and use the counter to assign the next cultivar
	}
}
	open (OUTFILE, ">>/home/dkuhn/CoreSNPs/diversity_SNP_table.txt")|| die "Can't open outfile";
	foreach $key (sort keys %hash3)
{
	$header .= "$hash3{$key}\t";  #append all the cultivar names to make the header row and print to file.
}
	print OUTFILE "$header\n";
	foreach $key (sort keys %hash1)
{	
	
	print OUTFILE "$key\t$hash1{$key}\n";  #Print to the same file the genotypes for each SNP with a tab between the genotype for each cultivar.
}	
