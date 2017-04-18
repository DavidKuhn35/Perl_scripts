#!/usr/bin/perl -w
#This program makes a table out of the genotype data from the SNP chip.
my $SNP = ();	#Initialize variables as undefined.
my $cultivar = ();
my $allele = ();
my @CGD = ();
#my $v2SNP = ();

my $header = "DNA_name\t";  #Initialize header variable with DNA name
#my $lastSNP = ARGV[0];
%hash1 = ();
%hash2 = ();
%hash3 = ();

open (INFILE, "/home/dkuhn/FrostyPod_2011_working/FPQTL_haplotypes/CATIE_Type1_mapped_SNP_genotypes.txt")|| die "Can't open infile";  #open the file to be read.
	while (<INFILE>) #read in from file, line by line
{
	chomp;	#remove end of line character
#	if (m/^[P U]/) {	#skip any header information, start at genotype data
#print "$_\n";

		($cultivar, $SNP, $allele) = split /\s+/, $_; #assign variables from genotype data
		($CGD) = split /\|/, $SNP;
#print  "Cultivar is $cultivar\t SNP is $SNP\t CGD is $CGD\n";
		open (DNA_NAMES, "/home/dkuhn/FrostyPod_2011_working/FPQTL_haplotypes/Frostypod_DNA_names.txt")|| die "Can't open DNA names file";
		while (<DNA_NAMES>)  #read in DNA names to choose among all data
{
		chomp;
#print "$_\n";
		if (m/$cultivar/)  {  #if cultivar matches one in the list
			open	(SNP_NAMES,"/home/dkuhn/FrostyPod_2011_working/FPQTL_haplotypes/FPQTL3_CGD_names.txt")|| die "Can't open SNP names file"; 
			while (<SNP_NAMES>)
{
#print "$_\n";
		if (m/$CGD/)  {  #if SNP name matches one in the list
		$hash1{$cultivar} .= "$allele\t";  #append genotype for first SNP for this individual
		$hash2{$cultivar} += 1; #set a counter to see when this cultivar is read in again
		$hash3{$hash2{$cultivar}} = $SNP;  #assign the SNP name for this set of SNPs and use the counter to assign the next cultivar
			}
			}
		}
		}
}

	open (OUTFILE, ">>/home/dkuhn/FrostyPod_2011_working/FPQTL_haplotypes/FPQTL3_CATIE_haplotype.txt")|| die "Can't open outfile";
	foreach $key (sort keys %hash3)
{
	$header .= "$hash3{$key}\t";  #append all the cultivar names to make the header row and print to file.
}
	print OUTFILE "$header\n";
	foreach $key (sort keys %hash1)
{	
	
	print OUTFILE "$key\t$hash1{$key}\n";  #Print to the same file the genotypes for each SNP with a tab between the genotype for each cultivar.
}	
