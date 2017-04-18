#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
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
my %hash1 = ();
my %hash2 = ();
my %SNPlist = ();
my %cultivarlist = ();
my %cultivar_lookup = ();
my %SNP_rev_lookup = ();
my %cultivar_rev_lookup = ();
my %cultquery = ();
my %SNPquery = ();
my $header = "Tree_name\tWell_number\t";
my $v2_header = "Tree_name\tWell_number\t";
#my $genotype = ();
#Read in the SNP conversion file for the entire SNP chip and make a hash
open (SNP, "48kSNPs_Matv1.1_Name_Conversion.txt") || die "Can't open 48kSNP name conversion infile";
	while (<SNP>)
{
	chomp;
	if (m/^Tcm/)  	{	#If line has SNP name with Tcm
	(my $Tcm_name,my $v2_name) = split (/\t/);
	$SNP_lookup{$v2_name} = $Tcm_name; # new SNP name with genome position => old SNP name
	$SNP_rev_lookup{$Tcm_name} = $v2_name;
			}
}
open (MANIFEST, "DNA_Manifest_renamed_complete.txt") || die "Can't open manifest infile";
	while (<MANIFEST>)
{
	chomp;
	if (m/^WG/)	{	#If line has cultivar name in it
	(my $wellnumber, my $cultivarname) = split (/\t/);
	$cultivar_lookup{$wellnumber} = $cultivarname; #cultivar name => well number
	$cultivar_rev_lookup{$cultivarname} = $wellnumber;
			}
} 
close MANIFEST;
#Read in the SNPchip data and convert SNP names to SNP position names using hash and wellnumbers to cultivar names using hash
open (INFILE, "6k_valid_data.txt")|| die "Can't open infile";  
	while (<INFILE>) #read in from file, line by line
{
	chomp;	#remove end of line character
	if (m/^contig/) {	#skip any header information, start at genotype data
	(my $SNP, my $wellnumber, my $allele1,my $allele2) = split (/\t/); #assign variables from genotype data
	$hoh{$cultivar_lookup{$wellnumber}}{$SNP_lookup{$SNP}} = "$allele1$allele2"; #hoh with cultivar name as first key and SNP name as second key, genotype as value
			}
}
	
	foreach my $key (sort keys %cultivar_lookup)	{	#look through all well numbers
#	if ($cultivar_lookup{$key} =~ /^CATIE/	{	#Gets all CATIE
	if ($cultivar_lookup{$key} =~ /PNG_8/)  {  #Gets all trees from PNG population 1
		if (exists $hoh{$cultivar_lookup{$key}})	{
		my $cultname = $cultivar_lookup{$key};	#gives cultname variable the cultivar name, not the well number
		$cultquery{$cultname} = $cultname;	#makes a hash of cultivar names that can be sorted

			foreach my $key (sort keys %SNP_lookup)	{	#look through all SNP v2 names
			if ($SNP_lookup{$key} =~ /^Tcm/)	{	#Would get all SNPs
			if (exists $hoh{$cultname}{$SNP_lookup{$key}})	{ #hoh is cultivar name, Tcm => genotype
			$SNPquery{$SNP_lookup{$key}} = $SNP_lookup{$key}; #make hash of Tcm names that can be sorted on Tcm name
			}
}
}
}
}
}
for $cultivar (sort keys %cultquery) {
		if (exists $hoh{$cultivar})  {
			for my $key (sort keys %SNPquery) {
			if (exists $hoh{$cultivar}{$key}) {
			$hoh2{$cultivar}{$key} = $hoh{$cultivar}{$key};	#Make a new hash of hashes that only contains the information you are looking for.  This will make it easier to print out as a table. Cultivar key 1 Tcm name key 2 => genotype
			}
			}
		}
		}

open (OUTFILE, ">>PNG_8_SNP_out.txt") || die "Can't open outfile";  
#Step through the new hoh, make two new hashes, one for the SNP names and one for the appended genotypes.
for my $x (sort keys %hoh2) {
		for my $y (sort keys %{ $hoh2{$x} }) {
		$hash1{$x} .= "$hoh2{$x}{$y}\t";
		$hash2{$y} = $y;		
        }
	}
	for my $key (sort keys %hash2) {
		$header .= "$key\t";
		$v2_header .= "$SNP_rev_lookup{$key}\t";
		}
	print OUTFILE "$v2_header\n$header\n";
	for my $key (sort keys %hash1) {	
	print OUTFILE "$key\t$cultivar_rev_lookup{$key}\t$hash1{$key}\n";
}
