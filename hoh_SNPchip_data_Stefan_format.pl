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
#Read in the SNP conversion file for the entire SNP chip and make a hash
open (SNP, "/home/dkuhn/6kSNPchip/6kSNPs_sortable_Aug11.txt") || die "Can't open SNP infile";
	while (<SNP>)
{
	chomp;
	if (m/^Tcm/)  	{	#If line has SNP name with Tcm
	($dummy1,$SNPpos1, $SNP) = split (/\t/);
	$SNP_lookup{$SNP} = $SNPpos1; #old SNP name => new SNP name with genome position
			}
}
open (MANIFEST, "/home/dkuhn/6kSNPchip/DNA_Manifest_renamed_complete.txt") || die "Can't open manifest infile";
	while (<MANIFEST>)
{
	chomp;
	if (m/^WG/)	{	#If line has cultivar name in it
	($wellnumber, $cultivarname) = split (/\t/);
	$cultivar_lookup{$wellnumber} = $cultivarname; #well number => cultivar name
			}
} 
close MANIFEST;
#Read in the SNPchip data and convert SNP names to SNP position names using hash and wellnumbers to cultivar names using hash
open (INFILE, "/home/dkuhn/6kSNPchip/6k_valid_data.txt")|| die "Can't open infile";  
	while (<INFILE>) #read in from file, line by line
{
	chomp;	#remove end of line character
	if (m/^contig/) {	#skip any header information, start at genotype data
	($SNP, $wellnumber, $allele1,$allele2) = split (/\t/); #assign variables from genotype data
	$hoh{$cultivar_lookup{$wellnumber}}{$SNP_lookup{$SNP}} = "$allele1$allele2";
			}
}
	
#	open (CULTIVAR_LIST, "/home/dkuhn/6kSNPchip/DNA_Manifest_renamed_complete.txt") || die "Can't open wellnumber_list_test";  #Put cultivar list in here!  Should change to ARGV to make it more convenient.  Determine format for correct splitting.
#	while (<CULTIVAR_LIST>)
#{		
#	chomp;
#	if (m/^WG/)	{
#	($cultivar) = split /\t/ ;
#	$cultname = $cultivar_lookup{$cultivar};
#	($wellnumber, $cultname) = split /\t/;
	foreach $key (sort keys %cultivar_lookup)	{
	if ($cultivar_lookup{$key} =~ /^(CATIE_P|CATIE_T1)/)  {
		if (exists $hoh{$cultivar_lookup{$key}})	{
		$cultname = $cultivar_lookup{$key};
		$cultquery{$cultname} = $cultname;
#	$cultquery{$cultname} = $cultname;
#		open (SNP_LIST, "/home/dkuhn/6kSNPchip/6kSNPs_sortable_Aug11.txt") || 			die "Can't open SNP list";  #Put SNPlist in here, if necessary.  Can change to ARGV to make more convenient.  Determine format for correct splitting.
#		while (<SNP_LIST>)
#{
#		chomp;
#		($oldSNP, $newSNP) = split /\t/ ;  #Depends on SNP list format.  Really need to get that standardized.  Perhaps do lookup for all SNPs.
#		$SNPlist{$SNP_lookup{$oldSNP}} = $SNP_lookup{$oldSNP};
#		$newSNP = $SNP_lookup{$oldSNP};
	foreach $key (sort keys %SNP_lookup)	{
		if ($SNP_lookup{$key} =~ /001s/)	{	#Would get all LG1 SNPs
		if (exists $hoh{$cultname}{$SNP_lookup{$key}})	{
		$SNPquery{$SNP_lookup{$key}} = $SNP_lookup{$key};
		}
}
}
}
}
}
for $cultivar (sort keys %cultquery) {
		if (exists $hoh{$cultivar})  {
			for $key (sort keys %SNPquery) {
			if (exists $hoh{$cultivar}{$key}) {
			$hoh2{$cultivar}{$key} = $hoh{$cultivar}{$key};	#Make a new hash of hashes that only contains the information you are looking for.  This will make it easier to print out as a table.
			}
			}
		}
		}

open (OUTFILE, ">>CATIE_T1_All_SNPs_genotypes_out.txt") || die "Can't open outfile";  #Rename outfile.  Should be ARGV based on Cultivar infile
#Step through the new hoh, make two new hashes, one for the SNP names and one for the appended genotypes.
for $cultivar (sort keys %hoh2) {
		for $key (sort keys %{ $hoh2{$cultivar} }) {
		$hash1{$cultivar} .= "$hoh2{$cultivar}{$key}\t";
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
