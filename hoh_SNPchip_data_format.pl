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
	
	open (CULTIVAR_LIST, "/home/dkuhn/6kSNPchip/Nanay_4_test.txt") || die "Can't open wellnumber_list_test";
	while (<CULTIVAR_LIST>)
{		
	chomp;
#	($cultivar) = split /\t/ ;
#	$cultivarlist{$cultivar_lookup{$cultivar}} = $cultivar_lookup{$cultivar};
	($wellnumber, $cultname) = split /\t/;
#	if ($cultname =~ /^DIV/)  {
#		$cultquery{$cultname} = $cultname;
	if (exists $hoh{$cultname})	{
	$cultquery{$cultname} = $cultname;
		open (SNP_LIST, "/home/dkuhn/6kSNPchip/SNPlist_test10.txt") || 			die "Can't open SNP list";
		while (<SNP_LIST>)
{
		chomp;
		($oldSNP, $newSNP) = split /\t/ ;
#		$SNPlist{$SNP_lookup{$oldSNP}} = $SNP_lookup{$oldSNP};
		if (exists $hoh{$cultname}{$newSNP})	{
		$SNPquery{$newSNP} = $newSNP;
		}
}
}
}
for $cultivar (sort keys %cultquery) {
		if (exists $hoh{$cultivar})  {
			for $key (sort keys %SNPquery) {
			if (exists $hoh{$cultivar}{$key}) {
			$hoh2{$cultivar}{$key} = $hoh{$cultivar}{$key};	
			}
			}
		}
		}

open (OUTFILE, ">>test_out.txt") || die "Can't open outfile";

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
