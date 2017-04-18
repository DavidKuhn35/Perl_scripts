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
#my %hoh2 = ();
my %SNPlist = ();
my %cultivarlist = ();
#my $header = "Tree name\t";
%same_genotype = ();
%different_genotype = ();
%homo_SNP = ();
%het_SNP = ();
%missing_data_SNP = ();
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
open (INFILE, "/home/dkuhn/6kSNPchip/6k_valid_data_recalled.txt")|| die "Can't open infile";  
	while (<INFILE>) #read in from file, line by line
{
	chomp;	#remove end of line character
	if (m/^contig/) {	#skip any header information, start at genotype data
	($SNP, $wellnumber, $allele1,$allele2) = split (/\t/); #assign variables from genotype data
	$hoh{$cultivar_lookup{$wellnumber}}{$SNP_lookup{$SNP}} = "$allele1$allele2";
			}
}



	open (LIST, "/home/dkuhn/6kSNPchip/$ARGV[0]")|| die "Can't open infile";  #open the file to be read.
	while (<LIST>)
{
	chomp;
	(undef,$cultivar1) = split (/\t/); 
	for $key (sort keys %{$hoh{$cultivar1}}) {
		if ($hoh{$cultivar1}{$key} =~ /AA|CC|GG|TT/) {
		$homo_SNP{$cultivar1} += 1; #make a counter to see how many times SNP is homozygous
		}elsif ($hoh{$cultivar1}{$key} =~ /AC|CA|AG|GA|CT|TC|GT|TG/) {
		$het_SNP{$cultivar1} += 1; #make a counter to see how many times SNP is heterozygous		
		}else {
		$missing_data_SNP{$cultivar1} += 1;
	}
}
}
	open (OUTFILE, ">>$ARGV[0]_heterozygosity_out.txt")|| die "Can't open outfile";
	for $key (sort keys %homo_SNP)
{
	printf OUTFILE "$key\t$homo_SNP{$key}\t%5.2f\t$het_SNP{$key}\t%5.2f\t$missing_data_SNP{$key}\t%5.2f\n", $homo_SNP{$key}/($homo_SNP{$key}+$het_SNP{$key}+$missing_data_SNP{$key}), $het_SNP{$key}/($homo_SNP{$key}+$het_SNP{$key}+$missing_data_SNP{$key}), $missing_data_SNP{$key}/($homo_SNP{$key}+$het_SNP{$key}+$missing_data_SNP{$key}) ;
}
