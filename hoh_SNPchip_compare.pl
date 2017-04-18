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
	$ref = $ARGV[0];	
		open (LIST, "/home/dkuhn/6kSNPchip/$ARGV[1]")|| die "Can't open infile";  #open the file to be read.
	while (<LIST>)
{
	chomp;
	(undef,$cultivar1) = split (/\t/); 
	for $key (sort keys %{$hoh{$ref}}) {
		if ($hoh{$ref}{$key} eq $hoh{$cultivar1}{$key}) {
		$same_genotype{$key} = $hoh{$ref}{$key};
#	print "The value of same_genotype is $same_genotype{$key}\n";
		$same_SNP{$key} += 1; #make a counter to see how many times SNP is same genotype
#		print "SNP $key is the same SNP $same_SNP{$key} times\n";
		$same_SNP_cultivar{$key} .= "$cultivar1\t$hoh{$cultivar1}{$key}\t"; #says which cultivars have the same genotype at the SNP
		}else {
		$different_genotype{$key} = $hoh{$cultivar1}{$key};
		$different_SNP{$key} += 1; #make a counter to see how many times SNP is a different genotype
		$different_SNP_cultivar{$key} .= "$cultivar1\t$hoh{$cultivar1}{$key}\t";  #says which cultivars have a different genotype at the SNP and what the genotype is.
#		print "The SNP where the genotypes differ is $key.\n";
		}
	}
	$same_count = keys %same_genotype;
	$diff_count = keys %different_genotype; 
open (OUTFILE, ">>$ARGV[0]_vs_$ARGV[1]")|| die "Can't open outfile";
	printf OUTFILE "$ARGV[0] vs $cultivar1\tsame is ($same_count / ($same_count + $diff_count))\t %5.2f \t different is ($diff_count / ($same_count + $diff_count)) %5.2f \n", $same_count/($same_count + $diff_count), $diff_count/($same_count + $diff_count);
	%same_genotype = ();
	%different_genotype = (); 
}
#open (OUTFILE, ">>compare_$ARGV[0]_vs_$ARGV[1]_out.txt")||die "Can't open compare outfile";
	for $key (sort keys %different_SNP_cultivar)
{
#		print OUTFILE "SNP $key is the same SNP in these cultivars $same_SNP_cultivar{$key}\n";		
	print OUTFILE "$key\t$hoh{$ARGV[0]}{$key}\n$different_SNP_cultivar{$key}\n";
}		
	for $key (sort keys %different_SNP)
{
	print OUTFILE "$key\t$different_SNP{$key}\n";
}
