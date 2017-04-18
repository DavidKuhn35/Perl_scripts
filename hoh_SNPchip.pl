#!/usr/bin/perl -w
#This program makes a hash of hashes from SNPchip data read in from a file.
my $allele1 = ();
my $allele2 = ();
my $SNP = ();
my $cultivar = ();
my %hoh = ();
my %same_genotype = ();
my %different_genotype = ();
my %same_SNP = ();
my %different_SNP = ();
my %same_SNP_cultivar = ();
my %different_SNP_cultivar = ();
my $key = ();

my $SNPkey = ();
open (INFILE, "/home/dkuhn/CoreSNPs/coresnp_diversity_data.txt")|| die "Can't open infile";  #open the file to be read.
	while (<INFILE>) #read in from file, line by line
{
	chomp;	#remove end of line character
	if (m/^contig/) {	#skip any header information, start at genotype data
		($SNP, $cultivar, $allele1,$allele2) = split (/\t/); #assign variables from genotype data
	$hoh{$cultivar}{$SNP} = "$allele1$allele2";
#	print "For $cultivar, the $SNP genotype is $hoh{$cultivar}{$SNP}\n";
	}
}
	open (LIST, "/home/dkuhn/CoreSNPs/$ARGV[1]")|| die "Can't open infile";  #open the file to be read.
	while (<LIST>)
{
	chomp;
	$cultivar1 = $_; 
#	print "The value of cultivar $ARGV[0] and SNP $SNP is $hoh{$ARGV[0]}{$SNP}\n";
	for $key (sort keys %{ $hoh{$ARGV[0]} }) {
		if ($hoh{$ARGV[0]}{$key} eq $hoh{$cultivar1}{$key}) {
		$same_genotype{$key} = $hoh{$ARGV[0]}{$key};
#	print "The value of same_genotype is $same_genotype{$key}\n";
		$same_SNP{$key} += 1; #make a counter to see how many times SNP is same genotype
		print "SNP $key is the same SNP $same_SNP{$key} times\n";
		$same_SNP_cultivar{$key} .= "$cultivar1\t"; #says which cultivars have the same genotype at the SNP
		print "SNP $key is the same SNP in these cultivars $same_SNP_cultivar{$key}\n";
		}else {
		$different_genotype{$key} = $hoh{$ARGV[0]}{$key};
		$different_SNP{$key} += 1; #make a counter to see how many times SNP is a different genotype
		print "SNP $key is a different SNP $different_SNP{$key} times\n";
		$different_SNP_cultivar{$key} .= "$cultivar1\t";  #says which cultivars have a different genotype at the SNP
#		print "The SNP where the genotypes differ is $key.\n";
		print "SNP $key is a different SNP in these cultivars $different_SNP_cultivar{$key}\n";
	}
}
	$same_count = keys %same_genotype;
	$diff_count = keys %different_genotype; 
open (OUTFILE, ">>/home/dkuhn/CoreSNPs/$ARGV[0]_vs_$ARGV[1]")|| die "Can't open outfile";
	printf OUTFILE "$ARGV[0] vs $cultivar1\tsame is ($same_count / 5388)\t %5.2f \t different is ($diff_count / 5388) %5.2f \n", $same_count/5388, $diff_count/5388;
	%same_genotype = ();
	%different_genotype = (); 
}
		
open (OUTFILE, ">>/home/dkuhn/CoreSNPs/$ARGV[0]_vs$ARGV[1]_SNP_out.txt")|| die "Can't open outfile";
	foreach $key (sort keys %same_SNP)
{
	print OUTFILE "$key\t$same_SNP{$key}\n";
	print OUTFILE "$key\t$same_SNP_cultivar{$key}\n";
	print OUTFILE "$key\t$different_SNP{$key}\n";
	print OUTFILE "$key\t$different_SNP_cultivar{$key}\n";
}	 
