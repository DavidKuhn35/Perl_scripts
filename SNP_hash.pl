#!C:\strawberry\perl\bin\perl

# Loading data into hash arrays
%hash1_SNP_to_genotype = ();
#%hash2_probe_to_genename = ();
open(fp, "SNP_little_test.txt") || die "Cant open file\n";
while(<fp>)
{
	($SNP, $cultivar, $allele1, $allele2) = split(/\t/);
#	$concat=join("::", $genbank, $locuslink, $genename);
#	$hash_concat{$probe} = $concat;
	$hash1_SNP_to_genotype{$SNP} = "$allele1 \t $allele2";
#	$hash2_probe_to_genename{$probe} = $genename;
}
close(fp);

	foreach $SNP (%hash1_SNP_to_genotype)
{
	print "%hash1_SNP_to_genotype{$SNP} \n";
}

# Using the data in hash arrays to create a report
#@someprobes = ("1020_s_at", "1058_at", "37902_at", "999_at");
#	@keys = keys %hash1_probe_to_genbank;
#foreach $probe (@someprobes)
#{
#	print "For probe $probe, the genbank id is $hash1_probe_to_genbank{$probe}";
#	print " and the gene name is $hash2_probe_to_genename{$probe}\n";
#	print "For probe $probe, the concatenated annotation is $hash_concat{$probe}\n";
#}

