#!c:\perl\bin

open(fp1, "gene_info.txt");
print "Reading gene information file.....";
while(<fp1>)
{
	$current_line = $_;
	chomp($current_line);
	@entries = split(/\t/,$current_line);
	# $entries[0] = gene id
	# $entries[1] = gene symbol
	# $entries[2] = chromosome number
	# $entries[3] = genbank id
	$hash_gene_to_symbol{$entries[0]} = $entries[1];
	$hash_gene_to_chr{$entries[0]} = $entries[2];
	$hash_gene_to_genbank{$entries[0]} = $entries[3];
}
close(fp1);
print "Done!\n\n";


open(fp1, "gene_data.txt");
print "Reading gene data file.....";
while(<fp1>)
{
	$current_line = $_;
	@entries = split(/\t/,$current_line);
	# $entries[0] = gene id
	# $entries[1] = Score for sample 1
	# $entries[2] = Call for sample 1
	# $entries[3] = Score for sample 2
	# $entries[4] = Call for sample 2
	$hash_gene_to_s1{$entries[0]} = $entries[1];
	$hash_gene_to_s2{$entries[0]} = $entries[3];
}
close(fp1);
print "Done!\n\n";

print "Creating merged report....";
open(fp1, ">perlproject_report.txt");
print fp1 "Gene Name\tScore1\tScore2\tAverageScore\tGene Symbol\tChromosome\tGenbank\n";
@genes = sort {$a cmp $b} keys %hash_gene_to_s1;
foreach $gene (@genes)
{
	$s1 = $hash_gene_to_s1{$gene};
	$s2 = $hash_gene_to_s2{$gene};
	$avg = mean($s1,$s2);
	$symb = $hash_gene_to_symbol{$gene};
	$chr = $hash_gene_to_chr{$gene};
	$chr =~ s/chr//;
	$genbank = $hash_gene_to_genbank{$gene};
	print fp1 "$gene\t$s1\t$s2\t$avg\t$symb\t$chr\t$genbank\n";
}
print "Done\n";
close(fp1);

print "Creating statistics for Sample 1....\n";
@s1_scores = values %hash_gene_to_s1;
print "The min value of sample 1 is ".min(@s1_scores)."\n";
print "The max value of sample 1 is ".max(@s1_scores)."\n";
print "The mean value of sample 1 is ".mean(@s1_scores)."\n";

print "Creating statistics for Sample 2....\n";
@s2_scores = values %hash_gene_to_s2;
print "The min value of sample 2 is ".min(@s2_scores)."\n";
print "The max value of sample 2 is ".max(@s2_scores)."\n";
print "The mean value of sample 2 is ".mean(@s2_scores)."\n";

#########FUNCTIONS
sub mean
{
	@input_nums = @_;
	$avg1 = 0;
	for($i=0; $i<scalar(@input_nums); $i=$i+1) 
	{
		$avg1 = $avg1 + $input_nums[$i];
	}
	$avg1 = $avg1/scalar(@input_nums);
	return $avg1;
}

sub min
{
	@input_nums = @_;
	@sorted = sort {$a <=> $b} @input_nums;
	return $sorted[0];
}

sub max
{
	@input_nums = @_;
	@sorted = sort {$a <=> $b} @input_nums;
	return $sorted[scalar(@input_nums)-1];
}

