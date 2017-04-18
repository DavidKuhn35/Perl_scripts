#!/usr/bin/perl -w
%SNP_lookup = ();
$SNPpos1 = ();
$SNP = ();
%hash1 = ();
%hash2 = ();
%hash3 = ();
%hash4 = ();
%hash5 = ();
%hash6 = ();
%hash7 = ();
$SNPnocontig = ();
$SNPcontig = ();
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

open (INFILE, "/home/dkuhn/6kSNPchip/CATIE_Type1_SNP_Chip_data_JM_coded.txt") || die "Can't open CATIE coded SNPchip data";
while (<INFILE>)	
{
chomp;
	if (m/^\d/)	{
	($SNPnocontig, $SNPcode) = split (/\t/);
	$SNPcontig = $SNPnocontig;
	$SNPcontig =~ s/$SNPnocontig/contig$SNPnocontig/;
	print "$SNPcontig\n";
	$hash1{$SNP_lookup{$SNPcontig}} = $SNP_lookup{$SNPcontig};
	$hash2{$SNP_lookup{$SNPcontig}} = $SNPcode;
	print "$SNP_lookup{$SNPcontig}\t$hash2{$SNP_lookup{$SNPcontig}}\n";
	$newSNPtonocontig{$SNP_lookup{$SNPcontig}} = $SNPnocontig;
	}
}
#Sort mapped SNPs for LG2 and nnxnp code
	foreach $key (sort keys %hash2)
{
		if ($key =~ /Tcm005s/)	{
			if ($hash2{$key} =~ /<nnxnp>|<lmxll>|<hkxhk>/)	{
			$hash3{$key} = $hash2{$key};
			}
		}
}
open (INFILE, "/home/dkuhn/6kSNPchip/CATIE_SNP_map.txt") || die "Can't open SNP map infile";
	while (<INFILE>)	
{
	chomp;
	if (m/LG5/)	{
	($LG, $SNP_nocontig, $cMpos) = split /\t/;
	$SNP_contig = $SNP_nocontig;
	$SNP_contig =~ s/$SNP_nocontig/contig$SNP_nocontig/;
	$hash4{$SNP_lookup{$SNP_contig}} = $SNP_lookup{$SNP_contig};
	$cMpos =~ s/(\d+).\d+/$1/;
	$hash5{$SNP_lookup{$SNP_contig}} = $cMpos;
	%hash6 = reverse %hash5;
	%hash7 = reverse %hash6;
	foreach $key (sort keys %hash7)	{
		if (exists $hash3{$key})	{
		$hash8{$key} = $newSNPtonocontig{$key};
		}
		}
	}
}
	
open (OUTFILE, ">>/home/dkuhn/6kSNPchip/CATIE_LG5_spacedSNPs_out.txt") || die "Can't open CATIE spaced SNPs outfile";
	foreach $key	(sort keys %hash8)
{
	print OUTFILE "$key\t$hash8{$key}\n";
}
open (OUTFILE, ">>/home/dkuhn/6kSNPchip/CATIE_LG5_spaced_SNPs_all_out.txt") || die "Can't open CATIE coded table outfile";
open (INFILE, "/home/dkuhn/6kSNPchip/CATIE_Type1_SNP_Chip_data_JM_coded.txt") || die "Can't open CATIE coded SNPchip data";
while (<INFILE>)	
{
chomp;
	if (m/SNP name/)	{
	print OUTFILE "$_\n";
}	else	{
		foreach $key (sort keys %hash8)
		{
			if (m/^$hash8{$key}/)	{
			print OUTFILE "$_\n";
			}
		}
}		
}



	
