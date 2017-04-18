#!/usr/bin/perl -w
my $dummy1 = ();
my $SNP = ();
my $SNPpos1 = ();
#my $newSNPpos1 = ();
#my $cultivar = ();
#my $cultivar1 = ();
#my $allele1 = ();
#my $allele2 = ();
my %SNP_lookup = ();
my %hoh = ();
my $contig = ();
my $phase = ();
my $code = ();
#my %hoh2 = ();
#my %SNPlist = ();
#my %cultivarlist = ();
#my $header = "Tree name\t";
#%same_genotype = ();
#%different_genotype = ();
#%homo_SNP = ();
#%het_SNP = ();
#%missing_data_SNP = ();
#my $genotype = ();
#Read in the SNP conversion file for the entire SNP chip and make a hash
open (SNP, "/home/dkuhn/6kSNPchip/6kSNPs_sortable_Aug11.txt") || die "Can't open SNP infile";
	while (<SNP>)
{
	chomp;
	if (m/^Tcm/)  	{	#If line has SNP name with Tcm
	($dummy1,$SNPpos1, $contig) = split (/\t/);
	$SNP = $contig;
	$SNP =~ s/^contig// ;
	$SNP_lookup{$SNP} = $SNPpos1; #old SNP name => new SNP name with genome position
			}
}

open (LOC, "/home/dkuhn/map_loc/Maps_for_Assembly/3_CATIE_Type1_SNP_map/CATIE_MAP_LOC_files/CATIE_T1_LG02_13_inverted_truncated.loc") || die "Can't open LOC file.";
	while (<LOC>)
{
	chomp;
	if (m/^;/)	{
	$comment = $_;
	}elsif(m/^name/){
	$name = $_;
	}elsif (m/^popt/){	
	$popt = $_;
	}elsif(m/^nloc/)  {
	$nloc = $_;
	}elsif (m/^nind/)  {
	$nindo = $_;
	}elsif	(m/^\d+_S_/)	{
	($contig, $code, $phase) = split (/\s+/);
	print "$contig\t$code\t$phase\n";
	$hoh{$SNP_lookup{$contig}}{code} = $code;
	$hoh{$SNP_lookup{$contig}}{phase} = $phase;
	}elsif	(m/^\s+h|^\s+l|^\s+n|^\s+k/)	{
	$hoh{$SNP_lookup{$contig}}{genotype} .= $_;
	open (LOG, ">>map_loc_log.txt") || die "Can't open log file";
	print LOG "$contig\t$SNP_lookup{$contig}\t$hoh{$SNP_lookup{$contig}}{genotype}\n";
	close LOG;
	}
}

open (OUTFILE, ">>CATIE_T1_loc_out.txt") || die "Can't open LOC outfile.";
	print OUTFILE "$comment\n\n$name\n$popt\n$nloc\n$nindo\n\n";

	for $key (sort keys %hoh)
{ 
	print OUTFILE	"$key\t$hoh{$key}{code}\t$hoh{$key}{phase}\n$hoh{$key}{genotype}\n";
}
