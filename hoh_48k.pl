#!/usr/bin/perl -w
#Read in the 48k SNP file and make a hash of hashes
open (SNP, "/home/david/coreSNPs/6kSNPs_Matv11_Name_Conversion.txt") || die "Can't open 6k SNP infile";
	while (<SNP>)
{
	chomp;
	if (m/^Tcm/)  	{	#If line has SNP name with Tcm
	($matv11, $v2contig, $lg, $pos, $chip, $seq, $v09cgd, $matv10, $matv09) = split (/\t/);
	$hoh {$matv11} {matv11} = $matv11;
	$hoh {$matv11} {v2contig} = $v2contig;
	$hoh {$matv11} {lg} = $lg;
	$hoh {$matv11} {pos} = $pos;
	$hoh {$matv11} {chip} = $chip;
	$hoh {$matv11} {seq} = $seq;
	$hoh {$matv11} {v09cgd} = $v09cgd;
	$hoh {$matv11} {matv10} = $matv10;
	$hoh {$matv11} {matv09} = $matv09;
			}
}
#Print out fasta file with v11 name and 121mer sequence
open (OUTFILE, ">>6kSNPs.fasta") || die "Can't open outfile";

for $v11_id (sort keys %hoh) {
		
	print OUTFILE "\>$v11_id\n $hoh{$v11_id}{seq}\n";
}
