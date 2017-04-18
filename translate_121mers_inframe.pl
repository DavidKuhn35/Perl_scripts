#!/usr/bin/perl
#Takes a parsed blastx file with the name of a 121mer that had an X in its translation and the translation frame and creates a fasta file with the translated versions of both possible sequences for the SNP
use Bio::Seq; 
use Bio::SeqIO;

my %hoh = ();
my @SNP = ();

open (INFILE, "$ARGV[0]") || die "Can't open infile.\n";
	while (<INFILE>) #read in from file, line by line
{
	chomp;	#remove end of line character
	if (m/^Tcm/) {	#skip any header information, start at genotype data
		@SNP = split (/\t/); #assign variables from genotype data
		$SNPframe = shift(@SNP);
		print "$SNPframe\n";
		$frame = pop(@SNP);
		print "$frame\n";
		$hoh{$SNPframe}{frame} = $frame;
		}
}
	

open (OUTFILE, ">>$ARGV[0]\_name_frame.txt") || die "Can't open outfile";  
#Step through the new hoh, make two new hashes, one for the SNP names and one for the appended genotypes.
for $key (sort keys %hoh) {
	print OUTFILE "$key\_1\t$hoh{$key}{frame}\n$key\_2\t$hoh{$key}{frame}\n";	
        }
	
	
