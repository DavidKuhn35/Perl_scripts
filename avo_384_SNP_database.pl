#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			:  April 10, 2017
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	:  creates a 384 SNP genotype file for a list of cultivars
# <input>		:  SNPchip data, DNA manifest, list of 384 SNPs
# <output>		:  genotype file for a particular cultivar with SNP names
#
# Dependencies	: N/A
# Usage			: joinmap_rawdata.pl
#
##################################################################################
use strict;
use warnings;
use diagnostics;
use Bio::Seq;
use Bio::SeqIO;
use Bio::SeqUtils;
use Bio::Tools::OddCodes;
use Bio::PrimarySeq;

#####################################################################################################################
#
my %hoh         = ();
my %hoh2        = ();
my %snp_hash    = ();
my $snp         = 0;
my $dna         = 0;
my $geno1       = 0;
my $geno2       = 0;
my %WG          = ();
my %CV			= ();
my %snp384		=  ();
my %pop1        = ();
my $popname     = ();
my %mapsnp_hash = ();
my $mapsnp      = 0;
my $seg_type    = 0;
my $parent1     = 0;
my $parent2     = 0;
my %p1mapsnp    = ();
my %p2mapsnp    = ();
my %pop1_raw    = ();
my $out_header  = "SNP name\t";
my $out_geno    = ();

#####################################################################################################################
#	main	Read in SNP chip data.
open( SNPDATA, $ARGV[0] ) || die "Can't open avo SNP chip data $!\n";
while (<SNPDATA>) {
    chomp;
    if (m/^SHRSPa/)
    {    #if line has SNP name in first column, avocado SHRSPa  cacao contig
        ( $snp, $dna, $geno1, $geno2 )
            = split /,/;    #split on tab for cacao, comma for avocado
        my @genotype = ( $geno1, $geno2 );    #put genotypes in an array
        my @genosort
            = sort (@genotype);    #sort genotypes ASCIIbetically A,C,G,T
        $hoh{$snp}{$dna} = "$genosort[0]$genosort[1]"
            ;    #hoh with SNP key and dna key and genotype value
        $snp_hash{$snp} = $snp;
        $hoh2{$dna}{$snp} = "$genosort[0]$genosort[1]"
            ;    #hoh2 with dna key and snp key and genotype value
    }
}

######################################################################################################################
# Read in manifest data and make hashes of parents
open( DNADATA, $ARGV[1] ) || die "Can't open avo manifest data $!\n";
#open (CVNAMES_TEST, ">>cvnamestest.txt") || die "Can't open cvnamestest $!\n";
while (<DNADATA>) {
    chomp;
    if (m/^\d/) {    #if line has DNA name data  (WG cacao)  (\d avocado)
        my @dnadata = split /\t/;
			if ($dnadata[6] =~ (m/Control|GERM/) )	{
			$WG{ $dnadata[3] } = $dnadata[6]  ; #make hash with well numbers as keys and common names as values ([0] and [1] for cacao) ([3] [6] for avocado)
#			print CVNAMES_TEST "$dnadata[3]  $WG{$dnadata[3]} \n";
            $CV{$dnadata[6]} = $dnadata[3]; #make hash with cultivar name as key and DNA name as value, will get last of control genotypes when multiple names are identical
#			$out_header.= "$key\t"; #put names of cultivars as column headers
		}
	}
}
#######################################################################################################################
#  Get all members of a population and make hashes of WG names
#$popname = $ARGV[2];
#foreach my $key ( sort keys %CV ) {
#    if ( $WG{$key} =~ /$popname/ )
#     if ( $key =~ (/Control/ || /GERM/))     {         #get all controls and germplasm cultivars
#        $pop1{$key} = $key  ; #makes hash of well numbers that had popname as value in manifest
            
			foreach my $key (sort keys %WG)	{
            $out_header.= "$WG{$key}\t"; #put names of cultivars as column headers
            
    }

######################################################################################################################
#	Read in mappable snps, segregation patterns and parent 1 and parent 2 genotypes at snp.  Create hashes with seg_type and parent1 genotype and parent2 genotype
#open (MAPSNP, $ARGV[2]) || die "Can't open mappable snps file $!\n";
#	while (<MAPSNP>)		{
#		if (m/^SHRSPa/)	{	# (cacao contig) (avocado SHRSPa)
#			($mapsnp, $seg_type, $parent1, $parent2) = split /\t/;	#parent1 x parent2 cross
#			$mapsnp_hash{$mapsnp}= $seg_type;
#			$p1mapsnp{$mapsnp} = $parent1;
#			$p2mapsnp{$mapsnp} = $parent2;
#		}
#	}
#
#########################################################################################################################
#	Go through mappable SNPs for each genotype in a mapping population and make a file with genotype data.  snps are rows, genotypes are columns
#oreach my $key ( sort keys %mapsnp_hash )
#{    #go through mappable SNPs one at a time
###########################################################################################################################
#  Read in 384 SNPs
open( SNP384, $ARGV[2] ) || die "Can't open avo SNP 384 data $!\n";
while (<SNP384>) {
   chomp;
    if (m/SHRSPa/)	{
		 my $snppy = $_;
#		 print "$snppy\n";
    $snp384{$snppy} = $snppy;
}
    
}
	
#############################################################################################################################




	foreach my $s (sort keys %snp384)	{
		
			foreach my $x ( sort keys %WG ) 	{    # go through all pop1 genotypes
			$pop1_raw{$s}.= "$hoh{$s}{$x}\t"; #creates a raw data file for pop1 wit snp rows and well number columns

		}
		}		       


#
##############################################################################################################################
#  Get all SNP data for a particular cultivar

##########################################################################################################################
#	Print outfile with all genotype data for a particular cultivar
#my $date = `date`;
#$date =~ s/\s//g;
#open( RAW, ">>$ARGV[2]\_raw\_$date\_.txt" )   || die "Can't open outfile $!\n";   # ARGV[2] is popname
 open( RAW, ">>avo_germplasm_384.txt" )   || die "Can't open outfile $!\n" ;   #    
    
   print RAW "$out_header\n";
foreach my $x ( sort keys %pop1_raw ) {
      print RAW "$x\t$pop1_raw{$x}\n";
}
close RAW;

######################################################################################################################
#

#######################################################################################################################
#

#########################################################################################################################
#

##########################################################################################################################
#

