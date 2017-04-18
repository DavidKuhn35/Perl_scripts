#!/usr/bin/perl
##################################################################################
# Author        : David Kuhn
# Date          :March 6, 2014
# Version       : ver. 1.0
# History       : ver. 1.0
# Description   :  creates a SNP genotype file for each germplasm accession and calculates heterozygosity
# <input>       :  SNPchip data,
# <output>      :  genotype file for a particular cultivar with SNP names
#
# Dependencies  : N/A
# Usage         : joinmap_rawdata.pl
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
my $het_header
    = "Well number\thomo\thomopercent\thet\thetpercent\tnosnp\tnosnppercent\ttotal snps\n";
my $out_geno = ();
my %nosnp    = ();
my %homosnp  = ();
my %hetsnp   = ();

#####################################################################################################################
#   main    Read in SNP chip data.
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
        $hoh2{$dna}{$snp} = "$genosort[0]\t$genosort[1]"
            ;    #hoh2 with dna key and snp key and genotype value}

    }
}
#################################################################################################################################################
#Determine heterozygosity  of each germplasm accession
#
foreach my $key ( sort keys %hoh2 ) {
    foreach my $x ( sort keys %{ $hoh2{$key} } ) {
        ( my $z1, my $z2 ) = split /\t/, $hoh2{$key}{$x};

        #           print "for $key, z1 is $z1, z2 is $z2\n";
        if ( ( $z1 eq "-" ) || ( $z1 eq " " ) ) {
            $hoh2{$key}{$x}  = "XX";
            $nosnp{$key}{$x} = "XX";

        }

        elsif ( $z1 eq $z2 ) {
            $homosnp{$key}{$x} = "$z1$z2";
        }
        else {
            $hetsnp{$key}{$x} = "$z1$z2";

        }
    }
}
######################################################################################################################
# Read in manifest data and make hashes of parents
open( DNADATA, $ARGV[1] ) || die "Can't open avo manifest data $!\n";
while (<DNADATA>) {
    chomp;
    if (m/^\d/) {    #if line has DNA name data  (WG cacao)  (\d avocado)
        my @dnadata = split /\t/;
        $WG{ $dnadata[3] } = $dnadata[6]
            ; #make hash with well numbers as keys and common names as values ([0] and [1] for cacao) ([3] [6] for avocado)
    }
}
#######################################################################################################################
#  Get all members of a population and make hashes of WG names
$popname = $ARGV[2];
foreach my $key ( sort keys %WG ) {
    if ( $WG{$key} =~ /$popname/ )
    {         #put in selector for population name (e.g. PNG_8)
        $pop1{$key} = $key
            ; #makes hash of well numbers that had popname as value in manifest
        $out_header .= "$key\t";
    }
}
######################################################################################################################
#   Read in mappable snps, segregation patterns and parent 1 and parent 2 genotypes at snp.  Create hashes with seg_type and parent1 genotype and parent2 genotype
#open (MAPSNP, $ARGV[2]) || die "Can't open mappable snps file $!\n";
#   while (<MAPSNP>)        {
#       if (m/^SHRSPa/) {   # (cacao contig) (avocado SHRSPa)
#           ($mapsnp, $seg_type, $parent1, $parent2) = split /\t/;  #parent1 x parent2 cross
#           $mapsnp_hash{$mapsnp}= $seg_type;
#           $p1mapsnp{$mapsnp} = $parent1;
#           $p2mapsnp{$mapsnp} = $parent2;
#       }
#   }
#
#########################################################################################################################
#   Go through mappable SNPs for each genotype in a mapping population and make a file with genotype data.  snps are rows, genotypes are columns
#oreach my $key ( sort keys %mapsnp_hash )
#{    #go through mappable SNPs one at a time
#foreach my $key ( sort keys %hoh ) {

foreach my $x ( sort keys %pop1 ) {    # go through all pop1 genotypes
    my $homo    = keys %{ $homosnp{$x} };
    my $het     = keys %{ $hetsnp{$x} };
    my $nosnp   = keys %{ $nosnp{$x} };

    my $total   = $homo + $het + $nosnp;
    unless ($total == 0){
    my $perhomo = $homo / $total;
    my $perhet  = $het / $total;
    my $perno   = $nosnp / $total;
       $pop1_raw{$x}
            = "$x\t$homo\t$perhomo\t$het\t$perhet\t$nosnp\t$perno\t$total\n";

        #creates a raw data file for pop1 wit snp rows and well number columns
    }
}


#}
#
##############################################################################################################################
#  Get all SNP data for a particular cultivar

##########################################################################################################################
#   Print outfile with all genotype data for a particular cultivar
my $date = `date`;
$date =~ s/\s//g;
open( RAW, ">>$ARGV[2]\_het\_raw\_$date\_.txt" )
    || die "Can't open outfile $!\n";    # ARGV[2] is popname
print RAW "$het_header\n";
foreach my $x ( sort keys %pop1_raw ) {
    print RAW "$pop1_raw{$x}\n";
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

