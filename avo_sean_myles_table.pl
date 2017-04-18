#!/usr/bin/perl
##################################################################################
# Author        : David Kuhn
# Date          : July 22, 2014
# Version       : ver. 1.0
# History       : ver. 1.0
# Description   : compares chip and sequence data and calculates data for Sean Myles congruence table
# <input>       : avo1_14f_6kSNP_list_for_Illumina.txt, 2bacon_seq_geno_out (same for HST) 2BBB_b_chip_geno_out (same for HST)
# <output>      :  data for Sean Myles congruence table
#
# Dependencies  : N/A
# Usage         :
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

############################################################################################
#
my %hoh         = ();
my %hoh2        = ();
my %snp_hash    = ();
my %snp         = ();
my %refbase     = ();
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
my $out_geno    = ();
my %s1          = ();
my %s2          = ();
my %c1          = ();
my %c2          = ();
my %allele1     = ();
my %allele2     = ();
my %SHRS        = ();
my %sm          = ();

#####################################################################################
#   main    Read in SNP chip data.
open( SNPDATA, $ARGV[0] ) || die "Can't open avo Illumina data $!\n";
while (<SNPDATA>) {
    chomp;
    if (m/^SHRSPa/)
    {    #if line has SNP name in first column, avocado SHRSPa  cacao contig
        my @data    = split(/\t/);  #split on tab for cacao, comma for avocado
        my $SHRS    = $data[0];     #SHRSPa snp name
        my $seq     = $data[1];     #121mer with SNP in brackets
        my $refbase = $data[6];     #Hass ref base
        if ( $seq =~ m!\[(\w)/(\w)\]! ) {    #get SNP from brackets
            my @snp     = ( $1, $2 );
            my @sortsnp = sort (@snp);       #sort SNP
            $allele1{$SHRS} = $sortsnp[0];   #store allele 1 in hash
            $allele2{$SHRS} = $sortsnp[1];   #store allele 2 in hash
            $refbase{$SHRS} = $refbase;      #store refbase in hash
            $SHRS{$SHRS} = $SHRS;  #make hash of SNP names to sort for alleles
        }    #have SNP and refbase for all PaSHRS names

    }
}

#my $a = keys(%SHRS);
#print "$a\n";
###########################################################################################
# Read in sequence data for one parent
open( SEQDATA, $ARGV[1] ) || die "Can't open bseq data $!\n";
while (<SEQDATA>) {
    chomp;
    if (m/^SHRSPa/) {    #if line has avo SNP name
        my @seq  = split /\t/;
        my $name = $seq[0];
        my $z1   = &code( $seq[1], $name );
        my $z2   = &code( $seq[2], $name ); #alleles are sorted but not codes!
        my @sortcode = sort ( $z1, $z2 );
        $s1{$name} = $sortcode[0];
        $s2{$name} = $sortcode[1];
        open SEQOUT, (">>seqout.log");
        print SEQOUT "$name\t$s1{$name}\t$s2{$name}\n";
    }

}

#my $a = keys %s1;
#print "s1 has $a keys\n";
############################################################################################
# Read in chip data for one parent
open( CHIPDATA, $ARGV[2] ) || die "Can't open bseq data $!\n";
while (<CHIPDATA>) {
    chomp;
    if (m/^SHRSPa/) {    #if line has avo SNP name
        my @chip     = split /\t/;
        my $name     = $chip[0];
        my $x1       = &code( $chip[1], $name );
        my $x2       = &code( $chip[2], $name );
        my @sortchip = sort ( $x1, $x2 );
        $c1{$name} = $sortchip[0];
        $c2{$name} = $sortchip[1];
        open CHIPOUT, (">>chipout.log");
        print CHIPOUT"$name\t$c1{$name}\t$c2{$name}\n";

    }

}

#my $b = keys %c1;
#print "c1 has $b keys\n";
####################################################################################
#sort seq and chip data and calculate sean myles categories as separate hashes
open CAT, ">>cat.log" || die "can't open catlog $!\n";
foreach my $key ( sort keys %SHRS ) {
    if ( ( exists $s1{$key} ) && ( exists $c1{$key} ) ) {
        my $cat = &sean( $s1{$key}, $s2{$key}, $c1{$key}, $c2{$key} );

        print CAT "$key\t$cat\n";
        $sm{$cat}{$key} = $key;

    }

}
open SHRSOUT, ">>$ARGV[3]_category_keys.out" || die "can't open SHRSOUT $!\n";
foreach my $y ( sort keys %sm ) {
    my $smkeys = keys %{ $sm{$y} };
    print "number of keys for $ARGV[3] in  category $y is $smkeys\n";
    foreach my $p ( sort keys %{ $sm{$y} } ) {
        print SHRSOUT "$ARGV[3]\t$y\t$sm{$y}{$p}\n";

    }

}
###########################################################################################
# sub to determine sean myles table)
sub sean {
    my $ret = "bufu";
    my ( $s1, $s2, $c1, $c2 ) = @_;
    if ( ( $s1 == 1 ) && ( $s2 == 1 ) && ( $c1 == 1 ) && ( $c2 == 1 ) ) {
        $ret = 1;
    }
    elsif ( ( $s1 == 1 ) && ( $s2 == 1 ) && ( $c1 == 1 ) && ( $c2 == 2 ) ) {
        $ret = 2;
    }
    elsif ( ( $s1 == 1 ) && ( $s2 == 1 ) && ( $c1 == 2 ) && ( $c2 == 2 ) ) {
        $ret = 3;
    }
    elsif ( ( $s1 == 1 ) && ( $s2 == 2 ) && ( $c1 == 1 ) && ( $c2 == 1 ) ) {
        $ret = 4;
    }
    elsif ( ( $s1 == 1 ) && ( $s2 == 2 ) && ( $c1 == 1 ) && ( $c2 == 2 ) ) {
        $ret = 5;
    }
    elsif ( ( $s1 == 1 ) && ( $s2 == 2 ) && ( $c1 == 2 ) && ( $c2 == 2 ) ) {
        $ret = 6;
    }
    elsif ( ( $s1 == 2 ) && ( $s2 == 2 ) && ( $c1 == 1 ) && ( $c2 == 1 ) ) {
        $ret = 7;
    }
    elsif ( ( $s1 == 2 ) && ( $s2 == 2 ) && ( $c1 == 1 ) && ( $c2 == 2 ) ) {
        $ret = 8;
    }
    elsif ( ( $s1 == 2 ) && ( $s2 == 2 ) && ( $c1 == 2 ) && ( $c2 == 2 ) ) {
        $ret = 9;
    }
    else {
        $ret = 666;
    }

}

##########################################################################################
# subroutine to give ref allele 1 or alt allele 2 coding
sub code {
    my ( $s, $ref ) = @_;
    if ( $s eq $refbase{$ref} ) {
        $s = 1;
    }
    else {
        $s = 2;
    }
}
########################################################################################################
# Printing out summary of sean myles table and list of snps in each cat
open SMOUT, ">>sm_summary_$ARGV[3].txt"
    || die "can't print sm_summary_$ARGV[3].log $! \n";
foreach my $cat ( sort keys %sm ) {
    my $x = keys %{ $sm{$cat} };    #gets number of snps in each category
    print SMOUT "$cat\t$x\t$ARGV[3]\n";
}
