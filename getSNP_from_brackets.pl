#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: July 22, 2014
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: takes SNP out of brackets from avo1_14f_6kSNP_list_for_Illumina.txt and Hass ref allele, puts it with BHST chip and seq genos to calculate Sean Myles table
# <input>		:  	Avo1_14f...., 2bacon_seq_geno_out, (same for HST), 2BBB_b_chip_geno_out (same for HST)
# <output>		:  Sean Myles table data of seq chip geno congruence
#
# Dependencies	: N/A
# Usage			:
#
##################################################################################
use strict;
use warnings;
use diagnostics;

#####################################################################################################################
#	main
open (SNP, "$ARGV[0]") || die "Can't open SNP infile";
	while (<SNP>)
{
	chomp;
	if (m/^SEQID/)  	{	#If line has SNP name for avocado
	my @snp = split (/\t/);
	$snpname = $snp[1];	#snpname is snp name
	if (($snp[2] eq "N")||($snp[3] eq "N"))	{	# if most frequent nucleotide or reference nucleotide is N remove from further filtering
		$allN{$snpname} = $snpname;
	}	elsif (($snp[6] eq "-")||($snp[9] eq "-")||($snp[14] eq "-")|| ($snp[17] eq "-") || ($snp[22] eq "-") || ($snp[25] eq "-") || ($snp[30] eq "-") || ($snp[33] eq "-")) 	{
			$allindel{$snpname} = $snpname; #if any nucleotide is an indel "-" remove the variant from further filtering
			$isotig = $snp[36];	#isotig name
			$hoh{$isotig}{$snpname} = $snpname; #make hash of hash with isotig as first key and snpname as second key
	} else {
		 	 @bacon = ($snp[1], $snp[2], $snp[3], $snp[4], $snp[5], $snp[6], $snp[7], $snp[8], $snp[9], $snp[10], $snp[11]);
open (SNPlist, "avo1_14f_6kSNP_list_for_Illumina.txt") || die "Can't open avo1_14f... SNP list\n";
	while (<SNPlist>)
{
	chomp;
	
{ 
    unless (
        (      ( $seq =~ m/-/ )
            || ( $seq =~ m!\[\w/\w/\w\]! )
            || ( $seq =~ m/N/ )
            || ( $seq =~ m!\[\w/\w/\w/\w\]!)
        )
        )
    {
        if ( ( $seq =~ m!\[A/G\]! ) || ( $seq =~ m!\[G/A\]! ) ) {
            $seq =~ s!\[\w/\w\]!R!;

            #print "The converted sequence is $seq.\n";

            #	print OUTFILE ">$fasta_name\n$seq\n";
        }
        elsif ( ( $seq =~ m!\[C/G\]! ) || ( $seq =~ m!\[G/C\]! ) ) {
            $seq =~ s!\[\w/\w\]!S!;

#	open (OUTFILE, ">>IUPAC_outfile$ARGV[0].fasta") || die "Can't open outfile\n";
#	print OUTFILE ">$fasta_name\n$seq\n";
#	print "The converted sequence is $seq.\n";
        }
        elsif ( ( $seq =~ m!\[T/G\]! ) || ( $seq =~ m!\[G/T\]! ) ) {
            $seq =~ s!\[\w/\w\]!K!;

#	open (OUTFILE, ">>IUPAC_outfile$ARGV[0].fasta") || die "Can't open outfile\n";
#	print OUTFILE ">$fasta_name\n$seq\n";
#	print "The converted sequence is $seq.\n";
        }
        elsif ( ( $seq =~ m!\[A/C\]! ) || ( $seq =~ m!\[C/A\]! ) ) {
            $seq =~ s!\[\w/\w\]!M!;

#	open (OUTFILE, ">>IUPAC_outfile$ARGV[0].fasta") || die "Can't open outfile\n";
#	print OUTFILE ">$fasta_name\n$seq\n";
#	print "The converted sequence is $seq.\n";
        }
        elsif ( ( $seq =~ m!\[T/C\]! ) || ( $seq =~ m!\[C/T\]! ) ) {
            $seq =~ s!\[\w/\w\]!Y!;

#	open (OUTFILE, ">>IUPAC_outfile$ARGV[0].fasta") || die "Can't open outfile\n";
#	print OUTFILE ">$fasta_name\n$seq\n";
#	print "The converted sequence is $seq.\n";
        }
        elsif ( ( $seq =~ m!\[A/T\]! ) || ( $seq =~ m!\[T/A\]! ) ) {
            $seq =~ s!\[\w/\w\]!W!;

            #	print "The converted sequence is $seq.\n";
        }
        $fname{$fasta_name}       = $fasta_name;
        $convert_seq{$fasta_name} = $seq;
    }
}

######################################################################################################################
#
#######################################################################################################################
#make fasta file with 201 mers through filter 6
#
#my $xmerobj = Bio::SeqIO->new( -file => "$ARGV[0]", -format => "fasta" )
#   ;    #get fasta objects
#while ( my $xmer = $xmerobj->next_seq ) {
#   $header = $xmer->display_id;
#  my $xmerseq  = $xmer->seq;
#    my $xmerdesc = $xmer->desc;
#    if ( exists $filter6{$header} ) {
foreach my $x ( sort keys %fname ) {
    my $mer_out = Bio::Seq->new(
        -seq        => $convert_seq{$x},
        -display_id => $x,

        #            -desc       => $xmerdesc
    );
    my $mer_outIO = Bio::SeqIO->new(
        -file   => ">>avo60K_convert_bracket_IUPAC.fasta",
        -format => "fasta"
        )
        ; #creates new fasta file with IUPAC substitutions for Blast against mango CDS
    $mer_outIO->write_seq($mer_out);
}

#######################################################################################################################
#

#########################################################################################################################
#

##########################################################################################################################
#

