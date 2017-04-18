#!/usr/bin/perl
##################################################################################
# Author		: David Kuhn
# Date			: March 26, 2014
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: takes fasta of SNP sequences and converts bracket format to IUPAC
# <input>		:  	fasta of SNP sequences with brackets
# <output>		:   fasta of SNP sequences with IUPAC
#
# Dependencies	: N/A
# Usage			:
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
my $seqio_obj   = ();
my $seq_obj     = ();
my $fasta_name  = ();
my %fname       = ();
my %convert_seq = ();
my $seq         = ();
#####################################################################################################################
#	main
use Bio::Seq;
use Bio::SeqIO;
$seqio_obj = Bio::SeqIO->new( -file => $ARGV[0], -format => "fasta" );

#open( OUTFILE, ">>IUPAC_outfile_$ARGV[0]" ) || die "Can't open outfile\n";
while ( $seq_obj = $seqio_obj->next_seq ) {
    $fasta_name = $seq_obj->display_id;

    #	$seq_desc = $seq_obj-> desc;
    $seq = $seq_obj->seq;

    #print "I made it this far.The fasta id is $fasta_name. \n";
    #print "The sequence is $seq.\n";
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

