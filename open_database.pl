#!/usr/bin/perl -w
##################################################################################
# Author		: David Kuhn
# Date			: 4/20/2012
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: Opens database, prints outfile of keys and selected data in database
# <input>		:  (database name) (outfile header)			
# <output>		:  (outfile with key and selected values from database)  
# 			
# Dependencies	: N/A
# Usage			: open_database.pl (database name) (
#                    
##################################################################################
use strict;
use warnings;
use diagnostics;
use DB_File;
use Fcntl;
use Statistics::Descriptive;
my @snp = ();
my %hoh = ();
#my $badsnp = 0;
my %isotig = ();
my @annot = ();
my $annot = ();

tie (%hoh, "DB_File", "$ARGV[0]") || die ("Cannot open database\n");	#open a database
tie (%isotig, "DB_File", "avo_isotig.db", O_CREAT, 0644, $DB_HASH) || die ("Cannot create or open avo_isotig.db");	#create database for isotig and annotation data
for my $key(sort keys %hoh)	{
	@snp = split (/:/,$hoh{$key});
	my $annot = join (":", $snp[36], $snp[51], $snp[53], $snp[54], $snp[55]);
	$isotig{$snp[36]} = $annot;
	}

open (OUTFILE, ">>ZZ$ARGV[1]\_db.out") || die ("Cannot open outfile\n");

for my $key(sort keys %isotig)	{
	@annot = split (/:/, $isotig{$key});
	print OUTFILE "$key\t$annot[1]\t$annot[2]\n";
#	print OUTFILE "$key\t$isotig{$key}\n";
	}
	 
