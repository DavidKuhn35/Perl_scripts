#!/usr/bin/perl -w
##################################################################################
# Author		: David Kuhn
# Date			: 4/23/2012
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: Read in datafile, make database.  Commented out portion was for avo_blastx.db.  Then same program was used to create the APVO db (although it looks like it was still
#	 under development as APVO was $ARGV[1] and $ARGV[0] was commented out.
# <input>		:  tab delimited datafile			
# <output>		:   database 
# 			
# Dependencies	: N/A
# Usage			: create_avoblastx_apvo_database.pl (avo_blastx datafile) (apvo datafile)
#                    
##################################################################################
use strict;
use warnings;
use diagnostics;
use DB_File;
use Fcntl;
use Statistics::Descriptive;
my @snp = ();
my %avo_blastx = ();
my %apvo = ();
my @apvo_val = ();
my @atg = ();
my @apvo_header = ();
my @blastx_header = ();

#tie (%avo_blastx, "DB_File", "avo_blastx.db", O_CREAT|O_RDWR, 0644, $DB_HASH) || die "Cannot create or open avo_blastx.db: $!";	#create a database from the hash from the read in data
#Read in the avo blastx file and get rid of rows with missing data or Ns
#open (SNP, "$ARGV[0]") || die "Can't open SNP infile: $!";
#	while (<SNP>)
#{
#	chomp;
#	if (m/^g\d+t/)	{
#	@snp = split (/\t/);
#	$avo_blastx{$snp[0]}= join ":", @snp;
#	}	else	{
#	@blastx_header = split (/\t/);
#	$avo_blastx{header} = join ":", @blastx_header;
#	}
#}

tie (%apvo, "DB_File", "apvo.db", O_CREAT|O_RDWR, 0644, $DB_HASH) || die "Cannot create or open apvo.db: $!";

open (APVO, "$ARGV[1]") || die "Cannot open apvo infile: $!";
	while (<APVO>)
{
	chomp;
	if (m/^SUP/)	{
	@apvo_header = split(/\t/);
	$apvo{header} = join ":", @apvo_header;
	}	else	{
	@apvo_val = split(/\t/);
	@atg = split (/\|/, $apvo_val[17]);
	$apvo{$atg[2]} = join ":", @apvo_val;
	}
}

