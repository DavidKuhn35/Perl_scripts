#!/usr/bin/perl -w
##################################################################################
# Author		: David Kuhn
# Date			: 4/24/2012
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: opens avo_blastx.db, determines isotig count and unique number of isotigs by reversing and re-reversing the hash.  Good subroutines for keycount, inverthash and
#			hashprint with correct dereferencing.
# <input>		:  access to avo_blastx.db			
# <output>		: file of hash with summary of number of keys   
# 			
# Dependencies	: N/A
# Usage			: use_avo_db.pl (outfile_name_header)
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
my %isotig = ();
my %uniq_iso = ();


tie (%avo_blastx, "DB_File", "avo_blastx.db", O_CREAT|O_RDWR, 0644, $DB_HASH) || die "Cannot create or open avo_blastx.db: $!";	#create a database from the hash from the read in data
	for my $isotig (sort keys %avo_blastx)	{	
		if ($isotig =~ m/header/)	{
	my @header = split (":",$avo_blastx{$isotig});
	}	elsif ($isotig =~ m/^g/)	{
	@snp = split (":", $avo_blastx{$isotig});		
	$isotig{$snp[0]}= $snp[1];
	}
}

	my $isotig_count = &keycount (\%isotig);
	&hashprint ("isotig", \%isotig);

	%uniq_iso = &inverthash (\%isotig);
	my $uniq_iso_count = &keycount (\%uniq_iso);

	&hashprint ($ARGV[0], \%uniq_iso);  
	

#tie (%apvo, "DB_File", "apvo.db", O_CREAT|O_RDWR, 0644, $DB_HASH) || die "Cannot create or open apvo.db: $!";

#open (APVO, "$ARGV[1]") || die "Cannot open apvo infile: $!";
#	while (<APVO>)
#{
#	chomp;
#	if (m/^SUP/)	{
#	@apvo_header = split(/\t/);
#	$apvo{header} = join ":", @apvo_header;
#	}	else	{
#	@apvo_val = split(/\t/);
#	@atg = split (/\|/, $apvo_val[17]);
#	$apvo{$atg[2]} = join ":", @apvo_val;
#	}
#}

sub keycount	{
	my $keycount = keys (%{$_[0]});
}

sub inverthash	{
	my %z = reverse %{$_[0]};
	my %x = reverse %z;
	}

sub hashprint {
	open (HASHOUT, ">>$_[0]\_out.txt") || die "Cannot open hashout: $!";
	print HASHOUT "name of hash is $_[0]\n";
	my $x = keys (%{$_[1]});
	print HASHOUT "hash $_[0] has $x keys\n";
	for my $key (sort keys %{$_[1]})	{
	print HASHOUT "$key\t${$_[1]}{$key}\n";
	}
}


