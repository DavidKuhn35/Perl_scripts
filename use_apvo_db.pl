#!/usr/bin/perl -w
##################################################################################
# Author		: David Kuhn
# Date			: 4/24/2012
# Version		: ver. 1.0
# History		: ver. 1.0
# Description	: opens APVO database and creates a hash and database of the 959 sequences found in APVO.  Has some good subroutines.
# <input>		:  access to apvo.db			
# <output>		: apvo959.db   
# 			
# Dependencies	: N/A
# Usage			: use_apvo_db.pl
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
my %hoh = ();
my %apvo959 = ();


#tie (%avo_blastx, "DB_File", "avo_blastx.db", O_CREAT|O_RDWR, 0644, $DB_HASH) || die "Cannot create or open avo_blastx.db: $!";	#create a database from the hash from the read in data
#	for my $isotig (sort keys %avo_blastx)	{	
#		if ($isotig =~ m/header/)	{
#	my @header = split (":",$avo_blastx{$isotig});
#	}	elsif ($isotig =~ m/^g/)	{
#	@snp = split (":", $avo_blastx{$isotig});		
#	$isotig{$snp[0]}= $snp[1];
#	}
#}

#	my $isotig_count = &keycount (\%isotig);
#	&hashprint ("isotig", \%isotig);

#	%uniq_iso = &inverthash (\%isotig);
#	my $uniq_iso_count = &keycount (\%uniq_iso);

#	&hashprint ($ARGV[0], \%uniq_iso);  
	

tie (%apvo, "DB_File", "apvo.db", O_CREAT|O_RDWR, 0644, $DB_HASH) || die "Cannot create or open apvo.db: $!";

#open (APVO, "$ARGV[1]") || die "Cannot open apvo infile: $!";
#	while (<APVO>)
#{
#	chomp;
	for my $x (sort keys %apvo)	{
		if ($x =~ m/header/)	{
		my @header = split (":", $apvo{$x});
		}	elsif	($x =~ m/^AT/)	{
		my @y = split (":", $apvo{$x});
		 $hoh{$x}{array}= join ":",@y;
		 $hoh{$x}{Acount} = $y[9];
		 $hoh{$x}{Pcount} = $y[7];
		 $hoh{$x}{Vcount} = $y[6];
		 $hoh{$x}{Ocount} = $y[5];
		}
}

tie (%apvo959, "DB_File", "apvo959.db", O_CREAT|O_RDWR, 0644, $DB_HASH) || die "Cannot create or open apvo959.db : $!";		
	for my $key (sort keys %hoh)	{
		if (($hoh{$key}{Acount}==1) && ($hoh{$key}{Pcount}==1) && ($hoh{$key}{Vcount}==1) && ($hoh{$key}{Ocount}==1))	{
			$apvo959{$key} = $hoh{$key}{array};
		}
}
#&hashprint ("apvo959", \%apvo959);


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


