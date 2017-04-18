#!/usr/bin/perl -w
##################################################################################
# Author		: Ping Zheng
# Target		: Main Lab
# Date			: 04-27-10
# Version		: ver. 1.0.1
# History		: ver. 1.0.0, 12-18-09
# Description	: Parse the default blast result,print Top Hit only add the description and Best Evalue 
# <input>		: BLAST-report-file			
# <output>		: Three tab delimited text files 
# Dependencies	: N/A
# Usage			: %>perl parse_blastReportFile_separateOut.pl 
#                 <BLAST-report-file> 
##################################################################################

use strict;
use Bio::SearchIO;

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                     PARAMETERS                          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                   GLOBAL VARIABLES                      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
my $file 			= $ARGV[0];
my $timestamp 		= get_timestamp();
my $hit_result 		= $file."synonymous_TopHit_".$timestamp."_parsed.txt";
my $noHit_result 	= $file."synonymous_noTopHit_".$timestamp."_parsed.txt";
my $summary 		= $file."synonymous_".$timestamp."_TopHit_summary.txt";
my %count_total		= ();
my %count_noHit		= ();
my %count_e_hit		= ();
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                         MAIN                             #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
main();
exit;
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                      SUBROUTINES                         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
sub main {
	
	my $blast_report = new Bio::SearchIO(-format => 'blast',-file => $file );
	my ($db_name,$db_entries)= ('','');
	open(HIT,">$hit_result")		or die "Can't open $hit_result: $!\n";
	open(NOHIT,">$noHit_result")	or die "Can't open $noHit_result: $!\n";
	open(SUMMARY,">$summary")		or die "Can't open $summary: $!\n";	
	print HIT "Query_name\tQuery_length\tSubject_name\tSubject_source\tSubject_description\tSubject_length\tBest_exp\tScore\tPID\tAlign_length\tStrand\tQuery_start\tQuery_stop\tSubject_start\tSubject_stop\tQuery_sequence\tQuery_frame\tHit_Sequence\tMatch_Sequence\n";
	print NOHIT "Query_name\tQuery_length\n";
	while ( my $result = $blast_report->next_result ) {
		$db_name 	= $result->database_name;
		$db_entries = $result->database_entries;
		$count_total{$result->query_name}++;
		if ($result->num_hits > 0){
			print_blast_results($result);
		}else{
			$count_noHit{$result->query_name}++;
			print NOHIT  join("\t",$result->query_name,$result->query_length). "\n";
		}
	}
	my $count_noHit = scalar keys %count_noHit;
	my $total 		= scalar keys %count_total;
	my $count_hit 	= $total-$count_noHit;
	print SUMMARY "Database name    :\t$db_name\n";
   	print SUMMARY "Database entries :\t$db_entries\n";
	print SUMMARY "Total queries    :\t$total\n";
	printf SUMMARY "=================================================\n";
	print SUMMARY "\tQuery_hit\tQuery_noHit\n";
	print SUMMARY "#query  :\t$count_hit\t$count_noHit\n"; 
	printf SUMMARY ("Coverage:\t%.2f%%\t%.2f%%\n",100*$count_hit/$total,100*$count_noHit/$total);
	close(HIT);
	close(NOHIT);
	close(SUMMARY);
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# print_blast_results :print out parsed result            #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
sub print_blast_results{
	my $result = shift;
   	
   	my $best_hit;
   	my @hits = $result->hits;
 	for my $hit ( sort { $a-> significance <=> $b->significance } @hits ) {
   		$best_hit = $hits[0];
	}
	my $best_hsp;
	my @hsps = $best_hit->hsps;
	for my $hsp ( sort { $a-> evalue <=> $b->evalue } @hsps ) {
   		$best_hsp = $hsps[0];
	 }
	if ($best_hsp->percent_identity >= 95) {
		unless ($best_hsp->query_string =~/X/) {
	my $PID = sprintf ("%.2f",$best_hsp->percent_identity);
	my $sub_name = $best_hit->name;
	my @cols = split(/\|/,$sub_name);
	my $db = trim($cols[0]);
	my $acc = trim($cols[1]);

	if ($best_hsp->strand('hit')< 0){
		my $strand = "Plus/Minus";
	 	print HIT  join("\t",$result->query_name,$result->query_length,$acc,$db,$best_hit->description,$best_hit->length,$best_hsp->evalue,$best_hsp->bits,$PID,$best_hsp->length('total'),$strand,$best_hsp->start('query'),$best_hsp->end('query'),$best_hsp->end('hit'),$best_hsp->start('hit'),$best_hsp->query_string,$best_hsp->query->frame, $best_hsp->hit_string, $best_hsp->homology_string). "\n";
	}else{
		my $strand = "Plus/Plus";
	 	print HIT  
	 join("\t",$result->query_name,$result->query_length,$acc,$db,$best_hit->description,$best_hit->length,$best_hsp->evalue,$best_hsp->bits,$PID,$best_hsp->length('total'),$strand,$best_hsp->start('query'),$best_hsp->end('query'),$best_hsp->start('hit'),$best_hsp->end('hit'),$best_hsp->query_string, $best_hsp->query->frame, $best_hsp->hit_string, $best_hsp->homology_string). "\n";
		}
		}
	}
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# get_timestamp for the file         					  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
sub get_timestamp {

	my @months = qw(01 02 03 04 05 06 07 08 09 10 11 12);
	my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
	my $year = 1900 + $yearOffset;
	if ($dayOfMonth < 10 ){
		$dayOfMonth = "0$dayOfMonth";
	}
	my $thetime = "$year$months[$month]$dayOfMonth";
	return $thetime;
}
#---------------------------------------------------------#
# trim : remove white spaces at front and end
#---------------------------------------------------------#
sub trim($) {
	my $string = shift;
	if (!defined $string){ 
		return '';
	}
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

__END__
