#!/usr/bin/perl -w
#
# Fetch_seq script usage:
#
#    - markers.txt must exist in the same directory, and each line should have the following format:
#         <marker name>,LG,QTL
#
#    - markers that don't exist in GenBank or fail to download will be logged in errors.txt
#
#    - for markers with more than one GenBank entry, the last entry will be downloaded.
#


use Bio::DB::GenBank;


open(FILE, "<libraries.txt") or die "file not found";
open(LOGFILE, ">errors.txt") or die "error creating log file";

while (<FILE>)
{ 
    my($library) = $_;
    chomp($library);
#    @marker = split(/,/, $line);
   
    print "\nProcessing marker: " , $library , "\n"; 

    $gb = new Bio::DB::GenBank;

    # get a stream via a query string
    my $query = Bio::DB::Query::GenBank->new
        (-query   =>'Theobroma cacao[Organism] AND ' . $library,
	 -db      => 'est');
    my $seqio = $gb->get_Stream_by_query($query);

    $cnt = 0;
    while( my $seq =  $seqio->next_seq ) 
    {
      print "Accession: ", $seq->accession_number, " length: ", $seq->length,"\n";
      $lastacc = $seq->accession_number;
      $cnt++;
    }

#    if ($cnt > 1) { die "multiple query results, refine search"; }
    if (!$cnt) 
    { 
	print LOGFILE $library . ": no results for marker\n\n"; 
	next;
    }

    $gb = new Bio::DB::GenBank(-retrievaltype => 'tempfile' , 
			       -format => 'Fasta');
    my $seqio = $gb->get_Stream_by_id($lastacc);
    while( my $clone =  $seqio->next_seq ) 
    {

#      print "Cloneid is ", $clone->display_id, " ", 
#             $clone->accession_number, "\n";

    $fname = $library . ".fas";

    print "Writing output to file: " . $fname . "\n";

    open(MYOUTFILE, ">>", $fname);

    print MYOUTFILE ">>" ;
    print MYOUTFILE $clone->seq;

    close(MYOUTFILE);

#      print "sequence:\n" , $clone->seq , "\n";

    sleep(5);

    }
}

close(FILE);
close(LOGFILE);
