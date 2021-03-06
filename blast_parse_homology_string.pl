 #!/usr/bin/perl -w
 
 use strict;
 use Bio::SeqIO;
 use Bio::Tools::Run::StandAloneBlast;
 my $query1_in  = Bio::SeqIO->newFh ( -file   => "mus-betaglobin-bh0.fas",
  				    -format => 'fasta' );
 my $query1 = <$query1_in>; 
 my $query2_in  = Bio::SeqIO->newFh ( -file   => "mus-betaglobin-bh3.fas", 
  				      -format => 'fasta' );
 my $query2 = <$query2_in>; 
 
 $factory = Bio::Tools::Run::StandAloneBlast->new('program'  => 'blastn');
 $report = $factory->bl2seq($query1, $query2);
 
 while (my $result = $report->next_result) {
    print "Query: ".$result->query_name."\n";
    while (my $hit = $result->next_hit) {
 	while ($hsp = $hit->next_hsp) {
 	    print $hsp->algorithm, ": identity ", 100*$hsp->frac_identical, "\%, rank ", $hsp->rank, " (E:", $hsp->evalue, ")\n";
 	    printf("%7s: %s\n", "subj", $hsp->query_string);
 	    printf("%7s: %s\n", "", $hsp->homology_string);
  	    printf("%7s: %s\n", "hom", $hsp->hit_string);
 	    print "\n";
 	}
    print "\n";
    }
 }
