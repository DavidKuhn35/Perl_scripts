  use Bio::SearchIO;
   # format can be 'fasta', 'blast', 'exonerate', ...
  my $hsp = ();
  my $result = ();
  my $hit = ();
  my $searchio = ();
 my $file = $ARGV[0];
   my $searchio = Bio::SearchIO->new( -format => 'blast',
                                     -file   => $file );
   while ( my $result = $searchio->next_result() ) {
       while( my $hit = $result->next_hit ) {
        # process the Bio::Search::Hit::HitI object
           while( my $hsp = $hit->next_hsp ) { 
            # process the Bio::Search::HSP::HSPI object
		print $hsp->algorithm, ": identity ", 100*$hsp->frac_identical, "\%, rank ", $hsp->rank, " (E:", $hsp->evalue, ")\n";
 	    printf("%7s: %s\n", "subj", $hsp->query_string);
 	    printf("%7s: %s\n", "", $hsp->homology_string);
  	    printf("%7s: %s\n", "hom", $hsp->hit_string);
 	    print "\n";
           }
       }
   }
