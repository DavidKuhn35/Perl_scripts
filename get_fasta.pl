#!C:blueberry\bin\perl

#use Bio::Perl;
use Bio::Seq; 
use Bio::SeqIO;
 
$seqio_obj = Bio::SeqIO->new(-file => "C:\\Research_Sync\\Genome Sequencing Project\\GEX\\Theobroma_cacao.main_genome.scaffolds.fasta", -format => "fasta" );
while ($seq_obj = $seqio_obj->next_seq){   
    # print the sequence   
    print $seq_obj->seq,"\n";
	$fasta_name = $seq_obj->display_id;
	$seqio_obj = Bio::SeqIO->new(-file => '>$fasta_name.fasta', -format => 'fasta' );
 
$seqio_obj->write_seq($seq_obj);
}
