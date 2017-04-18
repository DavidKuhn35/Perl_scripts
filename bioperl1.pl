#!C:\perl\bin

use Bio::Perl;
 
# this script will only work if you have an internet connection on the
# computer you're using, the databases you can get sequences from 
# are 'swiss', 'genbank', 'genpept', 'embl', and 'refseq'
 
$seq_object = get_sequence('genbank',"ROA1_HUMAN");
write_sequence(">roa1.fasta",'fasta',$seq_object);

$blast_result = blast_sequence($seq_object);
write_blast(">roa1.blast",$blast_result);
