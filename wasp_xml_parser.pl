#This program parses the output of the standalone allele specific primer program asap.py
# Section A
use XML::Parser;
@AoA = ();
@snpname = ();
@PrimerList = ();
$primernum = 0;
$snpcount = 0;
@text = ();

# Section B
$parser = new XML::Parser(Style => 'Stream');
$parser->parsefile($ARGV[0]);	#$ARGV[0] is file to be parsed.  Must have template as opening and closing <>

# Section C
open (OUTFILE, ">>$ARGV[0]\_parse") || die "Can't open outfile";
my $count = keys %SNPwprimer;
print OUTFILE "The total number of SNPs that were analyzed was $snpcount\n  The total number of SNPs where primers could be designed was $count\n";
#for($i = 0; $i < $snpcount; $i+=1) {
#  print  "snp name  " . ($i+1) . "  is $snpname[$i]\t The sequence is $snpseq[$i]\n",
# }
print OUTFILE "The SNPs where primers could be designed are:\n";
for $key (sort keys %SNPwprimer)	{
	print OUTFILE "$key\n";
}
close OUTFILE;
open (XOUT, ">>XX$ARGV[0]\_primer_test") || die "Can't open outfile";
	for $i (0 .. $#AoA) {
		$row = $AoA[$i];
		for $j (0 .. $#$row ) {
		$AoA[$i][$j]= "\t$AoA[$i][$j]";
		}
	}

	for $row (@AoA)	{
	$a = shift (@$row);
	$a =~ s/^\t//;
	unshift (@$row, $a);
	print XOUT "@$row\n";
	}


	open (PRIMEOUT, ">>YY$ARGV[0]\_primer_order") || die "Can't open outfile";
	for my $i (0 .. $#AoA) {
	$row = $AoA[$i];
		$x = $#$row;	#gets last cell in array
		my $primer_sets = ($x-3)/26;	#determines how many sets of primers were designed
		my @snp = ($AoA[$i][0], $AoA[$i][1]);	#gives SNP name and length of original sequence
		print PRIMEOUT "@snp\n";
		for (my $k = 0; $k <= $primer_sets; $k++)	{ 
				for (my $q = 0; $q <=2 ; $q++ )	{
					my @oligo1 = ($AoA[$i][3+($k*26)], $AoA[$i][4+($k*26)]);	#gets direction and size
					my @mfp = ($AoA[$i][5+($q*8)+ $k*26], $AoA[$i][6+($q*8)+ $k*26],$AoA[$i][7+($q*8)+ $k*26],$AoA[$i][8+($q*8)+ $k*26],$AoA[$i][9+($q*8)+ $k*26],$AoA[$i][10+($q*8)+ $k*26],$AoA[$i][11+($q*8)+ $k*26],$AoA[$i][12+($q*8)+ $k*26]);	#for a particular oligo, print out sets of primers on separate lines to make ordering easier
					print PRIMEOUT "@oligo1@mfp\n";
		}	
		}
}

exit;

# Section D
sub StartTag {
  my $Handler = shift(@_);
  my $Name = shift(@_);
  my %Attr = %_;
	if ($Name eq 'template') {	#Do nothing
  	} elsif($Name eq 'snp') {
	$snplabel = $Attr{'id'};
	push @{$AoA[$snpcount]} , $Attr{'id'};	#SNP id is [0]
 	} elsif($Name eq 'sequence') {
#	my $fixsnp = ($snpcount-1);
#	push @{$AoA[$fixsnp]}, $sequence;
	$seqlength = $Attr{'length'};
	$seqlength =~ s/^\n//;
	push @{$AoA[$snpcount]}, $seqlength;	#seqlength is [1]
	}elsif	($Name eq 'properties') {  #Do nothing
  	}elsif ($Name eq 'item')	{  #Do nothing
  	}elsif ($Name eq 'oligo')	{
	$SNPwprimer{$snplabel} = {$snplabel};
	$direction = $Attr{'direction'};
	$size = $Attr{'size'};
	push @{$AoA[$snpcount]}, $direction;	#direction is [3]
	push @{$AoA[$snpcount]}, $size;	#amplicon size is [4]
	} elsif ($Name eq 'primer')	{
#	@text = split /,/, $fix;
#	$y = pop (@text);
#	$y =~ s/\n$//;
#	push @text, $y;
#	push @{$AoA[$snpcount]}, @text;	
	$name = $Attr{'name'};
	$name =~ s/^\n//;
	push @{$AoA[$snpcount]} , $name;	#primer name is [5]
	push @{$AoA[$snpcount]} , $Attr{'position'};	#position is [6]
	push @{$AoA[$snpcount]} , $Attr{'length'};	#length is [7]
	push @{$AoA[$snpcount]} , $Attr{'seq'};		#sequence is [8]
}
	
 
  # Section D6
  else { # Here if there was an unknown start tag
    die "Found an unknown start tag:\n$_\n";
  }


# Section E
sub Text { $TheText = $_ };

# Section F
sub EndTag {
  my $Handler = shift(@_);
  my $Name = shift(@_);
	if ($Name eq 'sequence')	{
  	$sequence  = $TheText;
	$sequence =~ s/^\n//;
	push@{$AoA[$snpcount]}, $sequence;	#sequence is [2]
	}elsif	($Name eq 'primer')	{
	$fix = $TheText;
	$fix =~ s/^\n//;
	@text = split /,/, $fix;
	$y = pop (@text);
	$y =~ s/\n$//;
	push @text, $y;
	push @{$AoA[$snpcount]}, @text;
	}elsif ($Name eq 'snp')	{	
	$snpcount++;
}}
}
