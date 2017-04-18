#!/usr/local/bin/perl

# Linda Ballard, 9 June 2008, USDA-ARS-MSA Genomics Lab
#   Written for David Kuhn
#   reads fasta sequence file
#   returns sequences with linkers masked
#      and headers modified with linker information

$answer = 1;

do {
	print "Mask Linker & Add info to header\n\n";
	print "FASTA file name:> ";
	chomp ($FileName = <STDIN>);

	my ($ct, $ct01, $ct02, $nocount, $flag) = (0,0,0,0,0);
	
	open(DATA2, "$FileName") or die "Couldn't open $FileName $!\n"; 
	$FileName =~ s/\.\w+$//g;
	$OutFile = $FileName."\_maskedlinker.fasta";
	open(OUT2, ">$OutFile") || die "Couldn't open $OutFile $!\n";
	
	$/ = ">"; #adjust record separator for fasta input
	
	while (<DATA2>) {
		($header, $seq) = split (/\n/,$_,2);
		next unless ($seq && $header);
		$seq =~ s/\s//g;
		$seq =~ s/>//g;
		$flag = 0;
	
		# masking & header modification
		if ($seq =~ /ACGTT{4}/ ) { # linker 1   ACGTTTTT
			$seq =~ s/ACGTTTTT/NNNNNNNN/g;  #g means match all occurences
			$header =~ s/$/ L01/;  #put L01 on end of header line
			$ct01++;
			$flag++;		
		}
		if ($seq =~ /TTAAGGCC/ ) { # linker 2
			$seq =~ s/TTAAGGCC/NNNNNNNN/g;
			$header =~ s/$/ L02/;  #put L02 on end of header line
			$ct02++;
			$flag++;
		}
		#insert additional if statements as needed	
				
   		printf OUT2 "$header\n$seq\n";
   		if ($flag == 0) { $nocount++; }
   		$flag = 0;
       	$count++;
  	}
	close (OUT2);

	# what happened in while loop	
	print "$count sequences\n";
	print "$nocount had no linkers\n";
	print "$ct01 linker 1\n";
	print "$ct02 linker 2\n";  
	#insert additional print statements for other linkers 
	
	$/ = "\n";  #reset record separator to default
	
	print "\n Do another (yes = 1)?";
	chomp ($an = <STDIN>);
	if ($an != 1) { $answer = 0; }
		
}while ($answer == 1);

exit;