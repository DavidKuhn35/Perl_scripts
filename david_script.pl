	#!c:\perl\bin

	$dir_contents = `dir *.exe`;
	@lines = split(/\n/,$dir_contents);
	@file_list = @lines[5..($#lines-2)];

	foreach (@file_list)
	{
		@x = split(/\s/);
		$file = $x[$#x];	# last element is the filename
		$newfile = $file;	# make duplicate	
		$newfile =~ s/.abl/_pcr.f.abl/;
		$dos_command = "ren $file $newfile";
		print "$dos_command\n";
		`$dos_command`;		# executes the command string
	}