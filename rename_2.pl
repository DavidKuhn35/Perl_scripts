#!c:\perl\bin

$dir_contents = `dir *.ab1` ;
print "These are the directory contents  " . $dir_contents . "\n";

@lines = split(/\n/,$dir_contents);
print "These are the line contents  " . "@lines" . "\n";
@file_list = @lines[5..($#lines-2)];
print "These are the file_list contents  " . " @file_list" . "\n";

foreach (@file_list)
{
	@x = split(/\s/);
	$file = $x[$#x];
	$newfile = $file;
	$newfile =~ s/\.\.\.\.f\.ab1/\.f\.ab1/;
	$dos_command = "rename $file $newfile" ;
	print "$dos_command \n" ;
	`$dos_command`;
}
foreach (@file_list)
{
	print "This is the value of newfile  $newfile";
	@y = split(/_/), $newfile;
	print "This is the value of array y  @y";
	$bufu = $y[$#y];
	$noofile = $bufu;
	$noofile = join "_", $y[0], $y[2], $y[3];
	print "This is the value of noofile  $noofile \n";
	#$noofile = $file;
	$dos_command_2 = "rename $newfile $noofile[$#noofile]";
	print "$dos_command_2 \n";
	`$dos_command_2`;
}