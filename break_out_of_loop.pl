#!c:\perl\bin

#This program illustrates the concept of breaking out of a loop

$i=0;
while(1)	# Using it like this causes infinite loop
{
	$i = $i + 1;
	print "$i\r";
	if($i>1000000)
	{
		last;
	}
}
print "Done!!!\n";