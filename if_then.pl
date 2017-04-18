#!c:\perl\bin

@scores = (10, 100, 100, 100.1, 200, 500);
$noise = 100;
print "Number of elements is ".scalar(@scores)."\n";

for($i = 0; $i < scalar(@scores); $i = $i + 1)
{
	print "Original value: $scores[$i]";
	if($scores[$i] < $noise)
	{
		$scores[$i] = $noise;	
	}
	print "\tNew value: $scores[$i]\n";
}