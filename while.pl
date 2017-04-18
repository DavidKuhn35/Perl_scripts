#!c:\perl\bin

@array = (3,5,7,9);

print "First way: Using for loops\n";
for($i =0; $i<= $#array; $i=$i+1)
{
	print "$array[$i]\t";
}

print "\nSecond way: Using foreach loops\n";
foreach $array_element (@array)
{
	print "$array_element\t";
}

print "\nThird way: Using while loops\n";
$i = 0;
while($i<=$#array)
{
	print "$array[$i]\t";
	$i = $i + 1;
}

print "\nDone!\n";