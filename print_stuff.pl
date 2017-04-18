#!C:\perl\bin


print "Displaying data in DOS window...\n";
print "Creating variables...\n";
$a = 2;
@b = (3, 4, 5, 10);

print "The value of number A is $a\n";
print "The values in array B are @b\n";
for ($i=0; $i<3; $i=$i+1)
{
	print "$b[$i]\t";
}
print "\n\n";

@array = (3, 5, 7, 9);
print "The first element of array is $array[0]\n\n";


print "Performing hundreds of iterations.... please wait!\n";
$i=0;
while($i < 100000)
{
	$i = $i + 1;
	print "$i\r";  #Carriage return
}
print "\nDone!\n";