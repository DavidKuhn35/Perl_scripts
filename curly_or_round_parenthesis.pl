#!C:\perl\bin -w

use strict;

@array1 = ("Tom", "John", "Harry");
foreach (@array1)
{
	print "$_\t";
}
print "\n";

@array2 = {"Tom", "John", "Harry","Potter"};
foreach (@array2)
{
	print "$_\t";
}
print "\n";