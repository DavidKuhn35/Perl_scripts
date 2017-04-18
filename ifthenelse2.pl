#!c:\perl\bin


# Multiple comparison criteria
$a = 5;
$b = 6;
$x = "a";
$y = "A";

if (($a > $b) && ($x gt $y))
{
	print "\$a is greater and \$x is greater\n";	
}
elsif (($a > $b) && ($x lt $y))
{
	print "\$a is greater and \$x is lesser\n";
}
elsif (($a < $b) && ($x gt $y))
{
	print "\$a is lesser and \$x is greater\n";
}
else
{	
	print "\$a is lesser and \$x is lesser\n";
}