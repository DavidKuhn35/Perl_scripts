#c:\perl\bin\perl

$a = 1;
print "Enter the value for \$b : ";
$b = <STDIN>;
chomp $b; 	# This removes last newline character, if there is one
$c = $a + $b;

print "The value of variable \$a is $a, \$b is $b, and the sum of \$a and \$b is $c\n";