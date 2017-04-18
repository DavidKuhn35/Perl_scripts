#c:\perl\bin\perl

$a = "Perl ";

print "Enter a number : ";
$b = <STDIN>;
chomp $b;

print "Concatenating \$a to your input \$b...\n";
print "The result is ".$a.$b."\n";