#!C:\perl\bin

print "I will make a function call now, this will print stars below this line\n";
&printstars;
&printstars;
&printstars;
&printstars;
&printstars;
print "Printed, didn't it?\n";


sub printstars
{
	print "*********************************************************\n";
}
