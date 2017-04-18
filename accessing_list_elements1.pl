#c:\perl\bin\perl

#This program creates a report of the list of variables on the DOS prompt

$no_of_variables = 5;
$var_names = {'Height', 'Weight', 'Age', 'Score', 'Group'}


print "The variables are:\n";
for $i=1:$num_of_variables
	print "$var_names($i)\n";
end