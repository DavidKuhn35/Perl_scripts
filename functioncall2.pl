#!C:\perl\bin


@nums = (3,5,6,1,2,3,4,5,7,4,5,8,9,3,2,2,4,5,6,6,7);


print "I will call a function that produces the mean of some numbers...\n";
&printstars;
$avg = &mean(@nums);
print "The means of numbers is $avg\n";
&printstars;





sub printstars
{
	print "**********************************************\n";
}


sub mean
{
	@input_nums = @_;
	$avg1 = 0;
	for($i=0; $i<scalar(@input_nums); $i=$i+1) 
	{
		$avg1 = $avg1 + $input_nums[$i];
	}
	$avg1 = $avg1/scalar(@input_nums);
	return $avg1;
}