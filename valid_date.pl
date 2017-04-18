#Program validates the date entered by a user

print "Enter Numeric : month day year \n";

@numofdays = (31,28,31,30,31,30,31,31,30,31,30,31);
$check = 0;
$_ = <STDIN>;
($month,$day,$year)=split;
if($year > 9999 || $year < 0 || $year == NULL ){
	print "Year $year is not a valid year\n";
	$check = 1;
	exit;
}
else{
	if($year%4 == 0){
	#Case of Leap year - we modify the array
	$numofdays[1]=29;
	}
}

if(1 > $month || $month > 12){
	print "There is no month $month. Invalid Month Input \n";	
	$check = 1;
	exit;
}
if($day > $numofdays[$month-1]){
	print "Month $month doesnt have $day days. Invalid Day Input \n";
	$check = 1;
	exit;
}

if($check == 0){
	print "\nValid Date : $month/$day/$year \n\n";
}




