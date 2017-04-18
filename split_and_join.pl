#!C:\perl\bin

$old_date_format = "11-7-2006";

($month, $date, $year) = split(/-/,$old_date_format);
$new_date_format = join(":", $month, $date, $year);

print "Old format is $old_date_format\n";
print "New format is $new_date_format\n";
