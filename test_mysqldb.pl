use strict;
use DBI;
my $url = "DBI:mysql:database=test;host=localhost";
my $user     = "root";
my $password = "mysqlclass";
DBI->connect( $url, $user, $password ) or die "Failure!\n";
print "Success!\n";