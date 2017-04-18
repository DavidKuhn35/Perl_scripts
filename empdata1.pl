#!/usr/bin/perl -w

# empdata1.pl: a very simple employee database application.
#  This program reads a text file that is an employee
#  database and lets you query it by ID number.
#
#  Each record is on a single line. Each field in the
#  record is separated by a single tab character
#  ("\t"). The database has four fields:
#  - Last name
#  - First name
#  - ID number
#  - Telephone extension

# The name of the database file
$TheDB = 'edata.txt';

# Open the database file but quit if it doesn't exist
open(INDB, $TheDB) or die "The database $TheDB could " .
  "not be found.\n";

while(1) { # Loop forever
  print "\nDo you want to search by employee ID (I), " .
    " or quit (Q): ";
  $DoSearch = <STDIN>;
  chomp($DoSearch);
  $DoSearch =~ tr/A-Z/a-z/;
  # Check if they want to quit
  if($DoSearch eq 'q') { last }
  # Check if they did *not* say i or I
  unless($DoSearch eq 'i') {
    print "You must enter either I or Q.\n";
    next; # Go out to the while loop
  }

  # Ask them what ID they want to search for
  print "Search for ID number: ";
  $SearchFor = <STDIN>;
  chomp($SearchFor);
  # Go to the top of the database in case this isn't
  #  the first time they are searching
  seek(INDB, 0, 0);
  # Reset the count of successes
  $SuccessCount = 0;
  # Loop through the records in the file
  while(<INDB>) {
    $TheRec = $_;
    chomp($TheRec);
    ($LastName, $FirstName, $ID, $Tel) =
      split(/\t/, $TheRec);
    if($ID eq $SearchFor) {
      $SuccessCount = $SuccessCount + 1;
      print "$ID: $FirstName $LastName, ext. ".
        "$Tel\n";
    } # End of if
  } # End of while(<INDB>)
  if($SuccessCount == 0) { print "No records found.\n" }
  else { print "$SuccessCount records found.\n" }
} # End of while(1)

print "Program finished.\n";
