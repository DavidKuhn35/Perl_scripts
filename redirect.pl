$WebPerson = 'webmaster@yourcompany.com';

# Look at the string they handed to us in the URL
$TheSearch = $ENV{QUERY_STRING};

# Look in the search string
unless($TheSearch eq '') {
  $RefURL = ''; # The URL we will send them to
  # Open the database
  open(IN, 'AdURLs') or die "Could not open the URL " .
    "database. Please send mail to $WebPerson so " .
    "we can fix this.\n";
  # Read the database line-by-line
  while (<IN>) {
    $TheString = $_; chomp($TheString);
    # Split the line into two parts, based on
    #  whitespace between them
    ($TheName, $TheURL) = split(/\s+/,$TheString, 2);
    # Check if the part we just got is the same as
    #  the part in the search string
    if ($TheName eq $TheSearch) {
      # If so, set some variables and leave
      $RefName = $TheName;
      $RefURL = $TheURL;
      last;
    }
  }
  close(IN);
  unless ($RefURL eq '') {  # We found a match
    &PrintFound;
  }
  else {  # We didn't find a match
    &PrintNotFound;
    # Change the $RefName for the log
    $RefName = "ErrorNotInDatabase:$TheSearch";
  }
}
else { #The search string was empty
  print "HTTP/1.0 400 Bad Request\n";
  print "Content-Type: text/html\n\n";
  print "<head><title>Incorrect Request</title>" .
    "</head><body>\n";
  print "There was no query.\n";
  print "</body>\n";
  $RefName = 'ErrorNoRequest'; # This is for the log
}

# We displayed what we wanted to, now make a log entry

# Get the time right now into a bunch of variables
($Second, $Minute, $Hour, $DayOfMonth, $Month, $Year,
  $WeekDay, $DayOfYear, $IsDST) = localtime(time);
$RealMonth = $Month + 1; # Since $Month starts at 0
# Format the variables to be two digits long
$DateAndTime = sprintf("%02d/%02d/%02d\t%02d:%02d:%02d",
  $RealMonth, $DayOfMonth, $Year, $Hour, $Minute,
  $Second);

# Get the host name or IP address for the log
unless($ENV{REMOTE_HOST} eq '')
  { $TheWhere = $ENV{REMOTE_HOST} }
else { $TheWhere = $ENV{REMOTE_ADDR} }
# Open the log and write the line
open(ADLOG, '>>AdLog') or die "Could not append to the " .
  "URL log. Please send mail to $WebPerson so we can " .
  "fix this.\n";
print ADLOG "$DateAndTime\t$TheWhere\t$RefName\n";
close(ADLOG);

sub PrintFound {
    print <<"EndFound";
HTTP/1.0 302 Found
Location: $RefURL
Content-Type: text/html

<head><title>Redirection Information</title></head><body>
If you are reading this, your Web browser does not support
redirection. That's OK: simply select <a href=\"$RefURL\">
this</a> and you will be taken to the advertiser you
selected.
</body>
EndFound
}
sub PrintNotFound {
    print <<"EndNotFound";
HTTP/1.0 404 Not Found
Content-Type: text/html

<head><title>Incorrect Request</title></head><body>
There is no entry in the URL database for '$TheSearch'.
This may indicate an error in our database of URLs.
If you are sure you got this message by simply choosing
one of the links on our page, please send email to
<a href=\"mailto:$WebPerson\">$WebPerson</a>.
</body>
EndNotFound
}
