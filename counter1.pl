#!/usr/bin/perl -w

# counter1.pl: one way to count the characters, words, and
#  lines in a text file.

# The name of the file that will be counted
$TheFile = "sample.txt";

# Open the file but quit if it doesn't exist
open(INFILE, $TheFile) or die "The file $TheFile could " .
  "not be found.\n";

# Initialize the counters
$CharCount = 0;
$WordCount = 0;
$LineCount = 0;

while(<INFILE>) {
  $TheLine = $_; # Save the line's contents
  chomp($TheLine); # Get rid of the line break
  $LineCount = $LineCount + 1; # This is the easy part!
  $LineLen = length($TheLine);
  # We now know how many characters (minus the line
  #  break) are on this line
  $CharCount = $CharCount + $LineLen;

  # Do the tricky stuff for counting words on the line

  # The line is empty, so we're done
  if($TheLine eq "") { next };
  # The line has at least one word on it
  $WordCount = $WordCount + 1;
  # Now loop through each character on this line
  #  to look for words
  $CharPos = 0; # Position we are in the line
  # Check for line end; if not, check for a space
  until($CharPos == $LineLen) {
    if(substr($TheLine, $CharPos, 1) eq " ")
      { $WordCount = $WordCount + 1 }
    $CharPos = $CharPos + 1;
  } # End of until
} # End of while(<INFILE>) loop

# All finished, so print out the results

print "For the file $TheFile:\n";
print "Number of characters $CharCount\n";
print "Number of words    $WordCount\n";
print "Number of lines    $LineCount\n";
