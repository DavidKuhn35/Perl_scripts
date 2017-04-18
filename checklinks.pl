$TheURL = $ARGV[0]; # Get the URL from the command line

# Section A
use LWP::Simple;
$TheMain = get($TheURL); # Get the main page
if($TheMain eq '')
  { die "The starting URL was not valid.\n" }

# Section B
$TheMain =~ tr/A-Z/a-z/; # Lowercase everything
$TheMain =~ s/\"/'/g;   # Double quotes to single
$TheMain =~ s/\n/ /g;   # Get rid of the linefeeds
$TheMain =~ s/\s+/ /g;  # Compress spaces/tabs

# Section C
# First find any "base href='URL'" for relative URLs.
$BasePos = index($TheMain, '<base href');
if($BasePos > -1) {
  $OpenQuote = index($TheMain, "'", $BasePos);
  $BaseStart = $OpenQuote + 1;
  $CloseQuote = index($TheMain, "'", $BaseStart);
  $BaseLen = $CloseQuote - $OpenQuote - 1;
  $URLBase = substr($TheMain, $BaseStart, $BaseLen);
}
else { $URLBase = '' }

# Section D
# Keep a list of the URLs we find
@URLList = ();
while(1) { # Loop forever, rely on "last"
  $TheA = index($TheMain, '<a href'); # The next link
  if($TheA == -1) { last };
  # Find the URL in the quotes
  $OpenQuote = index($TheMain, "'", $TheA);
  $URLStart = $OpenQuote + 1;
  $CloseQuote = index($TheMain, "'", $URLStart);
  $URLLen = $CloseQuote - $OpenQuote - 1;
  $TheURL = substr($TheMain, $URLStart, $URLLen);
  # Put the URL into @URLList
  push(@URLList, $TheURL);
  # Chop off beginning of $TheMain to our position
  $TheMain = substr($TheMain, $CloseQuote);
}

# Section E
@CheckURLs = ();
foreach $TheURL (@URLList) {
  # Determine which URLs are relative. Absolute URLs
  #  start with 'http:' or 'ftp:' or 'mailto:';
  #  others are relative to the base. Relative URLs have
  #  $URLBase added to the beginning.
  unless($TheURL =~ /^(http:|ftp:|mailto:)/)
    { $TheURL = $URLBase . $TheURL }
  # Don't check URLs starting with 'mailto:'
  unless($TheURL =~ /^mailto:/)
    { push(@CheckURLs, $TheURL) }
}

# Section F
# Now check each URL in the list for error messages
foreach $TheURL (@CheckURLs) {
  $TheCheck = get($TheURL);
  if($TheCheck ne "") { print "OK:  $TheURL\n" }
  else { print "Fail: $TheURL\n" }
}
