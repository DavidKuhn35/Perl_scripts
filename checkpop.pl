$Name = 'larry';
$Pass = 'AouerGyvvt';
$Serv = 'mail.yourisp.com';
$Subj = 'Subject: Daily logging report';

use Mail::POP3Client;
$Client = new Mail::POP3Client($Name, $Pass, $Serv);
$TheState = $Client->State;
if($TheState eq 'AUTHORIZATION')
  { die "Bad user name or password.\n" }
elsif($TheState eq 'DEAD')
  { die "Mail server unreachable or unavailable.\n" }
# Find out how many messages there are
$NumMsg = $Client->Count;
#Loop through the messages (starting at 1)
for($i = 1; $i<=$NumMsg; $i +=1) {
  $Headers = $Client->Head($i);
  @HeadList = split(/\n/, $Headers);
  foreach $Line (@HeadList) {
    if($Line =~ /^$Subj/) {
      # Found the message; get the body, then delete
      $Body = $Client->Body($i);
      $Client->Delete($i);
      # Process the report and leave
      &ProcessReport($Body);
      last;
    }
  }
}
# Close the connection so the delete happens
$Client->Close;
sub ProcessReport {
  my($Report) = pop(@_);
  # Do something here that is processing the report
  return;
}
