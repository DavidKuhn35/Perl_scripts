# Section A
use XML::Parser;

$ItemCount = 0;
@ItemList = ();

# Section B
$parser = new XML::Parser(Style => 'Stream');
$parser->parsefile('theorder.xml');

# Section C
print "The total number of items ordered was $ItemCount\n";
for($i = 1; $i <= $ItemCount; $i+=1) {
  print "Item $i is $ItemList[$i][3]\n",
    "  Quantity: $ItemList[$i][0]\n",
    "  Color: $ItemList[$i][1]\n",
    "  Size: $ItemList[$i][2]\n\n";
}
exit;

# Section D
sub StartTag {
  my $Handler = shift(@_);
#	print "$Handler\n";
  my $Name = shift(@_);
#	print "$Name\n";
  my %Attr = %_;
  # Section D1
  if($Name eq 'order') { } # Nothing to do
  # Section D2
  elsif($Name eq 'purchaser') {
    $PurchName = $Attr{'name'};
	print "$PurchName\n";
    $PurchID = $Attr{'id'};
	print "$PurchID\n";
    if($PurchID eq '')
      { die "There was no purchaser ID\n" }
  }
  # Section D3
  elsif($Name eq 'shipping') {
    $ShipMethod = $Attr{'method'};
	print "$ShipMethod\n";
  }
  # Section D4
  elsif($Name eq 'item') {
    # Add 1 to the count of items ordered
    $ItemCount += 1;
  }
  # Section D5
  elsif($Name eq 'itemid') {
    $ItemList[$ItemCount][0] = $Attr{'quantity'};
    $ItemList[$ItemCount][1] = $Attr{'color'};
    $ItemList[$ItemCount][2] = $Attr{'size'};
  }
  # Section D6
  else { # Here if there was an unknown start tag
    die "Found an unknown start tag:\n$_\n";
  }
}

# Section E
sub Text { $TheText = $_  };

# Section F
sub EndTag {
  my $Handler = shift(@_);
  my $Name = shift(@_);
  if($Name eq 'itemid') {
    $ItemList[$ItemCount][3] = $TheText;
  }}
