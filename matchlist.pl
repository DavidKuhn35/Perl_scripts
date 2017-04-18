#!c:\perl\bin -w
#making lists of singlecopy matches

#open SINGLECOPY_LIST, ARGV[1];
#open MATCH_LIST, >>matchlist.out;
#select MATCH_LIST;

while (<>) {
chomp;
@parse = split /\t/, $_;
print "$parse[0] \t $parse[2] \n";
$match{"$parse[0]} = "$parse[2]";
}
while ( ($key, $value) = each %match ) {
	print "$key => $value\n;
}