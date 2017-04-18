#!C:\perl/bin

@array = (1,3,7,5,9);

print "Original array is @array\n";

print "Popping one element...\n";
$popped_item = pop(@array);
print "Array after popping is @array\n";

print "Adding one element...\n";
push (@array, 6);
print "Array after pushing is @array\n";

print "Shifting one element...\n";
shift (@array);
print "Array after shifting is @array\n";

print "Unshifting one element...\n";
unshift (@array,11);
print "Array after unshifting is @array\n";


