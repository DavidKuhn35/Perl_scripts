#!c:\perl\bin

@array = (1,2,11,23,5,7,9,99);
# @sorted = sort(@array);
@sorted = sort{$a <=> $b} @array;

print "@array\n";
print "@sorted\n";