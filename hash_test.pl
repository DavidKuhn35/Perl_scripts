#!/usr/bin/perl -w
my %hash = ("cgd1" => "lib1", "cgd2" => "lib2", "cgd3" => "lib3");
print "Original value is $hash{cgd1}\n";
$hash{cgd1} = "$hash{cgd1}\tlib4";
print  "New value is $hash{cgd1}\n";
