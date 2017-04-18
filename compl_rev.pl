#Opens the pure sequence file "forward.txt" and reverse compliments it
#Output is the file "compl_rev.txt"
#compl_rev.pl


open(FP,"<forward.txt");
open(FP2,">compl_rev.txt");
while(<FP>) {
chomp;
s/a/x/g;
s/t/a/g;
s/x/t/g;
s/g/y/g;
s/c/g/g;
s/y/c/g;
@temp=split(//,$_);
@temp=reverse(@temp);
print FP2 @temp;
}

close(FP);
close(FP2);
