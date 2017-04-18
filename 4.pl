print STDOUT "Tell me something:";
while($input=<STDIN>){
	print STDOUT "You said, quote:$input endquote\n";
	chop $input;
	print STDOUT "Without the newline: $input endquote\n";
	if($input eq '') {print STDERR "Null input!\n"};
	print STDOUT "Tell me more, or ^z to end:\n";
	}
print STDOUT "That's all!\n";