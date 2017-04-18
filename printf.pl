$Str = 'Dummies';
$Num = 3.141593;
$Int = 4142;
@StrFmts = ('s', '4s', '.4s', '4.4s', '10.4s',
  '10.9s', '-10.9s' );
print "The original string is \"$Str\"\n";
&OutList(@StrFmts, $Str);
@NumFmts = (
  'd', '4d', '.4d', '4.4d', '10.4d', '10.10d',
  '-4.4d', '+4.4d', ' 4.4d', '04.4d',
  'f', '4f', '.4f', '4.4f', '10.4f', '10.10f',
  '-10.4f', '+10.4f', ' 10.4f', '010.4f',
  'e', '4e', '.4e', '4.4e', '10.4e', '10.10e',
  '-10.4e', '+10.4e', ' 10.4e', '010.4e',
  'g', '4g', '.4g', '4.4g', '10.4g', '10.10g',
'-10.4g', '+10.4g', ' 10.4g', '010.4g',
);
print "\nThe original number is $Num\n";
&OutList(@NumFmts, $Num);
@HOFmts = ('x', 'X', '3x', '.3x', '3.3x', '10.3x',
  '10.10x', '#x', 'o', '#o' );
print "\nAnd a few hexa and octal examples on $Int.\n";
&OutList(@HOFmts, $Int);

sub OutList {
  local($TestStr) = pop(@_);
  local(@TestFmts) = @_;
  local(@OutFmts, $Fmt);
  foreach $Fmt (@TestFmts) {
    unshift(@OutFmts, sprintf("%s%$Fmt%s", ("%$Fmt=|",
      $TestStr, "|")));
  }
  until(@OutFmts == 0) {
    printf('%-29s%-29s%s',
    pop(@OutFmts), pop(@OutFmts), "\n");
  }
}
