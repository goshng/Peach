my $a = "abC";

print "Do not use == for string comparison!\n";
if ($a == "abc")
{
  print "$a == abc\n";
}
else
{
  print "Not that $a == abc\n";
}

if ($a eq "abc")
{
  print "$a eq abc\n";
}
else
{
  print "Not that $a eq abc\n";
}

