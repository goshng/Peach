my @array1 = qw( a b c d );
my @array2 = qw( c d e f );
@union = @intersection = @difference = ();
%count = ();
foreach $element (@array1, @array2) { $count{$element}++ }
foreach $element (keys %count) {
  push @union, $element;
  push @{ $count{$element} > 1 ? \@intersection : \@difference }, $element;
}

print "array1: ";
print join(" ",@array1);
print "\n";

print "array2: ";
print join(" ",@array2);
print "\n";

print "union: ";
print join(" ",@union);
print "\n";

print "intersection: ";
print join(" ",@intersection);
print "\n";

print "difference: ";
print join(" ",@difference);
print "\n";
