my @array = qw( your array here );
my $search_for = "ere";
my ($index) = grep { $array[$_] eq $search_for } 0..$#array;
print "Array: ";
print join(" ",@array);
print "\n";
if (defined $index)
{
print "here is found at ";
print $index, "\n";
}
