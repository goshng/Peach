# print $ARGV[0], "\n";

$ARGV[0] =~ /(IL.*)\.csv/;
my $ilname = $1;
$ilname =~ /IL(\d+)-/;
my $chrN = $1;

open IN, "$ARGV[0]";
my $line = <IN>;
my @a = split /\s+/, $line;
my $start = $a[0] - 1000000;
if ($start < 0)
{
  $start = 1;
}
my $end = $start + 10;
print "sl$chrN $start $end $ilname\n";
close (IN);
