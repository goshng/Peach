# print $ARGV[0], "\n";

$ARGV[0] =~ /(IL.*)\.csv/;
my $ilname = $1;
$ilname =~ /IL(\d+)-/;
my $chrN = $1;

open IN, "$ARGV[0]";
my $line = <IN>;
my @a = split /\s+/, $line;
my $end = $a[0] + 10;
print "sl$chrN $a[0] $end $ilname\n";
close (IN);
