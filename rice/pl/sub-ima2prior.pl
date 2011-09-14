
sub ima2prior_makeDefault ($)
{
  my ($f) = @_;
  open FILE, '>', $f or die "Could not open $f $!";
  print FILE "# Default\n";
  print "What is the number populations [e.g., 4]? ";
  chomp (my $npop = <>);
  my $nsplit = $npop - 1;
  print "------------------------------\n\n";

  print "What is the species tree [e.g., ((0,1):5,(2,3):4):6]? ";
  my $tree = <>;
  print FILE $tree;
  print "------------------------------\n\n";

  print "What is the species tree with maximum split time?\n";
  print <<PRIOR;
For example, with the population string from ((0,1):5,(2,3):4):6 and supposing
that the upper bound on the most recent splitting event is 0.5 and that it is
5.0 for the older two events, the string would be :
((0:5.0,1:5.0):5:5.0,(2:0.5,3:0.5):4:5.0):6.
PRIOR
  print "So, what is your splitt time prior? ";
  $tree = <>;
  print FILE $tree;
  print "------------------------------\n\n";

  print "What is the the species tree with maximum population size?\n";
  print <<PRIOR;
For example, if populations 0, 1 and 4 have an upper bound of 5.0 and the
remainder have an upper bound of 10.0 the string would be:
((0:5.0,1:5.0):5:10.0,(2:10.0,3:10.0):4:5.0):6:10.0
PRIOR
  print "So, what is your population size prior? ";
  $tree = <>;
  print FILE $tree;
  print "------------------------------\n\n";

  print "What is the maximum value of migration prior? ";
  chomp (my $m = <>);
  print FILE <<PRIOR;
0 $m $m 0 0
$m 0 $m $m 0
$m $m 0 0 0
0 $m 0 0 0
0 0 0 0 0
PRIOR

  # I could have a different prior for one of three species tree.
  # I think that this appears to be more direct.

  close FILE;
}

1;
