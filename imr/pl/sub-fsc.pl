use strict;

sub fsc_convert2gal {
  my ($in,$out) = @_;
  my $l;
  open OUT, ">$out" or die "Could not open $out $!";
  open IN, $in or die die "Could not open $in $!"; 

  # Read lines until [Data].
  while ($l = <IN>)
  {
    last if $l =~ /^\[Data\]/;
  }

  # Read lines until SampleName.
  my $SampleName;
  my $SampleSize;
  while ($l = <IN>)
  {
    if ($l =~ /SampleName=\"(.+)\"/)
    {
      $SampleName = $1;
      last;
    }
  }

  # Read SampleSize.
  $l = <IN>;
  if ($l =~ /SampleSize=(\d+)/)
  {
    $SampleSize = $1;
  }
  else
  {
    die "Not found SampleSize";
  }

  # Read SampleData
  $l = <IN>;
  die "Not found SampleSize" unless ($l =~ /SampleData=/);

  my @msa;
  my $ColumnSize;
  for my $i (1..$SampleSize)
  {
    $l = <IN>;
    chomp $l;
    my @e = split /\t/, $l;
    my @sequence = split //, $e[2];
    $ColumnSize = scalar @sequence;
    push @msa, @sequence;
  }

  for (my $c = 0; $c < $ColumnSize; $c++)
  {
    my @column;
    for (my $r = 0; $r < $SampleSize; $r++)
    {
      push @column, $msa[$c + $r * $ColumnSize];
    }
    my %seen = (); my @uniqu = grep { ! $seen{$_} ++ } @column;
    
    unless (scalar (@uniqu) == 2)
    {
      die "Not a dimorphic site at $c @column"
    }
    for (my $r = 0; $r < $SampleSize; $r++)
    {
      if ($column[$r] eq $uniqu[0])
      {
        $column[$r] = 0;
      }
      elsif ($column[$r] eq $uniqu[1])
      {
        $column[$r] = 1;
      }
    }

    #for (my $r = 0; $r < $SampleSize; $r++)
    #{
      #print $column[$r];
    #}
    #print "\n";
     
    for (my $r = 0; $r < $SampleSize; $r++)
    {
      $msa[$c + $r * $ColumnSize] = $column[$r];
    }
  }

  print OUT "$in\n";
  for (my $r = 0; $r < $SampleSize; $r++)
  {
    for (my $c = 0; $c < $ColumnSize; $c++)
    {
      print OUT $msa[$c + $r * $ColumnSize];
    }
    print OUT "\n";
  }

  close IN;
  close OUT;
}

1;
