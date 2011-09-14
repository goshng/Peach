sub ima2data_findNumberLoci ($)
{
  my ($f) = @_;
  my $nloci = -1;
  my @sampleSizes;

  my $line;
  open FILE, $f or die "Could not open $f $!";
  $line = <FILE>;
  while ($line = <FILE>)
  {
    chomp $line;
    last if $line !~ /^#/; 
  }
  my $npop = $line;
  $line = <FILE>; chomp $line;
  my @popName = split /\s+/, $line;
  $line = <FILE>; chomp $line;
  my $tree = $line;
  $line = <FILE>; chomp $line;
  $nloci = $line;

  # Ready to read sequence data for each locus.
  $nloci = 0;
  while ($line = <FILE>)
  {
    chomp $line;
    my @e = split /\s+/, $line;
    my $locusName = $e[0];
    my $sampleSize = 0;
    for my $i (1..$npop)
    {
      $sampleSize += $e[$i];
    }
    for my $i (1..$sampleSize)
    {
      $line = <FILE>; chomp $line;
      my $seqName = substr $line, 0, 10;
      my $seq = substr $line, 10;
      my @e = split //, $seq;
      for my $nucleotide (@e)
      {
        if ($nucleotide ne 'A' 
            and $nucleotide ne 'C' 
            and $nucleotide ne 'G' 
            and $nucleotide ne 'T')
        {
          die "$f $locusName, $i-th sample $seqName has $nucleotide";
        }
      }
    }
    $nloci++;
    push @sampleSizes, $sampleSize;
  }  
  close FILE;
  return ($nloci, \@sampleSizes);
}

1;
