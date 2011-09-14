sub fasta_getNumberSequences ($) {
  my ($f) = @_;
  my $c = 0;
  open FILE, $f or die "Could not open $f $!";
  while (<FILE>)
  {
    if (/^>/)
    {
      $c++;
    }
  }
  close FILE;
  return $c;
}

sub fasta_getIndividual ($) {
  my ($f) = @_;
  my @name; 
  open FILE, $f or die "Could not open $f $!";
  while (<FILE>)
  {
    if (/^>.+organism=([\w\s]+).+specimen_voucher=(\w+)/)
    {
      push @name, "$1:$2";
    }
  }
  close FILE;
  return @name;
}

sub fasta_getPopulation ($) {
  my ($f) = @_;
  my @name; 
  open FILE, $f or die "Could not open $f $!";
  while (<FILE>)
  {
    if (/^>.+organism=([\w\s]+)/)
    {
      push @name, $1;
    }
  }
  close FILE;
  return @name;
}

# Produces a file in IMa2 input format.
# $f is a FASTA formated file that contains DNA sequences.
# $pops is a reference of an array of population names.
sub fasta_getPolymorphicSites ($$$) {
  my ($f,$pops,$assignment) = @_;
  my $population;
  my $individual;
  my %sampleSizePopulation;
  my $previousPopulation = "";
  my $numberPopulationWithIndividual = 0;
  my @individualName; 
  my @sequence;
  my %individualPopulation; # population of an indiviudal

  my $useTheSequence = 0;
  open FILE, $f or die "Could not open $f $!";
  while (<FILE>)
  {
    # organism is the name of a population to which a sequence belongs.
    # specimen_voucher is the name of an individual from which the sequence was
    # sampled.
    if (/^>.+organism=([\w\s]+).+specimen_voucher=(\w+)/)
    {
      $population = $1;
      $individual = $2;

      # Indivduals that do not belong to one of populations in the input
      # argument are excluded. We simply skip the sequence as well by using
      # $useTheSequence.
      $useTheSequence = 0;
      for my $p (@$pops)
      {
        if ($p eq $population)
        {
          $useTheSequence = 1;
        }
      }
      next if $useTheSequence == 0;
      
      # Individual name is the name of the sequence.
      # We keep track of the population. The format of IMa2 input file requires
      # that sequences are grouped by populations in order. This may be the case
      # for a raw sequence data. I check if the number of populations and the
      # number of changes population groups are equal.
      $individual =~ /(\d+)/;
      $individual = $1;
      push @individualName, $individual;
      if (exists $individualPopulation{$individual})
      {
        die "Duplicate individual $individual in $f";
      }
      $individualPopulation{$individual} = $population;
      unless ($previousPopulation eq $population)
      {
        $numberPopulationWithIndividual++; 
      }
      if (exists $sampleSizePopulation{$population})
      {
        $sampleSizePopulation{$population}++; 
      }
      else
      {
        $sampleSizePopulation{$population} = 1;
      }
      $previousPopulation = $population;
    }
    elsif ($useTheSequence == 1)
    {
      # If the sequence is a part of the data, i.e., I use the sequence, then I
      # split the sequence by characters and concatenate the characters to
      # sequences. I then use the cacatenated characters to extract columns to
      # decide whether a site is polymorphic. 
      chomp;
      my @e = split //;
      push @sequence, @e;
    }
  }
  close FILE;

  
  my $numberSequence = scalar @individualName;
  my $numberCharacter = scalar @sequence;
  my $numberSite = $numberCharacter/$numberSequence;

=cut
  ################################################
  # For DEBUG
  for (my $i = 0; $i < $numberSite; $i++)
  {
    for (my $j = 0; $j < $numberSequence; $j++)
    {
      print $sequence[$i + $j * $numberSite]; 
    }
    print "\n";
  }
  print "\n";
  #
  ################################################
=cut

  die "$numberCharacter is not divisible by $numberSequence" unless $numberCharacter % $numberSequence == 0;
  my @polymorphicColumn;
  my @numberAllele;
  for (my $i = 0; $i < $numberSite; $i++)
  {
    my @column;
    for (my $j = 0; $j < $numberSequence; $j++)
    {
      push @column, $sequence[$i + $j * $numberSite]; 
    }
    my %seen = (); my @uniqu = grep { ! $seen{$_} ++ } @column;

    # Remove gaps
    my $containGap = 0;
    for my $c (@uniqu)
    {
      if ($c eq '-')
      {
        $containGap = 1;
      }
    }

    # If N and another character are only characters,
    # we remove the column.
    my $nDimorphicSite = 0;
    if (scalar (@uniqu) == 2)
    {
      for my $c (@uniqu)
      {
        if ($c eq 'N' or $c eq 'n')
        {
          $nDimorphicSite = 1;
        }
      }
    }

    # Remove if columns contain characters not
    # A, C, G, T.
    my $notACGT = 0;
    for my $c (@uniqu)
    {
      if ($c eq 'W' 
          or $c eq 'Y' 
          or $c eq 'S' 
          or $c eq 'M' 
          or $c eq 'K' 
          or $c eq 'R')
      {
        $notACGT = 1;
      }
    }
   
    if ($notACGT == 0
        and $nDimorphicSite == 0
        and $containGap == 0 
        and scalar (@uniqu) > 1)
    {
      # This is a polymorphic site.
      push @polymorphicColumn, [ @column ];
      push @numberAllele, scalar (@uniqu);
    }
  }

  my $numberPopulation = scalar (keys %sampleSizePopulation);
  die "Individuals are not consecutive in order of population" unless $numberPopulationWithIndividual == $numberPopulation;

  # Remove individuals with 'N'
  my %exculdeIndividual;
  my $numberPolymorphicSite = scalar (@polymorphicColumn);
  for (my $i = 0; $i < $numberSequence; $i++)
  {
  #  printf "%-15s", $individualName[$i]; 
    my $withN = 0; 
    for (my $j = 0; $j < $numberPolymorphicSite; $j++)
    {
      if ($polymorphicColumn[$j][$i] eq 'N')
      {
        $withN = 1;
      }
    }
    if ($withN == 1)
    {
      $exculdeIndividual{$i} = 0;
    }
  }
  # Remove sites with no segregating sites.
  @numberAllele = ();
  for (my $j = 0; $j < $numberPolymorphicSite; $j++)
  {
    my @s;
    for (my $i = 0; $i < $numberSequence; $i++)
    {
      unless (exists $exculdeIndividual{$i})
      {
        push @s, $polymorphicColumn[$j][$i];
      }
    }
    my %seen = (); my @uniqu = grep { ! $seen{$_} ++ } @s;
    push @numberAllele, scalar(@uniqu);
  }

  # Adjusts the sample size to reflect the removal of some individuals.
  for (my $i = 0; $i < $numberSequence; $i++)
  {
    if (exists $exculdeIndividual{$i})
    {
      my $pop = $individualPopulation{$individualName[$i]}; 
      $sampleSizePopulation{$pop}--;
    }
  }

  # Find the number of filtered polymorphic sites. 
  my $lengthPolymophicSite = 0;
  for my $j (@numberAllele)
  {
    if ($j == 2)
    {
      $lengthPolymophicSite++;
    }
  }

  # Do not use loci with 1 site segregating site or more than 12 segregating
  # sites.
  if ($lengthPolymophicSite <= 1
      or $lengthPolymophicSite > 12)
  {
    return 0;
  }

  ###############################################################
  # Starts to print the IMa2 input file.
  my $locusName;
  $f =~ /(\w+)\.fas/;
  $locusName = $1;
  print $locusName;

  # Prints the sample size.
  if ($assignment == 0)
  {
    for my $key (@$pops)
    {
      if (exists $sampleSizePopulation{$key})
      {
        print " $sampleSizePopulation{$key}";
      }
      else
      {
        print " 0";
      }
    }
    print " $lengthPolymophicSite";
    # Mutation model and inheritance scalar
    print " I 1";
    print "\n";
  }
  else
  {
    my $samplesize = 0;
    for my $key (@$pops)
    {
      if (exists $sampleSizePopulation{$key})
      {
        $samplesize += $sampleSizePopulation{$key};
      }
      print " 0";
    }
    print " $lengthPolymophicSite";
    print " I 1 A$samplesize";
    print "\n";
  }

  #print "A ";
  my $i = 0;
  my $countSite = 0;
  for my $a (@numberAllele)
  {
    $i++;
    if ($a == 1)
    {
      #print " X";
    }
    elsif ($a > 2)
    {
      #print " ZZZZZZZZZZZZZZ($a:$i)";
    }
    elsif ($a == 2)
    {
      #print " 2";
      $countSite++;
    }
    else
    {
      die "Impossible $a";
    }
  }
  # printf ("%d", $countSite);
  # print "\n";

  for (my $i = 0; $i < $numberSequence; $i++)
  {
    unless (exists $exculdeIndividual{$i})
    {
      printf "%-10s", $individualName[$i]; 
      for (my $j = 0; $j < $numberPolymorphicSite; $j++)
      {
        if ($numberAllele[$j] == 2)
        {
          print $polymorphicColumn[$j][$i];
        }
      }
      print "\n";
    }
  }
  return 1;
}

sub fasta_getSampleSize ($$) {
  my ($f,$pops) = @_;
  my $population;
  my $individual;
  my %sampleSizePopulation;
  my $previousPopulation = "";
  my $numberPopulationWithIndividual = 0;
  open FILE, $f or die "Could not open $f $!";
  while (<FILE>)
  {
    if (/^>.+organism=([\w\s]+).+specimen_voucher=(\w+)/)
    {
      $population = $1;
      $individual = $2;
      unless ($previousPopulation eq $population)
      {
        $numberPopulationWithIndividual++; 
      }
      if (exists $sampleSizePopulation{$population})
      {
        $sampleSizePopulation{$population}++; 
      }
      else
      {
        $sampleSizePopulation{$population} = 1;
      }
      $previousPopulation = $population;
    }
  }
  close FILE;
  my $numberPopulation = scalar (keys %sampleSizePopulation);
  print "\n";
  print "$f ($numberPopulationWithIndividual == $numberPopulation)\n";
  if ($numberPopulationWithIndividual == $numberPopulation)
  {
    for my $key (@$pops)
    {
      if (exists $sampleSizePopulation{$key})
      {
        print "$key:$sampleSizePopulation{$key}\n";
      }
      else
      {
        print "$key:0\n";
      }
    }
  }
}

1;
