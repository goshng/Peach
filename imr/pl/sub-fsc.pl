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

sub imr_header {
  my ($f) = @_;
my $header = <<HEADER;
<?xml version=\"1.0\" ?>
<ImrOutput>
    <Seed Value=\"12423467\" />
    <GenerationTime year=\"25\" />
    <Chain Length=\"100000\" Burnin=\"50000\" Thin=\"100\" />
    <MC3 Number=\"25\" A=\"0.995\" B=\"0.5\" />
    <Populations Number=\"3\">
        <Population Name=\"East\" ID=\"1\" />
        <Population Name=\"West\" ID=\"2\" />
        <Population Name=\"Central\" ID=\"3\" />
    </Populations>
    <Individuals Number=\"2\">
        <Individual Name=\"Alice\" ID=\"1\" Population=\"1\" />
        <Individual Name=\"Bob\" ID=\"2\" Population=\"2\" />
    </Individuals>
    <DemographicModel Model=\"Tree\" Tree=\"((1,2)4,3)5\">
        <PriorSplitTime>((1:1,2:1)4:2,3:2)5</PriorSplitTime>
        <PriorPopulationSize>((1:5,2:5)4:5,3:5)5:2</PriorPopulationSize>
        <PriorMigration>0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0</PriorMigration>
    </DemographicModel>
HEADER
  print $f $header;
}

sub imr_footer {
  my ($f) = @_;
my $footer = <<FOOTER;
</ImrOutput>
FOOTER
  print $f $footer;

}

sub fsc_convert2imr {
  my ($in, $out, $length) = @_;
  my $l;
  open OUT, ">$out" or die "Could not open $out $!";
  open IN, $in or die die "Could not open $in $!"; 

  # Read lines until [Data].
  while ($l = <IN>)
  {
    last if $l =~ /^\[Data\]/;
  }

  # Read lines until "polymorphic positions"
  while ($l = <IN>)
  {
    last if $l =~ /polymorphic positions on chromosome/;
  }
  $l = <IN>;
  $l =~ s/#//g;
  $l =~ s/,//g;
  my @position = split /\s+/, $l;

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

  imr_header (*OUT);
  print OUT "    <Data NumberLoci=\"1\">\n";
  print OUT "        <Locus Name=\"locus1\" NumberSequence=\"2\" NumberSite=\"6\">\n";
  print OUT "            <Position Length=\"$length\">";
  print OUT join (" ", @position);
  print OUT "</Position>\n";

  for my $i (1..$SampleSize)
  {
    $l = <IN>;
    chomp $l;
    my @e = split /\t/, $l;
    my $sequence = $e[2];
    # print OUT "$e[0]\t$sequence\n";
    $e[0] =~ /\d+\_(\d+)/;
    my $id = $1;
    print OUT "            <Sequence IndividualID=\"$id\">$sequence</Sequence>\n";
  }

  print OUT "        </Locus>\n";
  print OUT "    </Data>\n";
  imr_footer (*OUT);

  close IN;
  close OUT;
}

sub fsc_convert2imrWithTree {
  my ($in, $inTree, $out, $length, $repetition) = @_;
  my $l;
  open OUT, ">$out" or die "Could not open $out $!";
  open IN, $in or die die "Could not open $in $!"; 

  # Read lines until [Data].
  while ($l = <IN>)
  {
    last if $l =~ /^\[Data\]/;
  }

  # Read lines until "polymorphic positions"
  while ($l = <IN>)
  {
    last if $l =~ /polymorphic positions on chromosome/;
  }
  $l = <IN>;
  $l =~ s/#//g;
  $l =~ s/,//g;
  my @position = split /\s+/, $l;

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

  imr_header (*OUT);
  print OUT "    <Data NumberLoci=\"1\">\n";
  print OUT "        <Locus Name=\"locus1\" NumberSequence=\"2\" NumberSite=\"6\">\n";
  print OUT "            <Position Length=\"$length\">";
  print OUT join (" ", @position);
  print OUT "</Position>\n";

  for my $i (1..$SampleSize)
  {
    $l = <IN>;
    chomp $l;
    my @e = split /\t/, $l;
    my $sequence = $e[2];
    # print OUT "$e[0]\t$sequence\n";
    $e[0] =~ /\d+\_(\d+)/;
    my $id = $1;
    print OUT "            <Sequence IndividualID=\"$id\">$sequence</Sequence>\n";
  }
  print OUT "        </Locus>\n";
  print OUT "    </Data>\n";
  close IN;

  open IN, $inTree or die die "Could not open $in $!"; 
  print OUT "    <RecombinantTrees Number=\"1\">\n";
  print OUT "        <RecombinantTree Length=\"$length\">\n";
  while ($l = <IN>)
  {
    # tree 1_1_pos_0 = [&U] ((1.1:8466, (5.1:6673, 3.1:6673):1793):3301, (2.1:7015, 4.1:7015):4752);
    chomp $l;
    if ($l =~ /\s+tree\s+($repetition)\_(\d+)\_pos\_(\d+) = \[\&U\] (.+)/)
    {
      print OUT "            <LocalTree Position=\"$3\">$4</LocalTree>\n";
    }
  }
  print OUT "        </RecombinantTree>\n";
  print OUT "    </RecombinantTrees>\n";
  imr_footer (*OUT);
  close IN;
  close OUT;
}

1;
