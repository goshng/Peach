#line 45085 "noweb/imamp.nw"
package Structure2IM;
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(read write);

sub new 
{
  my $type = shift;
  my %parm = @_;
  my $this = {};
  $this->{'filename'} = $parm{'filename'};
  $this->{'comment'} = $parm{'comment'};
  $this->{'npops'} = $parm{'npops'};
  $this->{'first'} = $parm{'first'};
  $this->{'second'} = $parm{'second'};
  $this->{'im'} = $parm{'im'};
  $this->{'poptree'} = $parm{'poptree'};
  bless $this, $type;
}

sub read
{
  my $class = shift @_;
  my $line;
  my $isData;
  my $isDimension;
  my $isLoci;
  my $ii;
  my $li;
  my $nitems;
  my @items;
  my @items2;
  my $i;
  my @inds;
  my @genotypes;
  my @genotype_ind;
  my @numgenes;
  my @numsamps;
  my @numsampsperpop;
  my $ngenes;
  my $pi;

  my %samppop;
  # pop0
  for ($i = 1; $i <= 16; $i++)
    {
      $samppop{$i} = 0;
    }
  # pop1
  for ($i = 17; $i <= 22; $i++)
    {
      $samppop{$i} = 1;
    }
  # pop2
  for ($i = 23; $i <= 64; $i++)
    {
      $samppop{$i} = 2;
    }
  # pop3
  for ($i = 65; $i <= 75; $i++)
    {
      $samppop{$i} = 3;
    }
  # pop4
  for ($i = 76; $i <= 78; $i++)
    {
      $samppop{$i} = 4;
    }
  # pop5
  for ($i = 79; $i <= 84; $i++)
    {
      $samppop{$i} = 5;
    }
  for ($i = 0; $i < 6; $i++) 
    {
      push @numsampsperpop, 0;
    }
  for ($i = 0; $i < 84; $i++)
    {
      push @numsamps, [ @numsampsperpop ];
    } 
  # per locus

  $isData = 0;
  $isDimension = 0;
  $isLoci = 0;
  open STRIN, $class->{'filename'};

  # the only comment 
  $line = <STRIN>;
  @items = split /\t/, $line;
  $class->{'nloci'} = $#items;
  $ii = 0;
  while (<STRIN>)
    {
      $ii++;
    }
  $ii = int ($ii/3);
  die "$ii must be 84" unless $ii == 84;
  $class->{'nind'} = $ii;
  
  close STRIN;
  $class->{'locusnames'} = \@items;

  open STRIN, $class->{'filename'};
  $line = <STRIN>;
  for ($ii = 0; $ii < $class->{'nind'}; $ii++)
    {
      next unless exists $samppop{$ii}; 
      $pi = $samppop{$ii}; 
       
      @genotype_ind = ();
      $line = <STRIN>;
      @items = split /\s+/, $line;
      $line = <STRIN>;
      @items2 = split /\s+/, $line;

      $nitems = scalar (@items);
      if (length $items[0] > 10)
        {
          push @inds, substr ($items[0], -10, 10);
        }
      else
        {
          push @inds, $items[0];
        }
      for ($i = 1; $i < $nitems; $i++)
        {
          if ($items[$i] > 0)
            {
              push @genotype_ind, $items[$i];
              $numsamps[$i-1][$pi]++;
            }
          else
            {
              push @genotype_ind, -1;
            }
          if ($items2[$i] > 0)
            {
              push @genotype_ind, $items2[$i];
              $numsamps[$i-1][$pi]++;
            }
          else
            {
              push @genotype_ind, -1;
            }
        }
      $genotypes[$ii] = [ @genotype_ind ];
      $line = <STRIN>;
    }
  close STRIN;

  for ($li = 0; $li < $class->{'nloci'}; $li++)
    {
      $ngenes = 0;
      for ($ii = 0; $ii < $class->{'nind'}; $ii++)
        {
          if ($genotypes[$ii][2*$li] > 0)
            {
              $ngenes++;
            }
          if ($genotypes[$ii][2*$li + 1] > 0)
            {
              $ngenes++;
            }
        }
      push @numgenes, $ngenes;
    }

  $class->{'individuals'} = \@inds;
  $class->{'genotypes'} = \@genotypes;
  $class->{'numgenes'} = \@numgenes;
  $class->{'numsamps'} = \@numsamps;

  close STRIN;
}

sub write
{
  my $class = shift @_;
  my $line;
  my $i;
  my $li;
  my $ii;
  my @individuals = @{$class->{'individuals'}};
  my @genotypes = @{$class->{'genotypes'}};
  my @numgenes = @{$class->{'numgenes'}};
  my @numsamps = @{$class->{'numsamps'}};

  print $class->{'comment'}, "\n";
  print "#\n";
  print $class->{'npops'}, "\n";
  for ($i = 0; $i < $class->{'npops'}; $i++)
    {
      print "pop$i ";
    }
  print "\n";
  print $class->{'poptree'}, "\n";
  print $class->{'nloci'}, "\n";

  for ($li = 0; $li < $class->{'nloci'}; $li++)
    {
      print "locus$li ";
      if ($class->{'im'} eq 'yes')
        {
          for ($i = 0; $i < $class->{'npops'}; $i++)
            {
              print $numsamps[$li][$i];
              print " ";
            }
        }
      else
        {
          for ($i = 0; $i < $class->{'npops'}; $i++)
            {
              print "0 ";
            }
        }
      print "1 S1 ";
      if ($class->{'im'} eq 'yes')
        {
        }
      else
        {
          print "A";
          print $numgenes[$li], "\n";
        }
      for ($ii = 0; $ii < $class->{'nind'}; $ii++)
        {
          if ($genotypes[$ii][2*$li] > 0)
            {
              printf ("%-10s", $individuals[$ii]);
              print $genotypes[$ii][2*$li], "\n";
            }
          if ($genotypes[$ii][2*$li + 1] > 0)
            {
              printf ("%-10s", $individuals[$ii]);
              print $genotypes[$ii][2*$li + 1], "\n";
            }
        }
    }
}

