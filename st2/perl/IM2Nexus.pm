package IM2Nexus;
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(read write writemathematica);


# Declare the subroutines
sub trim($);
sub ltrim($);
sub rtrim($);

# Create a test string
# my $string = "  \t  Hello world!   ";

# Here is how to output the trimmed text "Hello world!"
# print trim($string)."\n";
# print ltrim($string)."\n";
# print rtrim($string)."\n";

# Perl trim function to remove whitespace from the start and end of the string
sub trim($)
{
        my $string = shift;
        $string =~ s/^\s+//;
        $string =~ s/\s+$//;
        return $string;
}
# Left trim function to remove leading whitespace
sub ltrim($)
{
        my $string = shift;
        $string =~ s/^\s+//;
        return $string;
}
# Right trim function to remove trailing whitespace
sub rtrim($)
{
        my $string = shift;
        $string =~ s/\s+$//;
        return $string;
}

sub new 
{
  my $type = shift;
  my %parm = @_;
  my $this = {};
  $this->{'inputfile'} = $parm{'inputfile'};
  $this->{'isdiploid'} = $parm{'isdiploid'};
  $this->{'isasn'} = $parm{'isasn'};
  bless $this, $type;
}

sub read
{
  my $class = shift @_;
  my $line;
  my %genenames;
  my %genes;
  my $i;
  my $j;
  my @individuals;
  my @genotypes;
  my %avoiddiploid;

  open IMAIN, $class->{'inputfile'};

  # the only comment 
  $line = <IMAIN>;
  chomp $line;
  $class->{'comment'} = $line;

  $line = <IMAIN>;
  while ($line =~ /^#/)
    {
      $line = <IMAIN>;
    }
  # We assume that we are at the start line, or number of populations.

  # number of populations
  $line =~ /^(\d+)/;
  $class->{'npops'} = $1;

  # names of populations
  $line = <IMAIN>;
  my @elements = split /\s+/, $line; 
  $class->{'names'} = \@elements;

  # population tree
  $line = <IMAIN>;
  chomp $line;
  $class->{'tree'} = $line;

  # number of loci
  $line = <IMAIN>;
  $line =~ /^(\d+)/;
  $class->{'nloci'} = $1;

  # iterate the loci
  my $ngenes = 0;
  my $nind = 0;
  for ($i = 0; $i < $class->{'nloci'}; $i++)
    {
      $line = <IMAIN>;

      # Find the number of genes
      my $n;
      if ($class->{'isasn'} == 1)
        {
          $line =~ /\s+A(\d+)/;
          $n = $1;
        }
      else
        {
          my @es = split /\s+/, $line;
          $n = 0;
          for ($j = 0; $j < $class->{'npops'}; $j++)
            {
              $n += $es[$j+1];
            }
        }

      my @items = split /\s+/, $line;

      for ($j = 0; $j < $n; $j++)
        {
          $line = <IMAIN>;
          my @items = split /\s+/, $line; 
          my $genename = substr $line, 0, 10;
          $genename = trim ($genename);
          unless (exists $genenames{$genename})
            {
              $genenames{$genename} = $nind;
              $nind++;
            }
          my $gene = substr $line, 10;
          $gene = trim ($gene);
          unless (exists $genes{$gene})
            {
              $genes{$gene} = $ngenes;
              $ngenes++;
            }
        }
    }
  close IMAIN;

  my @nullgenes = ();
  for ($i = 0; $i < $class->{'nloci'}; $i++)
    {
      push @nullgenes, "?";
    }

  @individuals = ();
  for ($i = 0; $i < $nind; $i++)
    {
      push @individuals, "x";
      push @genotypes, [ @nullgenes ];
    }
  foreach (keys %genenames)
    {
      $individuals[$genenames{$_}] = $_;
    }

  open IMAIN, $class->{'inputfile'};
  $line = <IMAIN>;
  $line = <IMAIN>;
  while ($line =~ /^#/)
    {
      $line = <IMAIN>;
    }
  $line = <IMAIN>;
  $line = <IMAIN>;
  $line = <IMAIN>;

  # iterate the loci
  $ngenes = 0;
  for ($i = 0; $i < $class->{'nloci'}; $i++)
    {
      %avoiddiploid = ();
      $line = <IMAIN>;

      my $n;
      if ($class->{'isasn'} == 1)
        {
          $line =~ /\s+A(\d+)/;
          $n = $1;
        }
      else
        {
          my @es = split /\s+/, $line;
          $n = 0;
          for ($j = 0; $j < $class->{'npops'}; $j++)
            {
              $n += $es[$j+1];
            }
        }

      my @items = split /\s+/, $line; 
      for ($j = 0; $j < $n; $j++)
        {
          $line = <IMAIN>;
          my @items = split /\s+/, $line; 
          my $genename = substr $line, 0, 10;
          $genename = trim ($genename);
          my $gene = substr $line, 10;
          $gene = trim ($gene);
          my $ii =$genenames{$genename}; 
          unless (exists $avoiddiploid{$genename})
            {
              $genotypes[$ii][$i] = $genes{$gene};
            }
          $avoiddiploid{$genename} = $ii;
        }
    }
  close IMAIN;

  $class->{'individuals'} = \@individuals;
  $class->{'genotypes'} = \@genotypes;
  $class->{'nind'} = $nind;
}

sub write
{
  my $class = shift @_;
  my $line;
  my $i;
  my $j;
  my $li;
  my $ii;
  my @individuals = @{$class->{'individuals'}};
  my @genotypes = @{$class->{'genotypes'}};

  print "#NEXUS\n\n";
  print "[".$class->{'comment'}, "]\n\n";
  print "begin data;\n";
  print "   dimensions nind=";
  print $class->{'nind'};
  print " nloci=";
  print $class->{'nloci'};
  print ";\n";
  print "   info\n";
  for ($i = 0; $i < $class->{'nind'}; $i++)
    {
      print "   ".$individuals[$i]." ";
      for ($j = 0; $j < $class->{'nloci'}; $j++)
        {
          #print "(".$genotypes[$i][$j].", ?) ";
          print "(".$genotypes[$i][$j].") ";
        }
      if ($i < $class->{'nind'} - 1)
        {
          print ",";
        }   
      print "\n";
    }
  print "  ;\n";
  print "end;\n";

#  print "\n";
#  print "\nbegin structurama;\n";
#  print "  model numpops=2;\n";
#  print "  mcmc printfreq=100000;\n";
#  print "end;\n";

}

sub writemathematica
{
  my $class = shift @_;
  my $line;
  my $i;
  my $j;
  my $li;
  my $ii;
  my @individuals = @{$class->{'individuals'}};
  my @genotypes = @{$class->{'genotypes'}};

  print "{";
  for ($j = 0; $j < $class->{'nloci'}; $j++)
    {
      print "{";
      for ($i = 0; $i < $class->{'nind'}; $i++)
        {
          print "{".$genotypes[$i][$j]."} ";
          if ($i < $class->{'nind'} - 1)
            {
              print ",";
            }   
        }
      print "}";
      if ($j < $class->{'nloci'} - 1)
        {
          print ",";
        }
    }
  print "}";
}
1;

