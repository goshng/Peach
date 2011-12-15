package S2M;
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(read write);


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
  my $numberIndividual;

  open IMAIN, $class->{'inputfile'};

  # the first line is a list of locus
  $line = <IMAIN>;
  chomp $line;
  my @loci = split /\s+/, $line; 
  $class->{'nloci'} = scalar @loci;

  $numberIndividual = 0;
  while ($line = <IMAIN>)
  {
    # Diploid Individuals are assumed.
    chomp $line;
    my @genotype1 = split /\s+/, $line; 
    $line = <IMAIN>;
    chomp $line;
    my @genotype2 = split /\s+/, $line; 
    #print @genotype1, "\n";

    $numberIndividual++; 
    push @individuals, $genotype1[0];


    my @individualGenotype;
    for ($j = 0; $j < $class->{'nloci'}; $j++)
    {
      if ($genotype1[5+$j] eq "-9")
      {
        $genotype1[5+$j] = "?";
      }
      if ($genotype2[5+$j] eq "-9")
      {
        $genotype2[5+$j] = "?";
      }
      my $g = sprintf ("(%s,%s)", $genotype1[5+$j], $genotype2[5+$j]);
      push @individualGenotype, $g;
    }
    push @genotypes, \@individualGenotype;
  }
  close IMAIN;

  $class->{'individuals'} = \@individuals;
  $class->{'nind'} = scalar @individuals;
  $class->{'comment'} = "comment";
  $class->{'genotypes'} = \@genotypes;
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
    print "   ";
    print $individuals[$i];
    for ($j = 0; $j < $class->{'nloci'}; $j++)
    {
      print " ";
      print $genotypes[$i][$j];
    }
    if ($i == $class->{'nind'} - 1)
    {
      print "\n";
    }
    else
    {
      print ",\n";
    }
  }

  print "   ;\n";
  print "end;\n";

#  print "\n";
#  print "\nbegin structurama;\n";
#  print "  model numpops=2;\n";
#  print "  mcmc printfreq=100000;\n";
#  print "end;\n";

}

1;

