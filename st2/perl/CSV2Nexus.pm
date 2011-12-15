#line 45081 "noweb/imamp.nw"
package CSV2Nexus;
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(read write);
#line 44999 "noweb/imamp.nw"
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


#line 45087 "noweb/imamp.nw"
sub new 
{
  my $type = shift;
  my %parm = @_;
  my $this = {};
  $this->{'filename'} = $parm{'filename'};
  $this->{'comment'} = $parm{'comment'};
  $this->{'npops'} = $parm{'npops'};
  $this->{'popns'} = $parm{'popns'};
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
  my $j;
  my @inds;
  my @genotypes;
  my @genotype_ind;
  my @numgenes;
  my @numsamps;
  my @numsampsperpop;
  my $ngenes;
  my $pi;

  my %samppop;
  $pi = 0;
  $ii = 0;
  # per locus

  $isData = 0;
  $isDimension = 0;
  $isLoci = 0;

  open IN, $class->{'filename'};
  $line = <IN>;
  $line = trim ($line);
  @items = split /\s+/, $line;
  $class->{'nloci'} = $#items + 1;
  $class->{'locusnames'} = \@items;
  $ii = 0;
  while ($line = <IN>)
    {
      @genotype_ind = ();
      @items = split /\s+/, $line;
      $items[0] =~ s/"(\w+)"/$1/; 
      push @inds, $items[0];
      for ($i = 1; $i < $#items + 1; $i++)
        {
          push @genotype_ind, $items[$i];
        }
      push @genotypes, [ @genotype_ind ];
      $ii++;
    }
  close IN;

  $class->{'nind'} = $ii;
  $class->{'individuals'} = \@inds;
  $class->{'genotypes'} = \@genotypes;

  close STRIN;
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
          print "(";
          if ($genotypes[$i][$j*2] > 0)
            {
              print $genotypes[$i][$j*2];
            }
          else
            {
              print " ? ";
            }
          print ",";
          if ($genotypes[$i][$j*2 + 1] > 0)
            {
              print $genotypes[$i][$j*2 + 1];
            }
          else
            {
              print " ? ";
            }
          print ") ";
        }
      if ($i < $class->{'nind'} - 1)
        {
          print ",";
        }   
      print "\n";
    }
  print "  ;\n";
  print "end;\n";

  print "\n";
  print "\nbegin structurama;\n";
  print "  model numpops=2;\n";
  print "  mcmc printfreq=100000;\n";
  print "end;\n";

}


