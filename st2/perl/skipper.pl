#line 79 "noweb/imamp.nw"
#!/usr/bin/perl -w

use strict;
use Bio::SeqIO;

sub skipper_populations ();
sub skipper_nseqs ();
sub skipper_seq ($$);
sub skipper_seqname ($$);
sub skipper_seqlen ();

#line 42918 "noweb/imamp.nw"
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


#line 92 "noweb/imamp.nw"
my @populations = ();
my %nseqs;
my $i;

print "Hebert et al., Ten species in one: DNA barcoding reveals cryptic species in the neotropical skipper butterfly Astraptes fulgerator. PNAS 101:14812-14817. 2004.\n";
print "# Hebert et al. (2004)\n";
print "12\n";
%nseqs= skipper_nseqs ();
foreach (keys %nseqs)
  {
    print $_, " ";
  }
print "\n";
print "0\n";
print "1\n";
print "COI ";
foreach (keys %nseqs)
  {
    print $nseqs{$_}, " ";
  }
my $len = skipper_seqlen ();
print "$len H 1.0\n";

foreach (keys %nseqs)
  {
    my $nseq = $nseqs{$_};
    for ($i = 0; $i < $nseq; $i++)
      {
        my $strseq = skipper_seq ($_, $i); 
        my $name = skipper_seqname ($_, $i); 
        print "$name$strseq\n";
      }
  }


# functions 
#line 131 "noweb/imamp.nw"
sub skipper_populations ()
{
  my %populations = ();
  my $infile = "noweb/skipper440.fa";
  my $inseq = Bio::SeqIO->new(-format => 'fasta', -file => $infile);
  while (my $seq = $inseq->next_seq)
    {
      my @elements = split /[>|\s]/, $seq->desc;
      unless (exists $populations{$elements[0]})
        {
          $populations{$elements[0]} = 1;
        }
    }
  $inseq->close ();
  return keys %populations;
}

sub skipper_nseqs ()
{
  my %populations = ();
  my $infile = "noweb/skipper440.fa";
  my $inseq = Bio::SeqIO->new(-format => 'fasta', -file => $infile);
  while (my $seq = $inseq->next_seq)
    {
      my @elements = split /[>|\s]/, $seq->desc;
      if (exists $populations{$elements[0]})
        {
          $populations{$elements[0]}++;
        }
      else
        {
          $populations{$elements[0]} = 1;
        }
    }
  $inseq->close ();
  return %populations;
}

sub skipper_seq ($$)
{
  my ($popnname, $ei) = @_;

  my $strseq;
  my $i = 0;
  my %populations = ();
  my $infile = "noweb/skipper440.fa";
  my $inseq = Bio::SeqIO->new(-format => 'fasta', -file => $infile);
  while (my $seq = $inseq->next_seq)
    {
      my @elements = split /[>|\s]/, $seq->desc;
      if ($elements[0] eq $popnname)
        {
          if ($i == $ei)
            {
              $strseq = $seq->seq;
              last;
            }
          $i++;
        }
    }
  $inseq->close ();
  return $strseq;
}

sub skipper_seqname ($$)
{
  my ($popnname, $ei) = @_;

  my $name;
  my $i = 0;
  my %populations = ();
  my $infile = "noweb/skipper440.fa";
  my $inseq = Bio::SeqIO->new(-format => 'fasta', -file => $infile);
  while (my $seq = $inseq->next_seq)
    {
      my @elements = split /[>|\s]/, $seq->desc;
      if ($elements[0] eq $popnname)
        {
          if ($i == $ei)
            {
              my @elements = split /[>|\s]/, $seq->display_id;
              $name = $elements[0];
              last;
            }
          $i++;
        }
    }
  $inseq->close ();
  return $name;
}

sub skipper_seqlen ()
{
  my $len;
  my $infile = "noweb/skipper440.fa";
  my $inseq = Bio::SeqIO->new(-format => 'fasta', -file => $infile);
  my $seq = $inseq->next_seq;
  $len = $seq->length;
  while ($seq = $inseq->next_seq)
    {
      die "Sequence lengths are different!" unless $len == $seq->length; 
    }
  $inseq->close ();
  return $len;
}




