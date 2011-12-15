#line 41963 "noweb/imamp.nw"
#!/usr/bin/perl -w

use strict;
use Bio::SeqIO;

#line 42055 "noweb/imamp.nw"
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


#line 41970 "noweb/imamp.nw"
sub getsequence ($$);

my $li;
my $ii;
my $line;
my $strain;
my $allele;
my @alleles = ();
my @individuals = ();
my @allelelens = ();
my @allelenames = ();
my $allelename;

my $nind = 64;
my $nloci = 8;
my $inputfile = "noweb/burgdorferi.txt";
open BURG, $inputfile or die;

open BURG, "noweb/burgdorferi.txt" or die;
$line = <BURG>;
$allelename = substr $line, 60;
@allelenames = split /\s+/, $allelename; 

while ($line = <BURG>)
  {
    $strain = substr $line, 0, 55;
    $allele = substr $line, 60;
    $allele = trim ($allele);
    my @strains = split /,\s+/, $strain;
    my $nstrain = scalar @strains;
    my @allelesperind = split /\s+/, $allele;
    for ($ii = 0; $ii < $nstrain; $ii++)
      {
        my $trimmedname = trim ($strains[$ii]);
        push @individuals, $trimmedname;
        push @alleles, [ @allelesperind ];
      }
  }
close BURG;

for ($li = 0; $li < $nloci; $li++)
  {
    my $infilename = "noweb/".$allelenames[$li].".fa";
    my $str_seq = getsequence ($infilename, 1);
    push @allelelens, length ($str_seq);
  }

open BURG, $inputfile or die;
print "B. burgdorferi\n";
print "# Margos et al., MLST of housekeeping genes captures geographic\n";
print "# population structure and suggests a European origin of\n";
print "# Borrelia burgdorferi, PNAS 2008, 105, 8730-8735.\n";
print "2\n";
print "pop0 pop1\n";
print "(0,1):2\n";
print "$nloci\n";
for ($li = 0; $li < $nloci; $li++)
  {
    printf ("%s 0 0 %d %d H 1.0\n", $allelenames[$li], $nind, $allelelens[$li]);
    for ($ii = 0; $ii < $nind; $ii++)
      {
        my $infilename = "noweb/".$allelenames[$li].".fa";
        my $str_seq = getsequence ($infilename, $alleles[$ii][$li]);
        printf ("%-10s%s\n", $individuals[$ii], $str_seq);
      }
  }

sub getsequence ($$)
{
  my ($infile, $n) = @_;

  my $i;
  my $inseq; 
  my $seqin = Bio::SeqIO->new(-format => 'fasta', -file => $infile);

  die "$n is not positive" unless $n > 0;
  for ($i = 0; $i < $n; $i++)
    {
      $inseq = $seqin->next_seq ();
    }
  $seqin->close();
  return $inseq->seq;
}

