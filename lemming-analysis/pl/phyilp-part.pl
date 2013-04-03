#!/usr/bin/perl
###############################################################################
# Copyright (C) 2013 Sang Chul Choi
#
# This file is part of Lemming Analysis.
#
# Lemming Analysis is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Lemming Analysis is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Lemming Analysis.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

###############################################################################
# COMMAND LINE
###############################################################################
$| = 1; # Do not buffer output
my $VERSION = 'phylip-part.pl 1.0';

my $cmd = "";
sub process {
  my ($a) = @_;
  $cmd = $a;
}

my $man = 0;
my $help = 0;
my %params = ('help' => \$help, 'h' => \$help, 'man' => \$man);
GetOptions( \%params,
            'help|h',
            'man',
            'verbose',
            'version' => sub { print $VERSION."\n"; exit; },
            'start=i',
            'end=i',
            'in=s',
            'out=s',
            '<>' => \&process
            ) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

my $infile;
my $outfile;

if (exists $params{in})
{
  open ($infile, "<", $params{in}) or die "cannot open > $params{in}: $!";
}
else
{
  $outfile = *STDOUT;   
}
if (exists $params{out})
{
  open ($outfile, ">", $params{out}) or die "cannot open > $params{out}: $!";
}
else
{
  $outfile = *STDOUT;   
}

###############################################################################
# DATA PROCESSING
###############################################################################
if ($cmd eq "phylip")
{
  # print $params{in}, "\n";
  # print $params{start}, "\n";
  # print $params{end}, "\n";
 
  my $line = <$infile>;
  my ($numberOfSequences, $numberOfColumns) = split /\s+/, $line;
  # print $numberOfSequences, "\n";
  # print $numberOfColumns, "\n";

  my @sequenceNames;
  my @sequences;
  
  # Read the first section.
  my $numberOfReadSites = 0;
  $line = <$infile>;
  for (my $i = 0; $i < $numberOfSequences; $i++)
  {
    $line = <$infile>;
    chomp $line;
    my $nameOfSequence = substr $line, 0, 10;
    my $sequence = substr $line, 10;
    $sequence =~ s/\s+//g;
    if ($i == 0)
    {
      $numberOfReadSites += length $sequence;
    }
    push @sequenceNames, $nameOfSequence;
    push @sequences, $sequence;
  }

  # Read the rest of sections.
  while ($numberOfReadSites < $numberOfColumns)
  {
    $line = <$infile>;
    for (my $i = 0; $i < $numberOfSequences; $i++)
    {
      $line = <$infile>;
      chomp $line;
      my $sequence = substr $line, 10;
      $sequence =~ s/\s+//g;
      if ($i == 0)
      {
        $numberOfReadSites += length $sequence;
      }
      $sequences[$i] .= $sequence;
    }
  }

  # Extract a part using start and end positions.
  print $outfile "$numberOfSequences $numberOfColumns\n";
  my @extractedSequences;
  foreach (@sequences)
  {
    my $offset = $params{start} - 1;
    my $len = $params{end} - $params{start} + 1;
    push @extractedSequences, substr($_, $offset, $len);
    #print $offset, " ", $len, "\n";
  }

  # Print the part.
  print $outfile "\n";
  for (my $i = 0; $i < $numberOfSequences; $i++)
  {
    print $outfile $sequenceNames[$i];
    print $outfile $extractedSequences[$i];
    print $outfile "\n";
  }
}

if (exists $params{in})
{
  close $infile;
}
if (exists $params{out})
{
  close $outfile;
}
__END__
=head1 NAME

phylip-part - extract a part of a phylip file

=head1 VERSION

phylip-part 1.0

=head1 SYNOPSIS

perl phylip-part.pl phylip -in in.phylip -out out.phylip -start 10 -end 15

perl pl/phyilp-part.pl phylip -in /Volumes/Elements/Documents/Projects/Peach/lemming-analysis/data/Mega_mtGenome_23112011_Fin.txt -start 10 -end 15

=head1 DESCRIPTION

phylip-part will help you to extract something from or extract parts of a PHYLIP
file.

Command:
  phylip - count short reads mapped on genes using shortreadfile and geneposfile.

=head1 OPTIONS

=over 8

=item B<-in> <file>

A file should be an interleaved PHYLIP file.  Standard input is used
unless an input file name is specified.

=item B<-out> <file>

An output file should be an interleaved PHYLIP file.  Standard output is used
unless an output file name is specified.

=item B<-start> <number>

A start position of a part in the PHYLIP alignment file.

=item B<-end> <number>

A end position of a part in the PHYLIP alignment file.

=back

=head1 AUTHOR

Sang Chul Choi, C<< <goshng_at_yahoo_dot_co_dot_kr> >>

=head1 BUGS

If you find a bug please post a message 
so that I can make Lemming Analysis better.

=head1 COPYRIGHT

Copyright (C) 2013  Sang Chul Choi

=head1 LICENSE

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see <http://www.gnu.org/licenses/>.

=cut
