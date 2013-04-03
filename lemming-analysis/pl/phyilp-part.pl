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
            'format=s',
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
  unless (exists $params{format})
  {
    $params{format} = $cmd;
  }
  $params{format} = lc $params{format};
  # print $params{in}, "\n";
  # print $params{start}, "\n";
  # print $params{end}, "\n";
 
  my $line = <$infile>;
  my ($numberOfSequences, $numberOfColumns) = split /\s+/, $line;
  unless (exists $params{start})
  {
    $params{start} = 1;
  }
  unless (exists $params{end})
  {
    $params{end} = $numberOfColumns;
  }
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
  my @extractedSequences;
  foreach (@sequences)
  {
    my $offset = $params{start} - 1;
    my $len = $params{end} - $params{start} + 1;
    push @extractedSequences, substr($_, $offset, $len);
  }

  if ($params{format} eq "phylip")
  {
    # Print the part out in PHYLIP format.
    # 
    # Segment 1  Segment 2  Segment 3
    # AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA -> Chunk 1
    # AAAAAAAAAA AAAAAAAAAA AAAA       -> Chunk 2
    # Note that numberOfColumns is now the length of a part.
    $numberOfColumns = $params{end} - $params{start} + 1;
    print $outfile "$numberOfSequences $numberOfColumns\n";
    my $numberOfChunks = int($numberOfColumns / 60);
    if ($numberOfColumns % 60 != 0)
    {
      $numberOfChunks++;
    }

    for (my $i = 0; $i < $numberOfChunks; $i++)
    {
      print $outfile "\n";
      for (my $j = 0; $j < $numberOfSequences; $j++)
      {
        if ($i == 0)
        {
          print $outfile $sequenceNames[$j];
        }
        else
        {
          print $outfile "          ";
        }

        my $chunkOfSequence = substr $extractedSequences[$j], 10*6*$i, 60;
        my $lengthOfChunk = length($chunkOfSequence);
        my $numberOfSegments = int($lengthOfChunk / 10);
        if ($lengthOfChunk % 10 != 0)
        {
          $numberOfSegments++;
        }

        my $chunk = "";
        for (my $k = 0; $k < $numberOfSegments; $k++)
        {
          $chunk .= substr $chunkOfSequence, 10*$k, 10;
          if ($k < $numberOfSegments) 
          {
            $chunk .= " ";
          }
        }
        print $outfile "$chunk\n";
      }
    }
  }
  elsif ($params{format} eq "fasta")
  {
    for (my $i = 0; $i < $numberOfSequences; $i++)
    {
      $sequenceNames[$i] =~ s/\s+//g;
      print $outfile ">$sequenceNames[$i]\n";

      $numberOfColumns = $params{end} - $params{start} + 1;
      my $numberOfChunks = int($numberOfColumns / 60);
      if ($numberOfColumns % 60 != 0)
      {
        $numberOfChunks++;
      }
      for (my $j = 0; $j < $numberOfChunks; $j++)
      {
        my $chunkOfSequence = substr $extractedSequences[$i], 10*6*$j, 60;
        print $outfile "$chunkOfSequence\n";
      }
    }
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

perl pl/phyilp-part.pl phylip -in Mega_mtGenome_23112011_Fin.txt -start 10 -end 15

perl pl/phyilp-part.pl phylip -in Mega_mtGenome_23112011_Fin.txt -start 10 -end 15 -format fasta

=head1 DESCRIPTION

phylip-part will help you to extract something from or extract parts of a PHYLIP
file.

Command:
  phylip - convert a PHYLIP alignment file to another file.

=head1 OPTIONS

=over 8

=item B<-in> <file>

A file should be an interleaved PHYLIP file.  Standard input is used
unless an input file name is specified.

=item B<-out> <file>

An output file should be an interleaved PHYLIP file.  Standard output is used
unless an output file name is specified.

=item B<-start> <number>

A start position of a part in the PHYLIP alignment file. The default is 1 unless
the start is specified.

=item B<-end> <number>

A end position of a part in the PHYLIP alignment file. The default is the length
of the alignment if the end is not specified.

=item B<-format> <string>

The input foramt is given by the command of this PERL script. The output format
is given using this option.

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
