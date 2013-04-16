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
my @subsetOfSequences;
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
            'subsetOfSequences=s' => \@subsetOfSequences,
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

if (length @subsetOfSequences > 0)
{
  @subsetOfSequences = split(/,/,join(',',@subsetOfSequences));
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
    $nameOfSequence =~ s/\s+//g;
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

  # Subset sequences
  if (scalar @subsetOfSequences > 0)
  {
    my %countHash = ();
    foreach (@sequenceNames, @subsetOfSequences) { $countHash{$_}++ }
    my @count = grep { $_ > 1 } values(%countHash);
    unless (scalar(@count) == scalar(@subsetOfSequences))
    {
      die "Error: check -subsetOfSequences option if all sequences names exist in the input sequence alignment file\n";
    }

    my @tempSequences;
    for (my $i = 0; $i < scalar @subsetOfSequences; $i++)
    {
      my ($index) = grep { $sequenceNames[$_] eq $subsetOfSequences[$i] } 0..$#sequenceNames;
      die "Assert: index must be defined" unless defined $index;
      push @tempSequences, $sequences[$index];
    }
    @sequenceNames = @subsetOfSequences;
    @sequences = @tempSequences;
    $numberOfSequences = scalar @sequenceNames;
  }

  # Extract a part using start and end positions.
  my @extractedSequences;
  foreach (@sequences)
  {
    my $offset = $params{start} - 1;
    my $len = $params{end} - $params{start} + 1;
    push @extractedSequences, substr($_, $offset, $len);
  }

  #############################################################################
  # Print out sequences
  $numberOfColumns = $params{end} - $params{start} + 1;
  if ($params{format} eq "phylip")
  {
    # Print the part out in PHYLIP format.
    # 
    # Segment 1  Segment 2  Segment 3
    # AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA -> Chunk 1
    # AAAAAAAAAA AAAAAAAAAA AAAA       -> Chunk 2
    # Note that numberOfColumns is now the length of a part.
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
          printf $outfile "%-10s", $sequenceNames[$j];
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
  elsif ($params{format} eq "migrate-n")
  {
    print $outfile " 1 1\n";
    print $outfile "$numberOfColumns\n";
    print $outfile "$numberOfSequences Population1\n";
    for (my $i = 0; $i < $numberOfSequences; $i++)
    {
      $sequenceNames[$i] =~ s/\s+//g;
      printf $outfile "%-10s%s\n", $sequenceNames[$i], $extractedSequences[$i];
    }
  }
  elsif ($params{format} eq "genepop")
  {
    print $outfile "Automatically generated by phylip-part.pl\n";
    print $outfile "locus S_M_1 <M>\n";
    print $outfile "Pop\n";
    for (my $i = 0; $i < $numberOfSequences; $i++)
    {
      $sequenceNames[$i] =~ s/\s+//g;
      printf $outfile "%-10s , <[%s]>\n", $sequenceNames[$i], $extractedSequences[$i];
    }
  }
  elsif ($params{format} eq "fasta")
  {
    for (my $i = 0; $i < $numberOfSequences; $i++)
    {
      $sequenceNames[$i] =~ s/\s+//g;
      print $outfile ">$sequenceNames[$i]\n";

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
  elsif ($params{format} eq "nexus")
  {
    # Print the part out in NEXUS format.
    # 
    # Segment 1  Segment 2  Segment 3
    # AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA -> Chunk 1
    # AAAAAAAAAA AAAAAAAAAA AAAA       -> Chunk 2
    # Note that numberOfColumns is now the length of a part.
    $numberOfColumns = $params{end} - $params{start} + 1;
    my $numberOfChunks = int($numberOfColumns / 60);
    if ($numberOfColumns % 60 != 0)
    {
      $numberOfChunks++;
    }

    # Print the taxa block.
    print $outfile "#NEXUS\n";
    print $outfile "begin taxa;\n";
    print $outfile "  dimensions ntax=$numberOfSequences;\n";
    print $outfile "  taxlabels\n";
    for (my $i = 0; $i < $numberOfSequences; $i++)
    {
      $sequenceNames[$i] =~ s/\s+//g;
      print $outfile "    $sequenceNames[$i]\n";
    }
    print $outfile ";\nend;\n";
    print $outfile "begin characters;\n";
    print $outfile "  dimensions nchar=$numberOfColumns;\n"; 
    print $outfile "  format missing=? gap=- matchar=. datatype=nucleotide interleave=yes;\n";
    print $outfile "  matrix\n";

    for (my $i = 0; $i < $numberOfChunks; $i++)
    {
      print $outfile "\n";
      for (my $j = 0; $j < $numberOfSequences; $j++)
      {
        printf $outfile "%-10s", $sequenceNames[$j];

        my $chunkOfSequence = substr $extractedSequences[$j], 10*6*$i, 60;
        print $outfile "$chunkOfSequence\n";
      }
    }
    print $outfile ";\nend;\n";
  }
  else
  {
    die "No such output format: $params{format}";
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

perl pl/phyilp-part.pl phylip -in Mega_mtGenome_23112011_Fin.txt -start 10 -end 15 -format nexus

perl pl/phyilp-part.pl phylip -in Mega_mtGenome_23112011_Fin.txt -start 10 -end 15 -format genepop

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
is given using this option: phylip, migrate-n, fasta, nexus, or genepop. 

=item B<-subsetOfSequences> <string,string,...>

A comma-separated string of sequence identifiers: A,B,C would be extracted from
the sequence alignment of five sequences A, B, C, D, E. 

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
