#!/opt/local/bin/perl -w
#===============================================================================
#   Author: Sang Chul Choi, BSCB @ Cornell University, NY
#
#   File: summarize-sequence-data.pl
#   Date: Fri Jun 10 15:37:12 EDT 2011
#===============================================================================

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

$| = 1; # Do not buffer output

my $VERSION = 'summarize-sequence-data.pl 1.0';

my $man = 0;
my $help = 0;
my %params = ('help' => \$help, 'h' => \$help, 'man' => \$man);        
GetOptions( \%params,
            'help|h',
            'man',
            'verbose',
            'version' => sub { print $VERSION."\n"; exit; },
            'data=s',
            'out=s',
            'command=s',
            'assignment'
            ) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

=head1 NAME

summarize-sequence-data.pl - Summarize FASTA files.

=head1 VERSION

summarize-sequence-data.pl 1.0

=head1 SYNOPSIS

perl summarize-sequence-data.pl [-h] [-help] [-version] 
  [-data directory] 
  [-out file] 
  [-command string] 
  [-assignment] 

=head1 DESCRIPTION

We summarize FASTA files in a directory.

=head1 OPTIONS

=over 8

=item B<-help> | B<-h>

Print the help message; ignore other arguments.

=item B<-man>

Print the full documentation; ignore other arguments.

=item B<-version>

Print program version; ignore other arguments.

=item B<-verbose>

Prints status and info messages during processing.

=item B<***** INPUT OPTIONS *****>

=item B<-data> <directory>

A directory contains FASTA formatted files.

=item B<-out> <file>

An output file.

=item B<-command> <string>

1. number-sequences: List files and the number of sequences
2. individual: List individuals
3. population: List populations
4. polymorphicsite: Make a file in IM2 input format
5. samplesize: List populations with sample sizes for each locus

=item B<-assignment>

With option -command polymorphicsite this produces an IM2 input file for
population assignment and tree inference.

=back

=head1 AUTHOR

Sang Chul Choi, C<< <goshng_at_yahoo_dot_co_dot_kr> >>

=head1 BUGS

If you find a bug please post a message rnaseq_analysis project at codaset dot
com repository so that I can make summarize-sequence-data.pl better.

=head1 COPYRIGHT

Copyright (C) 2011  Sang Chul Choi

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

require "pl/sub-error.pl";
require "pl/sub-fasta.pl";

my $data;
my $out;
my $verbose;
my $command = "number-sequences";
my $assignment = 0;

if (exists $params{data})
{
  $data = $params{data};
}
else
{
  &printError("you did not specify a data directory of FASTA files");
}

if (exists $params{out})
{
  $out = $params{out};
}

if (exists $params{command})
{
  $command = $params{command};
}


if (exists $params{verbose})
{
  $verbose = 1;
}

if (exists $params{assignment})
{
  $assignment = 1;
}


################################################################################
## DATA PROCESSING
################################################################################

if ($command eq "individual")
{
  my @individual;
  opendir DIR, $data;
  my @FILES= readdir(DIR); 
  for my $file (@FILES)
  {
    if ($file =~ /fas$/)
    {
      my @s = fasta_getIndividual ("$data/$file");
      push @individual, @s;
    }
  }
  my %seen = (); my @uniqu = grep { ! $seen{$_} ++ } @individual;
  print join ("\n", @uniqu);
}
elsif ($command eq "number-sequences")
{
  # List files
  opendir DIR, $data;
  my @FILES= readdir(DIR); 
  for my $file (@FILES)
  {
    if ($file =~ /fas$/)
    {
      my $numberSequence = fasta_getNumberSequences ("$data/$file");
      print "$file\t$numberSequence\n";
    }
  }
}
elsif ($command eq "population")
{
  my @population;
  opendir DIR, $data;
  my @FILES= readdir(DIR); 
  for my $file (@FILES)
  {
    if ($file =~ /fas$/)
    {
      my @s = fasta_getPopulation ("$data/$file");
      push @population, @s;
    }
  }
  my %seen = (); my @uniqu = grep { ! $seen{$_} ++ } @population;
  print join ("\n", @uniqu);
}
elsif ($command eq "polymorphicsite")
{
  my @population = ("Oryza sativa indica group",
                    "Oryza rufipogon",
                    "Oryza sativa japonica group");

  ###############################################################
  # Starts to print a header of the IMa2 input file.
  print <<HEADIMA2INPUT;
Molina et al. (2011) of PNAS
# This file is in IMa2 input format
#
# The original sequence data set was cited by:
# Molina, J.; Sikora, M.; Garud, N.; Flowers, J. M.; Rubinstein, S.;
# Reynolds, A.; Huang, P.; Jackson, S.; Schaal, B. A.; Bustamante, C.  D.; Boyko,
# A. R. & Purugganan, M. D., 2011 Molecular evidence for a single evolutionary
# origin of domesticated rice. Proc Natl Acad Sci U S A 108: 8351-8356.
# It was downloaded from 
# http://puruggananlab.bio.nyu.edu/Rice_data/RiceScan.zip.
#
# The preferred tree is
# ((0,2):3,1):4
# and the other two possible trees are
# (0,(1,2):3):4 and 
# ((0,1):3,2):4.
#
# Replace NLOCI with a positive integer.
3
indica rufipogon japonica
((0,2):3,1):4
NLOCI
HEADIMA2INPUT

  opendir DIR, $data;
  my @FILES= readdir(DIR); 
  my $nloci = 0;
  for my $file (@FILES)
  {
    
    if ($file =~ /fas$/)
    {
      $nloci += fasta_getPolymorphicSites ("$data/$file", \@population, $assignment);
    }
  }
  print "TOTAL NUMBER OF LOCI IS $nloci\n";
}
elsif ($command eq "samplesize")
{
  my @population = ("Oryza barthii",
                    "Oryza sativa indica group",
                    "Oryza meridionalis ",
                    "Oryza nivara",
                    "Oryza rufipogon",
                    "Oryza sativa japonica group");
  opendir DIR, $data;
  my @FILES= readdir(DIR); 
  for my $file (@FILES)
  {
    if ($file =~ /fas$/)
    {
      fasta_getSampleSize("$data/$file", \@population);
    }
  }
}

exit;

