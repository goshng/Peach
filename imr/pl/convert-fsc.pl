#!/usr/bin/perl
#===============================================================================
#   Author: Sang Chul Choi, BSCB @ Cornell University, NY
#
#   File: convert-fsc.pl
#   Date: Thu Jun 23 22:27:18 EDT 2011
#   Version: 1.0
#===============================================================================

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;

$| = 1; # Do not buffer output

my $VERSION = 'convert-fsc.pl 1.0';

my $man = 0;
my $help = 0;
my %params = ('help' => \$help, 'h' => \$help);        
GetOptions( \%params,
            'help|h',
            'verbose',
            'version' => sub { print $VERSION."\n"; exit; },
            'outputformat=s',
            'in=s',
            'inTree=s',
            'out=s',
            'length=i',
            'repetition=i'
            ) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

=head1 NAME

convert-fsc.pl - Convert fsc simulated files.

=head1 VERSION

convert-fsc.pl 1.0

=head1 SYNOPSIS

perl convert-fsc.pl.pl [-h] [-help] [-version] [-verbose]
  [-outputformat string] 
  [-in file] 
  [-inTree file] 
  [-out file] 
  [-length number] 
  [-repetition number] 

=head1 DESCRIPTION

convert-fsc.pl converts a FSC simulated file.

=head1 OPTIONS

=over 8

=item B<-help> | B<-h>

Print the help message; ignore other arguments.

=item B<-version>

Print program version; ignore other arguments.

=item B<-verbose>

Prints status and info messages during processing.

=item B<***** INPUT OPTIONS *****>

=item B<-outputformat> <string>

A few output formats can be specified by users.

1. galledtree
2. imr

=item B<-out> <file>

An output file.

=item B<-in> <file>

An input file.

=item B<-inTree> <file>

The true tree file created by Fastsimcoal along with the sequence data.

=item B<-length> <number>

The length of a block.

=item B<-repetition> <number>

The index of repetition in Fastsimcoal.

=back

=head1 AUTHOR

Sang Chul Choi, C<< <goshng_at_yahoo_dot_co_dot_kr> >>

=head1 BUGS

If you find a bug please post a message rnaseq_analysis project at codaset dot
com repository so that I can make convert-fsc.pl better.

=head1 COPYRIGHT

Copyright (C) 2011 Sang Chul Choi

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

################################################################################
## COMMANDLINE OPTION PROCESSING
################################################################################

require "pl/sub-error.pl";
require "pl/sub-fsc.pl";

my $outputformat = "imr";
my $in;
my $inTree = "";
my $out;
my $length;
my $repetition;

if (exists $params{outputformat}) {
  $outputformat = $params{outputformat};
}

if (exists $params{in}) {
  $in = $params{in};
} else {
  &printError("you did not specify an in file name");
}

if (exists $params{inTree}) {
  $inTree = $params{inTree};
} 

if (exists $params{repetition}) {
  $repetition = $params{repetition};
} else {
  &printError("you did not specify a repetition");
}

if (exists $params{out}) {
  $out = $params{out};
} else {
  &printError("you did not specify an out file name");
}

if (exists $params{length}) {
  $length = $params{length};
} else {
  &printError("you did not specify a length file name");
}

################################################################################
## DATA PROCESSING
################################################################################

if ($outputformat eq "galledtree")
{
  fsc_convert2gal ($in, $out);
}
elsif ($outputformat eq "imr")
{
  if ($inTree eq "")
  {
    fsc_convert2imr ($in, $out, $length);
  }
  else
  {
    fsc_convert2imrWithTree ($in, $inTree, $out, $length, $repetition);
  }
}
else
{
  die "$outputformat is not available";
}
