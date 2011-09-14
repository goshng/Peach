#!/opt/local/bin/perl -w
#===============================================================================
#   Author: Sang Chul Choi, BSCB @ Cornell University, NY
#
#   File: examine-ima2-data.pl
#   Date: Sat Jun 11 13:10:25 EDT 2011
#===============================================================================

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

$| = 1; # Do not buffer output

my $VERSION = 'examine-ima2-data.pl 1.0';

my $man = 0;
my $help = 0;
my %params = ('help' => \$help, 'h' => \$help, 'man' => \$man);        
GetOptions( \%params,
            'help|h',
            'man',
            'verbose',
            'version' => sub { print $VERSION."\n"; exit; },
            'ima2data=s',
            'command=s'
            ) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

=head1 NAME

examine-ima2-data.pl - Examine an IMa2 input file.

=head1 VERSION

examine-ima2-data.pl 1.0

=head1 SYNOPSIS

perl examine-ima2-data.pl [-h] [-help] [-version] 
  [-ima2data file] 
  [-command string] 

=head1 DESCRIPTION

We inspect an IMa2 input file.

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

=item B<-ima2data> <file>

An IMa2 input file.

=item B<-command> <string>

1. number-loci: Shows the number of loci

=back

=head1 AUTHOR

Sang Chul Choi, C<< <goshng_at_yahoo_dot_co_dot_kr> >>

=head1 BUGS

If you find a bug please post a message rnaseq_analysis project at codaset dot
com repository so that I can make examine-ima2-data.pl better.

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
require "pl/sub-ima2data.pl";

my $ima2data;
my $verbose;
my $command = "number-loci";

if (exists $params{ima2data})
{
  $ima2data = $params{ima2data};
}
else
{
  &printError("you did not specify an IMa2 input file");
}

if (exists $params{command})
{
  $command = $params{command};
}

if (exists $params{verbose})
{
  $verbose = 1;
}

################################################################################
## DATA PROCESSING
################################################################################

if ($command eq "number-loci")
{
  my @individual;

  my ($nloci, $sampleSizes) = ima2data_findNumberLoci ($ima2data);
  print "# The number of loci is $nloci of $ima2data\n";
  print "# The following is a list of sample sizes for the $nloci loci\n";
  for my $s (@$sampleSizes)
  {
    print "\t$s";
  }
  print "\n";
}

exit;

################################################################################
## END OF THE MAIN SCRIPT.
################################################################################
