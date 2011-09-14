#!/opt/local/bin/perl -w
#===============================================================================
#   Author: Sang Chul Choi, BSCB @ Cornell University, NY
#
#   File: prior-file.pl
#   Date: Sat Jun 11 16:11:37 EDT 2011
#===============================================================================

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

$| = 1; # Do not buffer output

my $VERSION = 'prior-file.pl 1.0';

my $man = 0;
my $help = 0;
my %params = ('help' => \$help, 'h' => \$help, 'man' => \$man);        
GetOptions( \%params,
            'help|h',
            'man',
            'verbose',
            'version' => sub { print $VERSION."\n"; exit; },
            'out=s',
            'command=s'
            ) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

=head1 NAME

prior-file.pl - Produces a prior file for IMa2 analysis

=head1 VERSION

prior-file.pl 1.0

=head1 SYNOPSIS

perl prior-file.pl [-h] [-help] [-version] 
  [-out file] 
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

=item B<-out> <file>

An IMa2 prior file.

=item B<-command> <string>

1. number-loci: Shows the number of loci

=back

=head1 AUTHOR

Sang Chul Choi, C<< <goshng_at_yahoo_dot_co_dot_kr> >>

=head1 BUGS

If you find a bug please post a message rnaseq_analysis project at codaset dot
com repository so that I can make prior-file.pl better.

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
require "pl/sub-ima2prior.pl";

my $out;
my $verbose;
my $command = "default";

if (exists $params{out})
{
  $out = $params{out};
}
else
{
  &printError("you did not specify an IMa2 prior file");
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

if ($command eq "default")
{
  my @individual;
  ima2prior_makeDefault ($out);
}

exit;

################################################################################
## END OF THE MAIN SCRIPT.
################################################################################
