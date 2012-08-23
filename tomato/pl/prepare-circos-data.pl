#!/usr/bin/perl
###############################################################################
# Copyright (C) 2012 Sang Chul Choi
#
# This file is part of RNAseq Analysis.
# 
# RNAseq Analysis is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# RNAseq Analysis is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with RNAseq Analysis.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
$| = 1; # Do not buffer output
my $VERSION = 'prepare-circos-data.pl 1.0';

my $cmd = ""; 
sub process {
  my ($a) = @_; 
  $cmd = $a; 
}

my $man = 0;
my $help = 0;
my %params = ('help' => \$help, 'h' => \$help);        
GetOptions( \%params,
            'help|h',
            'verbose',
            'version' => sub { print $VERSION."\n"; exit; },
            'chr=i',
            'nchr=i',
            'il=s',
            'ilpositionfile=s',
            'ilpositiondir=s',
            'in=s',
            'infile=s',
            'outdir=s',
            'out=s',
            '<>' => \&process
            ) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

################################################################################
## COMMANDLINE OPTION PROCESSING
################################################################################

my $in;
my $infile;
my $out;
my $outfile;

if (exists $params{in})
{
  $in = $params{in};
  open ($infile, "<", $in) or die "cannot open < $in: $!";
}
else
{
  $infile = *STDIN;   
}

if (exists $params{out})
{
  $out = "$params{out}";
  open ($outfile, ">", $out) or die "cannot open > $out: $!";
}
else
{
  $outfile = *STDOUT;   
}

if ($cmd eq "")
{
  print STDERR "ERROR: You need a command.\n";
	exit(0);
}
elsif ($cmd eq "links")
{
  unless (defined $params{il}
          and defined $params{ilpositiondir}
          and defined $params{chr})
  {
    print STDERR "ERROR: You need options for $cmd command.\n";
    exit(0);
  }
}
elsif ($cmd eq "exp")
{
  unless (defined $params{il}
          and defined $params{ilpositiondir}
          and defined $params{chr})
  {
    print STDERR "ERROR: You need options for $cmd command.\n";
    exit(0);
  }
}
elsif ($cmd eq "inter")
{
  unless (defined $params{il}
          and defined $params{ilpositiondir}
          and defined $params{chr})
  {
    print STDERR "ERROR: You need options for $cmd command.\n";
    exit(0);
  }
}
elsif ($cmd eq "exppen")
{
  unless (defined $params{infile}
          and defined $params{nchr}
          and defined $params{outdir})
  {
    print STDERR "ERROR: You need options for $cmd command.\n";
    exit(0);
  }
}
elsif ($cmd eq "snp")
{
  unless (defined $params{chr})
  {
    print STDERR "ERROR: You need options for $cmd command.\n";
    exit(0);
  }
}
elsif ($cmd eq "othersnp")
{
  unless (defined $params{chr})
  {
    print STDERR "ERROR: You need options for $cmd command.\n";
    exit(0);
  }
}
elsif ($cmd eq "pen")
{
  unless (defined $params{infile}
          and defined $params{nchr}
          and defined $params{outdir})
  {
    print STDERR "ERROR: You need options for $cmd command.\n";
    exit(0);
  }
}

################################################################################
## DATA PROCESSING
################################################################################
if ($cmd eq "links")
{
  # Find the IL position.
  my $sourceStartPosition = 0;
  open IL, "$params{ilpositiondir}/$params{il}.csv"
    or die "cannot open $params{ilpositiondir}/$params{il}.csv";
  my $line = <IL>;
  my @a = split /\s+/, $line;
  $sourceStartPosition = $a[2];
  if ($sourceStartPosition == 0)
  {
    print STDERR "ERROR: No IL position available in $params{ilpositionfile}\n";
    exit(0);
  }
  my $sourceEndPosition = $sourceStartPosition + 10;

  # Print the links file.
  my $i = 0;
  while (<$infile>)
  {
    chomp;
    $i++;
    next if $i == 1;
    my @a = split /\t/;
    # Solyc00g009120.2.1	2.98296807960194	8831059
    if ($a[0] =~ /Solyc(\d\d)g/)
    {
      my $chrN = int($1);
      my $endPosition = $a[2] + 10;
      print $outfile "sl$params{chr} $sourceStartPosition $sourceEndPosition sl$chrN $a[2] $endPosition\n";
    }
    else
    {
      die "Unknown mRNA: $_\n";
    }
  }
  close($infile);
}
elsif ($cmd eq "exp" 
       or $cmd eq "inter")
{
  my $sourceStartPosition = 0;
  open IL, "$params{ilpositiondir}/$params{il}.csv"
    or die "cannot open $params{ilpositiondir}/$params{il}.csv";
  my $line = <IL>;
  my @a = split /\s+/, $line;
  $sourceStartPosition = $a[0];
  my $sourceEndPosition = $a[1];
  if ($sourceStartPosition == 0)
  {
    print STDERR "ERROR: No IL position available in $params{ilpositiondir}\n";
    exit(0);
  }

  # Print the links file.
  my $i = 0;
  while (<$infile>)
  {
    chomp;
    $i++;
    next if $i == 1;
    my @a = split /\t/;
    # Solyc00g009120.2.1	2.98296807960194	8831059
    if ($a[0] =~ /Solyc(\d\d)g/)
    {
      my $chrN = int($1);
      my $endPosition = $a[2] + 10;

      if ($cmd eq "exp")
      {
        if ($chrN eq $params{chr})
        {
          if ($sourceStartPosition <= $a[2] and $a[2] <= $sourceEndPosition)
          {
            print $outfile "sl$chrN $a[2] $endPosition $a[1]\n";
          }
        }
      }
      else
      {
        if ($chrN eq $params{chr})
        {
          unless ($sourceStartPosition <= $a[2] and $a[2] <= $sourceEndPosition)
          {
            print $outfile "sl$chrN $a[2] $endPosition $a[1]\n";
          }
        }
        else
        {
          print $outfile "sl$chrN $a[2] $endPosition $a[1]\n";
        }
      }
    }
    else
    {
      die "Unknown mRNA: $_\n";
    }
  }
  close($infile);
}
elsif ($cmd eq "exppen")
{
  for (my $chrI = 0; $chrI <= $params{nchr}; $chrI++)
  {
    open OUT, ">$params{outdir}/chr$chrI"
      or die "cannot open $params{outdir}/chr$chrI"; 
    open IN, "$params{infile}" 
      or die "cannot open $params{infile}";
    my $i = 0;
    while (<IN>)
    {
      chomp;
      $i++;
      next if $i == 1;
      my @a = split /\t/;
      # Solyc00g009120.2.1	2.98296807960194	8831059
      if ($a[0] =~ /Solyc(\d\d)g/)
      {
        my $chrN = int($1);
        my $endPosition = $a[2] + 10;

        if ($chrN eq $chrI)
        {
          print OUT "sl$chrN $a[2] $endPosition $a[1]\n";
        }
      }
      else
      {
        die "Unknown mRNA: $_\n";
      }
    }
    close(IN);
    close(OUT);
  }
}
elsif ($cmd eq "snp" or $cmd eq "othersnp")
{
  # Print the links file.
  my $count = 0;
  my $countM = 0;
  my $countI = 0;
  my $countD = 0;
  my $i = 0;
  while (<$infile>)
  {
    chomp;
    $i++;
    next if $i == 1;
    my @a = split /\t/;
    if ($a[2] =~ /SL2.40ch(\d\d)/)
    {
      my $chrN = int($1);
      my $startPosition = $a[3];
      my $endPosition = $a[3] + 10;
      if (($cmd eq "snp" and $chrN eq $params{chr})
          or ($cmd eq "othersnp" and $chrN ne $params{chr}))
      {
        print $outfile "sl$chrN $startPosition $endPosition color=white\n";
        $count++;
        if ($a[0] eq 'M')
        {
          $countM++;
        } 
        elsif ($a[0] eq 'I')
        {
          $countI++;
        } 
        elsif ($a[0] eq 'D')
        {
          $countD++;
        } 
      }
    }
    else
    {
      die "Unknown mRNA: $_\n";
    }
  }
  close($infile);
  print STDERR "sl$params{chr}\t$params{in}\t$countM\t$countI\t$countD\t$count\n";
}
elsif ($cmd eq "pen")
{
  for (my $chrI = 0; $chrI <= $params{nchr}; $chrI++)
  {
    open OUT, ">$params{outdir}/chr$chrI"
      or die "cannot open $params{outdir}/chr$chrI"; 
    open IN, "$params{infile}" 
      or die "cannot open $params{infile}";
    my $i = 0;
    my $count = 0;
    while (<IN>)
    {
      chomp;
      $i++;
      next if $i == 1;
      my @a = split /\t/;
      if ($a[2] =~ /SL2.40ch(\d\d)/)
      {
        my $chrN = int($1);
        my $startPosition = $a[3];
        my $endPosition = $a[3] + 10;
        if ($chrN eq $chrI)
        {
          print OUT "sl$chrN $startPosition $endPosition color=white\n";
          $count++;
        }
      }
      else
      {
        die "Unknown mRNA: $_\n";
      }
    }
    close(IN);
    close(OUT);
    print $outfile "sl$chrI\t$count\n";  
  }
}
elsif ($cmd eq "circos")
{
  my @colors = ("black","red","blue","purple","orange","grey","magenta","brown","violet","pink","tomato");
  for (my $chrI = 1; $chrI <= $params{nchr}; $chrI++)
  {
#  my $chrI = 10;
  my %il;

    my $colori = 0;
    open IN, "$params{infile}" 
      or die "cannot open $params{infile}";
    while (<IN>)
    {
      chomp;
      if (/^IL$chrI-.+/)
      {
        $colori++;
        my @a = split /,/;
# print "$_\t$a[0]\t$a[1]\n";
        # my @a = split /\s+/;
        my $rec = {};
        $il{$a[0]} = $rec;
        if (1 < scalar(@a))
        {
          $rec->{order} = $a[1];
        }
        else
        {
          $rec->{order} = $colori;
        }
        if (2 < scalar(@a))
        {
          $rec->{color} = $a[2];
        }
        else
        {
          $rec->{color} = $colors[$colori];
        }
      }
    }
    close(IN);
 
foreach my $key (keys %il)
{
print STDERR "$key,$il{$key}{order},$il{$key}{color}\n";
}

    open OUT, ">$params{outdir}/chr$chrI.conf"
      or die "cannot open $params{outdir}/chr$chrI.conf"; 

print OUT <<OUT 
# The chromosome sizes
karyotype = data/karyotype/karyotype.sol.txt
chromosomes_units = 1000000

# Enlarge 
chromosomes_scale = sl$chrI:5

<image>                                                                         
auto_alpha_colors = yes                                                         
auto_alpha_steps  = 5
</image>

<colors>
<<include etc/colors.conf>>
<<include etc/colors.unix.txt>> 
</colors>
OUT
;

print OUT "<links>\n";

foreach my $key (keys %il)
{
print OUT <<OUT
<link>
file          = output/links/$key.csv
radius        = 0.4r
bezier_radius = 0r
color         = $il{$key}{color}_a2
thickness     = 2

# <rules>
# <rule>
# condition     = var(intrachr)
# show          = no
# </rule>
# </rules>

</link>
OUT
;
}

print OUT "</links>\n";


print OUT "<plots>\n";

# PEN SNPs
my $radiusBase = 0.70;
print OUT <<OUT
<plot>                                                                          
type = tile
file = output/pen/chr$chrI
r0   = 1.07r
r1   = 1.07r

layers = 1
margin      = 0.02u                                                             
thickness   = 50
padding     = 8 

layers_overflow       = collapse
stroke_thickness = 1
stroke_color     = green_a2
</plot>
OUT
;

# SNPs
foreach my $key (keys %il)
{
my $radius = $radiusBase - 0.06 * $il{$key}{order};
print OUT <<OUT
<plot>
type = tile
file = output/snp/${key}_SNP
r0   = ${radius}r
r1   = ${radius}r

layers = 1
margin      = 0.02u
thickness   = 50
padding     = 8 

layers_overflow       = collapse
stroke_thickness = 1
stroke_color     = $il{$key}{color}_a2
</plot>
OUT
;

=cut
print OUT <<OUT
<plot>
type = tile
file = output/othersnp/${key}_SNP
r0   = ${radius}r
r1   = ${radius}r

layers = 1
margin      = 0.02u
thickness   = 50
padding     = 8 

layers_overflow       = collapse
stroke_thickness = 1
stroke_color     = $il{$key}{color}_a2
</plot>
OUT
;
=cut
}

print OUT <<OUT
<plot>
type = scatter
file             = output/exppen/chr$chrI
fill_color       = green_a2
stroke_color     = green_a2
glyph            = circle
glyph_size       = 15

max   = 15
min   = -15
r1    = 1.35r
r0    = 1.15r

<axes>
<axis>
thickness = 2
position  = 0
color     = red
</axis>
<axis>
thickness = 1
color     = lgrey
spacing   = 5
</axis>
</axes>

</plot>
OUT
;

foreach my $key (keys %il)
{
# Intra-chromosomal eQTL
print OUT <<OUT
<plot>
type = scatter
file             = output/exp/${key}.csv
fill_color       = $il{$key}{color}_a2
stroke_color     = $il{$key}{color}_a2
glyph            = circle
glyph_size       = 15

max   = 15
min   = -15
r1    = 0.97r
r0    = 0.77r

<axes>
<axis>
thickness = 2
position  = 0
color     = red
</axis>
<axis>
thickness = 1
color     = lgrey
spacing   = 5
</axis>
</axes>

</plot>
OUT
;

# Inter-chromosomal eQTL
print OUT <<OUT
<plot>
type = scatter
file             = output/inter/${key}.csv
fill_color       = black_a2
stroke_color     = black_a2
glyph            = circle
glyph_size       = 15

max   = 15
min   = -15
r1    = 0.97r
r0    = 0.77r

</plot>
OUT
;
}
print OUT "</plots>\n";


print OUT <<OUT
<<include ideogram.conf>>

<<include ticks.conf>>

<image>
<<include etc/image.conf>>                
</image>

<<include etc/colors_fonts_patterns.conf>> 

<<include etc/housekeeping.conf>> 
OUT
;

    close(OUT);
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

prepare-circos-data.pl - Create data files for circos visualization.

=head1 VERSION

prepare-circos-data.pl 1.0

=head1 SYNOPSIS

perl prepare-circos-data.pl pen -infile input/raw/IL_SNP_gene/pen_SNP
  -nchr 12 -outdir output/pen

perl prepare-circos-data.pl links -in IL10-1-1.csv -chr 10 
  -ilpositionfile output/10.il -out output/links/IL10-1-1.csv

=head1 DESCRIPTION

prepare-circos-data.pl creates data files for circos visualization.

=head1 OPTIONS

  command: 

  pen - Create circos files for SNPs.
  snp - 
  links - Find core genes in a BED file using a mcl dump file.

pen: PEN snp. A file called, input/raw/IL_SNP_gene/pen_SNP, contains lines such
as "M  A->G  SL2.40ch00  552376". For each chromosome ranging from 0 to 12, for
example, we grab the position or the 4th position of lines that are part of the
chromosome by recognizing two digits in the 3rd column e.g., SL2.40ch00, and
create a line, sl0 552376 552386 color=white, in output/pen/chr0. This is used
to plot SNPs in a circos plot.

snp: SNPs.  Files, e.g., input/raw/IL_SNP_gene/IL9-3_SNP, contains lines of
positions for the SNPs (see command pen above). We create a file that contains
a line, sl9 552376 552386 color=white. The file is output/snp/IL9-3_SNP. 

othersnp: SNPs in other chromosomes. Files, e.g., 
input/raw/IL_SNP_gene/IL9-3_SNP, contains lines of
positions for the SNPs (see command pen above). We create a file that contains
a line, sl8 552376 552386 color=white. The file is output/snpother:b 1
/IL9-3_SNP. 



find-position-il.R: IL6-1_SNP no SNPs in chr06. 

links: links. 

exp: intra expression

inter: inter expression

circos: create circos files.


=over 8

=item B<-help> | B<-h>

Print the help message; ignore other arguments.

=item B<-version>

Print program version; ignore other arguments.

=item B<-verbose>

Prints status and info messages during processing.

=item B<***** INPUT OPTIONS *****>

=item B<-mcl> <file>

An mcl dump file. The file contains rows of strings separated by spaces. The
number of strings decreases as you go towards the end of the file.

=item B<-bed> <file>

A BED file for genes or some of the strings in the mcl dump file.

=item B<-out> <file>

An output file.

=item B<-coregenome> <string>

A string for a list of core genomes.

=item B<-orthomcl>

This allows the output file of OrthoMCL. The sequence names in the output file
from OrthoMCL consist of strain or species name and gene name separated by a
vertical line.

=item B<-orthomclgroupname>

The final output file of OrthoMCL can contain names for a group or cluster at
the first column until a colon. We remove the first column.

=back

=head1 AUTHOR

Sang Chul Choi, C<< <goshng_at_yahoo_dot_co_dot_kr> >>

=head1 BUGS

If you find a bug please post a message rnaseq_analysis project at codaset dot
com repository so that I can make prepare-circos-data.pl better.

=head1 COPYRIGHT

Copyright (C) 2012 Sang Chul Choi

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
