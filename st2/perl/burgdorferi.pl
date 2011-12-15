#line 41892 "noweb/imamp.nw"
#!/usr/bin/perl -w

use strict;

my $line;
my $strain;
my $allele;

my $nind = 0;
my $nloci = 0;

open BURG, "1" or die;
while ($line = <BURG>)
  {
    $strain = substr $line, 0, 55;
    $allele = substr $line, 60;
    my @strains = split /,\s+/, $strain;
    my @alleles = split /\s+/, $allele;
    foreach my $s (@strains)
      {
        $nind++;
        $nloci = 0;
        foreach my $a (@alleles)
          {
            $nloci++;
          }
      }
  }
close BURG;

open BURG, "1" or die;
print "#NEXUS\n\n";
print "[ burgdorferi ]\n\n";
print "begin data;\n";
printf ("   dimensions nind=%d nloci=%d;\n", $nind, $nloci);
print "   info\n";
while ($line = <BURG>)
  {
    $strain = substr $line, 0, 55;
    $allele = substr $line, 60;
    my @strains = split /,\s+/, $strain;
    my @alleles = split /\s+/, $allele;
    foreach my $s (@strains)
      {
        $s =~ /([\w-]+)/;
        printf ("   %-10s", $1);
        foreach my $a (@alleles)
          {
            if ($a =~ /\w/)
              {
                printf ("(%3d, \?) ", $a);
              }
          }
        print ",\n";
      }
  }
close BURG;
print "   ;\n";
print "end;\n";

