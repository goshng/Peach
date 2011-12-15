#!/usr/bin/perl -w

# perl s2m.pl --in inputfile-name --diploid 1 --asn 0
use strict;
use Getopt::Long;
use S2M;

my $inputfile;
my $isdiploid = 1;
my $isasn = 1;
my $result = GetOptions ("in=s" => \$inputfile,
                         "diploid:i" => \$isdiploid,
                         "asn:i" => \$isasn);
my $structure2structurama = new S2M ('inputfile' => $inputfile,
                        'isdiploid' => $isdiploid,
                        'isasn' => $isasn);
$structure2structurama->read();
$structure2structurama->write();


