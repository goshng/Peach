#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use IM2Nexus;

my $inputfile;
my $isdiploid = 0;
my $isasn = 1;
my $result = GetOptions ("in=s" => \$inputfile,
                         "diploid:i" => \$isdiploid,
                         "asn:i" => \$isasn);
my $im2nexus = new IM2Nexus ('inputfile' => $inputfile,
                             'isdiploid' => $isdiploid,
                             'isasn' => $isasn);
$im2nexus->read();
$im2nexus->write();
$im2nexus->writemathematica();


