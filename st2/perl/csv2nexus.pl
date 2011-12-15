#line 21160 "noweb/imamp.nw"
#!/usr/bin/perl -w
push (@INC, `pwd`);

#line 49281 "noweb/imamp.nw"
#!/usr/bin/perl -w
push (@INC, `pwd`);

#line 45074 "noweb/imamp.nw"
use strict;
use Getopt::Long;
use CSV2Nexus;

#line 45056 "noweb/imamp.nw"
my $filename;
my $comment;
my $npops;
my $poptree;
my $popns;
my $im;
my $result = GetOptions ("input=s"  => \$filename,
                         "comment=s"  => \$comment);
my $csv2nexus = new CSV2Nexus ('filename' => $filename,
                              'comment' => $comment,
                              'popns' => $popns,
                              'im' => $im,
                              'npops' => $npops,
                              'poptree' => $poptree);
$csv2nexus->read ();
$csv2nexus->write ();


