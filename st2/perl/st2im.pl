#line 21160 "noweb/imamp.nw"
#!/usr/bin/perl -w
push (@INC, `pwd`);

#line 49055 "noweb/imamp.nw"
#!/usr/bin/perl -w
push (@INC, `pwd`);

#line 45077 "noweb/imamp.nw"
use strict;
use Getopt::Long;
use Structure2IM;

#line 45052 "noweb/imamp.nw"
my $filename;
my $comment;
my $npops;
my $poptree;
my $first;
my $second;
my $im;
my $result = GetOptions ("input=s"  => \$filename,
                         "comment=s"  => \$comment,
                         "first=s" => \$first,
                         "second=s" => \$second,
                         "im=s" => \$im,
                         "npops=s"  => \$npops,
                         "poptree=s"  => \$poptree);
my $st2im = new Structure2IM ('filename' => $filename,
                              'comment' => $comment,
                              'first' => $first,
                              'second' => $second,
                              'im' => $im,
                              'npops' => $npops,
                              'poptree' => $poptree);
$st2im->read ();
$st2im->write ();


