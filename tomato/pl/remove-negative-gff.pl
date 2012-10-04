open IN, "ITAG2.3_gene_models.gff3.txt";

my $l;
my $i = 0;
while ($l = <IN>) {
  $i++;

  if ($l =~ /^#/) {
    print $l; next;
  }

  my @a = split /\t/, $l;
  unless ($a[4] - $a[3] < 0) {
    print "$l";
  }
}
close IN;
