#!/usr/bin/env perl

use strict;
use warnings;

my %counter;
while (<>) {
  chomp;
  my @x = split (/\t/);

  $x[2] = ($x[2] eq 'match') ? 'mRNA' : 'CDS';

  $x[8] =~ /Target \w+:(\S+)/ or die "Error, cannot parse target name from $x[8]";
  my $acc = $1;

  if ($x[2] eq 'mRNA') {
    $counter{$acc}++;
  }

  $acc .= "." . $counter{$acc};

  $x[8] = "GenePrediction $acc";
  
  print join ("\t", @x) . "\n";

}

exit(0);
