#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "\nusage: $0 asmbl_list_file gff_file\n\n";

my $asmbl_list_file = $ARGV[0] or die $usage;
my $gff_file = $ARGV[1] or die $usage;

my %asmbls;
{
  open (my $fh, $asmbl_list_file) or die "Error, cannot open file $asmbl_list_file";
  while (<$fh>) {
    chomp;
    s/\s//g;

    $asmbls{$_}++;
  }
  close $fh;
}

{
  open (my $fh, $gff_file) or die "Error, cannot open file $gff_file";
  while (<$fh>) {
    my @x = split (/\t/);
	unless ($asmbls{$x[0]}) {
      print;
    }
  }
  close $fh;
}


exit(0);

