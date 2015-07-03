#!/usr/bin/env perl

use strict;
use warnings;

use lib $ENV{EUK_MODULES};
use Fasta_reader;

my $fasta_file = $ARGV[0] || *STDIN{IO};

my $fasta_reader = new Fasta_reader($fasta_file);
while (my $seq_obj = $fasta_reader->next()) {
  
  my $acc = $seq_obj->get_accession();
  my $sequence = $seq_obj->get_sequence();

  my $seq_length = length($sequence);

  print join ("\t", $acc, "chromosome", "scaffold", 1, $seq_length, ".", ".", ".", "Scaffold $acc") . "\n";
}


exit(0);


