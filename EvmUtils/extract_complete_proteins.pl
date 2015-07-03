#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Fasta_reader;

my $usage = "usage: $0 fasta_file\n\n";

my $fasta_file = $ARGV[0] or die $usage;

my $fasta_reader = new Fasta_reader($fasta_file);

while (my $seq_obj = $fasta_reader->next()) {

    my $sequence = $seq_obj->get_sequence();
    if ($sequence =~ /^M/ && $sequence =~ /\*$/ ) {
        print $seq_obj->get_FASTA_format();
    }
}

exit(0);


