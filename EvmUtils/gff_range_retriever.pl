#!/usr/bin/env perl

use strict;
use warnings;
use Carp;

my $usage = "\nusage: $0 seq_id min_lend max_rend [adjust_min_lend_to_1=0] < gff_file > subset_gff_file\n\n";

my $seq_id = $ARGV[0] or die $usage;
my $min_lend = $ARGV[1] or die $usage;
my $max_rend = $ARGV[2] or die $usage;
my $adjust_to_1 = $ARGV[3] || 0; # default is false

my $adjust_coord = $min_lend - 1;

unless ($min_lend =~ /^\d+$/ && $max_rend =~ /^\d+$/) {
    die $usage;
}

my $got_spacer = 0;
while (<STDIN>) {
    unless (/\w/) { 
        print unless ($got_spacer);
        $got_spacer = 1;
        next;
    }

    if (/^\#/) {
        print;
        next;
    }

    my @x = split (/\t/);
    my ($contig, $lend, $rend) = ($x[0], $x[3], $x[4]);
    
    if ($contig eq $seq_id &&
        $lend >= $min_lend &&
        $rend <= $max_rend) {
        
        if ($adjust_to_1) {
            $x[3] -= $adjust_coord;
            $x[4] -= $adjust_coord;
        }
        
        print join ("\t", @x);
        $got_spacer = 0;
    }
}

exit(0);


                                            
