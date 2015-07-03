#!/usr/bin/env perl

use strict;
use warnings;

while (<STDIN>) {
    unless (/\w/) { next; }
    chomp;
    my @x = split (/\t/);
    print $x[0] . "_" . $x[5] . "_" . $x[6] . "_" . $x[3] . "\n";
    
}


exit(0);


