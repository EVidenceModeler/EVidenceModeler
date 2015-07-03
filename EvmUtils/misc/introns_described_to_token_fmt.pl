#!/usr/bin/env perl

use strict;
use warnings;

while (<STDIN>) {
    unless (/\w/) { next; }
    
    my @x = split (/\t/);
    print $x[0] . "_" . $x[3] . "_" . $x[4] . "_[" . $x[6] . "]\n";
    
}


exit(0);


