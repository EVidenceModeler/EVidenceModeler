#!/usr/bin/env perl

use strict;
use warnings;


my $prev_feat_type = "";
my $ofh;

while (<>) {
    if (/^\#/) { next; }
    unless (/\w/) { next; }
    
    my $line = $_;
    
    my @x = split(/\t/);
    my $feat_type = $x[1];
    unless ($feat_type) { next; }
    
    if ($feat_type ne $prev_feat_type) {
        $prev_feat_type = $feat_type;
        close $ofh if $ofh;
        
        print STDERR "writing to $feat_type.gff\n";
        open ($ofh, ">>$feat_type.gff") or die $!;
        
    }

    print $ofh $line if $ofh;

}

close $ofh if $ofh;
        
exit(0);

