#!/usr/bin/env perl

use strict;
use lib ($ENV{EUK_MODULES});
use Overlap_piler;

my $MAX_COVERAGE = 10;

my $usage = "\n\nusage: $0 promerCoordsFile\n\n";

my $coordsFile = $ARGV[0] or die $usage;

my @coordsData;
my @coordsIndex;
my $counter = 0;

my $piler = new Overlap_piler();


open (COORDS, "<$coordsFile") or die "ERROR, cannot open $coordsFile\n\n";
while (<COORDS>) {
    my $coordsLine = $_;
    my @x = split (/\t/, $coordsLine);
    my ($lend, $rend) = sort {$a<=>$b} ($x[0], $x[1]);
    $piler->add_coordSet($counter, $lend, $rend);
    
    $coordsData[$counter] = $coordsLine;
    $coordsIndex[$counter] = [$lend, $rend];
    
    $counter++;
    
}
close COORDS;

if ($counter) {
    my @clusters = $piler->build_clusters();
    foreach my $cluster (@clusters) {
        my ($lend, $rend) = &get_bounds($cluster);
        my $cluster_size = scalar (@$cluster);
        print STDERR "$lend-$rend\tsize: $cluster_size";
        if ($cluster_size > $MAX_COVERAGE) {
            print STDERR "\t*** ignoring cluster ***\n";
            next;
        }
        print STDERR "\n";
        
        foreach my $index (@$cluster) {
            print $coordsData[$index];
        }
    }
}

exit(0);



####
sub get_bounds {
    my ($cluster) = @_;
    my @eles = @$cluster;
    my $first_ele = shift @eles;
    my ($min_lend, $max_rend) = @{$coordsIndex[$first_ele]};
    foreach my $ele (@eles) {
        my ($lend, $rend) =  @{$coordsIndex[$ele]};
        if ($lend < $min_lend) {
            $min_lend = $lend;
        }
        if ($rend > $max_rend) {
            $max_rend = $rend;
        }
    }
    return ($min_lend, $max_rend);
}

