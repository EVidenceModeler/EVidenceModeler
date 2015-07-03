#!/usr/bin/env perl

use strict;

our $DEBUG = 0;

my $usage = "\n\nusage: $0 coordsFile\n\n";

my $coordsFile = $ARGV[0] or die $usage;

my $chainDist = 20000;

my %acc_to_coordsets;


open (COORDS, "<$coordsFile") or die "Error, cannot open $coordsFile\n";
while (<COORDS>) {
    chomp;
    my @x = split (/\t/);
    my ($end5, $end3, $genome_acc) = ($x[2], $x[3], $x[12]);
    
    my ($lend, $rend) = sort {$a<=>$b} ($end5, $end3);

    my $list_ref = $acc_to_coordsets{$genome_acc};
    unless (ref $list_ref) {
	$list_ref = $acc_to_coordsets{$genome_acc} = [];
    }
    push (@$list_ref, [$lend, $rend]);
}
close COORDS;


foreach my $acc (keys %acc_to_coordsets) {
    ## do chaining:
    
    my @coordsets = sort {$a->[0]<=>$b->[0]} @{$acc_to_coordsets{$acc}};
    
    if ($DEBUG) {
	print "\n\n## $acc, incoming coordsets:\n";
	foreach my $coordset (@coordsets) {
	    print "\t@$coordset\n";
	}
    }
    
    my @chains = shift @coordsets;
    
    while (@coordsets) {
	my ($chain_lend, $chain_rend) = @{$chains[$#chains]};
	my $next_coordset = shift @coordsets;
	my ($next_lend, $next_rend) = @$next_coordset;
	
	## see if within gap distance
	if ($next_lend - $chain_rend <= $chainDist ) { 
	    ## yes, within chaining distance
	    if ($next_rend > $chain_rend) {
		$chains[$#chains]->[1] = $next_rend; #update max rend coordinate for cluster entry
	    }
	} else {
	    ## no chaining done here, start new chain:
	    push (@chains, $next_coordset);
	}
    }
    
    foreach my $chain (@chains) {
	print "$acc\t" . join ("\t", @$chain) . "\n";
    }
    
}
    
