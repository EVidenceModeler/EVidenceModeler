#!/usr/bin/env perl

use strict;

my $usage = "\n\nusage: $0 sequenceFile btabFile\n\n";
my $sequenceFile = $ARGV[0];
my $btabFile = $ARGV[1];
my $min_per_ID = $ARGV[2] || 1;

unless ($sequenceFile && $btabFile) {
    die $usage;
}

unless (-s $sequenceFile) {
    die "Error, can't find $sequenceFile\n";
}

unless (-e $btabFile) {
    die "Error, can't find btab_file: $btabFile";
}

open (SEQUENCE, "<$sequenceFile") or die "Cannot read $sequenceFile";
my $seqLen = 0;
while (<SEQUENCE>) {
    if (/^>/) {
	next;
    }
    my $seq = $_;
    $seq =~ s/\s//g;
    $seqLen += length($seq);
}
close SEQUENCE;
print STDERR "sequence length: $seqLen\n";

## init conseqArray
my @conseqArray;
for (my $i=0; $i < $seqLen; $i++) {
    $conseqArray[$i] = 2;
}

if (-s $btabFile) {
    open (BTAB, "<$btabFile") or die "Error, cannot read $btabFile\n";
    while (<BTAB>) {
        my @x = split (/\t/);
        unless ($x[10] >= $min_per_ID) { next;} #invalid btab entry
        
        my ($lend, $rend) = sort {$a<=>$b} ($x[6], $x[7]);
        for (my $j=$lend; $j <= $rend; $j++) {
            $conseqArray[$j-1] = 1; #consider it a match
        }
    }
    close BTAB;
} else {
    print STDERR "Sorry, no btab file ($btabFile)\n";
}


print ">$sequenceFile.conseq\n" . join ("", @conseqArray) . "\n";

exit(0);

