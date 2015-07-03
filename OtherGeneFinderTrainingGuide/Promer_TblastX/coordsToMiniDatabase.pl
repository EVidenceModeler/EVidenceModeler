#!/usr/bin/env perl

use strict;

my $PAD_DIST = 5000; #amt sequence on each entry to retrieve

my $usage = "\n\nusage: $0 chainsFile FastaDatabase\n\n";

my $chainsFile = $ARGV[0] or die $usage;
my $fastaDB = $ARGV[1] or die $usage;

unless (-s $fastaDB) {
    die "Error, cannot find $fastaDB\n";
}

unless (-s "$fastaDB.cidx") {
    my $cmd = "cdbfasta $fastaDB";
    my $ret = system $cmd;
    die "Error, $cmd failed (ret: $ret)\n" if $ret;
}

my $counter = 0;
open (CHAINS, "<$chainsFile") or die "Error, cannot open $chainsFile ";
while (<CHAINS>) {
    $counter++;
    chomp;
    my ($acc, $start, $end) = split (/\t/);
    
    ($start, $end) = sort {$a<=>$b} ($start, $end);
    
    $start -= $PAD_DIST;
    if ($start < 1) { 
	$start = 1;
    }
    $end += $PAD_DIST;
    
    my $cmd = "cdbyank -a \"$acc $start $end\" -R $fastaDB.cidx";
    #print "CMD: $cmd\n";
    my $entry = `$cmd`;
    $entry =~ s/>/>${counter}_/;
    
    print $entry;
}
close CHAINS;

