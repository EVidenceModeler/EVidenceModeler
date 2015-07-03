#!/usr/bin/env perl

use strict;
use warnings;
use Carp;

my $usage = "usage: $0 command_list\n\n";

my $command_list = $ARGV[0] or die $usage;


my $failure_count = 0;
my $command_no = 0;
open (my $fh, $command_list) or die "Error, cannot open file $command_list";
while (<$fh>) {
    chomp;
    my $cmd = $_;
    $command_no++;
    
    my $ret = system $cmd;
    
    if ($ret) { $failure_count++; }

    print "#$command_no\t$ret\t$cmd\n";
}

if ($failure_count) {
    die "There were $failure_count of $command_no commands that failed execution.\n";
}
else {
    print STDERR "All commands executed successfully.\n";
    exit(0);
}

