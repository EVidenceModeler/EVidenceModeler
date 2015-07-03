#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "\n\nusage: $0 database query outputfilename blast_prog_name blast_option_list\n\n";

my $database = $ARGV[0] or die $usage;
my $query = $ARGV[1] or die $usage;
my $outputfile = $ARGV[2] or die $usage;
my $blast_prog = $ARGV[3] or die $usage;

shift @ARGV;
shift @ARGV;
shift @ARGV;
shift @ARGV;

my $blast_opts = join (" ", @ARGV);

if (! (-e $database && -s $query) ) {
    die "Error, cannot locate database ($database) or query ($query)\n";
}


if (! -s $database) {
    print "Nothing to search. Database is empty.\n";
    &process_cmd("touch $outputfile");
    &process_cmd("touch $outputfile.btab");
    exit(0);
}


my $cmd = "pressdb $database";
&process_cmd($cmd);

$cmd = "/usr/local/common/$blast_prog $database $query $blast_opts  > $outputfile";
&process_cmd($cmd);


$cmd = "btab2 $outputfile";
&process_cmd($cmd);

exit(0);
    
####
sub process_cmd {
    my $cmd = shift;
    my $ret = system $cmd;
    if ($ret) {
        die "ERROR CMD: ($cmd) died with ret ($ret)";
    } else {
        print "SUCCESS: $cmd ret:($ret) OK\n";
    }
}
