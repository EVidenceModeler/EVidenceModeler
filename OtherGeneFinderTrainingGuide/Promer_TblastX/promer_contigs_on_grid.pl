#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib qq ($FindBin::Bin);
use Run_Qsub64;
use Cwd;
use List::Util qw (shuffle);

my $usage = "\n\nusage: $0 target_genome_accs /full/path/to/other_genome_db\n\n";

my $target_genome_accs = $ARGV[0] or die $usage;
my $other_genome_db = $ARGV[1] or die $usage;

my $dir = cwd();

my @cmds;

open (my $fh, $target_genome_accs) or die $!;
while (<$fh>) {
    my $acc = $_;
    $acc =~ s/\s//g;

    my $target_dir = "$dir/genome_distributed/$acc";
    my $target_seq = "$target_dir/$acc.seq";
    
    my $cmd = "cd $target_dir; /usr/local/devel/ANNOTATION/EGC_utilities/MUMMER3/MUMmer3.18/promer -p $acc $target_seq $other_genome_db";

    push (@cmds, $cmd);
}

@cmds = shuffle(@cmds);

my $ret = &Run_Qsub64::run(@cmds);
if ($ret) {
    print "Sorry, at least one job failed.\n";
}
else {
    print " ***** Yippee, everything worked successfully *******\n";
}

exit(0);




