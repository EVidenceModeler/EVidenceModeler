#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib qq ($FindBin::Bin);
use Run_Qsub64;
use Cwd;
use List::Util qw (shuffle);
use File::Basename;

my $usage = "usage: $0 genome_accessions blast_prog_name blast_opts\n\n";

my $genome_accessions = $ARGV[0] or die $usage;
my $blast_prog_name = $ARGV[1] or die $usage;
shift @ARGV;
shift @ARGV;

my $blast_opts = join (" ", @ARGV);

my @accessions = `cat $genome_accessions`;

my $utilDir = $FindBin::Bin;
my $workdir = cwd();


my @cmds;

foreach my $genome_acc (@accessions) {
    
    $genome_acc =~ s/\s//g;

    my $genome_seq_dir = "$workdir/genome_distributed/$genome_acc";
    my $genome_seq = "$genome_seq_dir/$genome_acc.seq";
    my $miniDB = "$genome_seq_dir/$genome_acc.miniDB";

    my $outputfile = "$genome_seq_dir/$genome_acc.tblastx";
    
    my $cmd = "$utilDir/run_blast_search.pl $miniDB $genome_seq $outputfile $blast_prog_name $blast_opts";
    
    push (@cmds, $cmd);

}

@cmds = shuffle @cmds;

my $ret = &Run_Qsub64::run(@cmds);
if ($ret) {
    print "Sorry, at least one tblastx job failed.\n";
}
else {
    print " ***** Yippee, everything worked successfully *******\n";
}

exit(0);

