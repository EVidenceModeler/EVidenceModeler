#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use Cwd;

my $usage = "\n\nusage: $0 accesion_list /path/to/other_genome_fasta_file\n\n";

my $accession_list = $ARGV[0] or die $usage;
my $other_genome_fasta = $ARGV[1] or die $usage;

unless ($other_genome_fasta =~ /\//) {
    die "need full path to other genome fasta file.";
}

unless (-s $other_genome_fasta) {
    die "Sorry, can't locate file $other_genome_fasta";
}


my $utilDir = $FindBin::Bin;
my $workdir = cwd();

my @accessions = `cat $accession_list`;
foreach my $genome_acc (@accessions) {
    $genome_acc =~ s/\s//g;

    my $asmbl_dir = "$workdir/genome_distributed/$genome_acc";
    chdir $asmbl_dir or die "Error, couldn't cd to $asmbl_dir\n";
    
    my $cmd = "$utilDir/coordsToMiniDatabase.pl $genome_acc.coords.filtered.chains $other_genome_fasta > $genome_acc.miniDB";
    &process_cmd($cmd);
    
}


exit(0);

####
sub process_cmd {
    my $cmd = shift;
    my $ret = system $cmd;
    if ($ret) {
        die "CMD: ($cmd) died with ret ($ret)";
    } else {
        print "SUCCESS: $cmd: OK ret($ret)\n";
    }
}

