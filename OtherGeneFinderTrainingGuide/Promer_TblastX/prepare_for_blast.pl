#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use Cwd;

my $usage = "usage: $0 accession_list\n\n";

my $accession_list = $ARGV[0] or die $usage;

my $utilDir = $FindBin::Bin;

my $workdir = cwd();

my @asmbls = `cat $accession_list`;

foreach my $asmbl_id (@asmbls) {
    $asmbl_id =~ s/\s//g;
    
    warn "Processing $asmbl_id\n";
    
    eval {
	
        my $asmbl_dir = "$workdir/genome_distributed/$asmbl_id";
        chdir $asmbl_dir or die "Error, couldn't cd to $asmbl_dir\n";
        
        my $deltaFile = "$asmbl_dir/$asmbl_id.delta";
        
        my $cmd = "show-coords -T -H $deltaFile > $asmbl_id.coords";
        &process_cmd ($cmd);
        
        my $cmd = "$utilDir/promer_coords_overlap_filter.pl $asmbl_id.coords > $asmbl_id.coords.filtered";
        &process_cmd($cmd);
        
        my $cmd = "$utilDir/promer_coords_chainer.pl $asmbl_id.coords.filtered > $asmbl_id.coords.filtered.chains";
        &process_cmd($cmd);
        
    };
    
    if ($@) {
        print "Error processing $asmbl_id: $@\n";
    } else {
        print "Processed $asmbl_id OK\n";
    }
    
}

####
sub process_cmd {
    my $cmd = shift;
    my $ret = system $cmd;
    if ($ret) {
	die "CMD: ($cmd) died with ret ($ret)";
    }
}
