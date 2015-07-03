#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use Cwd;



my $usage = "usage: $0 accessions_list_file\n\n";

my $acc_list = $ARGV[0] or die $usage;
my $min_per_ID = $ARGV[1] || 1;

my @asmbls = `cat $acc_list`;

my $utildir = $FindBin::Bin;
my $working_dir = cwd();


foreach my $asmbl_id (@asmbls) {
    
    $asmbl_id =~ s/\s//g;
    
    print "Processing $asmbl_id\n";
    
    my $target_asmbl_dir = "$working_dir/genome_distributed/$asmbl_id";
        
    my $seq_file = "$target_asmbl_dir/$asmbl_id.seq";
    my $btab_file = "$target_asmbl_dir/$asmbl_id.tblastx.btab";
    
    my $cmd = "$utildir/conseq_btab.pl $seq_file $btab_file $min_per_ID > $target_asmbl_dir/$asmbl_id.conseq";
    my $ret = system $cmd;
    if ($ret) {
        die "$asmbl_id\tERROR\n";
        
    } else {
        print "$asmbl_id\tOK\n";
    }
}

exit(0);


