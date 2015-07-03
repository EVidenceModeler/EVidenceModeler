#!/usr/bin/env perl

use strict;
use warnings;


my %match_id_to_seg_list;
my %match_id_to_orient;
my %match_id_to_contig;

while (<STDIN>) {
    chomp;
    unless (/\w/) { next; }
    
    my @x = split (/\t/);
    my ($contig, $lend, $rend, $orient, $info) = ($x[0], $x[3], $x[4], $x[6], $x[8]);
    
    $info =~ /ID=([^;]+)/;
    my $match_id = $1 or die "Error, cannot extract match ID from $info";
    
    $match_id = "$contig/$match_id";
    
    $match_id_to_contig{$match_id} = $contig;
    $match_id_to_orient{$match_id} = $orient;

    my $segment_list_aref = $match_id_to_seg_list{$match_id};
    unless (ref $segment_list_aref) {
        $segment_list_aref = $match_id_to_seg_list{$match_id} = [];
    }

    push (@$segment_list_aref, [$lend, $rend]);
}

foreach my $match_id (keys %match_id_to_seg_list) {
    my $seg_list_aref = $match_id_to_seg_list{$match_id};
    
    my $num_segments = scalar (@$seg_list_aref);

   #  print "match_id: $match_id, $num_segments\n";
    
    if ($num_segments == 1) {
        # no candidate introns.
        next;
    }
    
    my $orient = $match_id_to_orient{$match_id};
    my $contig = $match_id_to_contig{$match_id};

    my @segments = @$seg_list_aref;
    @segments = sort {$a->[0]<=>$b->[0]} @segments;

    my $prev_segment = shift @segments;
    while (@segments) {
        my $next_segment = shift @segments;
        
        my ($prev_lend, $prev_rend) = @$prev_segment;
        my ($next_lend, $next_rend) = @$next_segment;
        
        my ($intron_lend, $intron_rend) = ($prev_rend + 1, $next_lend - 1);
        
        print "${contig}_${intron_lend}_${intron_rend}_[$orient]\n";
        
        $prev_segment = $next_segment;

    }
}



exit(0);

