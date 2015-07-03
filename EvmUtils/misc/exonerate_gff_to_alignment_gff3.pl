#!/usr/bin/env perl

use strict;
use warnings;


=EXONERATE_COMMAND

exonerate --model p2g --showvulgar no --showalignment no --showquerygff no --showtargetgff yes --percent 80 --ryo "AveragePercentIdentity: %pi\n" protein_db.pep target_genome.fasta

=cut


my $usage = "\n\nusage: $0 exonerate.gff [prot|est]\n\n";

my $alignment_gff = $ARGV[0] or die $usage;
my $MODE = $ARGV[1] || 'prot';

my $MATCH_COUNTER = 0;

main: {


    my $contig = "";
    my $target = "";
    my $orient = "";
    my @segments = ();
    my $per_id = "";
    my @rel_coords;

    open (my $fh, $alignment_gff) or die "Error, cannot open file $alignment_gff";
    ## get to first entry
    while (<$fh>) {
        if (/START OF GFF DUMP/) {
            last;
        }
    }
    
    while (<$fh>) {
        chomp;
        my $line = $_;
        
        if (/AveragePercentIdentity: (\S+)/) {
            $per_id = $1;
            
        }
        elsif (/^\#\#gff-version/) {
            &process_gff_entry($contig, $target, $orient, \@segments, $per_id, \@rel_coords);
            # re-init 
            $contig = "";
            $target = "";
            $orient = "";
            @segments = ();
            $per_id = ".";
            @rel_coords = ();
        }
        else {
            my @x = split(/\s+/);
            
            my $feat_type = $x[2];
            unless ($feat_type) { next; }

            if ($feat_type eq 'gene') {
                $line =~ /sequence (\S+)/ or die "Error, cannot parse target from $_";
                $target = $1;
            }
            elsif ($feat_type eq 'exon') {
                $contig = $x[0];
                $orient = $x[6];
                my ($lend, $rend) = sort {$a<=>$b} ($x[3], $x[4]);
                push (@segments, [$lend, $rend]);
            }
            
            elsif ($feat_type eq 'similarity') {
                while ($line =~ /Align \d+ (\d+) (\d+)/g) {
                    my ($pos, $len) = ($1, $2);
                    push (@rel_coords, [$pos, $len]);
                }
            }
                        
        }
    }
    close $fh;

    
    # get last one
    if (@segments) {
        &process_gff_entry($contig, $target, $orient, \@segments, $per_id, \@rel_coords);
    }

    exit(0);
}


####
sub process_gff_entry {
    my ($contig, $target, $orient, $segments_aref, $per_id, $rel_coords_aref) = @_;

    #use Data::Dumper;
    #print Dumper($segments_aref);
    #print Dumper($rel_coords_aref);


    $MATCH_COUNTER++;

    foreach my $segment (@$segments_aref) {
        my ($lend, $rend) = @$segment;
        
        my $rel_info = shift @$rel_coords_aref;
        my ($rel_pos, $rel_len) = @$rel_info;
        my $rel_lend = $rel_pos;
        my $rel_rend = ($MODE eq 'prot') ? $rel_pos + $rel_len/3 : $rel_pos + $rel_len;
        my $type = ($MODE eq 'prot') ? "nucleotide_to_protein_match" : "cDNA_match";
        

        print join("\t", $contig, "exonerate", $type,
                   $lend, $rend, $per_id, $orient,
                   ".", "ID=match.$$.$MATCH_COUNTER;Target=$target $rel_lend $rel_rend") . "\n";
    }
    print "\n"; # add spacer between alignments.

    return;
}


                       
            
            
            
