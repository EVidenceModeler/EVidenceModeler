package TrellisParser;

use strict;
use warnings;
use Carp;
use Trellis;

####
sub new {
    my $packagename = shift;
    my $self = {};

    bless ($self, $packagename);
    return ($self);
}

####
sub parse_trellis_file {
    my $self = shift;
    my $trellis_filename = shift;
    my $final_path_file = shift;
    my ($range_lend, $range_rend) = @_;

    my $trellis = new Trellis();
    
    open (my $fh, $trellis_filename) or die "Error, cannot open trellis file [$trellis_filename]";
    while (<$fh>) {
        unless (/\w/) { next; }
        if (/^\#/) { next; }
        
        if (/^(\w+)([\+\-])\~(\d)~(\d)\s(\d+),\s(\d+).*\sto\s(\w+)([\+\-])\~(\d)~(\d)\s(\d+),\s(\d+)/) {
            my ($typeA, $orientA, $phaseA_left, $phaseA_right, $lendA, $rendA, 
                $typeB, $orientB, $phaseB_left, $phaseB_right, $lendB, $rendB) = ($1, $2, $3, $4, $5, $6, 
                                                                                  $7, $8, $9, $10, $11, $12);
            
            if (defined $range_lend && defined $range_rend) {
                unless ($range_rend > $range_lend) { 
                    confess "error, range_rend !> $range_lend ($range_rend, $range_lend)";
                }
                unless (&_overlaps([$range_lend, $range_rend], [ $lendA, $rendA])
                        &&  
                        &_overlaps([$range_lend, $range_rend], [ $lendB, $rendB] ) ) { 
                    next; 
                }
                if ($lendA < $range_lend) { $lendA = $range_lend; }
                if ($rendA > $range_rend) { $rendA = $range_rend; }
                if ($lendB < $range_lend) { $lendB = $range_lend; }
                if ($rendB > $range_rend) { $rendB = $range_rend; }
                
            }
            
            my $phaseA = "$phaseA_left";#,$phaseA_right";
            my $phaseB = "$phaseB_left"; #,$phaseB_right";
            
            my $id_A = $trellis->add_exon($lendA, $rendA, $typeA, $orientA, $phaseA);
            my $id_B = $trellis->add_exon($lendB, $rendB, $typeB, $orientB, $phaseB);
            
            if (/compatible/) {
                $trellis->add_compatible_connection($id_A, $id_B);
            }
            #if (/BESTSOFAR/) {
            #    $trellis->add_best_connection($id_A, $id_B);
            #}
        }
        else {
            print STDERR "trellis: no match.\n";
        }
    }
    close $fh;

    open ($fh, $final_path_file) or die $!;
    my $prev_id = "";

    my %seen; #track those already added to the trellis.  Don't show them again.

    while (<$fh>) {
        if (/trellis: 2/) { 
            # last; #only showing first trellis right now. 
        }
        chomp;
        if (/^(\w+)([\+\-])\~(\d)\s(\d+),\s(\d+)/) {
            my ($typeA, $orientA, $phaseA, $lendA, $rendA) = ($1, $2, $3, $4, $5);
            if (defined $range_lend && defined $range_rend) {
                unless (&_overlaps([$range_lend, $range_rend], [ $lendA, $rendA])) {
                        next; 
                }
                if ($lendA < $range_lend) { $lendA = $range_lend; }
                if ($rendA > $range_rend) { $rendA = $range_rend; }
            }
            
            my $id_A = $trellis->add_exon($lendA, $rendA, $typeA, $orientA, $phaseA);
            if ($prev_id && !$seen{$prev_id}) {
                $trellis->add_best_connection($id_A, $prev_id);
                $seen{$prev_id} = 1;
            }
            $prev_id = $id_A;
        }
        else {
            print STDERR "$final_path_file, no match $_\n";
        }
    }
    close $fh;
    

    return ($trellis);

}

sub _overlaps {
    my ($coordsA_aref, $coordsB_aref) = @_;

    my ($lendA, $rendA) = @$coordsA_aref;
    my ($lendB, $rendB) = @$coordsB_aref;

    if ($lendA <= $rendB && $rendA >= $lendB) {
        return (1);
    }
    else {
        return (0);
    }
}


1; #EOM
