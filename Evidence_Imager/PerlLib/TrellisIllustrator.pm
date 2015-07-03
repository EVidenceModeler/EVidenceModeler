package TrellisIllustrator;

use strict;
use warnings;
use Carp;
use lib ($ENV{EUK_MODULES}); ## ToDo, include other mods in package
use Overlap_piler;
use EvidenceFeature;
use Data::Dumper;

####
sub new {
    my $packagename = shift;
    my $canvas_params = shift;
    
    unless (ref $canvas_params) { confess "need canvas params in constructor"; }
    
    my $self = {
        canvas_params => $canvas_params,
    };

    bless ($self, $packagename);

    return ($self);
}

####
sub draw_trellis {
    my $self = shift;
    my ($trellis, $image, $colors_href, $seq_start, $seq_end, $genome_seq_ref) = @_;

    my $genome_seq_length = length($$genome_seq_ref);

    my $canvas_params = $self->{canvas_params};

    my ($min_x, $max_x, $min_y, $max_y) = ($canvas_params->{x1}, $canvas_params->{x2}, 
                                           $canvas_params->{y1}, $canvas_params->{y2});    

    ## draw line separators:
    my $canvas_height = $max_y - $min_y + 1;
    my $canvas_width = $max_x - $min_x + 1;
    my $height_per_frame = int($canvas_height / 6);
    
    my $x_coord_range = $seq_end - $seq_start + 1;
    
    ####
    my $xcoord_transform = sub {
        my @coords = @_;
        my @ret_coords;
        foreach my $coord (@coords) {
            my $delta = $coord - $seq_start;
            my $ratio = $delta / $x_coord_range;
            
            my $x_coord = int ($ratio * $canvas_width) + $min_x;
            #print STDERR "coord $coord -> $x_coord, delta: $delta, $seq_start\n";
            
            push (@ret_coords, $x_coord);
        }
        return (@ret_coords);
    };


    #-------------------------------------
    ## ToDo: params to add to the constructor:
    my $ratio_feature_height = 0.75;
    my $ratio_start_height = 0.3;
    my $ratio_stop_height = 0.5;
    my $start_color = "green";
    my $stop_color = "red";
    my $exon_color = "blue";
    
    #------------------------------------
    
    ## get the exon list:
    my @exons = $trellis->get_exon_list();
    my %exon_ID_to_coords;
    
    &_add_frame_value_to_exons(\@exons, $genome_seq_length);
    ## draw each frame:
    foreach my $frame (1..6) {
        my @framed_rects;
        my $frame_y_pos_bottom = $frame * $height_per_frame + $min_y;
        my $frame_y_pos_top = ($frame-1) * $height_per_frame + $min_y;
        my @framed_exons = &_get_exons_in_frame($frame, \@exons);
        
        my @exon_clusters = &cluster_and_tier_exons(@framed_exons);
        
        foreach my $exon_cluster (@exon_clusters) {
            my $num_exons_in_cluster = scalar @$exon_cluster;
            my $height_per_exon = int($height_per_frame / $num_exons_in_cluster);
            
            my $subcluster_count = 0;
            foreach my $exon_subcluster (@$exon_cluster) {
                $subcluster_count++;
                my $y_pos_top = $frame_y_pos_top + ($subcluster_count  - 1) * $height_per_exon;
                my $y_pos_bottom = $y_pos_top + $height_per_exon;
                
                foreach my $exon (@$exon_subcluster) {
                    
                    

                    my $exon_ID = $exon->{ID};
                    my ($exon_lend, $exon_rend) = ($exon->{lend}, $exon->{rend});
                    if ($exon_lend > $exon_rend) {
                        confess "Error, exon coordinates in reverse lend-rend order: $exon_lend, $exon_rend\n";
                    }
                    ($exon_lend, $exon_rend) = &$xcoord_transform($exon_lend, $exon_rend);
                    
                    ## draw exon rect:
                    $image->filledRectangle($exon_lend, $y_pos_top, $exon_rend, $y_pos_bottom, $colors_href->{blue});
                    #$image->rectangle($exon_lend, $y_pos_top, $exon_rend, $y_pos_bottom, $colors_href->{black});
                    push (@framed_rects, [$exon_lend, $y_pos_top, $exon_rend, $y_pos_bottom, $colors_href->{black}]);
                    
                    $exon_ID_to_coords{$exon_ID} = "$exon_lend $y_pos_top $exon_rend $y_pos_bottom";
                }
            }
        }
        ## frame the rects:
        foreach my $rect (@framed_rects) {
            $image->rectangle(@$rect);
        }
        
        my @start_codons = &find_start_codons($frame, $seq_start, $seq_end, $genome_seq_ref);
        my @stop_codons = &find_stop_codons($frame, $seq_start, $seq_end, $genome_seq_ref);
        foreach my $stop_codon (@stop_codons) {
            my ($codon_pos) = &$xcoord_transform($stop_codon);
            $image->line($codon_pos, $frame_y_pos_bottom, $codon_pos, $frame_y_pos_bottom - int($ratio_stop_height * $height_per_frame), $colors_href->{red});
        }
        foreach my $start_codon (@start_codons) {
            my ($codon_pos) = &$xcoord_transform($start_codon);
            $image->line($codon_pos, $frame_y_pos_bottom, $codon_pos, $frame_y_pos_bottom - int($ratio_start_height * $height_per_frame), $colors_href->{green});
        }
        ## draw the baseline
        $image->line($min_x, $frame_y_pos_bottom, $max_x, $frame_y_pos_bottom, $colors_href->{black});
    }
    
    ## draws connections between the exons:
    my %compatible_connections = $trellis->get_compatible_connections();
    my %best_connections = $trellis->get_best_connections();

    foreach my $id_A (keys %compatible_connections) {
        
        my $exon_A = $trellis->get_exon_by_ID($id_A);
        my $exon_A_type = $exon_A->{type};
        #if ($exon_A_type eq 'bound') { next; }
        my $exon_A_coords = $exon_ID_to_coords{$id_A} or confess "error, no coords for ID: $id_A";
        my ($a_x1, $a_y1, $a_x2, $a_y2) = split (/ /, $exon_A_coords);

        my $midpt_a = int ( ($a_y1 + $a_y2) / 2);

        my $id_A_connections_href = $compatible_connections{$id_A};
        foreach my $id_B (keys %$id_A_connections_href) {
            my $exon_B = $trellis->get_exon_by_ID($id_B);
            my $exon_B_coords = $exon_ID_to_coords{$id_B} or confess "error, no coords for ID: $id_B";
            my ($b_x1, $b_y1, $b_x2, $b_y2) = split (/ /, $exon_B_coords);
            
            my $midpt_b = int ( ($b_y1 + $b_y2) / 2);
            if ($best_connections{$id_A}->{$id_B}) {
                $image->line($a_x1, $midpt_a, $b_x2, $midpt_b, $colors_href->{black});
            } else {
                #$image->line($a_x1, $midpt_a, $b_x2, $midpt_b, $colors_href->{cyan});
            }
        }
    }

    return;
}

####
sub _add_frame_value_to_exons {
    my ($exons_aref, $genome_seq_length) = @_;
    
    foreach my $exon (@$exons_aref) {
        my ($lend, $rend, $type, $orient, $phase) = ($exon->{lend},
                                                     $exon->{rend},
                                                     $exon->{type},
                                                     $exon->{orient},
                                                     $exon->{phase});
        if ($orient eq '+') {
            $phase--;
            $lend -= $phase;
            my $frame = $lend % 3;
            $frame = 3 if $frame == 0;
            $exon->{frame} = $frame;
        }
        else {
            # minus strand
            $phase--;
            $rend += $phase;
            my $frame = ($genome_seq_length - $rend + 1) % 3;
            $frame = 3 if $frame == 0;
            
            ## convert to 4,5,6
            if ($frame == 1) { $frame = 4; }
            elsif ($frame == 2) { $frame = 5; }
            elsif ($frame == 3) { $frame = 6; }
            $exon->{frame} = $frame;
        }
    }
    return;
}

####
sub _get_exons_in_frame {
    my ($frame, $exons_aref) = @_;
    
    my @exons;
    foreach my $exon (@$exons_aref) {
        if ($exon->{frame} == $frame) {
            push (@exons, $exon);
        }
    }
    return (@exons);
}

####
sub find_start_codons {
    my ($frame, $seq_start, $seq_end, $genome_seq_ref) = @_;
    
    my $forward_start_codon = "ATG";
    my $reverse_start_codon = "CAT";

    # array coords instead of seq coords
    $seq_start--;
    if ($seq_start < 0) { $seq_start = 0;}
    $seq_end--;

    my @start_codons;
    
    if ($frame <= 3) {
        ## forward strand:
        my $pos = index($$genome_seq_ref, $forward_start_codon, $seq_start);
        while ($pos != -1 && $pos <= $seq_end) {
            # print STDERR "F$frame: forward start search, found pos: $pos\n";
            my $pos_frame = ($pos + 1) % 3;
            $pos_frame = 3 if $pos_frame == 0;
            if ($pos_frame == $frame) {
                push (@start_codons, $pos);
            }
            $pos = index($$genome_seq_ref, $forward_start_codon, $pos + 1);
        }
    }
    else {
        # minus strand:
        my $genome_seq_length = length($$genome_seq_ref);
        my $pos = index($$genome_seq_ref, $reverse_start_codon, $seq_start);
        while ($pos != -1 && $pos <= $seq_end) {
            # print STDERR "F$frame: reverse start search, found pos: $pos\n";
            my $pos_frame = ($genome_seq_length - $pos + 2 + 1) % 3;
            if ($pos_frame == 0) { $pos_frame = 4; }
            elsif ($pos_frame == 1) { $pos_frame = 5; }
            elsif ($pos_frame == 2) { $pos_frame = 6; }
            if ($pos_frame == $frame) {
                push (@start_codons, $pos);
            }
            $pos = index($$genome_seq_ref, $reverse_start_codon, $pos+1);
        }
    }

    return (@start_codons);
}
            
####
sub find_stop_codons {
    my ($frame, $seq_start, $seq_end, $genome_seq_ref) = @_;
    
    # array coords instead of seq coords
    $seq_start--;
    if ($seq_start < 0) { $seq_start = 0;}
    $seq_end--; 

    my @forward_stop_codons = ("TAA", "TGA", "TAG");
    my @reverse_stop_codons = ("CTA", "TCA", "TTA");
    
    my @stop_codons;
    
    if ($frame <= 3) {
        ## forward strand:
        foreach my $forward_stop_codon (@forward_stop_codons) {
            my $pos = index($$genome_seq_ref, $forward_stop_codon, $seq_start);
            while ($pos != -1 && $pos <= $seq_end) {
                # print STDERR "F$frame: forward stop search, found pos: $pos\n";
                my $pos_frame = ($pos + 1) % 3;
                $pos_frame = 3 if $pos_frame == 0;
                if ($pos_frame == $frame) {
                    push (@stop_codons, $pos);
                }
                $pos = index($$genome_seq_ref, $forward_stop_codon, $pos + 1);
            }
        }
    } else {
        # minus strand:
        my $genome_seq_length = length($$genome_seq_ref);
        foreach my $stop_codon (@reverse_stop_codons) {
            my $pos = index($$genome_seq_ref, $stop_codon, $seq_start);
            while ($pos != -1 && $pos <= $seq_end) {
                # print STDERR "F$frame: reverse stop search, found pos: $pos\n";
                my $pos_frame = ($genome_seq_length - $pos + 2 + 1) % 3;
                if ($pos_frame == 0) { $pos_frame = 4; }
                elsif ($pos_frame == 1) { $pos_frame = 5; }
                elsif ($pos_frame == 2) { $pos_frame = 6; }
                if ($pos_frame == $frame) {
                    push (@stop_codons, $pos);
                }
                $pos = index($$genome_seq_ref, $stop_codon, $pos+1);
            }
        }
    }

    return (@stop_codons);
}


####
sub cluster_and_tier_exons {
    my @framed_exons = @_;
    
    unless (@framed_exons) { return (); }

    my $overlap_piler = new Overlap_piler();
    
    my %ID_to_exon;
    foreach my $exon (@framed_exons) {
        my ($ID, $lend, $rend) = ($exon->{ID}, $exon->{lend}, $exon->{rend});
        
        $ID_to_exon{$ID} = $exon;
        
        $overlap_piler->add_coordSet($ID, $lend, $rend);
    }
    
    my @exon_clusters;
    
    my @ID_clusters = $overlap_piler->build_clusters();
       
    foreach my $ID_aref (@ID_clusters) {
        my @IDs = @$ID_aref;
        
        if (scalar @IDs == 1) {
            my $id = shift @IDs;
            my $exon = $ID_to_exon{$id};
            push (@exon_clusters, [ [$exon] ]);
        }
        else {
            ## tier them
            my @features;
            foreach my $ID (@IDs) {
                my $exon = $ID_to_exon{$ID};
                my ($lend, $rend) = ($exon->{lend}, $exon->{rend});
                my $feature = new EvidenceFeature($ID, "nulltype", "nullsource");
                $feature->add_segment($lend, $rend, '+'); #keep all plus strand for ease
                push (@features, $feature);
                
            }
            my $evidenceTier = new EvidenceTier();
            $evidenceTier->tier_features(@features);
            my @tiers = $evidenceTier->get_top_strand_rows();
            
            my @exon_tiers;
            foreach my $row (@tiers) {
                my @features = @$row;
                my @exon_row;
                foreach my $feature (@features) {
                    my $ID = $feature->{ID};
                    my $exon = $ID_to_exon{$ID};
                    push (@exon_row, $exon);
                }
                push (@exon_tiers, [@exon_row]);
            }
            push (@exon_clusters, [ @exon_tiers ]);
        }

    }

    return (@exon_clusters);
}





1; #EOM
