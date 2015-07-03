package SequenceTickerIllustrator;

use strict;
use warnings;
use Carp;
use GD::SVG;

####
sub new {
    my $packagename = shift;
    my $constructor_params = shift;

    unless (ref $constructor_params) { 
        confess "Error, need constructor params href"; 
    }
    
    my $self = {
        seq_start => undef,
        seq_end => undef,
        canvas_position => undef,
        major_unit => undef,
        minor_unit => undef,
        label => "major",
        draw_label_flag => 1,
        color => "black",
        line_top_offset => 0.5, #default central
        minor_unit_height => 50, # height in pixels
        major_unit_height => 20,
        
    };


    foreach my $param (keys %$constructor_params) {
        if (exists $self->{$param}) {
            $self->{$param} = $constructor_params->{$param};
        }
        else {
            confess "error, $param is not an attribute";
        }
    }

    ## ensure minimial parameter settings:
    foreach my $required_att qw (seq_start  seq_end  canvas_position) {
        unless (exists $self->{$required_att}) {
            confess "Error, require setting of $required_att as an attribute in constructor";
        }
    }


    bless ($self, $packagename);
    return ($self);

}


####
sub draw {
    my $self = shift;
    my $image = shift;
    my $colors_href = shift;
    
    my $step = $self->{minor_unit};
    unless (defined $step) {
        $step = $self->{major_unit};
    }
    unless (defined $step) {
        confess "Error, step of minor or major unit is not defined for ticker";
    }

    my $start = $self->{seq_start};
    my $end = $self->{seq_end};
    unless ($start >= 0 && $end && $end > $start) {
        confess "Error, need start and end values and end must be greater than start: start: $start, end: $end";
    }

    ## unwrap the canvas rectangle:
    my $canvas_position = $self->{canvas_position};
    my ($x1, $y1, $x2, $y2) = ($canvas_position->{x1},
                               $canvas_position->{y1},
                               $canvas_position->{x2},
                               $canvas_position->{y2});
    
    unless ($x1 >= 0 && $y1 >= 0 && $x2 >= 0 && $y2 >= 0
            && $x2 >= $x1
            && $y2 >= $y1
            ) {
        confess "Error, rectangle is improper: ($x1, $y1) - ($x2, $y2)";
    }
    my $rect_height = $y2 - $y1;
    
    ## first, draw the central line:
    my $y_line = int ($rect_height * $self->{line_top_offset}) + $y1;
    
    my $ticker_color = $colors_href->{$self->{color}};
    
    $image->line($x1, $y_line, $x2, $y_line, $ticker_color);
    
    my $major_unit = $self->{major_unit};
    my $minor_unit = $self->{minor_unit};

    my $major_unit_height = $self->{major_unit_height};
    my $minor_unit_height = $self->{minor_unit_height};
    
    my $tick_label_type = $self->{label};

    ## draw ticks and add labels:
    for (my $i = $start; $i <= $end; $i += $step) {
        if ($major_unit && ($i - $start) % $major_unit == 0) {
            ## draw minor tick:
            my ($x_pos) = $self->_get_X_coords($i);
            $image->line($x_pos, $y_line, $x_pos, $y_line + $major_unit_height, $ticker_color);
            
            if ($i && $tick_label_type =~ /major/i) {
                # draw label:
                $image->string(gdSmallFont, $x_pos, $y_line + $major_unit_height, "$i", $ticker_color) unless $i == $end;
            }
            
        } elsif ( $minor_unit && ($i - $start) % $minor_unit == 0) {
            ## draw minor tick:
            my ($x_pos) = $self->_get_X_coords($i);
            $image->line($x_pos, $y_line, $x_pos, $y_line + $minor_unit_height, $ticker_color);
            if ($i && $tick_label_type =~ /minor/i) {
                # draw label:
                $image->string(gdSmallFont, $x_pos, $y_line + $minor_unit_height, "$i", $ticker_color) unless $i == $end;
            }
        }
    }

}


####
sub _get_X_coords {
    my $self = shift;
    my @coords = @_;

    my $canvas_position = $self->{canvas_position};
    my $seq_start = $self->{seq_start};
    my $seq_end = $self->{seq_end};
    
    my ($x1, $x2) = ($canvas_position->{x1}, $canvas_position->{x2});
    my $x_len = $x2 - $x1 + 1;
    
    my $seq_length = $seq_end - $seq_start + 1;
        
    my @ret_coords;
    foreach my $coord (@coords) {
        my $ratio_in_seq = ($coord - $seq_start) / $seq_length; 
        my $x_pos = int (($ratio_in_seq * $x_len)) + $x1;
        push (@ret_coords, $x_pos);
        # print STDERR "ticker: $coord -> $x_pos\n"; 
        
    }

    return (@ret_coords);
}



1; #EOM
        
    
    
