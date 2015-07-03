package VectorGraph;

use strict;
use warnings;
use Carp;
use GD::SVG;

####
sub new {
    my $packagename = shift;
    my $params_href = shift;
    
    unless (ref $params_href) {
        confess "need canvas params in constructor";
    }
    
    my $self = { 
        vectors => [], #list of data vectors
        min_x => undef,
        max_x => undef,
        min_y => undef,
        max_y => undef,

        canvas_params => undef,
        
        x_ticker_settings => undef,
        y_ticker_settings => undef,
        
    };
    
    if (ref $params_href eq "HASH") {
        foreach my $param (keys %$params_href) {
            if (exists $self->{$param}) {
                $self->{$param} = $params_href->{$param};
            }
            else {
                confess "Error, attribute $param does not exist in object";
            }
        }
    }
    
    my @required_atts = qw (canvas_params);
    foreach my $att (@required_atts) {
        unless ($self->{$att}) {
            confess "required attribute: $att is not set in constructor.";
        }
    }
    
    bless ($self, $packagename);
    
    return ($self);
}

####
sub add_vector {
    my $self = shift;
    my ($vector_name, $color, $points_aref) = @_;

    my $data_vector = VectorGraph::DataVector->new($vector_name, $color, $points_aref);
    
    if (! defined $self->{min_x}) {
        $self->{min_x} = $data_vector->{min_x};
        $self->{max_x} = $data_vector->{max_x};
        $self->{min_y} = $data_vector->{min_y};
        $self->{max_y} = $data_vector->{max_y};
    }
    else {
        # update min/max values:
        if ($data_vector->{min_x} < $self->{min_x}) { $self->{min_x} = $data_vector->{min_x};}
        if ($data_vector->{max_x} > $self->{max_x}) { $self->{max_x} = $data_vector->{max_x};}
        if ($data_vector->{min_y} < $self->{min_y}) { $self->{min_y} = $data_vector->{min_y};}
        if ($data_vector->{max_y} > $self->{max_y}) { $self->{max_y} = $data_vector->{max_y};}
    }

    push (@{$self->{vectors}}, $data_vector);
    
    return;
}

####
sub get_vectors {
    my $self = shift;
    return (@{$self->{vectors}});
}

####
sub draw_graph {
    my $self = shift;
    my $image = shift;
    my $colors_href = shift;
    my $axis_settings_href = shift;
    
    
    my $canvas_params = $self->{canvas_params};
    my ($x1, $y1, $x2, $y2) = ($canvas_params->{x1},
                               $canvas_params->{y1},
                               $canvas_params->{x2},
                               $canvas_params->{y2});
    
    
    my $x_pixels_len = $x2 - $x1 + 1;
    my $y_pixels_len = $y2 - $y1 + 1;
    
    my ($min_x, $max_x) = ($self->{min_x}, $self->{max_x});
    my ($min_y, $max_y) = ($self->{min_y}, $self->{max_y});

    if (ref $axis_settings_href) {
        if (defined ($axis_settings_href->{y_axis_max}) ) {
            $max_y = $axis_settings_href->{y_axis_max};
            if ($min_y > $max_y) { $min_y = $max_y; }
        }
        if (defined ($axis_settings_href->{y_axis_min})) {
            $min_y = $axis_settings_href->{y_axis_min};
            if ($max_y < $min_y) { $max_y = $min_y; }
        }
    }
    
    my $x_val_span = $max_x - $min_x;
    my $y_val_span = $max_y - $min_y;
    

    my $x_coord_func = sub { 
        my $x_coord = shift;
        my $offset_x = $x_coord - $min_x;
        my $ratio_x = $offset_x / $x_val_span;
        my $offset_pixels_x = int($ratio_x * $x_pixels_len);
        my $x_pixel_pos = $x1 + $offset_pixels_x;
        return ($x_pixel_pos);
    };

    my $y_coord_func = sub {
        my $y_coord = shift;
        
        if ($y_coord > $max_y) {
            $y_coord = $max_y; #positive peaked!
        }
        if ($y_coord < $min_y) {
            $y_coord = $min_y; # negative peaked!
        }

        my $offset_y = $y_coord - $min_y;
        my $ratio_y = $offset_y / $y_val_span;
        my $offset_pixels_y = int($ratio_y * $y_pixels_len);
        my $y_pixel_pos = $y2 - $offset_pixels_y;
        return ($y_pixel_pos);
    };
        

    ## connect the dots for each data set:
    print STDERR "connecting the dots.\n";
    foreach my $vector ($self->get_vectors()) {
        my $color_name = $vector->{color};
        
        my @data_points = $vector->get_data_points();
        
        my $prev_point = shift @data_points;
        while (@data_points) {
            my $next_point = shift @data_points;
            
            ## do coordinate transformations:
            my ($cx1, $cy1) = @$prev_point;
            my ($cx2, $cy2) = @$next_point;
            
            $cx1 = &$x_coord_func($cx1);
            $cx2 = &$x_coord_func($cx2);
            $cy1 = &$y_coord_func($cy1);
            $cy2 = &$y_coord_func($cy2);
            
            # print STDERR "graph line: ($cx1,$cy1) - ($cx2,$cy2)\n";
            
            $image->line($cx1, $cy1, $cx2, $cy2, $colors_href->{$color_name});
            

            $prev_point = $next_point;
        }
    }
    
    ## draw the X and Y axis

    ## x-axis as a central black line
    my $y_zero_pos = &$y_coord_func(0);
    $image->line($x1, $y_zero_pos, $x2, $y_zero_pos, $colors_href->{black});
        
    ## y-axis, add ticks for zero, min-y and max-y
    my $max_y_pos = &$y_coord_func($max_y);
    my $min_y_pos = &$y_coord_func($min_y);
    $image->line($x1, $max_y_pos, $x1, $min_y_pos, $colors_href->{black});
    
    my $tick_size = 4;
    foreach my $y_pos_info( [$max_y_pos, $max_y], 
                            [$y_zero_pos, 0],
                            [$min_y_pos, $min_y]) {
        my $x_pos = $x1 - $tick_size;
        my ($y_pos, $y_value) = @$y_pos_info;
        $image->line($x_pos, $y_pos, $x1, $y_pos, $colors_href->{black});
        my $string_length = length("$y_value");
        my $string_pixel_length = $string_length * 4;
        $image->string(gdTinyFont, $x_pos - $string_pixel_length - 3, $y_pos - 2, "$y_value", $colors_href->{black});
        
    }
    
    
    

    return;
}





##########################################################################
##########################################################################

package VectorGraph::DataVector;

use strict;
use warnings;
use Carp;

####
sub new {
    my $packagename = shift;
    
    my ($data_name, $color, $points_aref) = @_;

    my $self = {
        
        data_name => $data_name,
        color => $color,

        data_points => [ @$points_aref],
        min_x => undef,
        max_x => undef,
        min_y => undef,
        max_y => undef,
        
    };

    bless ($self, $packagename);
    
    $self->_set_min_max_vals();

    return ($self);
}

####
sub get_data_points {
    my $self = shift;
    return (@{$self->{data_points}});
}


####
sub _set_min_max_vals {
    my $self = shift;
    
    my $min_x = undef;
    my $max_x = undef;
    my $min_y = undef;
    my $max_y = undef;

    foreach my $point (@{$self->{data_points}}) {
        my ($x, $y) = @$point;
        if (!defined ($min_x)) {
            $min_x = $max_x = $x;
            $min_y = $max_y = $y;
        }
        else {
            if ($x < $min_x) { $min_x = $x; }
            if ($x > $max_x) { $max_x = $x; }
            if ($y < $min_y) { $min_y = $y; }
            if ($y > $max_y) { $max_y = $y; }
        }
    }

    $self->{min_x} = $min_x;
    $self->{max_x} = $max_x;
    $self->{min_y} = $min_y;
    $self->{max_y} = $max_y;

    return;
}


1; #EOM
    

