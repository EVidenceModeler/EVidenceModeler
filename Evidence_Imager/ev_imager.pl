#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/PerlLib");
use GD::SVG;
use Carp;
use EvidenceTierParser;
use EvidenceTier;
use SequenceTickerIllustrator;
use Config::IniFiles;
use VectorGraph;
use TrellisParser;
use Trellis;
use TrellisIllustrator;


my $usage = "usage: $0 conf_file start-end\n\n";
my $conf_file = $ARGV[0] or die $usage;
my $coords = $ARGV[1];
my ($range_start, $range_end);
if ($coords) {
    ($range_start, $range_end) = split (/-/, $coords);
}

my %conf = &parse_conf($conf_file);
my $cfg = new Config::IniFiles( -file => $conf_file);


my $evidence_tiers_list = $conf{EVIDENCE_TIERS} or confess "error, no evidence tiers listed in conf file";

my $genome_sequence_file = $conf{GENOME_SEQUENCE};
my $genome_sequence = &parse_genome_seq_file($genome_sequence_file);
my $genome_seq_length = length($genome_sequence);
unless (defined $range_start) {
    ($range_start, $range_end) = (0, $genome_seq_length);
}
unless (defined $range_start && defined $range_end) {
    confess "Error, range start and end is not set";
}

## Colors:
my %colors;

## canvas settings:
my $LEFT_MARGIN = $conf{LEFT_MARGIN} || 0;
my $RIGHT_MARGIN = $conf{RIGHT_MARGIN} || 0;
my $TOP_MARGIN = $conf{TOP_MARGIN} || 0;
my $BOTTOM_MARGIN = $conf{BOTTOM_MARGIN} || 0;

my $range_length = $range_end - $range_start + 1;

my %feature_type_to_RGB;

main: {

    my $total_tiers = 0;
    
    my @ev_tier_objs;
    
    if ($conf{DRAW_GENOME_VIEW} =~ /true/i) {
        foreach my $ev_tier_file (split (/,/, $evidence_tiers_list)) {
            
            my $ev_tier_parser = new EvidenceTierParser( {max_rows_per_tier => $conf{MAX_ROWS_PER_TIER},
                                                          range_start => $range_start,
                                                          range_end => $range_end,
                                                         } );
            
            my @tiers = $ev_tier_parser->parse_evidence_tiers($ev_tier_file);
            
            push (@ev_tier_objs, @tiers);
            
            #print $ev_tier_obj->toString();
            
        }
        
        ## count tiers:
        $total_tiers = &count_tiers(@ev_tier_objs);
    }
    

    ## Determine drawing dimensions:
    my $genome_view_height = $conf{PIXELS_ROW_HEIGHT} * $total_tiers;
    my $viewer_height = $genome_view_height + $TOP_MARGIN + $BOTTOM_MARGIN;
    my $viewer_width = int($conf{PIXELS_PER_KB_WIDTH} * $range_length / 1000) + $LEFT_MARGIN + $RIGHT_MARGIN;
    
    my @graphs = &parse_graphs();
    ## add graph to viewer height:
    if (@graphs) {
        my $graph_top_margin = $cfg->val("graphs", "GRAPH_TOP_MARGIN") or confess "error, no margin set for graph_top_margin";
        foreach my $graph (@graphs) {
            my $graph_height = $cfg->val("$graph", "GRAPH_HEIGHT") or confess "error, no graph height specified for $graph";
            $viewer_height += $graph_height + $graph_top_margin;
        }
    }
    
    if ($conf{SEQ_TICKER_DRAW} =~ /true/i) {
        $viewer_height += $conf{SEQ_TICKER_HEIGHT};
        $genome_view_height += $conf{SEQ_TICKER_HEIGHT};
    }
    
    if ($conf{DRAW_TRELLIS} =~ /true/i) {
        $viewer_height += $conf{TRELLIS_TOP_MARGIN} + $conf{TRELLIS_HEIGHT};
    }
    

    my $image = GD::SVG::Image->new($viewer_width, $viewer_height);
    &init_colors($image);

    
    ## draw from the top down
    
    my ($min_x, $max_x) = &convert_Xcoords($range_start, $range_end);
    # print STDERR "min_x: $min_x, max_x: $max_x\n";
    

    my $curr_y = $TOP_MARGIN;
    
    my $y_value = $curr_y;
    
    ## generate image:
    if ($conf{DRAW_GENOME_VIEW} =~ /true/i) {
        $y_value = &create_genome_view($image, \@ev_tier_objs, $min_x, $max_x, $curr_y, $genome_view_height);
    }
    
    foreach my $graph (@graphs) {
        ## determine rectangle for graph to be drawn in.  Ticker's drawn on the outer boundaries! Must have margins.
        my $graph_top_margin = $cfg->val("graphs", "GRAPH_TOP_MARGIN");
        $y_value += $graph_top_margin;
        my $graph_height = $cfg->val("$graph", "GRAPH_HEIGHT") or die "Error, no graph height for $graph in conf file";
        my $coord_struct = { x1 => $min_x, x2 => $max_x, y1 => $y_value, y2 => $y_value + $graph_height };
        
        &create_graph($image, $graph, $coord_struct);
        
        $y_value += $graph_height;
    }
    
    if ($conf{DRAW_TRELLIS} =~ /true/i) {
        $y_value += $conf{TRELLIS_TOP_MARGIN};
        
        my $canvas_params = { x1 => $min_x, x2 => $max_x, y1 => $y_value, y2 => $y_value + $conf{TRELLIS_HEIGHT} };
        &draw_trellis($image, $canvas_params);
        

    }




    my $svg_xml = $image->svg();
    
    $svg_xml =~ s/font=/font-face=/g;

    print $svg_xml;
    

}


exit(0);


####
sub count_tiers {
    my @tiers = @_;
    my $total_tiers = 0;
    
    foreach my $ev_tier_obj (@tiers) {
        my @top_rows = $ev_tier_obj->get_top_strand_rows();
        if (@top_rows) {
            $total_tiers += scalar(@top_rows);
        }
        my @bottom_rows = $ev_tier_obj->get_bottom_strand_rows();
        if (@bottom_rows) {
            $total_tiers += scalar (@bottom_rows);
        }
    }
    
    return ($total_tiers);
}
    

####
sub parse_graphs {
    my $graph_list = $conf{GRAPH_LIST};
    $graph_list =~ s/^\s+|\s+$//g; #trim surrounding whitespace
    my @graphs = split (/,/, $graph_list);
    foreach my $graph (@graphs) {
        $graph =~ s/\s//g;
    }

    return (@graphs);
}



####
sub parse_conf {
    my ($conf_file) = @_;
    
    my %conf;
    open (my $fh, $conf_file) or die $!;
    while (<$fh>) {
        chomp;
        unless (/\w/) { next;}
        if (/^\#/) { next; }
        unless (/=/) {
            print STDERR "line lacks key/value pair: $_\n";
            next;
        }
        s/^\s+|\s+$//g; # trim leading/trailing whitespace
        
        my ($key, $value) = split (/\s*=\s*/);
        $conf{$key} = $value;
    }
    close $fh;
    return (%conf);
}


####
sub create_genome_view {
    my ($image, $ev_tier_objs_aref, $min_x, $max_x, $curr_y, $viewer_height) = @_;
    
    ## make the genome view a black background:
    $image->rectangle($min_x, $curr_y, $max_x, $curr_y + $viewer_height, $colors{black});
    

    my @top_tiers;
    my @bottom_tiers;
    
    foreach my $ev_tier_obj (@$ev_tier_objs_aref) {
        my @top_rows = $ev_tier_obj->get_top_strand_rows();
        if (@top_rows) {
            unshift (@top_tiers, [reverse @top_rows]);
        }
        my @bottom_rows = $ev_tier_obj->get_bottom_strand_rows();
        if (@bottom_rows) {
            push (@bottom_tiers, [@bottom_rows]);
        }
    }

    my @labels_to_process;

    my $tier_num = 0;
    foreach my $tier (@top_tiers, "ticker", @bottom_tiers) {
        
        if ($tier eq "ticker") {
            if ($conf{SEQ_TICKER_DRAW} =~ /true/i) {
                my $seq_ticker_height = $conf{SEQ_TICKER_HEIGHT};
                &draw_ticker($image, $min_x, $curr_y, $max_x, $curr_y + $seq_ticker_height);
                $curr_y += $seq_ticker_height;
            };
            next;
            
        }
        
        $tier_num++;
        my $num_rows_in_tier = scalar @$tier;
        my $bottom_tier_y = $curr_y + $num_rows_in_tier * $conf{PIXELS_ROW_HEIGHT};
        
        ## add gray lighting to separate odd rows
        if ($tier_num % 2) {
            # $image->filledRectangle($min_x, $curr_y, $max_x, $bottom_tier_y, &get_color($image, (200,200,200)));
        }
        

        my $feature_source;
        my $top_tier_y = $curr_y;


        foreach my $row (@$tier) {
            ## determine y coordinate:
            my $top_y = $curr_y;
            my $bottom_y = $curr_y + $conf{PIXELS_ROW_HEIGHT};
            
            
            foreach my $feature (@$row) {
                &draw_feature($image, $feature, $top_y, $bottom_y);
                unless ($feature_source) {
                    $feature_source = $feature->{feature_source};
                }        
            }
                        
            
            $curr_y += $conf{PIXELS_ROW_HEIGHT};
            # on to the next row
        }

        
        ## add the feature type label to the left of the viewer
        my $source_color = &get_type_color($image, $feature_source);
        
        $feature_source = substr($feature_source, 0, 15);
        
        my $mid_y = int ( ($top_tier_y + $bottom_tier_y) / 2);
                
        push (@labels_to_process, { label => $feature_source,
                                    mid_y => $mid_y,
                                    color => $source_color, 
                                } );
        

    }
    
    ## draw the labels:
    my $max_string_length = 0;
    foreach my $label (@labels_to_process) {
        my $string_length = length($label->{label});
        if ($string_length > $max_string_length) {
            $max_string_length = $string_length;
        }
    }
    
    my $string_pixel_length = $max_string_length * 5 + 2;
    my $offset_x = $string_pixel_length * 3;
    foreach my $label (@labels_to_process) {
        my $label_text = $label->{label};
        my $mid_y = $label->{mid_y};
        my $source_color = $label->{color};
        $image->string(gdTinyFont, $min_x - $string_pixel_length, $mid_y, $label_text, $source_color);
    }
    
    return ($curr_y);
    
}


####
sub draw_feature {
    my ($image, $feature, $top_y, $bottom_y) = @_;
    
    ## can be as creative as needed here.
    
    my $feature_type = $feature->{feature_type};
    my $ratio_feature_height = $conf{RATIO_FEATURE_HEIGHT};
    my $feature_height = $bottom_y - $top_y;
    my $offset = int( ($feature_height - ($ratio_feature_height * $feature_height)) / 2);

    $top_y += $offset;
    $bottom_y -= $offset;
    

    if ($feature_type eq "gene") {
        &draw_gene($image, $feature, $top_y, $bottom_y);
    }
    else {
        &draw_generic_feature($image, $feature, $top_y, $bottom_y);
    }
 
    return;
}
    

####
sub draw_gene {
    my ($image, $feature, $top_y, $bottom_y) = @_;
    
    ## draw the central line:
    my ($lend, $rend) = ($feature->{lend}, $feature->{rend});
    my ($x_lend, $x_rend) = &convert_Xcoords($lend, $rend);
    
    my $feature_source = $feature->{feature_source};
    my $type_color = &get_type_color($image, $feature_source);

    ## draw a central line:
    my $midpt = int (($top_y + $bottom_y) / 2);
    
    $image->line($x_lend, $midpt, $x_rend, $midpt, $colors{black});
    

    my @exons = $feature->retrieve_descendants("exon");
    my @cds = $feature->retrieve_descendants("CDS");

    foreach my $structure_set ( [\@exons, $colors{black}, "exon"],
                                [\@cds, $type_color, "CDS"] ) {
        
        my ($features_aref, $color, $type) = @$structure_set;
        foreach my $feature (@$features_aref) {
            # print STDERR $feature->toString();
            foreach my $segment ($feature->get_segments()) {
                my ($lend, $rend) = @$segment;
                my ($x_lend, $x_rend) = &convert_Xcoords($lend, $rend);
                
                $image->filledRectangle($x_lend, $top_y, $x_rend, $bottom_y, $color);
            
                # outline the CDS regions:
                $image->rectangle($x_lend, $top_y, $x_rend, $bottom_y, $colors{black});
            }
        }
    }

    return;


}

####
sub draw_generic_feature {
    my ($image, $feature, $top_y, $bottom_y) = @_;

    my ($lend, $rend) = ($feature->{lend}, $feature->{rend});
    my ($x_lend, $x_rend) = &convert_Xcoords($lend, $rend);
    
    my $feature_source = $feature->{feature_source};
    my $type_color = &get_type_color($image, $feature_source);
    
    ## draw a central line:
    my $midpt = int (($top_y + $bottom_y) / 2);
    
    $image->line($x_lend, $midpt, $x_rend, $midpt, $colors{black});
    
    foreach my $segment ($feature->get_segments()) {
        my ($lend, $rend) = @$segment;
        my ($x_lend, $x_rend) = &convert_Xcoords($lend, $rend);
        
        ## filled rect for segments
        $image->filledRectangle($x_lend, $top_y, $x_rend, $bottom_y, $type_color);
        
        ## outline the segments:
        $image->rectangle($x_lend, $top_y, $x_rend, $bottom_y, $type_color);
        
    }
    
    return;
}


sub draw_ticker {
    my ($image, $x1, $y1, $x2, $y2) = @_;
    
    
    # color the canvas to separate the ticker from the rest of the view
    $image->filledRectangle($x1, $y1, $x2, $y2, &get_color($image, (200,200,200)));
                            


    my $ticker_illustrator = new SequenceTickerIllustrator ({ seq_start => $range_start,
                                                              seq_end => $range_end,
                                                              ## position in canvas
                                                              canvas_position => { x1 => $x1,
                                                                                   x2 => $x2,
                                                                                   y1 => $y1,
                                                                                   y2 => $y2},
                                                              
                                                              # tick positions
                                                              major_unit => $conf{SEQ_TICKER_MAJOR_UNIT},
                                                              minor_unit => $conf{SEQ_TICKER_MINOR_UNIT},

                                                              # tick heights
                                                              major_unit_height => $conf{SEQ_TICKER_MAJOR_UNIT_HEIGHT},
                                                              minor_unit_height => $conf{SEQ_TICKER_MINOR_UNIT_HEIGHT},
                                                              
                                                              label => $conf{SEQ_TICKER_LABEL_COORDINATE},
                                                              draw_label_flag => $conf{SEQ_TICKER_DRAW_LABEL},
                                                              line_top_offset => $conf{SEQ_TICKER_LINE_TOP_OFFSET_RATIO},
                                                              color => $conf{SEQ_TICKER_COLOR},
                                                              
                                                          });
    

    $ticker_illustrator->draw($image, \%colors);
    
    return;
}


####
sub convert_Xcoords {
    my @coords = @_;
    
    my @ret_coords;
    foreach my $coord (@coords) {
        
        if ($coord < $range_start || $coord > $range_end) {
            confess "Error, coordinate out of range: $coord, range: $range_start-$range_end";
        }
        
        my $x_coord = int (($coord - $range_start) * $conf{PIXELS_PER_KB_WIDTH} / 1000) + $LEFT_MARGIN;
        if ($x_coord < 0) {
            confess "Error, conversion of $coord => $x_coord negative";
        }
        # print STDERR "X: $coord -> $x_coord\n"; 
        push (@ret_coords, $x_coord);
    }

    return (@ret_coords);
}


####
sub init_colors {
    my $image = shift;
    $colors{blue} = $colors{"0:0:255"} = $image->colorAllocate(0,0,255);
    $colors{black} = $colors{"0:0:0"} = $image->colorAllocate(0,0,0);
    $colors{red} = $colors{"255:0:0"} = $image->colorAllocate(255,0,0);
    $colors{green} = $colors{"0:255:0"} = $image->colorAllocate(0,255,0);
    $colors{yellow} = $colors{"255:255:0"} = $image->colorAllocate(255,255,0);
    $colors{magenta} = $colors{"255:0:255"} = $image->colorAllocate(255,0,255);
    $colors{cyan} = $colors{"0:255:255"} = $image->colorAllocate(0,255,255);
    
    return;
}

####
sub get_color {
    my $image = shift;
    my ($r, $g, $b) = @_;
    
    my $color_token = "$r:$g:$b";
    if (my $color = $colors{$color_token}) {
        return ($color);
    }
    else {
        my $color = $colors{$color_token} = $image->colorAllocate($r, $g, $b);
        return ($color);
    }
}

####
sub parse_genome_seq_file {
    my $genome_seq_file = shift;
    
    my $genome_seq = "";
    open (my $fh, $genome_seq_file) or confess "error, cannot open file $genome_seq_file to read genome sequence";
    while (<$fh>) {
        if (/>/) { next; }
        s/\s+//g;
        $genome_seq .= uc $_;
    }
    close $fh;
    
    return ($genome_seq);
}

####
sub create_graph {
    my ($image, $graph_name, $canvas_rectangle) = @_;
    
    my $graph = new VectorGraph( { canvas_params => $canvas_rectangle,
                                                                      
                                 } );
    
    my $graph_input_files = $cfg->val("$graph_name", "GRAPH_INPUT_FILES") or confess "Error retrieving graph input files from cfg file";
    
    my $window_averaging_flag = $cfg->val("$graph_name", "GRAPH_WINDOW_AVERAGING");
    my $window_size = $cfg->val("$graph_name", "GRAPH_WINDOW_SIZE");

    my @inputs = split (/;/, $graph_input_files);
    foreach my $data_set (@inputs) {
        $data_set =~ s/\s+//g;
        my @params = split (/,/, $data_set);
        # expect the following keys: file, name, color
        my %data_params;
        foreach my $param (@params) {
            my ($token, $val) = split (/:/, $param);
            $data_params{$token} = $val;
        }
        my $data_file = $data_params{file} or confess "Error, data set $data_set lacks file param";
        my $data_color = $data_params{color} or confess "Error, no color specified for data set $data_set";
        my $data_name = $data_params{name} or confess "Error, no name specified for data set $data_set";

        my @points = &read_data_points($data_params{file}, $window_averaging_flag, $window_size);

        

        $graph->add_vector($data_name, $data_color, \@points);
    }
    
    ## get axis settings from config file
    my %axis_settings;
    

    my $y_axis_max = $cfg->val($graph_name, "Y_AXIS_MAX");
    if (defined($y_axis_max)) {
        $axis_settings{y_axis_max} = $y_axis_max;
    }
    my $y_axis_min = $cfg->val($graph_name, "Y_AXIS_MIN");
    if (defined($y_axis_min)) {
        $axis_settings{y_axis_min} = $y_axis_min;
    }
    
    
    $graph->draw_graph($image, \%colors, \%axis_settings);
    
    ## add graph title:
    my $graph_name_text = $graph_name;
    $graph_name_text =~ s/_/ /g;
    
    $image->string(gdSmallFont, $canvas_rectangle->{x1}, $canvas_rectangle->{y1} - 13, $graph_name_text, $colors{black});

    
    return;
    
}


####
sub read_data_points {
    my ($file, $window_averaging_flag, $window_size) = @_;
    
    open (my $fh, "$file") or die "Error, cannot read file [$file] containing data points.";
    my @points;
    while (<$fh>) {
        chomp;
        my ($x, $y) = split (/\s+/);
        
        ###############
        # todo: handle negative values
        #$y = abs($y);
        ###############
        
        unless ($x =~ /\d/ && $y =~ /\d/) {
            confess "Error, line $_ of file $file lacks an X Y value pair.";
        }
        if ($x >= $range_start && $x <= $range_end) {
            push (@points, [$x, $y]);
        }
    }
    close $file;
    
    @points = sort {$a->[0]<=>$b->[0]} @points; # make sure ordered by X-value

    ## average the data points:
    if ($window_averaging_flag =~ /true/i) {
        @points = &average_within_window(\@points, $window_size);
    }
    

    return (@points);
}

####
sub average_within_window {
    my ($points_aref, $window_size) = @_;
    
    print STDERR "Averaging points within window $window_size\n";

    ## Todo: don't trust points will be available for every base.  Check distance between data points.

    my @points;
    
    for (my $i = 0; $i <= $#$points_aref; $i+= $window_size) {
        my $num_pts = 0;
        my $sum_val = 0;
        for (my $j = $i; $j < $i + $window_size; $j++) {
            if ($j > $#$points_aref) { last; }
            $num_pts++;
            my $x_val = $points_aref->[$j]->[0]; # x-value
            my $y_val = $points_aref->[$j]->[1]; # y-value
            # print STDERR "incoming: ($x_val, $y_val)\n";
            $sum_val += $y_val;
        }
        my $x_val = $i + int($num_pts / 2);
        my $y_val = $sum_val / $num_pts;
        # print STDERR "\tavg -> ($x_val, $y_val)\n";
        
        push (@points, [$x_val, $y_val]);
    }

    return (@points);
}

####
sub draw_trellis {
    my ($image, $canvas_params) = @_;

    my $trellis_file = $conf{TRELLIS_FILE};
    my $trellis_final_path = $conf{TRELLIS_PATH};
    my $trellis_parser = new TrellisParser();
    my $trellis = $trellis_parser->parse_trellis_file($trellis_file, $trellis_final_path, $range_start, $range_end);
    
    my $trellis_illustrator = new TrellisIllustrator($canvas_params);
    $trellis_illustrator->draw_trellis($trellis, $image, \%colors, $range_start, $range_end, \$genome_sequence);

    ## add text header:
    $image->string(gdSmallFont, $canvas_params->{x1}, $canvas_params->{y1} - 13, "Highest Scoring Path Thru Candidate Exons", $colors{black});


    return;
    
}

####
sub get_RGB_from_type {
    my ($type) = @_;
    #pick random color
    my $sum_val = 0;
    my @chars = split (//, $type);
    foreach my $char (@chars) {
        $sum_val += ord($char);
    }
    $sum_val %= 156;
    #$sum_val += 100;
    
    my ($R, $G, $B) = ($sum_val, $sum_val, $sum_val);
    
    ## one of them should be maxxed for bright color
    my @colors = (\$R, \$G, \$B);
    my $rand_index = $sum_val % 3;
    ${$colors[$rand_index]} = 200;
    
    return ($R, $G, $B);
}

####
sub get_type_color {
    my ($image, $type) = @_;

    if (my $rgb_aref = $feature_type_to_RGB{$type}) {
        return (&get_color($image, @$rgb_aref));
    }
    else {
        my @rgb = &get_RGB_from_type($type);
        $feature_type_to_RGB{$type} = [@rgb];
        return (&get_color($image, @rgb));
    }
}

