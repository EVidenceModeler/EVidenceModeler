#!/usr/bin/env perl

package main;
our $SEE;

package Clone_viewer;

use Gene_obj;
use Tk;
use Tk::NoteBook;
use Tk::FileSelect;
use strict;
use TIGR_XML_parser;
use Gene_manager;
use Data::Dumper;
use GFF3_parser_simple;

## STATIC PACKAGE variables.
# @gene_container : static list of loaded gene objects.
# $currently_selected_gene : current gene selected (clone_viewer, feature_list)
use vars qw (@genes);


my $GENE_OUTLINE_COLOR = "purple";
my $SPLICE_BOUND_COLOR = "white";

my %COLOR = (annotated_genes=>'green');

my ($mainWindow, $canvas, $real_canvas, $viewer_tab); ##Main GUI components. 
my $scale_factor = 1.5; #used as the main scale factor.  scale_factor_X,Y simply track the scaling thru time.
my $scale_factor_X = 1;
my $scale_factor_Y = 1;
my $canvas_height = 700;
my $canvas_width = 1000;
my $initial_canvas_height = $canvas_height;
my $initial_canvas_width = $canvas_width;
my $feature_num = -1;
my $status_msg = "Status Bar";
my $status_msg_ref = \$status_msg;
my $fileselect_window = 0; # file selection dialogue
my $initial_zoom_factor = 5; ## needed to increase resolution
#$canvas_width *= $initial_zoom_factor;
#$canvas_height *= $initial_zoom_factor;
my $PIXELS_per_NT = 1;
my $controller; #object reference to the highest level object.
my %ANNOT_BOUNDARIES; #mark boundaries in evidence consistent with the curated models.
my $pixel_width_boundaries = 2;

sub new {
    my $package_name = shift (@_);
    
    $controller= shift (@_) or die "I need a mainWindow reference!";
    $mainWindow = $controller->{mainWindow};
    $viewer_tab = shift (@_) or die "I need a viewer widget reference!";
    
    ## Set up Menu Options
    my $top_menu_frame = $viewer_tab->Frame(-relief=>'groove')->pack(-side=>'top');
    my $status_frame = $viewer_tab->Frame(-relief=>"groove")->pack(-side=>'top', -fill=>'both'); 
   
    my $menu_button_FILE = $top_menu_frame->Menubutton(-relief => 'raised',-text => 'File');
    $menu_button_FILE->command(-label=>'Open',
                               -command=>\&open_file);
    $menu_button_FILE->command(-label=>'Close', -command=>\&close_file);
    $menu_button_FILE->command(-label=>'Save', -command=>\&save_file);
    $menu_button_FILE->command(-label=>'Exit', -command=>sub{$viewer_tab->destroy;});
    
    $menu_button_FILE->pack(-side => "left", -expand=>1, -fill=>'both');
    
    
    
    my $menu_button_VIEW = $top_menu_frame->Menubutton(-relief => 'raised',-text => 'View');
    $menu_button_VIEW->command(-label=>'get cDNA', -command=>\&get_cDNA_seq);
    $menu_button_VIEW->command(-label=>'get protein', -command=>\&get_protein);
    $menu_button_VIEW->pack(-side=>'left', -expand=>1, -fill=>'both');
    
    ## Set up ZOOM Buttons
    
    my $zoomX_btn = $top_menu_frame->Button(-text => "ZoomIN(X)", -command => \&zoomX_funct);
    $zoomX_btn->pack(-side=>'left', -expand=>1, -fill=>'both'); 
    my $unzoomX_btn = $top_menu_frame->Button(-text => "ZoomOUT(X)", -command => \&unzoomX_funct);
    $unzoomX_btn->pack(-side=>'left', -expand=>1, -fill=>'both'); 
    
    my $zoomY_btn = $top_menu_frame->Button(-text => "ZoomIn(Y)", -command => \&zoomY_funct);
    $zoomY_btn->pack(-side=>'left', -expand=>1, -fill=>'both');
    my $unzoomY_btn = $top_menu_frame->Button(-text => "ZoomOUT(Y)", -command => \&unzoomY_funct);
    $unzoomY_btn->pack(-side=>'left', -expand=>1, -fill=>'both');
    
    my $reset_btn = $top_menu_frame->Button(-text=>"reset", -command=>\&reset_display);
    $reset_btn->pack(-side=>'left', -expand=>1, -fill=>'both');
    $top_menu_frame->pack(-side=>'top', -expand=>1, -fill=>'both');
    
    
    ####
    ## Set up CANVAS
    ####
    
    $canvas = $viewer_tab->Scrolled("Canvas", -cursor=>"crosshair", -width=>$initial_canvas_width,-height=>$initial_canvas_height, -confine=>'true');
    $real_canvas = $canvas->Subwidget("canvas");
        
    #$real_canvas->Tk::bind('perimeter', '<Enter>', \&bind_annot_to_statusbar);
    $canvas->pack(-side=>'top', -fill=>'x');
    #$viewer_tab-
   
    $status_frame->Label(-textvariable=>$status_msg_ref, 
                       -borderwidth=>2, 
                       -relief=>'groove'
                       )->pack(
                               -fill=>'x', 
                               -side=>'bottom');
    
    ## return an instantiated object
    my $self = {};
    bless $self;
    return ($self);
}

####
sub display_annotation {
    my ($tigr_xml) = @_;
    ####
    # Get and display annotation
    ####
    &initialize_canvas();
    my $wait_message = &create_info_message("Loading $tigr_xml.  Only a moment, please.");
    my ($assembly_seq, $assembly_seq_length, @gene_container);
    my $TIGR_xml_parser = new TIGR_XML_parser();
    $TIGR_xml_parser->capture_genes_from_assembly_xml($tigr_xml);
    @gene_container = $TIGR_xml_parser->get_genes();
    $assembly_seq = $TIGR_xml_parser->get_assembly_sequence();
    &Gene_manager::set_assembly_seq($assembly_seq);
    $assembly_seq_length = length($assembly_seq);
    
    my (@top_section, @bottom_section);;
    $#top_section = 0;
    $#bottom_section = 0; #temp hack to keep central line from colliding with gene models.
    
    # position all elements within their respective rows.  Don't allow overlaps.
    my $MAX_X = &array_elements(\@gene_container, \@top_section, \@bottom_section);
    if ($assembly_seq_length > $MAX_X) {
        $MAX_X = $assembly_seq_length;
    }
   
    $PIXELS_per_NT = $canvas_width/$MAX_X;
    
    ## add evidence data
    my $assembly_evidence = $TIGR_xml_parser->{assembly_evidence};
    foreach my $evidence_type (keys %$assembly_evidence) {
        my $evidence_obj = $assembly_evidence->{$evidence_type};
        $evidence_obj->collapse_elements(2500);
        my @container = $evidence_obj->get_all_evidence();
       
        &array_elements(\@container, \@top_section, \@bottom_section);
       
    }
    
    
    ## draw elements on the canvas:
    
    &draw_elements(\@top_section, \@bottom_section);
    $wait_message->destroy();
    
    ## tell controller that annotation has finished loading
    $controller->update_other_components();
}



sub display_gff3_annots {
    my $self = shift;
    my ($files_listing_href) = @_;
    
    my $genome_seq_file = $files_listing_href->{genome_seq_file};
    my @gene_predictions_files = @{$files_listing_href->{gene_predictions_files}};
    my @search_evidence_files = @{$files_listing_href->{search_evidence_files}};

    
    ####
    # Get and display annotation
    ####
    &initialize_canvas();
    my $wait_message = &create_info_message("Loading Inputs.  Only a moment, please.");
    
    my @max_coords;
    
    ## set the genome sequence
    my $assembly_seq = &_get_seq_from_fasta_file($genome_seq_file);
    my $assembly_seq_length = length($assembly_seq);
    &Gene_manager::set_assembly_seq($assembly_seq);

    
    ## add the gene predictions:
    my (@top_section, @bottom_section);;
    $#top_section = 0;
    $#bottom_section = 0; #temp hack to keep central line from colliding with gene models.
    
    
    ## Add the gene prediction tiers:
    
    foreach my $gene_prediction_file (@gene_predictions_files) {
        
        if (my @predictions = &_get_genes_from_gff3($gene_prediction_file)) {
            foreach my $prediction (@predictions) {
                $prediction->{classification} = "genepred:$prediction->{source}";
            }
            my $max_coord = &array_elements(\@predictions, \@top_section, \@bottom_section);
            push (@max_coords, $max_coord);
        }
    }

        
    ## add search evidence data
    
    if (@search_evidence_files) {
        my $assembly_evidence_href = &_parse_search_evidence_tiers(@search_evidence_files);
        foreach my $evidence_type (keys %$assembly_evidence_href) {
            my $evidence_obj = $assembly_evidence_href->{$evidence_type};
            my @container = $evidence_obj->get_all_evidence();
            my $max_coord = &array_elements(\@container, \@top_section, \@bottom_section);
            push (@max_coords, $max_coord);
        }
    }
    
    @max_coords = reverse sort {$a<=>$b} @max_coords;
    my $MAX_X = shift @max_coords;
    

    if ($assembly_seq_length > $MAX_X) {
        $MAX_X = $assembly_seq_length;
    }
   
    $PIXELS_per_NT = $canvas_width/$MAX_X;

    ## draw elements on the canvas:
    
    &draw_elements(\@top_section, \@bottom_section);
    $wait_message->destroy();
    
    ## tell controller that annotation has finished loading
    $controller->update_other_components();
}



sub create_info_message {
    my ($text) = @_;
    my $message_frame = $mainWindow->Toplevel();
    my $message = $message_frame->Message(-text=>$text);
    $message->pack;
    $message->waitVisibility();
    return ($message_frame);
}



sub initialize_canvas {
    &Gene_manager::clear_contents();
    $real_canvas->delete("all");
    $canvas->createRectangle(0,0,$canvas_width, $canvas_height, -tags=>"perimeter", -fill=>'black');
    $real_canvas->bind("perimeter", "<Button-1>", &gene_deselected());
}

####
sub zoomX_funct {
    
    $canvas->scale("all", 0,0,$scale_factor,1);
    &reset_scroll_region();
    $scale_factor_X *= $scale_factor;
}


sub unzoomX_funct {
    if ($scale_factor_X >= $scale_factor) {
        $canvas->scale("all", 0,0,(1/$scale_factor),1);
        &reset_scroll_region();
        $scale_factor_X /= $scale_factor;
    }
}

sub zoomY_funct {
    $canvas->scale("all", 0,0,1,$scale_factor);
    &reset_scroll_region();
    $scale_factor_Y *= $scale_factor;
}

sub unzoomY_funct {
    if ($scale_factor_Y >= $scale_factor) {
        $canvas->scale("all", 0,0,1,(1/$scale_factor));
        &reset_scroll_region();
        $scale_factor_Y /= $scale_factor;
    }
}



####
sub add_rect_funct {
    $canvas->createRectangle(100,100,150,150, -fill=>'blue');
}

####
sub reset_display {
    my $scaleX_reduction = 1/$scale_factor_X;
    my $scaleY_reduction = 1/$scale_factor_Y;
    $canvas->scale("all", 0,0, $scaleX_reduction, $scaleY_reduction);
    $scale_factor_X = 1;
    $scale_factor_Y = 1;
    &reset_scroll_region();
}


####
sub reset_scroll_region {
    my @coords = $canvas->bbox("perimeter");
   
    $canvas->configure(-scrollregion=>[@coords]);
}


####
sub array_elements {
    my ($container_ref, $top_section_ref, $bottom_section_ref) = @_;
    my $max_x_coord = 0;
    my $start_top_row = $#{$top_section_ref} + 1; #> 0) ? $#{$top_section_ref} : 0;
    my $start_bottom_row = $#{$bottom_section_ref} + 1; #> 0) ? $#{$bottom_section_ref} : 0;
    
    
    foreach my $seq_element (@$container_ref) {
        my $type = $seq_element->{classification};
        my $strand = $seq_element->{strand};
        my @coords = sort {$a<=>$b} $seq_element->get_seq_span();
        my $local_max_x = $coords[1];
       
        if ($local_max_x > $max_x_coord) { $max_x_coord = $local_max_x;}
        if ($strand eq '+') {
            &position_element($seq_element, $top_section_ref, $start_top_row, @coords);
        } else {
            &position_element($seq_element, $bottom_section_ref, $start_bottom_row, @coords);
        }
    }
    return ($max_x_coord);
}


####
sub position_element {
    my ($seq_element, $section_ref, $start_row, $lend, $rend) = @_;
    my $found_position = 0;
    my $type = $seq_element->{classification};
    ## add to lookup table
    $feature_num++;
    my $feature_id = $type . '_' . $feature_num;
    
    
    
    ## store seq_element in a wrapper struct that includes type and lend, rend info.
    my $mid_pt =  (($lend + $rend)/2);
    my $current_struct_ref = {seq_element=>$seq_element,
                              type=>$type,
                              lend=>$lend,
                              rend=>$rend,
                              mip_pt => $mid_pt,
                              feature_id => $feature_id};
    
    &Gene_manager::add_gene_feature ($feature_id, $seq_element, $mid_pt);
    ## start at last type-occupied row. 
    my $row = $start_row;
 
    while (!$found_position) {
       
        my $row_set_ref = $section_ref->[$row];
       
        if (ref $row_set_ref) { #must be some elements here.
          
            my @existing_elements = 
                sort {$a->{lend} <=> $b->{lend}
                      ||
                          $b->{rend} <=> $a->{rend}
                  } @$row_set_ref;
            my $overlap = 0;
            foreach my $element (@existing_elements) {

                my $element_lend = $element->{lend};
                my $element_rend = $element->{rend};
                if ($element_lend > $rend) {
                    ## room in this row for current element
                    last;
                }
                ## check for overlap with element
                ## if overlap, must go to next row and start over.
                if ($lend <= $element_rend && $rend >= $element_lend) {
                    $row++;
                    $overlap = 1;
                    last;
                }
            }
            if (!$overlap && !$found_position) {
                ## room in this row for current element
                push (@$row_set_ref, $current_struct_ref);
                $found_position = 1;
            }
        } else {
            #create row, set to contain one element
            $section_ref->[$row] = [$current_struct_ref];
            $found_position = 1;
        }
    }
    
}



####
sub draw_elements {
    my ($top_sect_ref, $bottom_sect_ref) = @_;
    
    ## do some calculations for initial space allocations:
    my $num_top_rows = $#{$top_sect_ref} + 1;
    my $num_bottom_rows = $#{$bottom_sect_ref} + 1;
    my $max_row_sect = ($num_top_rows > $num_bottom_rows) ? $num_top_rows : $num_bottom_rows;
    my $total_num_rows = 2 * $max_row_sect + 1;
    my $max_element_height = $canvas_height/$total_num_rows;
    if ($max_element_height > 200) {
        $max_element_height = 100;
    }
    my $used_element_height = int(0.8 * $max_element_height);
    my $spacer_height = int (($max_element_height - $used_element_height)/2);
        
    ## draw strand separator
    my $central_pos_Y = int ($canvas_height/2);
    my $strand_sep_width = 10;
    #$canvas->createRectangle(0,$central_pos_Y+$strand_sep_width/2, $canvas_width, 
	#		     $central_pos_Y+$strand_sep_width/2, -fill=>'black', -tags=>'DNAstrand');
    
    &create_ticker (0, $central_pos_Y, $canvas_width);
    
    ## Draw each element
    ## Do top section first.
    for (my $i=0; $i < $num_top_rows; $i++) {
        my $row_num = $i + 1;
        my $top_Y = $central_pos_Y - ($row_num * $max_element_height);
        my $element_top_Y = $top_Y + $spacer_height;
        my $element_bottom_Y = $element_top_Y + $used_element_height;
        foreach my $element (@{$top_sect_ref->[$i]}) {
            &draw_single_seq_element($element, $element_top_Y, $element_bottom_Y);
        }
    }
    ## Now, do bottom section
    for (my $i=0; $i < $num_bottom_rows; $i++) {
        my $row_num = $i + 1;
        my $top_Y = $central_pos_Y + ($row_num * $max_element_height);
        my $element_bottom_Y = $top_Y - $spacer_height;
        my $element_top_Y =  $element_bottom_Y - $used_element_height;
        foreach my $element (@{$bottom_sect_ref->[$i]}) {
            &draw_single_seq_element($element, $element_top_Y, $element_bottom_Y);
        }
    }
}


####
sub draw_single_seq_element {
    my ($element, $element_top_Y, $element_bottom_Y) = @_;
    #unwrap struct
    my $seq_element = $element->{seq_element};
    my $type = $element->{type};
    my $feature_id = $element->{feature_id};
    print "Feature: $feature_id\n" if $SEE;
    ## draw line thru center of gene, connecting exons:
    my ($gene_lend, $gene_rend) = &convert_coords($seq_element->get_seq_span());
    my $gene_central_Y = (($element_top_Y +  $element_bottom_Y)/2);
    $canvas->createLine($gene_lend, $gene_central_Y, $gene_rend, $gene_central_Y, -fill=>'white');   
    ## Draw the exons
    my $element_color = &get_type_color($type);
    print "color $element_color for type $type\n" if $SEE;
    my @segments = $seq_element->get_segments();
    foreach my $segment (@segments) {
        my ($asmbl_end5, $asmbl_end3) = $segment->get_coords();
        
        my $seg_feat_id = "$feature_id$;$asmbl_end5-$asmbl_end3";
        
        my ($end5, $end3) = &convert_coords($asmbl_end5, $asmbl_end3);
        $canvas->createRectangle($end5, $element_top_Y, $end3, $element_bottom_Y, 
                                 -fill=>$element_color, 
                                 #-outline => 'black', 
                                 -tags=>[$feature_id, $seg_feat_id, "segment"]);
               
        #track asmbl coords for annotated genes.  Mark consistent evidence.
        if ($type =~ /genepred:/) {
            
            $ANNOT_BOUNDARIES{$asmbl_end5} .=  "$feature_id-splice ";
            $ANNOT_BOUNDARIES{$asmbl_end3} .= "$feature_id-splice ";
        } 
        
        if (my $pred_ids = $ANNOT_BOUNDARIES{$asmbl_end5}) {
            chomp $pred_ids;
            my @tags = split (/ /, $pred_ids);
            
            my ($end5_adj) = &convert_coords($asmbl_end5);
            $canvas->createRectangle($end5, $element_top_Y, $end5_adj, $element_bottom_Y, 
                                     -fill=>"black", 
                                     -outline => 'black', 
                                     -tags=>[@tags, "splicebounds"]
                                     );
            print "Adding tags: @tags to $asmbl_end5\n" if $SEE;
            
        }
        if (my $pred_ids = $ANNOT_BOUNDARIES{$asmbl_end3}) {
            chomp $pred_ids;
            my @tags = split (/ /, $pred_ids);
            print "TAGS: [" . join (";", @tags) . "\n" if $SEE;
            my ($end3_adj) = &convert_coords($asmbl_end3 - $pixel_width_boundaries);
            $canvas->createRectangle( $end3_adj, $element_top_Y, $end3, $element_bottom_Y, 
                                      -fill=>"black", 
                                      -outline => 'black', 
                                      -tags=>[@tags, "splicebounds"]
                                      );
        }
        
        
        
        $real_canvas->bind($seg_feat_id, '<Enter>', &enter_element_callback($seg_feat_id));
        $real_canvas->bind($seg_feat_id, '<Leave>', &leave_element_callback($seg_feat_id));
        
    }
    if ($element->{type} =~ /genepred:/) {
        $real_canvas->bind($feature_id, '<Button-1>', &gene_select_callback($feature_id));
    }

    
}

####
sub convert_coords {
    my @coords = @_;
    my @converted;
    foreach my $coord (@coords) {
        push (@converted, ($coord * $PIXELS_per_NT));
    }
    return (@converted);
}


####
sub get_type_color {
    my ($type) = @_;
    my $ret_color;
    if ($ret_color = $COLOR{$type}) {
        return ($ret_color);
    } else {
        #pick random color
        my $sum_val = 0;
        my @chars = split (//, $type);
        foreach my $char (@chars) {
            $sum_val += ord($char);
        }
        print "sum_val: $sum_val\n" if $SEE;
        $sum_val %= 156;
        $sum_val += 100;
        print "sum_val aft: $sum_val\n" if $SEE;

        my $R = sprintf("%02x", $sum_val);
        my $G = sprintf("%02x", $sum_val);
        my $B = sprintf("%02x", $sum_val);

        ## one of them should be maxxed for bright color
        my @colors = (\$R, \$G, \$B);
        my $rand_index = $sum_val % 3;
        ${$colors[$rand_index]} = sprintf ("%x", 255); # - ${$colors[$rand_index]});
        my $color = ("#" . $R . $G . $B); 
        print "color: $color\n" if $SEE;
        $COLOR{$type} = $color;
       
        return ($color);
    }
}

####
sub get_annot_text {
    my ($feature_id) = @_;
    my $seq_element = &Gene_manager::get_gene_via_feature_id($feature_id);
    my $text = $seq_element->get_annot_text();
    if (length($text) > 100) {
        $text = substr($text, 0, 100);
    }
    return ($text);
}

####
sub enter_element_callback {
    my ($feature_id) = @_;
    my ($feature_id_adj, $added_text) = split (/$;/, $feature_id);
    $feature_id = $feature_id_adj;
    my $annot_text = &get_annot_text($feature_id);
    if ($added_text) {
        $annot_text .= " $added_text";
    }
    my $subref = sub {
        #$real_canvas->itemconfigure($feature_id, -fill=>'blue');
        $$status_msg_ref = $annot_text;
        print "Entering $feature_id\n" if $SEE;
    };
    return ($subref);
}


####
sub leave_element_callback {
    my ($feature_id) = @_;
    my $subref = sub {
        $$status_msg_ref = "";
        print "Leaving $feature_id\n" if $SEE;
    };
    return ($subref);
}

####
sub gene_select_callback {
    my ($feature_id) = @_;
    my $subref = sub {
        $controller->update_with_new_gene_select($feature_id);
    };
    return ($subref);
}

####
sub gene_deselected {
    my $subref = sub {
        $controller->update_with_gene_deselected();
    };
    return ($subref);
}


####
sub update_with_new_gene_select {
    my $self = shift;
    my $feature_id = shift;
    &unhighlight_genes();
    &highlight_gene($feature_id);
}


####
sub update_with_gene_deselected {
    &unhighlight_genes();
}



####
sub highlight_gene {
    my ($feature_id) = @_;
    $real_canvas->itemconfigure($feature_id, -outline => $GENE_OUTLINE_COLOR);

    $real_canvas->itemconfigure("$feature_id-splice", -fill => $SPLICE_BOUND_COLOR);
    $real_canvas->itemconfigure("$feature_id-splice", -outline => $SPLICE_BOUND_COLOR);
    print "$feature_id-splice highlighted.\n" if $SEE;
    return;
}

####
sub unhighlight_genes {
    #$real_canvas->itemconfigure('annotated_genes', -fill => $COLOR{'annotated_genes'});
    $real_canvas->itemconfigure('segment', -outline => 'black');
    #$real_canvas->itemconfigure("splice_bounds", -fill => 'black');
   # $real_canvas->itemconfigure("splice_bounds", -outline => 'black');
    $real_canvas->itemconfigure("splicebounds", -fill => 'black');
    $real_canvas->itemconfigure("splicebounds", -outline => 'black');
    
    return;
}



##############################################################
##### Handlers for Menu Options ##############################
##############################################################

sub open_file{
    if (!Exists($fileselect_window)) {
        $fileselect_window = $viewer_tab->FileSelect();
        my $filename = $fileselect_window->Show();
        if ($filename) {
            &display_annotation($filename);
        }
        $fileselect_window = 0;
    }
}

sub close_file {
}

sub save_file {
}

sub get_cDNA_seq {
}

sub get_protein {
}

sub pick_random_hex {
    
    my $rand = rand(1);
    my $num = int ($rand * 256);
    #if ($num < 100) {
	#$num += 100;
    #}
    
    my $hex = sprintf ("%x", $num);
    print "color $num $rand $hex\t" if $SEE;
    while (length $hex > 2) {
        chop $hex;
    }
    while (length $hex < 2) {
        $hex = '0' . $hex;
    }
    print "$hex\n" if $SEE;
    return ($hex);
}


####
sub create_ticker {
    my ($xstart, $centraly, $xstop) = @_;
    #draw central line
    $canvas->createLine($xstart, $centraly, $xstop, $centraly, -fill=>'white');
    
    #draw ticks.
    my $short_ticker_height = 3;
    my $large_ticker_height = 10;
    my $xpos = 0;
    my $big_tick_at_nuc = 10000;
    my $short_tick_at_nuc = 1000; 
    for (my $x = 0; $xpos < $xstop; $x += $short_tick_at_nuc) {
        $xpos = $x * $PIXELS_per_NT;
        my $tick_height;
        if ($x % $big_tick_at_nuc == 0) {
            $tick_height = $large_ticker_height;
        } elsif ($x % $short_tick_at_nuc == 0) {
            $tick_height = $short_ticker_height;
        }
        if ($tick_height) {
            
            $canvas->createLine ($xpos, $centraly, $xpos, $centraly - $tick_height, -fill => 'white');
        }
    }
    
}


################################
# helper functions, private

####
sub _get_seq_from_fasta_file {
    my ($genome_seq_file) = @_;
    
    my $genome_seq = "";
    open (my $fh, $genome_seq_file) or die $!;
    while (<$fh>) {
        if (/>/) { next; }
        $genome_seq .= $_;
    }
    close $fh;

    $genome_seq =~ s/\s//g;

    return ($genome_seq);

}


####
sub _get_genes_from_gff3 {
    my ($filename) = @_;

    my $gff3_parser = new GFF3_parser_simple();
    
    $gff3_parser->parse_GFF3_file($filename);

    my @gene_objs = values %{$gff3_parser->{gene_objs}};

    foreach my $gene_obj (@gene_objs) {
        $gene_obj->trim_UTRs();
        $gene_obj->{com_name} .= " " . $gene_obj->{source} . " " . $gene_obj->{TU_feat_name} . " " . $gene_obj->{Model_feat_name};
        
    }

    

    return (@gene_objs);
}


####
sub _parse_search_evidence_tiers {
    my @search_evidence_files = @_;

    my %ev_type_to_ev_obj; ## maintain different ev obj for each ev type:

    foreach my $search_evidence_file (@search_evidence_files) {
        
        open (my $fh, $search_evidence_file) or die "Error, cannot open file $search_evidence_file";
        while (<$fh>) {
            chomp;
            unless (/\w/) { next; }
            if (/^\#/) {
                next; # comment line
            }

            my @x = split (/\t/);
            my ($asmbl_id, $ev_type, $feat_type, $lend, $rend, $percent_identity, $orient, $phase, $match_line) = @x;

            my ($end5, $end3) = ($orient eq '+') ? ($lend, $rend) : ($rend, $lend);

            ## parse the match line:
            $match_line =~ /ID=([^;]+)/;
            my $match_identifier = $1 or die "Error, no match_id for line $match_line";
            $match_line =~ /Target=([^\s;]+)/;
            my $acc = $1 or die "Error, on Target accession for line $match_line";
            
            
            my $match_length = $rend - $lend + 1;
            my ($match_lend, $match_rend) = (1, $match_length); # default
            if ($match_line =~ /\s+(\d+)\s+(\d+)\s*$/) {
                ($match_lend, $match_rend) = ($1, $2);
            }
            
            my $evidence_obj = $ev_type_to_ev_obj{$ev_type};
            unless (ref $evidence_obj) {
                $evidence_obj = $ev_type_to_ev_obj{$ev_type} = new Assembly::Evidence();
            }



            $evidence_obj->add_evidence($search_evidence_file,
                                        $ev_type,
                                        $acc,
                                        $match_identifier,
                                        $match_lend, $match_rend,
                                        $end5, $end3,
                                        "$ev_type, $acc");
            
        
        }
    }

    foreach my $ev_obj (values %ev_type_to_ev_obj) {
        
        $ev_obj->refine_evidence_coords();

        
    }
    
    return (\%ev_type_to_ev_obj);
    
}




1; #end of package.
