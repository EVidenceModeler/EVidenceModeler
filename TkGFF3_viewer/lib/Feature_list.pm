#!/usr/bin/env perl

package main;

our $SEE;


package Feature_list;

use Tk;
use strict;
use Gene_manager;

my $controller; ##caller
my $parent_window;
my $list_box;
my @feature_lookup;
my %index_lookup;

sub new {
    shift;
    $controller = shift;
    $parent_window = shift;
     #$list_box = $parent_window->Listbox(-width => 0, -height => 0);    
    $list_box = $parent_window->Scrolled("Listbox", -width =>100, -height => 30);
    $list_box->bind("<Button-1>", \&list_box_callback);
    $list_box->pack(-expand => 1, -fill => 'both');
    my $self ={};
    bless $self;
    return ($self);
}

sub populate_feature_list {
    &delete_list();
    my @feature_ids = &Gene_manager::get_ordered_feature_ids();
    my $feat_num = 0;
    foreach my $feature_id (@feature_ids) {
	my $gene_obj = &Gene_manager::get_gene_via_feature_id($feature_id);
	unless ($gene_obj->{classification} =~ /genepred:/) { next;}
	my $text = ($feat_num + 1) . ". ";
	my $pub_locus = $gene_obj->{pub_locus};
	my $locus = $gene_obj->{locus};
	my $com_name = $gene_obj->{com_name};
	my $pub_comment = $gene_obj->{pub_comment};
	if ($pub_locus) {
	    $text .= "$pub_locus ";
	}
	if ($locus) {
	    $text .= "$locus ";
	} 
	if ($com_name) {
	    $text .= ": $com_name";
	}
	if ($pub_comment) {
	    $text .= " : $pub_comment";
	}
	$list_box->insert("end", $text);
	$feature_lookup[$feat_num] = $feature_id;
	$index_lookup{$feature_id} = $feat_num;
	$feat_num++;
    }

  
}


####
sub list_box_callback {
    my ($call_ref) = @_;
    my  @selections = $call_ref->curselection();
    print "selections: @selections\n" if $SEE;
    my $selection = shift @selections;
    my $feature_id = $feature_lookup[$selection];
    $controller->update_with_new_gene_select($feature_id);
}

####
sub update_with_new_gene_select {
    my $self = shift;
    my $feature_id = shift;
    my $index = $index_lookup{$feature_id};
    &clear_list_box_selections();
    print "Activating index: $index\n" if $SEE;
    $list_box->selectionSet($index);
} 

####
sub update_with_gene_deselected {
    &clear_list_box_selections();
}


####
sub clear_list_box_selections {
    my $size = $list_box->size();
    if ($size) {
	$list_box->selectionClear(0, $size);
    }
}

####
sub delete_list {
    my $size = $list_box->size();
    if ($size) {
	$list_box->delete(0, $size);
    }
    @feature_lookup = ();
    %index_lookup = ();

}

1; #end of package.

