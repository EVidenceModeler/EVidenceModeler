#!/usr/bin/env perl

package main;
our $SEE;


package Fasta_seqs;
use strict;
use Gene_manager;
use Tk;

my $controller;
my $parent_window;
my $protein_seq;
my $cds_seq;
my $genomic_transcript;

## gui components.
my $protein_page;
my $cds_page;
my $genomic_trans_page;
my $cds_text;
my $protein_text;
my $genomic_trans_text;


sub new {
    shift;
    $controller = shift;
    $parent_window = shift;
    my $notebook = $parent_window->NoteBook(-dynamicgeometry=>'false');
    $protein_page = $notebook->add('protein_tab', -label => 'Protein');
    $cds_page = $notebook->add('CDS_page', -label => 'CDS');
    $genomic_trans_page = $notebook->add('genomic_transcript', -label =>'Genomic transcript');
    $notebook->pack(-expand=>1, -fill=>'both');
    &update_panel_content();
    
    ## return Fasta_seqs object.
    my $self = {};
    bless $self;
    return ($self);
}


####
sub update_panel_content {
    my $self = shift;
    &destroy_old_text();
    my $selected_gene_feature_id = &Gene_manager::get_selected_gene();
    print "Fasta_seqs: selected_gene_feature_id = $selected_gene_feature_id\n" if $SEE;
    my ($gene_seq, $cds_seq, $prot_seq);
    if ($selected_gene_feature_id) {
        my $selected_gene = &Gene_manager::get_gene_via_feature_id($selected_gene_feature_id);
        my $assembly_seq_ref = &Gene_manager::get_assembly_seq_ref();
       	$gene_seq = $selected_gene->create_gene_sequence($assembly_seq_ref, 1);
        $cds_seq = $selected_gene->create_CDS_sequence($assembly_seq_ref);
        $prot_seq = $selected_gene->get_protein_sequence();
    } else {
        $gene_seq = "No gene model is currently selected";
        $cds_seq = "No gene model is currently selected";
        $prot_seq = "No gene model is currently selected";
    }
    
    $protein_text = $protein_page->Text(-height => 10, -cursor => 'top_left_arrow')->pack(-expand=>1, -fill => 'both');
    $protein_text->insert('end', $prot_seq);
    $cds_text = $cds_page->Text(-height => 10, -cursor => 'top_left_arrow')->pack(-expand =>1, -fill=>'both');
    $cds_text->insert('end', $cds_seq);
    $genomic_trans_text = $genomic_trans_page->Text(-height => 10, -cursor => 'top_left_arrow')->pack(-expand=>1, -fill=>'both');
    $genomic_trans_text->insert('end', $gene_seq);
}

####
sub destroy_old_text {
    if ($protein_text) {
        $protein_text->destroy();
        $protein_text = 0;
    }
    if ($cds_text) {
        $cds_text->destroy();
        $cds_text = 0;
    }
    if ($genomic_trans_text) {
        $genomic_trans_text->destroy();
        $genomic_trans_text = 0;
    }
}


1; #end of module.




