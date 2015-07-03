#!/usr/bin/env perl

package main;
our $SEE;

package TkAnnotationViewer;
use strict;
use Tk;
use Clone_viewer;
use Feature_list;
use Gene_manager;
use Fasta_seqs;
use Gene_info;

####
sub new {
    my $package_name = shift;
    #instance members
    my $self = { mainWindow => 0,
		 notebook => 0,
		 clone_viewer => 0,
		 feature_list => 0,
		 gene_info => 0,
		 fasta_seqs => 0
		 };
    bless $self;
    return ($self);
}


####
sub set_up_gui {
    my $self = shift;
    my ($notebook); #gui components.
    my $mainWindow = new MainWindow();
    $self->{mainWindow} = $mainWindow;
    my $notebook = $mainWindow->NoteBook(-dynamicgeometry=>'false');
    $self->{notebook} = $notebook;
    $notebook->pack;

    ## set up clone_viewer
    my $viewer_tab = $notebook->add("viewer_tab", -label=>'Viewer');
    my $clone_viewer = new Clone_viewer($self, $viewer_tab);   
    $self->{clone_viewer} = $clone_viewer;

    ## set up feature_list
    my $feature_list_tab = $notebook->add("feature_list_tab", 
					  -label=>'Feature list',
					  -state=>'disabled');
    my $feature_list = new Feature_list($self,$feature_list_tab); 
    $self->{feature_list} = $feature_list;

    ## set up gene_info
    my $gene_info_tab = $notebook->add("gene_info_tab", 
				       -label=>'Gene info', 
				       -state => 'disabled', 
				       -raisecmd => sub {&update_gene_info($self)});
    my $gene_info = new Gene_info($self, $gene_info_tab);
    $self->{gene_info} = $gene_info;


    ## set up fasta_seqs tab
    my $fasta_seqs_tab = $notebook->add("fasta_seqs_tab", 
					-label=>'Fasta', 
					-state => 'disabled', 
					-raisecmd => sub {&update_fasta_seqs($self)});
    my $fasta_seqs = new Fasta_seqs ($self, $fasta_seqs_tab);
    $self->{fasta_seqs} = $fasta_seqs;

    ## set up live_connect
    # my $live_connect_tab = $notebook->add("live_connect_tab", -label=>'Live connect');
    
    
}

####
sub update_other_components {
    my $self = shift;
    my $feature_list = $self->{feature_list};
    $feature_list->populate_feature_list();
    my $notebook = $self->{notebook};
    $notebook->pageconfigure('feature_list_tab', -state => 'normal');
}

####
sub update_with_new_gene_select {
    my $self = shift;
    my $feature_id = shift;
    print "VIEWER: updating with new gene selected: $feature_id\n" if $SEE;
    $self->{clone_viewer}->update_with_new_gene_select($feature_id);
    $self->{feature_list}->update_with_new_gene_select($feature_id);
    #$self->{gene_info}->update_with_new_gene_select($feature_id);
    &Gene_manager::set_selected_gene($feature_id);
    ## enable gene-based pages
    my $notebook = $self->{notebook};
    $notebook->pageconfigure('gene_info_tab', -state => 'normal');
    $notebook->pageconfigure('fasta_seqs_tab', -state => 'normal');
}


####
sub update_with_gene_deselected {
    my $self = shift;
    $self->{clone_viewer}->update_with_gene_deselected();
    $self->{feature_list}->update_with_gene_deselected();
    #$self->{gene_info}->update_with_new_gene_deselected();
    &Gene_manager::set_selected_gene(0);
    ## disable gene-based pages.
    my $notebook = $self->{notebook};
    $notebook->pageconfigure('gene_info_tab', -state =>'disabled');
    $notebook->pageconfigure('fasta_seqs_tab', -state =>'disabled');
    
}

####
sub update_fasta_seqs {
    my $self = shift;
    my $fasta_seqs = $self->{fasta_seqs};
    print "updating Fasta Seqs $fasta_seqs\n" if $SEE;
    $fasta_seqs->update_panel_content();
}

####
sub update_gene_info {
    my $self = shift;
    my $gene_info = $self->{gene_info};
    $gene_info->update_gene_info();
}

1; #end of package.



