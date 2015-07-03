#!/usr/bin/env perl

package Gene_info;

use Tk;
use strict;
use Gene_manager;

my $controller; ##caller
my $parent_window;
my $list_box;

sub new {
    shift;
    $controller = shift;
    $parent_window = shift;
    $list_box = $parent_window->Scrolled("Listbox", -width =>100, -height => 30);
    #$list_box->bind("<Button-1>", \&list_box_callback);
    $list_box->pack(-expand => 1, -fill => 'both');
    my $self ={};
    bless $self;
    return ($self);
}

####
sub update_gene_info {
    &delete_list();
    my ($feature_id) = &Gene_manager::get_selected_gene();
    if ($feature_id) {
	my $gene_obj = &Gene_manager::get_gene_via_feature_id($feature_id);
	my @gene_infos = &get_gene_info($gene_obj);
	foreach my $gene_info (@gene_infos) {
	    $list_box->insert("end", $gene_info);
	}
    }
}


####
sub delete_list {
    my $size = $list_box->size();
    if ($size) {
	$list_box->delete(0, $size);
    }
}

#### 
sub get_gene_info {
    my ($gene_obj) = @_;
    my @gene_infos;
    ## get locus info
    my $locus_text;
    my $locus = $gene_obj->{locus};
    if ($locus) {
	$locus_text .= "[$locus] ";
    } 
    my $pub_locus = $gene_obj->{pub_locus};
    if ($pub_locus) {
	$locus_text .= "[$pub_locus] ";
    }
    if ($locus_text) {
	push (@gene_infos, "Locus: $locus_text");
    }
    ## get gene name.
    my $com_name = "Gene name: " . $gene_obj->{com_name};
    push (@gene_infos, $com_name);
    ## get pub comment
    my $pub_comment = $gene_obj->{pub_comment};
    if ($pub_comment) {
	push (@gene_infos, "Comment: $pub_comment");
    }
    ## get coordinate info
    my $orientation = $gene_obj->{strand};
    my @gene_span = $gene_obj->get_gene_span();
    push (@gene_infos, "Gene coordinates: $gene_span[0]-$gene_span[1]");
    my @exons = $gene_obj->get_exons();
    @exons = sort {$a->{end5} <=> $b->{end5}} @exons;
    if ($orientation eq '-') {
	@exons = reverse (@exons);
    }
    my $mRNA_exon_text = "mRNA exons: ";
    my $CDS_exon_text = "CDS exons: ";
    foreach my $exon (@exons) {
	my @mRNA_coords = $exon->get_mRNA_exon_end5_end3();
	$mRNA_exon_text .= "$mRNA_coords[0]-$mRNA_coords[1] ";
	my @CDS_coords = $exon->get_CDS_end5_end3();
	if (@CDS_coords) {
	    $CDS_exon_text .= "$CDS_coords[0]-$CDS_coords[1] ";
	}
    }
    push (@gene_infos, $mRNA_exon_text, $CDS_exon_text);
    return (@gene_infos);
}


1; #end of package.

