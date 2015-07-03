#!/usr/bin/env perl

package Gene_manager;
use strict;

## Do not instantiate this.  Use only as a namespace.

## private package variables:
my $selected_gene = 0; #identifies currently selected gene.
my %id_to_gene_obj = (); #keys gene objects to feature ids.
my $assembly_seq_ref = 0;

#### 
sub add_gene_feature {
    my ($feature_id, $gene_obj, $mid_pt) = @_;
    my $gene_struct = {feature_id => $feature_id,
		       gene_obj => $gene_obj,
		       mid_pt => $mid_pt};
    $id_to_gene_obj{$feature_id} = $gene_struct;

}

sub set_assembly_seq {
    my ($assembly_seq) = @_;
    $assembly_seq_ref = \$assembly_seq;
}

sub get_assembly_seq_ref {
    return ($assembly_seq_ref);
}

####
sub set_selected_gene {
    $selected_gene = shift;
}


####
sub get_selected_gene {
    return ($selected_gene);
}


####
sub get_feature_ids_to_genes {
    return (keys %id_to_gene_obj);
}


####
sub get_gene_via_feature_id {
    my $feature_id = shift;
    my $gene_struct = $id_to_gene_obj{$feature_id};
    return ($gene_struct->{gene_obj});
}

####
sub get_ordered_feature_ids {
    my @feature_ids = sort {$id_to_gene_obj{$a}->{mid_pt} <=> $id_to_gene_obj{$b}->{mid_pt}} keys %id_to_gene_obj;
    return (@feature_ids);
}

####
sub clear_contents {
    %id_to_gene_obj = ();
    $selected_gene = 0;
    $assembly_seq_ref = 0;
}

    
1; # end of package.
