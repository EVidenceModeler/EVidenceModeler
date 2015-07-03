package main;
our $SEE;


package Accuracy_Eval;

use strict;
use warnings;
use Gene_obj_comparator;

sub new {
    my $packagename = shift;

    my $self = {
        
        Nuc_TP => 0,
        Nuc_FP => 0,
        Nuc_TN => 0,
        Nuc_FN => 0,
        
        Exon_TP => 0,
        Exon_FP => 0,
        Exon_FN => 0,
        
        Transcript_TP => 0,
        Transcript_FP => 0,
        Transcript_FN => 0,
        
        Gene_TP => 0,
        Gene_FP => 0,
        Gene_FN => 0,
    
        ## Gene and Transcript Counts.
        Ref_Gene_Count => 0,
        Ref_Isoform_Count => 0,

        Pred_Gene_Count => 0,
        Pred_Isoform_Count => 0,
        
    };
    
    bless ($self, $packagename);
    
    return ($self);
    
}


####
sub add_entry {
    my $self = shift;
    my ($seq_start, $seq_end, $reference_gene_objs_aref, $other_gene_objs_aref) = @_;

    ################################################################################################
    ## Note Gene objects **MUST** be bundled into isoform clusters.  use $gene_obj->add_isoform() ##
    ################################################################################################
    my $num_input_ref_genes = scalar @$reference_gene_objs_aref;
    my $num_input_other_genes = scalar @$other_gene_objs_aref;
        
    $self->analyze_nucleotide_SnSp (@_);

    my @ref_objs_in_range = &retrieve_overlapping_features($seq_start, $seq_end, $reference_gene_objs_aref);
    my @other_gene_objs_in_range = &retrieve_overlapping_features($seq_start, $seq_end, $other_gene_objs_aref);
    
    my $num_ref_genes_in_range = scalar (@ref_objs_in_range);
    my $num_other_gene_objs_in_range = scalar (@other_gene_objs_in_range);
    
    if ($num_input_other_genes != $num_other_gene_objs_in_range) {
        die "Error, input $num_input_other_genes, filtered to $num_other_gene_objs_in_range!\n";
    }
    
    $self->analyze_exon_SnSp ($seq_start, $seq_end, \@ref_objs_in_range, \@other_gene_objs_in_range); # really just analyzing the coding segments (CDSs)
    
    $self->analyze_gene_and_transcript_SnSp ($seq_start, $seq_end, \@ref_objs_in_range, \@other_gene_objs_in_range);
    
    return;
}

sub get_summary_stats {
    my $self = shift;
    
    my ($Nuc_Sn, $Nuc_Sp) = $self->compute_nucleotide_SnSp();
    my ($Exon_Sn, $Exon_Sp) = $self->compute_exon_SnSp();
    my ($Transcript_Sn, $Transcript_Sp) = $self->compute_transcript_SnSp();
    my ($Gene_Sn, $Gene_Sp) = $self->compute_gene_SnSp();

    my $text = "Sn,Sp stats:\n-----\n";
    
    $text .= "Gene Sensitivity           " . sprintf ("%.2f", $Gene_Sn * 100) . "%\n";
    $text .= "Gene Specificity           " . sprintf ("%.2f", $Gene_Sp * 100) . "%\n\n";

    $text .= "Transcript Sensitivity     " . sprintf ("%.2f", $Transcript_Sn * 100) . "%\n";
    $text .= "Transcript Specificity     " . sprintf ("%.2f", $Transcript_Sp * 100) . "%\n\n";

    $text .= "Exon Sensitivity           " . sprintf ("%.2f", $Exon_Sn * 100) . "%\n";
    $text .= "Exon Specificity           " . sprintf ("%.2f", $Exon_Sp * 100) . "%\n\n";

    $text .= "Nucleotide Sensitivity     " . sprintf ("%.2f", $Nuc_Sn * 100) . "%\n";
    $text .= "Nucleotide Specificity     " . sprintf ("%.2f", $Nuc_Sp * 100) . "%\n\n";
    
    $text .= "\nFeature counts:\n------\n";
    
    
    my $R = "Reference";
    my $P = "Predicted";
    

    ## Note Gene and Isoform Counts are Actual: those w/ same CDS structures are still counted independently. (UTR diffs are ignored).
    ## gene counts.
    my $num_ref_genes = $self->{Ref_Gene_Count};
    my $num_pred_genes = $self->{Pred_Gene_Count};
    
    # genes.
    $text .= "$R genes: $num_ref_genes\t$P genes: $num_pred_genes\n";
    
    # isoforms:
    my $num_ref_isoforms = $self->{Ref_Isoform_Count};
    my $num_pred_isoforms = $self->{Pred_Isoform_Count};
    
    my $ref_isoforms_per_gene = sprintf ("%.2f", $num_ref_isoforms / $num_ref_genes);
    my $pred_isoforms_per_gene = sprintf ("%.2f", $num_pred_isoforms / $num_pred_genes);

    $text .= "$R transcripts: $num_ref_isoforms ($ref_isoforms_per_gene)\t$P transcripts: $num_pred_isoforms ($pred_isoforms_per_gene)\n";
    
    ## Note, Exon counts are distinct.  The same exon shared by multiple genes or isoforms are counted only once.
    # exons
    my $num_ref_exons = $self->{Exon_TP} + $self->{Exon_FN};
    my $num_pred_exons = $self->{Exon_TP} + $self->{Exon_FP};
    
    $text .= "$R exons: $num_ref_exons\t$P exons: $num_pred_exons\n";
    
    return ($text);
}


####
sub compute_nucleotide_SnSp {
    my $self = shift;
    
    my $TP = $self->{Nuc_TP};
    my $FP = $self->{Nuc_FP};
    my $FN = $self->{Nuc_FN};
    
    my ($Sn, $Sp) = &sensitivity_and_specificity($TP, $FP, $FN);
    
    return ($Sn, $Sp);
}


####
sub compute_exon_SnSp {
    my $self = shift;
    
    my $TP = $self->{Exon_TP};
    my $FP = $self->{Exon_FP};
    my $FN = $self->{Exon_FN};
    
    my ($Sn, $Sp) = &sensitivity_and_specificity($TP, $FP, $FN);
    
    return ($Sn, $Sp);
}


####
sub compute_transcript_SnSp {
    my $self = shift;
    
    my $TP = $self->{Transcript_TP};
    my $FP = $self->{Transcript_FP};
    my $FN = $self->{Transcript_FN};
    
    my ($Sn, $Sp) = &sensitivity_and_specificity($TP, $FP, $FN);
    
    return ($Sn, $Sp);
}

####
sub compute_gene_SnSp {
    my $self = shift;
    
    my $TP = $self->{Gene_TP};
    my $FP = $self->{Gene_FP};
    my $FN = $self->{Gene_FN};
    
    my ($Sn, $Sp) = &sensitivity_and_specificity($TP, $FP, $FN);
    
    return ($Sn, $Sp);
}




###################################################################################################################
#  pseudo-private methods below.  
###################################################################################################################

####
sub sensitivity_and_specificity {
    my ($TP, $FP, $FN) = @_;

    my $sensitivity = $TP / ($TP + $FN);
    
    my $specificity = $TP / ($TP + $FP);

    return ($sensitivity, $specificity);
}

####
sub analyze_gene_and_transcript_SnSp {
    my $self = shift;
    my ($seq_start, $seq_end, $ref_objs_in_range_aref, $other_gene_objs_in_range_aref) = @_;
    
    
    my $GENE_TP = 0;
    my $GENE_FP = 0;
    my $GENE_FN = 0;
    
    my $TRANSCRIPT_TP = 0;
    my $TRANSCRIPT_FP = 0;
    my $TRANSCRIPT_FN = 0;
    
    
    ## first do some counting
    my $num_ref_genes = 0;
    my $num_ref_isoforms = 0;


    ## assign each gene and isoform a temporary hidden unique identifier:
    foreach my $ref_gene (@$ref_objs_in_range_aref) {
        $num_ref_genes++;
        $ref_gene->{_tmp_gene_ID} = "Rg_$num_ref_genes";
        foreach my $ref_isoform ($ref_gene, $ref_gene->get_additional_isoforms()) {
            $num_ref_isoforms++;
            $ref_isoform->{_tmp_isoform_ID} = "Ri_$num_ref_isoforms";
        }
    }
    $self->{Ref_Gene_Count} += $num_ref_genes;
    $self->{Ref_Isoform_Count} += $num_ref_isoforms;
        
    my $num_other_genes = 0;
    my $num_other_isoforms = 0;
    foreach my $other_gene (@$other_gene_objs_in_range_aref) {
        $num_other_genes++;
        $other_gene->{_tmp_gene_ID} = "Og_$num_other_genes";
        foreach my $other_isoform ($other_gene, $other_gene->get_additional_isoforms()) {
            $num_other_isoforms++;
            $other_isoform->{_tmp_isoform_ID} = "Oi_$num_other_isoforms";
        }
    }
    $self->{Pred_Gene_Count} += $num_other_genes;
    $self->{Pred_Isoform_Count} += $num_other_isoforms;
    

    ## loop through and find matching transcripts:
    
    my %ref_matching_genes;
    my %ref_matching_isoforms;
    my %other_matching_genes;
    my %other_matching_isoforms;
    

    ## All-vs-all comparison of genes and isoforms:
    foreach my $ref_gene_obj (@$ref_objs_in_range_aref) {
        my $ref_tmp_gene_ID = $ref_gene_obj->{_tmp_gene_ID};

        foreach my $ref_isoform ($ref_gene_obj, $ref_gene_obj->get_additional_isoforms()) {
            my $ref_tmp_isoform_ID = $ref_isoform->{_tmp_isoform_ID};
            
            foreach my $other_gene_obj (@$other_gene_objs_in_range_aref) {
                my $other_tmp_gene_ID = $other_gene_obj->{_tmp_gene_ID};
                
                foreach my $other_isoform ($other_gene_obj, $other_gene_obj->get_additional_isoforms()) {
                    my $other_tmp_isoform_ID = $other_isoform->{_tmp_isoform_ID};
                    
                    compare_genes($ref_isoform, $other_isoform);
                    
                    if (are_CDS_same()) {
                        $ref_matching_genes{$ref_tmp_gene_ID} = 1;
                        $other_matching_genes{$other_tmp_gene_ID} = 1;
                        $ref_matching_isoforms{$ref_tmp_isoform_ID} = 1;
                        $other_matching_isoforms{$other_tmp_isoform_ID} = 1;
                    }
                }
                
            }
        }
    }
    
    my $num_matching_ref_genes = scalar (keys %ref_matching_genes);
    my $num_matching_ref_isoforms = scalar (keys %ref_matching_isoforms);
    my $num_matching_other_genes = scalar (keys %other_matching_genes);
    my $num_matching_other_isoforms = scalar (keys %other_matching_isoforms);
    
    if ($SEE) {
        print "\tnum_matching_ref_genes: $num_matching_ref_genes / $num_ref_genes\n"
            . "\tnum_matching_ref_isoforms: $num_matching_ref_isoforms / $num_ref_isoforms\n"
            . "\tnum_matching_other_genes: $num_matching_other_genes / $num_other_genes\n"
            . "\tnum_matching_other_isoforms: $num_matching_other_isoforms / $num_other_isoforms\n";
    }
    
    $GENE_TP = $num_matching_ref_genes;
    $GENE_FN = $num_ref_genes - $GENE_TP;
    $GENE_FP = $num_other_genes - $num_matching_other_genes;
    if ($GENE_FP < 0 || $GENE_FN < 0 || $GENE_FP < 0) { 
        die "Error, Gene TP ($GENE_TP) or FN ($GENE_FN) || FP (GENE_FP) are negative. ";
    }
    
    $TRANSCRIPT_TP = $num_matching_ref_isoforms;
    $TRANSCRIPT_FN = $num_ref_isoforms - $TRANSCRIPT_TP;
    $TRANSCRIPT_FP = $num_other_isoforms - $num_matching_other_isoforms;
    
    if ($TRANSCRIPT_TP < 0 || $TRANSCRIPT_FN < 0 || $TRANSCRIPT_FP < 0) {
        die "Error, Transcript TP ($TRANSCRIPT_TP) or FN ($TRANSCRIPT_FN) or FP ($TRANSCRIPT_FP) are negative ";
    }
    

    $self->{Transcript_TP} += $TRANSCRIPT_TP;
    $self->{Transcript_FP} += $TRANSCRIPT_FP;
    $self->{Transcript_FN} += $TRANSCRIPT_FN;

    $self->{Gene_TP} += $GENE_TP;
    $self->{Gene_FP} += $GENE_FP;
    $self->{Gene_FN} += $GENE_FN;
    
    return;
}


sub analyze_nucleotide_SnSp { 
    my $self = shift;
    my ($seq_start, $seq_end, $reference_gene_objs_aref, $other_gene_objs_aref) = @_;
        
    my $bp_delta = $seq_start;
    my $array_length = $seq_end - $bp_delta; # starting at zero.
    ## init array:
    my @ref_coding_array;
    my @other_coding_array;
    # prealloc
    $#ref_coding_array = $#other_coding_array = $array_length;
    # init
    for (my $i = 0; $i <= $array_length; $i++) {
        $ref_coding_array[$i] = $other_coding_array[$i] = 0;
    }
    
    $self->_populate_coding_array(\@ref_coding_array, $bp_delta, $array_length, $reference_gene_objs_aref);
    $self->_populate_coding_array(\@other_coding_array, $bp_delta, $array_length, $other_gene_objs_aref);


    my $TP = 0;
    my $TN = 0;
    my $FP = 0;
    my $FN = 0;

    for (my $i = 0; $i <= $array_length; $i++) {
        my $ref_val = $ref_coding_array[$i];
        my $other_val = $other_coding_array[$i];
        
        if ($ref_val && $other_val) {
            $TP++;
        }
        elsif ($ref_val && (! $other_val)) {
            $FN++;
        }
        elsif ( (! $ref_val) && $other_val) {
            $FP++;
        }
        elsif ( (! $ref_val) && (! $other_val)) {
            $TN++;
        }
        else {
            die "Error, should never end here.";
        }
    }

    $self->{Nuc_TP} += $TP;
    $self->{Nuc_TN} += $TN;
    $self->{Nuc_FP} += $FP;
    $self->{Nuc_FN} += $FN;

    return;
}

####
sub _populate_coding_array {
    my $self = shift;
    my ($coding_array_aref, $bp_delta, $array_length, $gene_objs_aref) = @_;

    foreach my $gene_obj (@$gene_objs_aref) {
        
        foreach my $isoform ($gene_obj, $gene_obj->get_additional_isoforms()) {

            foreach my $exon ($isoform->get_exons()) {

                if (my $cds = $exon->get_CDS_exon_obj()) {
                    
                    my ($lend, $rend) = sort {$a<=>$b} $cds->get_coords();
                    
                    for (my $i = $lend; $i <= $rend; $i++) {
                        
                        my $pos = $i - $bp_delta;
                        if ($pos >= 0 && $pos <= $array_length) {
                            $coding_array_aref->[$pos] = 1;
                        }
                    }
                }
            }
        }
    }

    return;
}


####
sub analyze_transcript_SnSp {
    my $self = shift;
    my ($seq_start, $seq_end, $reference_gene_objs_aref, $other_gene_objs_aref) = @_;
    
    


}


####
sub analyze_exon_SnSp {
    my $self = shift;
    my ($seq_start, $seq_end, $reference_gene_objs_in_range_aref, $other_gene_objs_in_range_aref) = @_;

    my @ref_objs_in_range = @$reference_gene_objs_in_range_aref;
    my @other_gene_objs_in_range = @$other_gene_objs_in_range_aref;
    
    my %ref_exons; # end5_end3 as key.
    foreach my $ref_gene (@ref_objs_in_range) {
        foreach my $ref_isoform ($ref_gene, $ref_gene->get_additional_isoforms()) {
            my @exons = $ref_isoform->get_exons();
            foreach my $exon (@exons) {
                if (my $cds = $exon->get_CDS_exon_obj()) {
                    my ($cds_lend, $cds_rend) = sort {$a<=>$b} $cds->get_coords();
                    if ($cds_lend < $seq_end && $cds_rend > $seq_start) {
                        
                        my $cds_key = join ("_", $cds->get_coords());
                        $ref_exons{$cds_key} = 1;
                    }
                }
            }
        }
    }

    my %other_exons; # end5_end3 as key.
    foreach my $other_gene (@other_gene_objs_in_range) {
        foreach my $other_isoform ($other_gene, $other_gene->get_additional_isoforms()) {
            my @exons = $other_isoform->get_exons();
            foreach my $exon (@exons) {
                if (my $cds = $exon->get_CDS_exon_obj()) {
                    my ($cds_lend, $cds_rend) = sort {$a<=>$b} $cds->get_coords();
                    if ($cds_lend < $seq_end && $cds_rend > $seq_start) {
                        
                        my $cds_key = join ("_", $cds->get_coords());
                        $other_exons{$cds_key} = 1;
                    }
                }
            }
        }
    }

    
    my $TP = 0;
    my $FP = 0;
    my $FN = 0;
    
    foreach my $exon_key (keys %ref_exons) {
        if ($other_exons{$exon_key}) {
            $TP++;
            delete $other_exons{$exon_key};
        }
        else {
            $FN++;
        }
    }
    foreach my $remaining_other_exon (keys %other_exons) {
        $FP++;
    }
    
    $self->{Exon_TP} += $TP;
    $self->{Exon_FP} += $FP;
    $self->{Exon_FN} += $FN;
    
    return;
}


####
sub retrieve_overlapping_features {
    my ($seq_start, $seq_end, $features_aref) = @_;
    ## Note feature must have method ->get_coords()

    my @feats_in_range;

    foreach my $feature (@$features_aref) {

        my ($feature_lend, $feature_rend) = sort {$a<=>$b} $feature->get_coords();

        if ($feature_lend < $seq_end && $feature_rend > $seq_start) {
            ## overlaps range > 1 bp.
            push (@feats_in_range, $feature);
        }
    }

    return (@feats_in_range);

}


    
1; #EOM
