
package SnSp::SnSp_analysis_entry;
use strict;
use warnings;
use base qw (SnSp::SnSp_analysis);
use Carp;
use Gene_obj_comparator;


sub new {
    my $packagename = shift;
    my $SnSp_analysis_manager = shift;
    unless (ref $SnSp_analysis_manager eq "SnSp::SnSp_analysis_manager") {
        confess "wrong args.";
    }
    
    my $template_gene_obj = shift;
    
    my $others_hashref = shift;

    my $self = {

        template_gene => $template_gene_obj,
        other_predictions_href => $others_hashref,
        
        manager => $SnSp_analysis_manager, 
    };

    bless ($self, $packagename);


    $self->SUPER::_init($SnSp_analysis_manager->get_prediction_names());
    $self->_init();
    
    return ($self);
}


sub get_gene_predictions {
    my $self = shift;
    my $ev_type = shift;

    my $list_ref = $self->{other_predictions_href}->{$ev_type};
    if (ref $list_ref eq "ARRAY") {
        return (@$list_ref);
    } 
    else {
        ## no predictions here
        confess "Error, no prediction for $ev_type.  If this is true, present an empty array ref.\n";
        return ();
    }
}


sub _init {
    my $self = shift;
    
    foreach my $ev_type ($self->{manager}->get_prediction_names()) {
        $self->{prediction_to_TFPN}->{$ev_type} = SnSp::TFPN->new();
    }
}

sub calc_SnSp {
    my $self = shift;

    $self->{number_genes_analyzed}++;
    
    my $template_gene_obj = $self->{template_gene};
    my $manager = $self->{manager};

    my $intergenic_included = $manager->get_intergenic_included();
    
    
    ## Process the Template Gene:
    my @gold_std_array;
    
    my ($model_lend, $model_rend) = sort {$a<=>$b} $template_gene_obj->get_model_span();
    
    my ($from_pos, $to_pos) = ($model_lend - $intergenic_included -1, $model_rend + $intergenic_included);
    unless ($from_pos > 0) { 
        # shorter upstream end included since close to seq boundary.
        $from_pos = 1;
    }
        
    ## init gold std array
    for (my $i = 0; $i <= $to_pos-$from_pos; $i++) {
        $gold_std_array[$i] = 0;
    }
    
    ## populate gold_std_array

    my %gold_exons;
    
    my @exons = $template_gene_obj->get_exons();
    foreach my $exon (@exons) {
        my @cds_coords = sort {$a<=>$b} $exon->get_CDS_end5_end3();
        if (@cds_coords) {
            $gold_exons { "$cds_coords[0]" . "_" . "$cds_coords[1]" } = 1; 
            $self->{number_exons_analyzed}++;
            for (my $i = $cds_coords[0]; $i <= $cds_coords[1]; $i++) {
                $gold_std_array[$i-$from_pos] = 1;
            }
        }
    }
    
    ## Score the actual positives and actual negatives:
    for (my $i = 0; $i <= $to_pos-$from_pos; $i++) {
        if ($gold_std_array[$i]) {
            $self->{AP}++;
        } else {
            $self->{AN}++;
        }
    }


    ## Process each of the gene structures:
    
    
    foreach my $ev_type ($self->{manager}->get_prediction_names()) {
        my $same_gene_structure_flag = 0;
        
        my %predicted_exons;
        
        my @prediction_array;
        ## init prediction array
        for (my $i=0; $i <= $to_pos-$from_pos; $i++) {
            $prediction_array[$i] = 0;
        }
        
        my @overlapping_predictions = $self->get_gene_predictions($ev_type);

        foreach my $predicted_gene_obj (@overlapping_predictions) {
            
            compare_genes($template_gene_obj, $predicted_gene_obj);
            if (are_CDS_same()) {
                $same_gene_structure_flag = 1;
            }
            
            my @exons = $predicted_gene_obj->get_exons();
            foreach my $exon (@exons) {
                my @cds_coords = sort {$a<=>$b} $exon->get_CDS_end5_end3();
                if (@cds_coords) {
                    
                    $predicted_exons{ "$cds_coords[0]" . "_" . "$cds_coords[1]" } = 1;
                    
                    for (my $i = $cds_coords[0]; $i <= $cds_coords[1]; $i++) {
                        my $pos = $i - $from_pos;
                        if ($pos >= 0 && $pos <= $to_pos - $from_pos) {
                            $prediction_array[$pos] = 1;
                        }
                    }
                }
            }
        }
		
        my $number_exons_predicted_correct = 0;
        foreach my $exon_key (keys %gold_exons) {
            if ($predicted_exons{$exon_key}) {
                $number_exons_predicted_correct++;
            }
        }

        my $TFPN_obj = $self->get_TFPN($ev_type);

        ## Score the Sensitivity and Specificity parameters:
        
        
        $TFPN_obj->{number_predicted_genes} = scalar (@overlapping_predictions);
        $TFPN_obj->{number_exons_correct} = $number_exons_predicted_correct;
        if ($same_gene_structure_flag) {
            #print "$ev_type, same structure.\n";
            $TFPN_obj->{genes_predicted_correct} = 1;
        } else {
            #print "$ev_type, diff structure.\n";
        }
        
        # score true positives, false positives, and negatives:
        
        for (my $i= 0; $i <= $to_pos-$from_pos; $i++) {
            
            ## TP, FP, TN, FN
            if ($gold_std_array[$i] != 0 && $prediction_array[$i] != 0) {
                $TFPN_obj->{TP}++;
            } elsif ($gold_std_array[$i] != 0 && $prediction_array[$i] == 0) {
                $TFPN_obj->{FN}++;
            } elsif ($gold_std_array[$i] == 0 && $prediction_array[$i] != 0) {
                $TFPN_obj->{FP}++;
            } elsif ($gold_std_array[$i] == 0 && $prediction_array[$i] == 0) {
                $TFPN_obj->{TN}++;
            }
            
            ## PP, PN (predicted positive, predicted negative)
            if ($prediction_array[$i]) {
                $TFPN_obj->{PP}++;
            } else {
                $TFPN_obj->{PN}++;
            }
            
        }
    }
}



1; #EOM
