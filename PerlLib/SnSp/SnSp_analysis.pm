#!/usr/local/bin/perl

=documentation


my $snsp_analyzer = new SnSp::SnSp_analysis_manager("genscan", "genemarkHMM", "fgenesh");


my %ev_type__to_gene_list = (    'genscan' => [ $geneA ],
                                 'genemarkHMM' => [ $geneB, $geneC], #predicted two Gene_obj's here where there's really one
                                 'fgenesh' => [ $geneD ]
                                 );


$snsp_analyzer->add_analysis_entry($gold_standard_gene_obj,  \%ev_type_to_gene_list, 1); # last param is for verbose descriptions



=cut




package main;
our $SEE;

package SnSp::SnSp_analysis;
use strict;
use warnings;
use Carp;

sub new {
    my $packagename = shift;
    

    confess "This class is abstract, never to be instantiated directly";
    
}

sub _init {
    my $self = shift;
    my @prediction_names = @_;

    unless (@prediction_names) {
        confess "Error, need prediction names for initialization.";
    }

    $self->{AP} = 0;
    $self->{AN} = 0;
    $self->{prediction_to_TFPN} = {};
    $self->{prediction_names} = [@prediction_names];
    $self->{number_genes_analyzed} = 0;
    $self->{number_exons_analyzed} = 0;
}


sub get_prediction_names {
    my $self = shift;
    return (@{$self->{prediction_names}});
}


sub get_TFPN {
    my $self = shift;
    my $ev_type = shift;

    return ($self->{prediction_to_TFPN}->{$ev_type});
}


sub summarize_SnSp {
    my $self = shift;
        
    my ($verbose_summary) = @_;
    
    my $number_genes_analyzed = $self->{number_genes_analyzed};
    my $number_exons_analyzed = $self->{number_exons_analyzed};

    my $AP = $self->{AP};
    my $AN = $self->{AN};

    my @prediction_names = $self->get_prediction_names();
    
    foreach my $ev_type (@prediction_names) {
        my $TFPN_obj = $self->get_TFPN($ev_type);
        
        my $TP = $TFPN_obj->{TP};
        my $TN = $TFPN_obj->{TN};
        my $FP = $TFPN_obj->{FP};
        my $FN = $TFPN_obj->{FN};
        
        my $PP = $TFPN_obj->{PP};
        my $PN = $TFPN_obj->{PN};
        
        my $number_predictions = $TFPN_obj->{number_predicted_genes};
        
        my $number_exons_correct = $TFPN_obj->{number_exons_correct};
        my $number_correct_genes = $TFPN_obj->{genes_predicted_correct};
        
        
        # Sn = TP/(TP+FN)
        my $sensitivity_val = $TP / ($TP + $FN);
        $TFPN_obj->{sensitivity} = $sensitivity_val;
        
        # Sp = TP/(TP+FP)
        my $specificity_val = 0;
        if ($TP || $FP) {
            $specificity_val = $TP / ($TP + $FP);
        }
       
        $TFPN_obj->{specificity} = $specificity_val;
        
        ## Correlation coeff
        my $correl_coeff = ($TP * $TN + $FP * $FN) * ( ($PP * $PN * $AP * $AN) ** (-0.5));
        $TFPN_obj->{correlation_coeff} = $correl_coeff;
        
        ## Fscore
        my $Fscore = 0;

        if ($sensitivity_val > 0 && $specificity_val > 0) {
            $Fscore = (2 * $sensitivity_val * $specificity_val) / ($sensitivity_val + $specificity_val);
            $TFPN_obj->{Fscore} = $Fscore;
        }
        
        if ($verbose_summary) {
            printf("EvType $ev_type ($number_predictions preds)\tNSn: %.2f\tNSp: %.2f\tCC: %.2f\tCexons: %d=%.2f\tCpred: %d=%.2f\tFscore: %.2f\n",  
                   $sensitivity_val, $specificity_val, $correl_coeff, 
                   $number_exons_correct, $number_exons_correct/$number_exons_analyzed*100, 
                   $number_correct_genes, $number_correct_genes/$number_genes_analyzed*100,
                   $Fscore*100);
        }
    } 
}



1; #EOM
