package SnSp::SnSp_analysis_manager;
use strict;
use warnings;
use Carp;
use SnSp::TFPN;
use SnSp::SnSp_analysis_entry;
use base qw (SnSp::SnSp_analysis);

sub new {
    my $packagename = shift;
    
    my @prediction_names = @_;
    unless (@prediction_names) {
        confess "wrong args";
    }

    my $self = { 
        intergenic_included => 500, #default
        
    };
    
    
    bless ($self, $packagename);
    

    $self->SUPER::_init(@prediction_names);
    $self->_init();

    return ($self);
}



sub _init {
    my $self = shift;
    
    foreach my $ev_type ($self->get_prediction_names()) {
        $self->{prediction_to_TFPN}->{$ev_type} = SnSp::TFPN->new();
    }
}


sub get_intergenic_included {
    my $self = shift;
    return ($self->{intergenic_included});
}


sub set_intergenic_included {
    my $self = shift;
    my $intergenic_included = shift;
    
    unless ($intergenic_included >= 0) {
        confess "Error, intergenic included value must be at least zero.";
    }
    
    $self->{intergenic_included} = $intergenic_included;
    
    return;
}

sub add_analysis_entry {
    my $self = shift;
    my ($template_gene_obj, $others_hashref, $verbose_flag) = @_;
    
    my $snsp_analysis_entry = SnSp::SnSp_analysis_entry->new($self, $template_gene_obj, $others_hashref);
    
    $snsp_analysis_entry->calc_SnSp();

    print "\nSingle comparison:\n" if $verbose_flag;
    $snsp_analysis_entry->summarize_SnSp($verbose_flag);
    
    
    ## sum current results to total results:
    $self->_append_analysis_results($snsp_analysis_entry);
    my $number_genes_analyzed = $self->{number_genes_analyzed};
    print "\nTotal so far ($number_genes_analyzed in training set):\n" if $verbose_flag;
    $self->summarize_SnSp($verbose_flag);
    

}


sub _append_analysis_results {
    my ($self, $analysis_ref) = @_;
    

    ## append results of individual genefinders.
    foreach my $pred_name ($self->get_prediction_names()) {
        my $global_TFPN_obj = $self->get_TFPN($pred_name);

        my $analysis_TFPN_obj = $analysis_ref->get_TFPN($pred_name);

        foreach my $att qw (TP FP TN FN PP PN 
                            number_exons_correct 
                            genes_predicted_correct 
                            number_predicted_genes) {
            
            $global_TFPN_obj->{$att} += $analysis_TFPN_obj->{$att};
        }
    }

    ## increment template tallies:
    foreach my $att qw (AP AN number_genes_analyzed number_exons_analyzed) {
        $self->{$att} += $analysis_ref->{$att};
    }
}

1; #EOM

