package SnSp::TFPN;
use strict;
use warnings;

sub new {
    my $packagename = shift;
 
    my $self = { TP => 0,  # true positive
                 FP => 0,  # false positive
                 TN => 0,  # true negative
                 FN => 0,  # false negative
                 PP => 0,  # predicted positive
                 PN => 0,  # predicted negative

                 number_exons_correct => 0,
                 genes_predicted_correct => 0,

                 number_predicted_genes => 0,

                 sensitivity => 0,
                 specificity => 0,
                 correlation_coeff => 0,
                 Fscore => 0, #  Fscore = 2SnSp/(Sn+Sp)
                 
             };
    
    bless ($self, $packagename);
    return ($self);
}

1; #EOM
