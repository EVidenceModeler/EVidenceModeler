package EVM_output_parser;

use strict;
use warnings;
use Carp;
use Gene_obj;

sub new {
    my $packagename = shift;

    my $evm_output_filename = shift;
    unless (defined $evm_output_filename) {
        confess "error, require evm_output filename as parameter";
    }
    
    unless (-e $evm_output_filename) {
        confess "Cannot find $evm_output_filename\n";
    }

    my $self = { evm_output_filename => $evm_output_filename };

    bless ($self, $packagename);

    return ($self);
}


sub get_genes {
    my $self = shift;

    my %model_coords;
    
    my $evm_file = $self->{evm_output_filename};
    my $model_id = 0;
    my %model_id_to_text;

    
    open (my $fh, $evm_file) or die "No such file $evm_file\n";
    while (<$fh>) {
        #print;
        chomp;
        if (/^\#/) {
            $model_id++;
            
            
        }
        elsif (/\w/)  {
            my @x = split (/\t/);
            $model_coords{$model_id}->{$x[0]} = $x[1];
        }
        
        if ($model_id && $_ =~ /\w/) {
            $model_id_to_text{$model_id} .= $_ . "\n";
        }
    }
    
    ## convert to gene objects
    my @gene_objs;
    foreach my $model (keys %model_coords) {
        my $coord_href = $model_coords{$model};
        my $gene_obj = new Gene_obj();
        $gene_obj->populate_gene_obj($coord_href, $coord_href);
        $gene_obj->{Model_feat_name} = $model;
        
        ## store original prediction text as the comment
        $gene_obj->{comment} = $model_id_to_text{$model};
        
        push (@gene_objs, $gene_obj);
    }
    
    close $fh;
    return (@gene_objs);
}


1; #EOM

