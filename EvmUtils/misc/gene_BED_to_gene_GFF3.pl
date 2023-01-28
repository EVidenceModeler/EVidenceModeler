#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;

use lib ("$FindBin::Bin/../../PerlLib");
use Gene_obj;



my $usage = "usage: $0 input.bed ev_type\n\n";

my $input_file = $ARGV[0] or die $usage;
my $model_type = $ARGV[1] or die $usage;

main: {

    open(my $fh, $input_file) or die "Error, cannot open file: $input_file";
    while(<$fh>) {
        chomp;
        if (/^\#/) { next; }

        my $bed_line = $_;

        my $gene_obj = Gene_obj::BED_line_to_gene_obj($bed_line);
        
        
        #$gene_obj->populate_gene_object($coords_href, $coords_href);
        #$gene_obj->{asmbl_id} = $scaffold;
        #$gene_obj->{TU_feat_name} = $gene_id;
        #$gene_obj->{Model_feat_name} = $trans_id;
        #$gene_obj->{com_name} = "$model_type prediction";
        
        #$gene_obj->join_adjacent_exons();
        
        print $gene_obj->to_GFF3_format(source => $model_type) . "\n";
        
        

    }


    exit(0);

}
        
        
