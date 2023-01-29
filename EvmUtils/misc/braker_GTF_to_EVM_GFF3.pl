#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;

use lib ("$FindBin::Bin/../../PerlLib");
use Gene_obj;



my $usage = "usage: $0 augustus.gff.output\n\n";

my $input_file = $ARGV[0] or die $usage;

main: {
    my %data;
    my %gene_to_source_type;
    
    my $scaffold = "";

    ## parse input file
    open (my $fh, $input_file) or die "Error, cannot open file $input_file";
    while (<$fh>) {
        if (/^\#/) { next; }
        chomp;
        unless (/\w/) { next; }
        
        my @x = split(/\t/);
        if ($x[2] eq 'CDS') { 
            my $scaffold = $x[0];
            my $source_type = $x[1];
            my $orient = $x[6];
            my $lend = $x[3];
            my $rend = $x[4];
            
            my ($end5, $end3) = ($orient eq '+') ? ($lend, $rend) : ($rend, $lend);

            my $info = $x[8];

            my $gene_id;
            my $trans_id;
            if ($info =~ /transcript_id \"([^\"]+)"/) {
                $trans_id = $1;
            }
            if ($info =~ /gene_id \"([^\"]+)"/) {
                $gene_id = $1;
            }
            if ($trans_id && $gene_id) {
                $trans_id = "$scaffold-$trans_id";
                $gene_id = "$scaffold-$gene_id";
            }
            else {
                die "Error, cannot parse gene_id and transcript_id from $_";
            }
            
            my $model = join("$;", $gene_id, $trans_id);
            
            $data{$scaffold}->{$model}->{$end5} = $end3;
            $gene_to_source_type{$model} = $source_type;
            
        }
    }
    
    close $fh;


    ## Generate gff3 output
    foreach my $scaffold (keys %data) {
        
        my $models_href = $data{$scaffold};

        my %gene_id_to_gene_obj;
        
        foreach my $model (keys %$models_href) {

            my $source_type = $gene_to_source_type{$model} or die "Error, no source type defined for $model";
            
            my ($gene_id, $trans_id) = split(/$;/, $model);

            my $coords_href = $models_href->{$model};

            my $gene_obj = new Gene_obj();
            $gene_obj->populate_gene_object($coords_href, $coords_href);
            $gene_obj->{asmbl_id} = $scaffold;
            $gene_obj->{TU_feat_name} = $gene_id;
            $gene_obj->{Model_feat_name} = $trans_id;
            $gene_obj->{com_name} = "$source_type prediction";
            $gene_obj->{source_type} = $source_type;
            
            $gene_obj->join_adjacent_exons();

            if (exists $gene_id_to_gene_obj{$gene_id}) {
                $gene_id_to_gene_obj{$gene_id}->add_isoform($gene_obj);
            }
            else {
                $gene_id_to_gene_obj{$gene_id} = $gene_obj;
            }

        }

        foreach my $gene_id (keys %gene_id_to_gene_obj) {

            my $gene_obj = $gene_id_to_gene_obj{$gene_id};

            $gene_obj->refine_gene_object();
            
            my $source_type = $gene_obj->{source_type};
            
            print $gene_obj->to_GFF3_format(source => $source_type) . "\n";

        }

    }
    

    exit(0);

}
        
        
