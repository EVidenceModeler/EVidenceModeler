#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Gene_obj;
use GFF3_utils;
use CdbTools;

my $usage = "usage: $0 gff3_file genome.fasta\n\n";

my $gff3_file = $ARGV[0] or die $usage;
my $genome_fasta = $ARGV[1] or die $usage;


my $gene_obj_indexer_href = {};
my $contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer_href);


open (my $genome_fh, ">genome.dna") or die $!;
open (my $coords_fh, ">genome.ann") or die $!; 


foreach my $asmbl_id (sort keys %$contig_to_gene_list_href) {
    
    my $genome_seq = &cdbyank($asmbl_id, $genome_fasta);
    
    ## write sequence entry
	print $genome_fh "$genome_seq\n";
    
    my @genes = @{$contig_to_gene_list_href->{$asmbl_id}};
    
    print $coords_fh ">$asmbl_id\n";

    foreach my $gene (@genes) {
                
        my $gene_obj = $gene_obj_indexer_href->{$gene};
		$gene_obj->trim_UTRs();
		
        my $model_feat_name = $gene_obj->{Model_feat_name};
        
        my @exons = $gene_obj->get_exons();
        my $num_exons = scalar (@exons);
        
        if ($num_exons == 1) {
            my ($end5, $end3) = $exons[0]->get_coords();
            print $coords_fh "Esngl\t$end5\t$end3\t$model_feat_name\n";
            

        } else {
            
            for (my $i=0; $i < $num_exons; $i++) {
                
                my $type = "Exon";
                
                my $exon = $exons[$i];
                my ($end5, $end3) = $exon->get_coords();
                
                if ($i==0) {
                    $type = "Einit";
                }
                if ($i == $num_exons - 1) { 
                    $type = "Eterm";
                }

                print $coords_fh "$type\t$end5\t$end3\t$model_feat_name\n";
                
            }
        }
        
        
    }

}


exit(0);


