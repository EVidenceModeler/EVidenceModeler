#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use Gene_obj;
use GFF3_utils;
use CdbTools;


my $usage = "usage: $0 gff3_file genome.fasta\n\n";

my $gff3_file = $ARGV[0] or die $usage;
my $genome_fasta = $ARGV[1] or die $usage;


my $gene_obj_indexer_href = {};
my $contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer_href);

my $core_filename = "glim.$$";

open (my $genome_fh, ">$core_filename.genome") or die $!;
open (my $coords_fh, ">$core_filename.coords") or die $!; 


foreach my $asmbl_id (keys %$contig_to_gene_list_href) {
    
    my $genome_seq = &cdbyank_linear($asmbl_id, $genome_fasta);
    
    print $genome_fh ">$asmbl_id\n$genome_seq\n";
    
	foreach my $gene_id (@{$contig_to_gene_list_href->{$asmbl_id}}) {
		        
        my $gene_obj = $gene_obj_indexer_href->{$gene_id};;
                
        $gene_obj->create_all_sequence_types(\$genome_seq);
        $gene_obj->trim_UTRs();

		
		my $protein = $gene_obj->get_protein_sequence();
        my ($is_5prime_partial, $is_3prime_partial) = (0,0);
        if ($protein !~ /^M/) {
            $is_5prime_partial = 0;
        }
        if ($protein !~ /\*$/) {
            $is_3prime_partial = 0;
        }

        my @exons = $gene_obj->get_exons();
        my $num_exons = scalar (@exons);
        
        for (my $i=0; $i < $num_exons; $i++) {
            
            my $exon = $exons[$i];
            my ($end5, $end3) = $exon->get_coords();
            
            if ($i==0 && $is_5prime_partial) {
                $end5 = "<$end5";
            }
            if ($i == $num_exons - 1 && $is_3prime_partial) {
                $end3 = ">$end3";
            }
            
            print $coords_fh "$asmbl_id\t$end5\t$end3\n";
        }
        print $coords_fh "\n";
        
    }

}


exit(0);

