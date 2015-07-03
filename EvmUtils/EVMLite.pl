#!/usr/bin/env perl

use strict;
use warnings;


use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use Fasta_reader;
use GFF3_utils;
use Overlap_piler;
use Data::Dumper;


my $usage = "\n\n\nusage: $0 weights_file fileA.gff3 [fileB.gff3 ...]\n\n";

my $weights_file = $ARGV[0] or die $usage;
shift @ARGV;
my @gff3_files = @ARGV;

unless (@gff3_files) {
    die $usage;
}

main: {


    
    my %weights = &parse_weights_file($weights_file);

    my %scaffold_to_gene_objs = &parse_genes(@gff3_files);
    
    foreach my $asmbl_id (keys %scaffold_to_gene_objs) {
        
        my @genes = @{$scaffold_to_gene_objs{$asmbl_id}};

        my @selected_genes = &DP_select_best_genes(\%weights, \@genes);
        
        foreach my $gene (@selected_genes) {
            
            print $gene->to_GFF3_format() . "\n";
            
        }
    }
    
    exit(0);
    
}

####
sub parse_weights_file {
    my ($weights_file) = @_;

    my %weights;
    
    open (my $fh, $weights_file) or die $!;
    while (<$fh>) {
        chomp;
        if (/^\#/) { 
            next; # comment line
        }
        unless (/\w/) { next; }
        my ($class, $source, $weight) = split(/\s+/);
        
        unless (defined ($source) && defined ($weight)) {
            die "Error, improperly formed weights file";
        }

        $weights{$source} = $weight;
    }
    close $fh;

    return(%weights);

}


####
sub parse_genes {
    my @gff3_files = @_;
    
    my %scaffold_to_gene_objs;
    

    foreach my $gff3_file (@gff3_files) {

        my $gene_obj_indexer_href = {};
        
        ## associate gene identifiers with contig id's.
        my $contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer_href);
        foreach my $asmbl_id (keys %$contig_to_gene_list_href) {
            my @gene_ids = @{$contig_to_gene_list_href->{$asmbl_id}};
            
            foreach my $gene_id (@gene_ids) {
                my $gene_obj = $gene_obj_indexer_href->{$gene_id};
                push (@{$scaffold_to_gene_objs{$asmbl_id}}, $gene_obj);
            }
        }

    }

    return(%scaffold_to_gene_objs);
}
    

####
sub DP_select_best_genes {
    my ($weights_href, $genes_aref) = @_;
    
    my @genes = @$genes_aref;
    
    my @selected_genes;

    ## compute piles of genes.
    my $overlap_piler = new Overlap_piler();
    my %acc_to_gene;
    foreach my $gene (@genes) {
        my ($lend, $rend) = $gene->get_coords();
        my $id = $gene->{Model_feat_name};
        
        $overlap_piler->add_coordSet($id, $lend, $rend);
        $acc_to_gene{$id} = $gene;
    }

    my @piles = $overlap_piler->build_clusters();

    foreach my $pile_aref (@piles) {
        my @accs = @$pile_aref;

        my @gene_pile;
        foreach my $acc (@accs) {
            my $gene = $acc_to_gene{$acc} or die "Error, cannot retrieve piled gene via acc: $acc";
            push (@gene_pile, $gene);
            
        }
        
        if (scalar(@gene_pile) == 1) {
            push (@selected_genes, @gene_pile);
        }
        else {
            my @dp_selected = &run_DP_scan($weights_href, \@gene_pile);
            push (@selected_genes, @dp_selected);
        }
    }
    
    return(@selected_genes);
}

####
sub run_DP_scan {
    my ($weights_href, $gene_pile_aref) = @_;

    if (scalar(@$gene_pile_aref) < 2) {
        die "Error, got pile less than 2: " . Dumper($gene_pile_aref);
    }
    
    foreach my $gene (@$gene_pile_aref) {

        my $source = $gene->{source};
        my $weight = $weights_href->{$source};
        if (! defined $weight) {
            die "Error, no weight defined for source type: $source";
        }
        
        my $cds_len = $gene->get_CDS_length();
        my $score = $cds_len * $weight;

        $gene->{base_score} = $score;
        $gene->{sum_score} = $score;
        $gene->{prev} = undef; # point to prev gene in trellis
        
        my ($lend, $rend) = sort {$a<=>$b} $gene->get_coords();
        $gene->{lend} = $lend;
        $gene->{rend} = $rend;
        
    }

    my @genes = @$gene_pile_aref;
    @genes = sort {$a->{lend}<=>$b->{lend}} @genes;

    ## build trellis
    for (my $i = 1; $i <= $#genes; $i++) {
        for (my $j = $i-1; $j >= 0; $j--) {

            my $gene_i = $genes[$i];
            my $gene_j = $genes[$j];
            
            my ($lend_i, $rend_i) = sort {$a<=>$b} $gene_i->get_coords();
            my ($lend_j, $rend_j) = sort {$a<=>$b} $gene_j->get_coords();

            if ($rend_j < $lend_i) {
                # don't overlap
                my $sum_score = $gene_i->{base_score} + $gene_j->{sum_score};
                if ($sum_score > $gene_i->{sum_score}) {
                    ## found better linkage
                    $gene_i->{sum_score} = $sum_score;
                    $gene_i->{prev} = $gene_j;
                }
            }
        }
    }

    ## find highest scoring entry
    my $highest_scoring_gene = shift @genes;
    foreach my $gene (@genes) {
        if ($gene->{sum_score} > $highest_scoring_gene->{sum_score}) {
            $highest_scoring_gene = $gene;
        }
    }

    my @final_genes = ($highest_scoring_gene);
    while ($highest_scoring_gene->{prev}) {
        $highest_scoring_gene = $highest_scoring_gene->{prev};
        if ($highest_scoring_gene) {
            push (@final_genes, $highest_scoring_gene);
        }
    }

    return(@final_genes);
}


    
    
