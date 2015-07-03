#!/usr/bin/env perl

use strict;
use warnings;


my %contig_to_gene_list;

main: {
    while (<STDIN>) {
        unless (/\w/) { next; }
        
        chomp;
        my @x = split (/\t/);
        
        unless ($x[2] eq 'gene') { next; }
        
        my ($contig, $lend, $rend, $orient, $id) = ($x[0], $x[3], $x[4], $x[6], $x[8]);
        
        if ($id =~ /ID=([^;]+)/) {
            $id = $1;
        }
        else {
            die "Error, couldn't extract the gene ID from $id";
        }
        
        my $gene_list_aref = $contig_to_gene_list{$contig};
        unless (ref $gene_list_aref) {
            $gene_list_aref = $contig_to_gene_list{$contig} = [];
        }

        my $gene_struct = { id => $id,
                            lend => $lend,
                            rend => $rend,
                            orient => $orient,
                        };
        
        push (@$gene_list_aref, $gene_struct);

    }
    

    foreach my $contig (keys %contig_to_gene_list) {
        
        my $gene_list_aref = $contig_to_gene_list{$contig};
        &find_nested_genes($gene_list_aref);
    }
}


####
sub find_nested_genes {
    
    my ($gene_list_aref) = @_;

    foreach my $gene_a (@$gene_list_aref) {

        my ($gene_a_id, $gene_a_lend, $gene_a_rend, $gene_a_orient) = ($gene_a->{id},
                                                                       $gene_a->{lend},
                                                                       $gene_a->{rend},
                                                                       $gene_a->{orient});
        

        foreach my $gene_b (@$gene_list_aref) {

            if ($gene_a eq $gene_b) { next; }
            
            
            my ($gene_b_id, $gene_b_lend, $gene_b_rend, $gene_b_orient) = ($gene_b->{id},
                                                                           $gene_b->{lend},
                                                                           $gene_b->{rend},
                                                                           $gene_b->{orient});

            
            if ($gene_a_id eq $gene_b_id) { next; }

            if ($gene_b_lend > $gene_a_lend && $gene_b_rend < $gene_a_rend) {
                
                print "$gene_b_id contained in $gene_a_id\n";
            }
        }



    }

    return;
}


            
    

        
