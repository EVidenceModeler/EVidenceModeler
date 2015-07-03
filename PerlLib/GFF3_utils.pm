#!/usr/local/bin/perl

package main;
our $SEE;


package GFF3_utils;

use strict;
use warnings;
use Gene_obj;
use Gene_obj_indexer;
use Carp;


####
sub index_GFF3_gene_objs {
    
    my ($gff_filename, $gene_obj_indexer, $contig_id) = @_;
    # contig_id is optional.
    

    my $hash_mode = 0;
    if (ref $gene_obj_indexer eq 'HASH') {
        $hash_mode = 1;
    }
    
    ## note can use either a gene_obj_indexer or a hash reference.
    
    my %gene_coords;
    my %asmbl_id_to_gene_id_list;
    my %transcript_to_gene;
    my %cds_phases;
    
    my %gene_names;

    open (my $fh, $gff_filename) or die $!;

    my %gene_id_to_source_type;
    
    my $counter = 0;
    # print STDERR "\n-parsing file $gff_filename\n";
    while (<$fh>) {
        chomp;
        
        unless (/\w/) { next;} # empty line
        
        if (/^\#/) { next; } # comment entry in gff3

        my @x = split (/\t/);
        my ($asmbl_id, $source, $feat_type, $lend, $rend, $orient, $cds_phase, $gene_info) = ($x[0], $x[1], $x[2], $x[3], $x[4], $x[6], $x[7], $x[8]);    
        
        if ($contig_id && $asmbl_id ne $contig_id) { next; }

        
        
        unless ($feat_type =~ /^(gene|mRNA|CDS|exon)$/) { next;} ## these are the only fields I care about right now.
        
        $gene_info =~ /ID=([^;]+);?/;
        my $id = $1 or die "Error, couldn't get the id field $_";
        
        $id = "$source$;$id";

        if ($feat_type eq 'gene') {
            my $gene_name = "";
            if ($gene_info =~ /Name=\"([^\"]+)/) {
                $gene_name = $1;
            }
            $gene_names{$id} = $gene_name;
        }
        
        if ($feat_type eq 'gene') { next;} ## beyond this pt, gene is not needed.
        
        $gene_info =~ /Parent=([^;]+);?/;
        my $parent = $1 or die "Error, couldn't get the parent info $_";
        $parent = "$source$;$parent";
        
        # print "id: $id, parent: $parent\n";
        
        if ($feat_type eq 'mRNA') {
            ## just get the identifier info
            $transcript_to_gene{$id} = $parent;
            next;
        }
        
        my $transcript_id = $parent;
        my $gene_id = $transcript_to_gene{$transcript_id} or die "Error, no gene_id for $parent";
        
        $gene_id_to_source_type{$gene_id} = $source;

        my ($end5, $end3) = ($orient eq '+') ? ($lend, $rend) : ($rend, $lend);
        
        $gene_coords{$asmbl_id}->{$gene_id}->{$transcript_id}->{$feat_type}->{$end5} = $end3;
        # print "$asmbl_id, $gene_id, $transcript_id, $feat_type, $end5, $end3\n";
        
        if ($cds_phase =~ /^\d+$/) {
            $cds_phases{$gene_id}->{$transcript_id}->{$end5} = $cds_phase;
        }
        
    }
    close $fh;
    
    ## 
    # print STDERR "\n-caching genes.\n";
    foreach my $asmbl_id (sort keys %gene_coords) {
        my $genes_href = $gene_coords{$asmbl_id};
        foreach my $gene_id (keys %$genes_href) {
                                          
            my ($gene_src, $real_gene_id) = split (/$;/, $gene_id);
            
            my $transcripts_href = $genes_href->{$gene_id};
            
            my @gene_objs;
            
            foreach my $transcript_id (keys %$transcripts_href) {
                
                my ($transcript_source, $real_transcript_id) = split (/$;/, $transcript_id);

                my $cds_coords_href = $transcripts_href->{$transcript_id}->{CDS};
                my $exon_coords_href = $transcripts_href->{$transcript_id}->{exon};
                
                unless (ref $cds_coords_href && ref $exon_coords_href) {
                    use Data::Dumper;
					print STDERR Dumper ($transcripts_href);
                    die "Error, missing cds or exon coords for $transcript_id, $gene_id\n";
                }
                
                my $gene_obj = new Gene_obj();
                
                $gene_obj->populate_gene_obj($cds_coords_href, $exon_coords_href);
                
                $gene_obj->{Model_feat_name} = $real_transcript_id;
                $gene_obj->{TU_feat_name} = $real_gene_id;
                $gene_obj->{asmbl_id} = $asmbl_id;
                
                $gene_obj->{com_name} = $gene_names{$gene_id} || $real_transcript_id;
        
                $gene_obj->{source} = $gene_id_to_source_type{$gene_id};
                
                ## set CDS phase info if available from the gff
                my $cds_phases_href = $cds_phases{$gene_id}->{$transcript_id};
                if (ref $cds_phases_href) {
                    ## set the cds phases
                    my @exons = $gene_obj->get_exons();
                    foreach my $exon (@exons) {
                        if (my $cds = $exon->get_CDS_obj()) {
                            my ($end5, $end3) = $cds->get_coords();
                            my $phase = $cds_phases_href->{$end5};
                            unless ($phase =~ /\d+/) {
                                confess "Error, should have phase set for cds $gene_id $transcript_id $end5, but I do not. ";
                            }
                            $cds->set_phase($phase);
                        }
                    }
                }
                        
                push (@gene_objs, $gene_obj);
            }
            
            ## want single gene that includes all alt splice variants here
            my $template_gene_obj = shift @gene_objs;
            foreach my $other_gene_obj (@gene_objs) {
                $template_gene_obj->add_isoform($other_gene_obj);
            }
            
            if ($hash_mode) {
                $gene_obj_indexer->{$gene_id} = $template_gene_obj;
            }
            else {
                $gene_obj_indexer->store_gene($gene_id, $template_gene_obj);
            }
            
            print "GFF3_utils: stored $gene_id\n" if $SEE;
            
            # add to gene list for asmbl_id
            my $gene_list_aref = $asmbl_id_to_gene_id_list{$asmbl_id};
            unless (ref $gene_list_aref) {
                $gene_list_aref = $asmbl_id_to_gene_id_list{$asmbl_id} = [];
            }
            push (@$gene_list_aref, $gene_id);
        }
    }
    return (\%asmbl_id_to_gene_id_list);
}


1; #EOM
